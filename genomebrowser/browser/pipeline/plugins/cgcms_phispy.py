import os
import shutil
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import Gene, Genome
"""
    This plugin runs PhiSpy for a set of genomes.
"""

def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'],
                               'phispy-plugin-temp'
                               )
    script_path = preprocess(annotator, genomes, working_dir)
    run(script_path)
    output_file = postprocess(annotator, genomes, working_dir)
    return(output_file)

def preprocess(annotator, genomes, working_dir):
    """
        Creates all directories and input files. 
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
        Output:
            path of the shell script running PhiSpy in conda environment
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    input_fasta_files = {}
    for genome in genomes:
        input_fasta_files[genome] = _export_proteins(genome, working_dir)

    phispy_script = os.path.join(working_dir, 'run_phispy.sh')
    
    with open(phispy_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['cgcms.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.phispy.conda_env'] + '\n'
                      )
        for genome in sorted(genomes.keys()):
            outfile.write(' '.join(['PhiSpy.py',
                                    '-o',
                                    '"' + os.path.join(output_dir, genome) + '"',
                                    '--threads',
                                    annotator.config['plugins.phispy.threads'],
                                    '--output_choice', '9', '-k', 'KEEP',
                                    '"' + genomes[genome] + '"']) + '\n')
        outfile.write('conda deactivate\n')

        # Run HMMSearch for each protein fasta file against pVOG HMMs
        for genome in input_fasta_files:
            faa_file = input_fasta_files[genome]
            domfile = os.path.join(output_dir, genome + '.domtblout.txt')
            outfile.write(' '.join(['hmmsearch',
                  '--domtblout',
                  '"' + domfile + '"',
                  '-o',
                  '/dev/null',
                  '-E',
                  annotator.config['plugins.phispy.evalue_threshold'],
                  '--cpu',
                  annotator.config['plugins.phispy.threads'],
                  '--noali',
                  '--notextw',
                  '"' + annotator.config['plugins.phispy.pvog_path'] + '"',
                  '"' + faa_file + '"'
                  ]) + '\n')
        
    return phispy_script
    
def run(script_path):
    """
    Runs PhiSpy script for GBK files of genomes.
    """
    cmd = ['/bin/bash', script_path]
    print(' '.join(cmd))
    # Close MySQL connection before starting external process because 
    # it may run for too long resulting in "MySQL server has gone away" error
    connection.close()
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        # Suppress false positive no-member error
        # (see https://github.com/PyCQA/pylint/issues/1860)
        # pylint: disable=no-member
        raise CalledProcessError(proc.returncode, proc.args)

def postprocess(annotator, genomes, working_dir):
    """
        Finds PhiSpy output files and creates file with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['cgcms.temp_dir'],
                               'phispy-plugin-output.txt'
                               )
    with open(output_file, 'w') as outfile:
        outfile.write('#Gene\tGenome\tSoure\tURL\tKey\tValue\tNote\n')
        for genome in sorted(genomes.keys()):
            phispy_regions = os.path.join(working_dir,
                                          'out',
                                          genome,
                                          'prophage_coordinates.tsv'
                                          )
            phispy_genes = os.path.join(working_dir,
                                        'out',
                                        genome,
                                        'prophage_information.tsv'
                                        )
            if not os.path.exists(phispy_regions):
                print('No prophages found in', genome)
                continue
            genes = {}
            with open(phispy_genes, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    if row[9] != '0':
                        genes[row[0]] = row[9]
            
            domains = parse_hmmsearch_output(
                    os.path.join(working_dir, 'out', genome + '.domtblout.txt'),
                    annotator.config['plugins.phispy.evalue_threshold']
                    )
            for gene in sorted(genes.keys()):
                if gene in domains:
                    for hmm in sorted(domains[gene].keys()):
                        outfile.write('\t'.join([
                            gene,
                            genome,
                            'pVOGs',
                            'http://dmk-brain.ecn.uiowa.edu/pVOGs/VOGtables/' +
                            hmm + '.html',
                            'Prophage gene',
                            hmm,
                            hmm + ' virus orthologous group. E-value: ' +
                            domains[gene][hmm]['evalue'] + '. Coordinates: ' +
                            ';'.join(domains[gene][hmm]['coords']) +
                            '. Predicted prophage region ' + genes[gene]
                            ]) + '\n')

    _cleanup(working_dir)
    return output_file

def _export_proteins(genome, output_dir):
    """Creates protein FASTA file"""
    try:
        genome_id = Genome.objects.filter(name=genome).values('id')[0]['id']
    except IndexError:
        print(genome, 'not found')
        raise
    genes = Gene.objects.filter(genome__id = genome_id).select_related('protein')
    out_path = os.path.join(output_dir, str(genome_id) + '.faa')
    with open(out_path, 'w') as outfile:
        for gene in genes:
            if gene.protein is not None:
                outfile.write('>' + gene.locus_tag + '\n')
                outfile.write(gene.protein.sequence + '\n')
    return out_path

def overlaps(start1, end1, coordinates):
    # coordinates is a list of strings "start..end"
    result = False
    overlap = 5
    for coord in coordinates:
        start2, end2 = coord.split('..')
        if int(end1) < int(start2) + overlap:
            continue
        elif int(start1) > int(end2) - overlap:
            continue
        else:
            result = True
            break
    return result
    
def parse_hmmsearch_output(hmmsearchfile, evalue_threshold):
    evalue_threshold = float(evalue_threshold)
    hits = defaultdict(dict)
    with open(hmmsearchfile, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n\r').split()
            evalue = float(row[6])
            if evalue > evalue_threshold:
                continue
            flag = True
            for hmm in hits[row[0]]:
                if overlaps(row[17], row[18], hits[row[0]][hmm]['coords']):
                    flag = False
            if flag:
                if row[3] in hits[row[0]]:
                    hits[row[0]][row[3]]['coords'].append(row[17] + '..' + row[18])
                else:
                    hits[row[0]][row[3]] = {
                           'evalue':row[6],
                           'coords':[row[17] + '..' + row[18]]
                           }
    return hits

def _cleanup(working_dir):
    shutil.rmtree(working_dir)
