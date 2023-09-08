import os
import shutil
import gzip
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import *
"""
    This plugin runs DefenseFinder for a set of genomes.
"""


def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'], 'macsyfinder-plugin-temp')
    script_path = preprocess(annotator, genomes, working_dir)
    run(script_path)
    output_file = postprocess(annotator, genomes, working_dir)
    return(output_file)


def preprocess(annotator, genomes, working_dir):
    """
        Creates all directories and input files. 
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
        Output:
            path of the shell script running DefenseFinder
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

    macsyfinder_script = os.path.join(working_dir, 'run_macsy_finder.sh')
    
    with open(macsyfinder_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source ' + annotator.config['cgcms.conda_path'] + '\n')
        outfile.write('conda activate ' + annotator.config['plugins.macsyfinder.macsyfinder_env'] + '\n')
        outfile.write('cd ' + working_dir + '\n\n')
        
        for genome in sorted(genomes.keys()):
            if os.path.getsize(input_fasta_files[genome]) == 0:
                continue
            genome_dir = os.path.join(output_dir, genome)
            os.mkdir(genome_dir)
            outfile.write(' '.join(['macsyfinder', '--models', annotator.config['plugins.macsyfinder.model'], 'all', '--models-dir', annotator.config['plugins.macsyfinder.models_dir'],
                                    '--out-dir', genome_dir, '--db-type', 'unordered', '--sequence-db', input_fasta_files[genome]]) + '\n')
        outfile.write('\nconda deactivate\n')
    return macsyfinder_script
    
def run(script_path):
    """
    Runs MacsyFinder script for FAA files.
    """
    cmd = ['/bin/bash', script_path]
    print(' '.join(cmd))
    # Close MySQL connection before starting external process because it may run for too long resulting in "MySQL server has gone away" error
    connection.close()
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
        # pylint: disable=no-member
        raise CalledProcessError(proc.returncode, proc.args)

def postprocess(annotator, genomes, working_dir):
    """
        Finds MacsyFinder output files and creates file with annotations for upload into DB
    """
    
    output_file = os.path.join(annotator.config['cgcms.temp_dir'], 'macsyfinder-plugin-output.txt')
    with open(output_file, 'w') as outfile:
        ref_data = defaultdict(dict)
        for genome in genomes:
            macsyfinder_outfile = os.path.join(working_dir, 'out', genome, 'all_systems.tsv')
            if not os.path.exists(macsyfinder_outfile):
                print('File does not exist:', macsyfinder_outfile)
                continue
            with open(macsyfinder_outfile, 'r') as infile:
                infile.readline()
                for line in infile:
                    line = line.rstrip('\n\r')
                    if line == '' or line.startswith('#') or line.startswith('replicon'):
                        continue
                    row = line.split('\t')
                    locus_tag = row[1]
                    df_type = row[4].split('/')[-1]
                    df_subtype = row[2]
                    outfile.write('\t'.join([locus_tag, genome, 'MacsyFinder',
                                  'https://github.com/gem-pasteur/macsyfinder',
                                  'Secretion System',
                                  df_type + '/' + df_subtype,
                                  'Type: ' + df_type + ', subtype: ' + df_subtype]) + '\n')
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

def _cleanup(working_dir):
    #shutil.rmtree(working_dir)
    pass
