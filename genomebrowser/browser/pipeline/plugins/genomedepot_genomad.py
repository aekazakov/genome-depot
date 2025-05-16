import os
import shutil
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.pipeline.util import export_nucl_bygenome
from browser.models import Gene, Genome
"""
    This plugin runs geNomad for a set of genomes.
"""


def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['core.temp_dir'],
                               'genomad-plugin-temp'
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
            path of the shell script running geNomad in conda environment
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    input_fasta_files =  export_nucl_bygenome(genomes, working_dir)

    tool_script = os.path.join(working_dir, 'run_genomad.sh')
    
    with open(tool_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['core.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.genomad.conda_env'] +
                      '\n'
                      )
        for genome in sorted(genomes.keys()):
            outfile.write(' '.join(['genomad end-to-end --cleanup',
                                    '"' + input_fasta_files[genome] + '"',
                                    '"' + os.path.join(output_dir, genome)  + '"',
                                    annotator.config['plugins.genomad.ref_db']
                                    ]) + '\n')
        outfile.write('conda deactivate\n')
    return tool_script

    
def run(script_path):
    """
    Runs geNomad script for FASTA files of genomes.
    """
    cmd = ['/bin/bash', script_path]
    print('Running ' + ' '.join(cmd))
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
        Finds geNomad output files and creates file
        with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['core.temp_dir'],
                               'genomad-plugin-output.txt'
                               )
    with open(output_file, 'w') as outfile:
        outfile.write('#Gene\tGenome\tSource\tURL\tKey\tValue\tNote\n')
        for genome in sorted(genomes.keys()):
            genome_id = str(
                Genome.objects.filter(name=genome).values_list('id', flat=True)[0]
            )
            proviruses = {}
            proviruses_output = os.path.join(working_dir, 'out',
                genome,
                genome_id + '_find_proviruses',
                genome_id + '_provirus.tsv',
            )
            with open(proviruses_output, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    proviruses[row[0]] = row[1]
            if not proviruses:
                continue
            genes_output = os.path.join(working_dir, 'out',
                genome,
                genome_id + '_summary',
                genome_id + '_virus_genes.tsv',
            )
            with open(genes_output, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    if row[8].endswith('VV') or row[8].endswith('VP') or \
                    row[8].endswith('VC'):
                        function = row[-1]
                        if function == 'NA':
                            function = 'Unknown function'
                        provirus_id = '_'.join(row[0].split('_')[:-1])
                        if provirus_id not in proviruses:
                            continue
                        contig_id = proviruses[provirus_id]
                        if row[4] == '-1':
                            end = int(row[1])
                            if end < 4:
                                end = 1
                            genes = Gene.objects.filter(
                                genome__name=genome,
                                contig__contig_id=contig_id,
                                start=end
                            )
                        else:
                            end = int(row[2])
                            genes = Gene.objects.filter(
                                genome__name=genome,
                                contig__contig_id=contig_id,
                                end=end
                            )
                        if len(genes) == 1:
                            outfile.write('\t'.join([genes[0].locus_tag, genome,
                                annotator.config['plugins.genomad.display_name'],
                                'https://portal.nersc.gov/genomad/',
                                'geNomad virus-specific marker',
                                row[8],
                                function]) + '\n'
                            )
    _cleanup(working_dir)
    return output_file


def _cleanup(working_dir):
    shutil.rmtree(working_dir)
