import os
import shutil
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.pipeline.util import export_proteins_bygenome
"""
    This plugin runs AMRFinder for a set of genomes.
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
                               'amrfinder-plugin-temp'
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
            path of the shell script running amrfinder in conda environment
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    input_fasta_files =  export_proteins_bygenome(genomes, working_dir)

    amrfinder_script = os.path.join(working_dir, 'run_amrfinder.sh')
    
    with open(amrfinder_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['core.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.amrfinder.conda_env'] +
                      '\n'
                      )
        for genome in sorted(genomes.keys()):
            outfile.write(' '.join(['amrfinder',
                                    '--threads',
                                    annotator.config['plugins.amrfinder.threads'],
                                    '-p',
                                    '"' + input_fasta_files[genome] + '"',
                                    '-o',
                                    '"' + os.path.join(output_dir,
                                                 genome + '.amrfinder.tsv"'
                                                 ),
                                    '--log',
                                    '"' + os.path.join(working_dir,
                                                 genome + '_log.txt"'
                                                 )
                                    ]) + '\n')
        outfile.write('conda deactivate\n')
    return amrfinder_script

    
def run(script_path):
    """
    Runs amrfinder script for GBK files of genomes.
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
        Finds amrfinder output files and creates file
        with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['core.temp_dir'],
                               'amrfinder-plugin-output.txt'
                               )
    with open(output_file, 'w') as outfile:
        outfile.write('#Gene\tGenome\tSoure\tURL\tKey\tValue\tNote\n')
        for genome in sorted(genomes.keys()):
            amrfinder_output = os.path.join(working_dir, 'out',
                                            genome + '.amrfinder.tsv'
                                            )
            with open(amrfinder_output, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    locus_tag = row[0]
                    function = row[4] + ': ' + row[7]
                    description = row[2] + ' [' + row[4] + '/' + row[5] + \
                                  '/' + row[6] + '/' + row[7] + ']'
                    outfile.write('\t'.join([locus_tag, genome, 'AMRFinderPlus',
                                  'https://www.ncbi.nlm.nih.gov/pathogens/' + 
                                  'antimicrobial-resistance/AMRFinder/',
                                  'Type/subclass',
                                  function,
                                  description]) + '\n')
    _cleanup(working_dir)
    return output_file


def _cleanup(working_dir):
    shutil.rmtree(working_dir)
