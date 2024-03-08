import os
import shutil
from subprocess import Popen, PIPE, CalledProcessError
"""
    This plugin runs abricate for a set of genomes.
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
                               'abricate-plugin-temp'
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
            path of the shell script running abricate in conda environment
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    abricate_script = os.path.join(working_dir, 'run_abricate.sh')
    
    with open(abricate_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['core.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.abricate.conda_env'] + '\n'
                      )
        for genome in sorted(genomes.keys()):
            outfile.write(' '.join(['abricate', '--threads',
                                    annotator.config['plugins.abricate.threads'],
                                    '"' + genomes[genome] + '"',
                                    '>',
                                    '"' + os.path.join(output_dir, genome + '.abricate.tsv') + '"']
                                    )
                          + '\n')
        outfile.write('conda deactivate\n')
    return abricate_script

    
def run(script_path):
    """
    Runs abricate script for GBK files of genomes.
    """
    cmd = ['/bin/bash', script_path]
    print(' '.join(cmd))
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
        # pylint: disable=no-member
        raise CalledProcessError(proc.returncode, proc.args)


def postprocess(annotator, genomes, working_dir):
    """
        Finds abricate output files and creates file with annotations for upload into DB
    """
    
    output_file = os.path.join(annotator.config['core.temp_dir'],
                               'abricate-plugin-output.txt'
                               )

    with open(output_file, 'w') as outfile:
        for genome in sorted(genomes.keys()):
            abricate_output = os.path.join(working_dir,
                                           'out',
                                           genome + '.abricate.tsv'
                                           )
            with open(abricate_output, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    if len(row) < 2:
                        continue
                    gbk_path = row[0]
                    if gbk_path != genomes[genome]:
                        raise ValueError('GBK path for genome ' + genome +
                                         'in abricate output differs from expected ' +
                                         genomes[genome]
                                         )
                    locus_tag = row[12] # it may report protein ID instead of locus tag
                    function = row[14]
                    description = row[5] + '(' + row[10] + '% identity, ' + \
                                  row[9] + '% coverage): ' + row[13]
                    outfile.write('\t'.join([locus_tag, genome, 'abricate',
                                  'https://github.com/tseemann/abricate',
                                  'Type',
                                  function,
                                  description]) + '\n')
    #TODO: uncomment
    #_cleanup(working_dir)
    return output_file


def _cleanup(working_dir):
    shutil.rmtree(working_dir)
