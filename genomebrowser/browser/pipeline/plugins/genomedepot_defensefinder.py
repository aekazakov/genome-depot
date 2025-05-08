import os
import shutil
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.pipeline.util import export_proteins_bygenome
from browser.models import Genome
"""
    This plugin runs DefenseFinder for a set of genomes.
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
                               'defensefinder-plugin-temp'
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
            path of the shell script running DefenseFinder
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    input_fasta_files =  export_proteins_bygenome(genomes, working_dir)

    defensefinder_script = os.path.join(working_dir, 'run_defense_finder.sh')
    
    with open(defensefinder_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['core.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.defensefinder.conda_env'] +
                      '\n'
                      )
        outfile.write('cd "' + working_dir + '"\n\n')
        
        for genome in sorted(genomes.keys()):
            if os.path.getsize(input_fasta_files[genome]) == 0:
                continue
            genome_dir = os.path.join(output_dir, genome)
            os.mkdir(genome_dir)
            outfile.write(' '.join(['defense-finder',
            'run',
            '--models-dir',
            '"' + annotator.config['plugins.defensefinder.defensefinder_models_dir'] + '"',
            '-o',
            '"' + genome_dir + '"',
            '"' + input_fasta_files[genome] + '"'])
            + '\n')
        outfile.write('\nconda deactivate\n')
    return defensefinder_script
    
def run(script_path):
    """
    Runs eCIS-screen script for GBK files of genomes.
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
        Finds DefenseFinder output files and creates file
        with annotations for upload into DB
    """
    
    output_file = os.path.join(annotator.config['core.temp_dir'],
                               'defensefinder-plugin-output.txt'
                               )
    with open(output_file, 'w') as outfile:
        for genome in genomes:
            genome_id = Genome.objects.get(name=genome).id
            defensefinder_outfile = os.path.join(working_dir,
                                                 'out',
                                                 genome,
                                                 str(genome_id) + '_defense_finder_genes.tsv'
                                                 )
            if not os.path.exists(defensefinder_outfile):
                print('File does not exist:', defensefinder_outfile)
                continue
            with open(defensefinder_outfile, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    locus_tag = row[1]
                    outfile.write('\t'.join([locus_tag, genome, 'DefenseFinder',
                                  'https://github.com/mdmparis/defense-finder',
                                  row[-1] + ' system',
                                  row[2],
                                  'Type: ' + row[-3] + ', subtype: ' + row[-2] + ', gene name: ' + row[2]
                                  ]) + '\n')
                        
    _cleanup(working_dir)
    return output_file

def _cleanup(working_dir):
    shutil.rmtree(working_dir)
