import os
import shutil
import gzip
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import *
"""
    This plugin runs eCIS-screen for a set of genomes.
"""


def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'], 'ecis-screen-plugin-temp')
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
            path of the shell script running eCIS-screen
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    input_dir = os.path.join(working_dir, 'genomes')
    os.mkdir(input_dir)
    # Create shell script
    ecis_screen_script = os.path.join(working_dir, 'run_ecis-screen.sh')
    for genome in genomes:
        genome_dir = os.path.join(input_dir, genome)
        os.mkdir(genome_dir)
        if genomes[genome].endswith('.gz'):
            file_handle = gzip.open(genomes[genome], 'rt')
        else:
            file_handle = open(genomes[genome], 'rt')
        with gzip.open(os.path.join(genome_dir, genome + '_genomic.gbff.gz'), 'wt') as outfile:
            for line in file_handle:
                outfile.write(line)
    with open(ecis_screen_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('cd ' + working_dir + '\n')
        outfile.write(' '.join([annotator.config['plugins.ecis_screen.ecis-screen_cmd'],
                                annotator.config['plugins.ecis_screen.ecis_hmm'],
                                input_dir,
                                os.path.join(working_dir, 'ecis-screen_summary.txt')]) + '\n')
    return ecis_screen_script

    
def run(script_path):
    """
    Runs eCIS-screen script for GBK files of genomes.
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
        Finds eCIS-screen output files and creates file with annotations for upload into DB
    """
    
    output_file = os.path.join(annotator.config['cgcms.temp_dir'], 'ecis-screen-plugin-output.txt')
    accessions = {}
    for genome in genomes:
        gbk_file = os.path.join(working_dir, 'genomes', genome, genome + '_genomic.gbff.gz')
        with gzip.open(gbk_file, 'rt') as infile:
            for line in infile:
                if line.startswith('ACCESSION   '):
                    accession = line[12:].rstrip('\n\r').split(' ')[0]
                    accessions[genome + '_' + accession] = genome

    with open(output_file, 'w') as outfile:
        with open(os.path.join(working_dir, 'ecis-screen_summary.txt'), 'r') as infile:
            header = infile.readline().rstrip('\n\r').split('\t')
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                try:
                    genome = accessions[row[0]]
                except KeyError:
                    print ('Unable to find genome name for sequence ID ' + row[0])
                    raise
                for ind, genes in enumerate(row):
                    if ind == 0:
                        continue
                    if genes == '-':
                        continue
                    family = header[ind]
                    description = family + ' subunit of extracellular contractile injection system (eCIS).'
                    for locus_tag in genes.split(','):
                        outfile.write('\t'.join([locus_tag, genome, 'eCIS-screen',
                                      'https://github.com/ipb-jianyang/eCIS-screen',
                                      'eCIS gene',
                                      family,
                                      description]) + '\n')
    _cleanup(working_dir)
    return output_file



def _cleanup(working_dir):
    shutil.rmtree(working_dir)