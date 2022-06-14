import os
import shutil
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from Bio import GenBank
from django.db import connection

"""
    This plugin runs antiSMASH for a set of genomes.
"""

def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'], 'antismash-plugin-temp')
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
            path of the shell script running abricate in conda environment
    """
    # Create directories
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    output_dir = os.path.join(working_dir, 'out')
    os.mkdir(output_dir)
    # Create shell script
    antismash_script = os.path.join(working_dir, 'run_antismash.sh')
    with open(antismash_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source ' + annotator.config['cgcms.conda_path'] + '\n')
        outfile.write('conda activate ' + annotator.config['plugins.antismash.conda_env'] + '\n')
        for genome in sorted(genomes.keys()):
            outfile.write(' '.join([annotator.config['plugins.antismash.antismash_cmd'],
                                    genomes[genome],
                                    '--output-dir', os.path.join(output_dir, genome),
                                    '--cpus',
                                    annotator.config['plugins.antismash.threads'],
                                    '--genefinding-tool', 'prodigal']) + '\n')
        outfile.write('conda deactivate\n')
    return antismash_script

    
def run(script_path):
    """
    Runs antiSMASH script for GBK files of genomes.
    """
    cmd = ['/bin/bash', script_path]
    print('Running ' + ' '.join(cmd))
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
        Finds antiSMASH output files and creates file with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['cgcms.temp_dir'], 'antismash-plugin-output.txt')
    antismash_ref = {}
    with open(annotator.config['plugins.antismash.antismash_ref'], 'r') as infile:
        for line in infile:
            key,value = line.rstrip('\n\r').split('\t')
            antismash_ref[key] = value

    with open(output_file, 'w') as outfile:
        outfile.write('#Gene\tGenome\tSoure\tURL\tKey\tValue\tNote\n')
        for genome in sorted(genomes.keys()):
            gbk_filename = genomes[genome].split('/')[-1]
            gbk_id = gbk_filename.split('.')[0]
            antismash_dir = os.path.join(working_dir, 'out', genome)
            for dirname in os.listdir(antismash_dir):
                subdir = os.path.join(antismash_dir, dirname)
                if not os.path.isdir(subdir):
                    continue
                for filename in os.listdir(subdir):
                    if 'region'  in filename:
                        outfile.write(get_genes(genome, os.path.join(subdir, filename), antismash_ref))

    _cleanup(working_dir)
    return output_file


def _cleanup(working_dir):
    shutil.rmtree(working_dir)

    
def get_genes(genome_id, gbk_file, antismash_ref):
    result = []
    protocores = {}
    # 1st pass: get protocores
    with open(gbk_file, 'r') as infile:
        parser = GenBank.parse(infile)
        for gbk_record in parser:
            contig_name = gbk_record.accession[0]
            for feature in gbk_record.features:
                if feature.key == 'proto_core':
                    start, end, strand = parse_location(feature.location)
                    for qualifier in feature.qualifiers:
                        if qualifier.key == '/protocluster_number=':
                            protocluster_number = qualifier.value[1:-1] + '_' + contig_name
                        elif qualifier.key == '/product=':
                            product = qualifier.value[1:-1]
                    protocores[protocluster_number] = {}
                    protocores[protocluster_number]['start'] = start
                    protocores[protocluster_number]['end'] = end
                    protocores[protocluster_number]['product'] = product

    # 2nd pass: select genes
    with open(gbk_file, 'r') as infile:
        parser = GenBank.parse(infile)
        for gbk_record in parser:
            contig_name = gbk_record.accession
            for feature in gbk_record.features:
                if feature.key == 'CDS':
                    start, end, strand = parse_location(feature.location)
                    protocluster = None
                    for protocluster_number in protocores:
                        if start >= protocores[protocluster_number]['start'] and end <= protocores[protocluster_number]['end']:
                            protocluster = protocluster_number
                            break
                    if protocluster is None:
                        continue
                    locus_tag = ''
                    source = 'antiSMASH'
                    url = 'https://antismash.secondarymetabolites.org/'
                    key = 'gene cluster'
                    value = protocores[protocluster]['product']
                    note = ''
                    functions = []
                    for qualifier in feature.qualifiers:
                        if qualifier.key == '/locus_tag=':
                            locus_tag = qualifier.value[1:-1]
                        elif qualifier.key == '/gene_functions=':
                            functions.append(parse_function(qualifier.value[1:-1], value))
                        elif qualifier.key == '/gene_kind=':
                            note = qualifier.value[1:-1]
                            if value in antismash_ref:
                                note = antismash_ref[value] + ' ' + note + ' gene.'
                            else:
                                note = note + ' gene from ' + value + ' cluster.'
                    if functions:
                        gene_line = '\t'.join([locus_tag, genome_id, source, url, key, value, note  + ' Functions: ' + '; '.join(functions)]) + '\n'
                        result.append(gene_line)
    return ''.join(result)


def parse_function(function, cluster_function):
    result = ''
    if '(rule-based-clusters)' in function:
        result = function.split('(rule-based-clusters)')[1].strip()
    elif '(smcogs)' in function:
        result = function.split('(smcogs)')[1].strip()
        if '(Score:' in result:
            result = result.split('(Score:')[0].strip()
    elif '(resist)' in function:
        result = function.split('(resist)')[1].strip()
        if '(Score:' in result:
            result = result.split('(Score:')[0].strip()
    elif '(t2pks)' in function:
        result = function.split('(t2pks)')[1].strip()
        if '(Score:' in result:
            result = result.split('(Score:')[0].strip()
    elif 'biosynthetic-additional' in function:
        result.replace('biosynthetic-additional', 'biosynthetic-additional gene')
        if '(Score:' in result:
            result = result.split('(Score:')[0].strip()
    else:
        result = function
    if cluster_function in result:
        result = result.replace(cluster_function, cluster_function + '-specific', 1)
    return result
    

def parse_location(location):
    # Copied from genome browser importer.py
    strand = 1
    if location.startswith('complement('):
        strand = -1
        location = location[11:-1]
    if location.startswith('join('):
        location = location[5:-1]
    location = location.replace('>','')
    location = location.replace('<','')
    location_coords = location.split('..')
    try:
        start = int(location_coords[0])
        end = int(location_coords[-1])
    except ValueError:
        print('Location', location)
        raise
    if start > end:
        # Feature that starts at the end of circular genome and ends at the beginning
        end = start
        start = 1
    return start, end, strand
