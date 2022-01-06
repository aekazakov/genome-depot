import os
import shutil
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import *
"""
    This plugin runs Fama for a set of genomes.
"""

collections = ['nitrogen_v11', 'cazy_v1', 'universal_v1.4']


def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'], 'fama-plugin-temp')
    project_files = preprocess(annotator, genomes, working_dir)
    
    for collection in project_files:
        run(annotator, project_files[collection])
    
    output_file = postprocess(annotator, genomes, project_files, working_dir)
    return output_file


def preprocess(annotator, genomes, working_dir):
    """
        Creates all directories and input files. 
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name as key and GBK path as value
    """
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    #input_dir = os.path.join(working_dir, 'input_files')
    #os.mkdir(input_dir)
    #fama_input_files = {}
    fama_input_file = os.path.join(working_dir, 'fama_input.faa')
    annotator.target_genes = defaultdict(list)
    proteins_written = set()
    with open(fama_input_file, 'w') as outfile:
        for genome in genomes:
            try:
                genome_id = Genome.objects.filter(name=genome).values('id')[0]['id']
            except IndexError:
                print(genome, 'not found')
                raise
            genes = Gene.objects.filter(genome__id = genome_id).select_related('protein')
            for gene in genes:
                if gene.protein is not None:
                    annotator.target_genes[gene.protein.protein_hash].append((genome, gene.locus_tag))
                    if gene.protein.protein_hash not in proteins_written:
                        outfile.write('>' + gene.protein.protein_hash + '\n')
                        outfile.write(gene.protein.sequence + '\n')
                        proteins_written.add(gene.protein.protein_hash)


#        fama_input_files[genome] = _export_proteins(genome, input_dir)
    project_files = {}
    seqlist_file = os.path.join(working_dir, 'sequences_list.txt')
    with open(seqlist_file, 'w') as outfile:
        outfile.write('\t'.join(['fama_input', fama_input_file]) + '\n')
#        for genome in sorted(fama_input_files.keys()):
#            outfile.write('\t'.join([genome, fama_input_files[genome]]) + '\n')
    for collection in collections:
        fama_project_file = os.path.join(working_dir, 'project_' + collection + '.ini')
        _run_prepare(annotator, seqlist_file, collection)
        project_files[collection]  = fama_project_file
    return project_files

    
def run(annotator, fama_project_file):
    """
    Runs fama.py for FASTA files of proteins.
    """
    cmd = ['python', os.path.join(annotator.config['plugins.fama.fama_dir'], 'fama.py'),
          '-c', annotator.config['plugins.fama.fama_config'],
          '-p', fama_project_file,
          '--prot'
          ]
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


def postprocess(annotator, genomes, project_files, working_dir):
    """
        Finds Fama output file and creates file with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['cgcms.temp_dir'], 'fama-plugin-output.txt')
    ref_data = defaultdict(dict)
    # Read Fama reference data
    with open(annotator.config['plugins.fama.fama_nitrate_lib'], 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            ref_data[row[0]]['description'] = row[1]
            ref_data[row[0]]['library'] = 'Nitrogen cycle genes'
    with open(annotator.config['plugins.fama.fama_universal_lib'], 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            ref_data[row[0]]['description'] = row[1]
            ref_data[row[0]]['library'] = 'Universal single-copy genes'
    with open(annotator.config['plugins.fama.fama_cazy_lib'], 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            ref_data[row[0]]['description'] = row[1]
            ref_data[row[0]]['library'] = 'Carbohydrate-active enzymes'

    with open(output_file, 'w') as outfile:
        for collection in project_files:
            fama_output = os.path.join(working_dir, collection, 'all_proteins.list.txt')
            with open(fama_output, 'r') as infile:
                infile.readline()
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    if len(row) < 2:
                        continue
                    #genome_name = row[0]
                    #locus_tag = row[1]
                    protein_hash = row[1]
                    functions = row[2].split(',')
                    for function in functions:
                        for item in annotator.target_genes[protein_hash]:
                            outfile.write('\t'.join([item[1], item[0], 'Fama',
                                          'https://iseq.lbl.gov/fama/reference',
                                          'Fama annotation',
                                          function,
                                          ref_data[function]['description'] +  ' (' + ref_data[function]['library'] + ')']) + '\n')
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


def _run_prepare(annotator, seqlist_file, collection):
    """
    Runs fama_prepare.py for FASTA files of proteins.
    """
    cmd = ['python', os.path.join(annotator.config['plugins.fama.fama_dir'], 'fama_prepare.py'),
          '-c', annotator.config['plugins.fama.fama_config'],
          '-p', collection,
          '-i', seqlist_file,
          '-r', collection,
          '--prot'
          ]
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
        # pylint: disable=no-member
        raise CalledProcessError(proc.returncode, proc.args)


def _cleanup(working_dir):
    shutil.rmtree(working_dir)