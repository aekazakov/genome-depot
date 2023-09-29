import os
import shutil
import logging
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import Gene, Genome

"""
    This plugin runs HMMSEARCH with Pfam HMM library for a set of genomes.
"""

logger = logging.getLogger("CGCMS")

def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'],
                               'hmmsearch-pfam-temp'
                               )
    script_path = preprocess(annotator, genomes, working_dir)
    run(script_path)
    output_file = postprocess(annotator, genomes, working_dir)
    _cleanup(working_dir)
    return output_file

def preprocess(annotator, genomes, working_dir):
    """
        Creates all directories and input files. 
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
    """
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)
    hmmsearch_outfile = os.path.join(working_dir,
                                     'hmmsearch_pfam.domtblout.txt'
                                     )
    input_file = os.path.join(working_dir, 'hmmsearch-pfam_input.faa')

    annotator.target_genes = defaultdict(list)
    proteins_written = set()
    with open(input_file, 'w') as outfile:
        for genome in genomes:
            try:
                genome_id = Genome.objects.filter(name=genome).values('id')[0]['id']
            except IndexError:
                logger.error('Genome %s not found', genome)
                raise
            for gene in Gene.objects.filter(genome__id = genome_id
                        ).select_related('protein'):
                if gene.protein is not None:
                    annotator.target_genes[gene.protein.protein_hash]\
                    .append((genome, gene.locus_tag))
                    if gene.protein.protein_hash not in proteins_written:
                        outfile.write('>' + gene.protein.protein_hash + '\n')
                        outfile.write(gene.protein.sequence + '\n')
                        proteins_written.add(gene.protein.protein_hash)
    hmmsearch_script = os.path.join(working_dir, 'run_hmmsearch-pfam.sh')
    
    with open(hmmsearch_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('export PATH=' + 
                      '/'.join(annotator.config['cgcms.conda_path'].split('/')[:-3]) +
                      '/envs/' +
                      annotator.config['plugins.hmmsearch_pfam.conda_env'] + 
                      '/bin:$PATH\n'
                      )
        outfile.write('source "' + annotator.config['cgcms.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.hmmsearch_pfam.conda_env'] +
                      '\n'
                      )
        outfile.write(' '.join([annotator.config[\
                                'plugins.hmmsearch_pfam.hmmsearch_command'],
                                '--domtblout', '"' + hmmsearch_outfile + '"',
                                '-o', '/dev/null', '--cut_tc',
                                '--cpu', annotator.config['cgcms.threads'],
                                '--noali', '--notextw',
                                '"' + annotator.config['plugins.hmmsearch_pfam.hmm_lib'] + '"',
                                '"' + input_file + '"'
                                ]) + '\n')
        outfile.write('conda deactivate\n')
    return hmmsearch_script

def run(script_path):
    """
    Runs fama.py for one or several project files.
    """
    cmd = ['/bin/bash', script_path]
    print('Running ' + ' '.join(cmd))
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
        logger.error('hmmsearch finished with error')
        logger.error('Args: %s', ' '.join(cmd))
        raise CalledProcessError(proc.returncode, proc.args)

def postprocess(annotator, genomes, working_dir):
    """
        Finds HMMSEARCH output file and creates file with annotations for upload into DB
    """
    hmmsearch_outfile = os.path.join(working_dir, 'hmmsearch_pfam.domtblout.txt')
    output_file = os.path.join(annotator.config['cgcms.temp_dir'],
                               'hmmsearch-pfam-plugin-output.txt'
                               )
    ref_hmm = {}
    with open(annotator.config['plugins.hmmsearch_pfam.ref_data'], 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n\r').split('\t')
            ref_hmm[row[0]] = {}
            ref_hmm[row[0]]['acc'] = row[1]
            ref_hmm[row[0]]['desc'] = row[-1]

    hits = {}
    with open(hmmsearch_outfile, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n\r').split()
            if (row[0], row[3]) in hits:
                hits[(row[0], row[3])]['coords'].append(row[17] + '..' + row[18])
            else:
                hits[(row[0], row[3])] = {'protein_hash':row[0],
                       'hmm_id':row[3],
                       'hmm_acc':row[4],
                       'evalue':row[6],
                       'coords':[row[17] + '..' + row[18]]
                       }
    with open(output_file, 'w') as outfile:
        for hit in hits.values():
            for item in annotator.target_genes[hit['protein_hash']]:
                outfile.write('\t'.join([item[1], item[0], 'Pfam database',
                              'https://www.ebi.ac.uk/interpro/entry/pfam/' +
                              ref_hmm[hit['hmm_id']]['acc'],
                              'Pfam domain',
                              hit['hmm_id'],
                              hit['hmm_id'] + ' (' + ref_hmm[hit['hmm_id']]['acc'] +
                              '): ' + ref_hmm[hit['hmm_id']]['desc'] +
                              '. E-value: ' +
                              hit['evalue'] + '. Coordinates: ' +
                              ';'.join(hit['coords'])]) + '\n')
    return output_file

def _cleanup(working_dir):
    shutil.rmtree(working_dir)


