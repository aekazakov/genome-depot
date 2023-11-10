import os
import shutil
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError
from django.db import connection
from browser.models import Gene, Genome
"""
    This plugin runs GapMind for a set of genomes.
"""

GAPMIND_REFERENCE = ['aa', 'carbon']

def application(annotator, genomes):
    """
        This function is an entry point of the plugin.
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name 
            as key and GBK path as value
        
    """
    working_dir = os.path.join(annotator.config['cgcms.temp_dir'],
                               'gapmind-plugin-temp'
                               )
    script_path = preprocess(annotator, genomes, working_dir)
    run(script_path)
    
    output_file = postprocess(annotator, genomes, working_dir)
    return output_file

def preprocess(annotator, genomes, working_dir):
    """
        Creates all directories and input files. 
        Input:
            annotator(Annotator): instance of Annotator class
            genomes(dict<str:str>): dictionary with genome name
            as key and GBK path as value
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

    gapmind_script = os.path.join(working_dir, 'run_gapmind.sh')
    
    with open(gapmind_script, 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('source "' + annotator.config['cgcms.conda_path'] + '"\n')
        outfile.write('conda activate ' +
                      annotator.config['plugins.gapmind.conda_env'] +
                      '\n'
                      )
        outfile.write('cd "' + annotator.config['plugins.gapmind.gapmind_dir'] + '"\n')
        
        for genome in sorted(genomes.keys()):
            if os.path.getsize(input_fasta_files[genome]) == 0:
                continue
            os.mkdir(os.path.join(output_dir, genome))
# perl bin/buildorgs.pl -out $TEMPDIR/orgs -orgs 'file:testinput.faa:Test genome'
            outfile.write('\n' + ' '.join(['perl',
                                           'bin/buildorgs.pl',
                                           '-out',
                                           '"' + os.path.join(output_dir, genome, 'orgs') + '"',
                                           '-orgs',
                                           "'file:" + input_fasta_files[genome] +
                                           ":" + genome + "'"
                                           ]) + '\n')
# diamond makedb --in $TEMPDIR/orgs.faa --db $TEMPDIR/orgs.faa
            outfile.write(' '.join(['diamond',
                                    'makedb',
                                    '--in',
                                    '"' + os.path.join(output_dir, genome, 'orgs.faa') + '"',
                                    '--db',
                                    '"' + os.path.join(output_dir, genome, 'orgs.faa') + '"'
                                    ]) + '\n')
                                    
            for collection in GAPMIND_REFERENCE:
# perl bin/gapsearch.pl -orgs $TEMPDIR/orgs -set aa -out $TEMPDIR/aa.hits -nCPU 8
                outfile.write(' '.join(['perl',
                                        'bin/gapsearch.pl',
                                        '-orgs',
                                        '"' + os.path.join(output_dir, genome, 'orgs') + '"',
                                        '-set',
                                        collection, 
                                        '-out',
                                        '"' + os.path.join(output_dir, 
                                                     genome,
                                                     collection + '.hits') + '"',
                                        '-nCPU',
                                        annotator.config['plugins.gapmind.threads']
                                        ]) + '\n')
# perl bin/gaprevsearch.pl -orgs $TEMPDIR/orgs -hits $TEMPDIR/aa.hits
# -curated tmp/path.aa/curated.faa.udb.dmnd -out $TEMPDIR/aa.revhits -nCPU 8
                outfile.write(' '.join(['perl',
                    'bin/gaprevsearch.pl',
                    '-orgs',
                    '"' + os.path.join(output_dir,
                                 genome,
                                 'orgs'
                                 ) + '"',
                    '-hits',
                    '"' + os.path.join(output_dir, 
                                 genome,
                                 collection + '.hits'
                                 ) + '"',
                    '-curated',
                    '"' + os.path.join(annotator.config['plugins.gapmind.gapmind_dir'],
                                 'tmp',
                                 'path.' + collection,
                                 'curated.faa.udb.dmnd') + '"',
                    '-out',
                    '"' + os.path.join(output_dir,
                                 genome,
                                 collection + '.revhits'
                                 ) + '"',
                    '-nCPU',
                    annotator.config['plugins.gapmind.threads']
                    ]) + '\n')
# perl bin/gapsummary.pl -orgs $TEMPDIR/orgs -set aa -hits $TEMPDIR/aa.hits 
# -rev $TEMPDIR/aa.revhits -out $TEMPDIR/aa.sum
                outfile.write(' '.join(['perl',
                                        'bin/gapsummary.pl',
                                        '-orgs',
                                        '"' + os.path.join(output_dir,
                                                     genome,
                                                     'orgs'
                                                     ) + '"',
                                        '-set',
                                        collection,
                                        '-hits',
                                        '"' + os.path.join(output_dir,
                                                     genome,
                                                     collection + '.hits'
                                                     ) + '"',
                                        '-rev',
                                        '"' + os.path.join(output_dir,
                                                             genome,
                                                             collection + '.revhits'
                                                             ) + '"',
                                        '-out',
                                        '"' + os.path.join(output_dir,
                                                             genome,
                                                             collection + '.sum'
                                                             ) + '"'
                                        ]) + '\n')
            
        outfile.write('conda deactivate\n')
    return gapmind_script

def run(script_path):
    """
    Runs GapMind script.
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
        raise CalledProcessError(proc.returncode, proc.args)

def postprocess(annotator, genomes, working_dir):
    """
        Finds GapMind output files and creates file with annotations for upload into DB
    """
    output_file = os.path.join(annotator.config['cgcms.temp_dir'],
                               'gapmind-plugin-output.txt'
                               )
    with open(output_file, 'w') as outfile:
        for collection in GAPMIND_REFERENCE:
            ref_data = defaultdict(dict)
            # Read GapMind reference data
            with open(os.path.join(annotator.config['plugins.gapmind.gapmind_dir'],
                                   collection + '.ref.txt'
                                   ), 'r') as infile:
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    ref_data[row[0]]['pathway'] = row[1]
                    ref_data[row[0]][row[2]] = row[3]
            for genome in genomes:
                cand_outfile = os.path.join(working_dir,
                                            'out',
                                            genome,
                                            collection + '.sum.cand'
                                            )
                if not os.path.exists(cand_outfile):
                    print('File does not exist:', cand_outfile)
                    continue
                rules_outfile = os.path.join(working_dir,
                                             'out',
                                             genome,
                                             collection + '.sum.rules'
                                             )
                if not os.path.exists(rules_outfile):
                    print('File does not exist:', rules_outfile)
                    continue
                rules = defaultdict(dict)
                with open(rules_outfile, 'r') as infile:
                    infile.readline()
                    for line in infile:
                        row = line.rstrip('\n\r').split('\t')
                        pathway = row[3]
                        score = float(row[8])
                        if score < 0:
                            continue
                        n_lo = int(row[7])
                        if n_lo > 0:
                            continue
                        for step in row[9].split():
                            rules[pathway][step] = 1
                with open(cand_outfile, 'r') as infile:
                    infile.readline()
                    for line in infile:
                        row = line.rstrip('\n\r').split('\t')
                        if len(row) < 6:
                            continue
                        pathway = row[3]
                        if pathway not in rules:
                            continue
                        step = row[4]
                        if step not in rules[pathway]:
                            continue
                        score = row[5]
                        if score == '0':
                            continue
                        locus_tag = row[6]
                        outfile.write('\t'.join([locus_tag,
                              genome,
                              'GapMind',
                              'https://papers.genomics.lbl.gov/cgi-bin/gapView.cgi',
                              'GapMind/' + collection,
                              pathway + '/' + step,
                              ref_data[pathway][step] +  ' (' +
                              ref_data[pathway]['pathway'] + ')'
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

def _cleanup(working_dir):
    shutil.rmtree(working_dir)
