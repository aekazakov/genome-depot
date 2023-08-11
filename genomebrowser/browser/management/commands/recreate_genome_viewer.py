import os
import gzip
import json
import shutil
from Bio import GenBank
from subprocess import Popen, PIPE, CalledProcessError
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer
from genomebrowser.settings import BASE_URL

class Command(BaseCommand):
    help = 'Deletes and re-creates static files for genome viewer'

    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        genome_id = options['g']
        # Check genome ID
        if genome_id == '':
            raise CommandError('Genome name required')
        print('Looking for genome', genome_id)
        genome_set = Genome.objects.filter(name=genome_id)
        if genome_set.count() == 0:
            print('Genome ' + genome_id + ' not found')
            raise CommandError()
        elif genome_set.count() > 1:
            print('Non-unique genome name: ' + genome_id)
            raise CommandError()
        print('Genome found:', genome_id)
        genome = genome_set[0]
        
        # Configure importer
        importer = Importer()
        importer.inputgenomes[genome_id]['gbk'] = genome.gbk_filepath
        importer.inputgenomes[genome_id]['url'] = genome.external_url
        importer.inputgenomes[genome_id]['external_id'] = genome.external_id
        if genome.strain is None:
            importer.inputgenomes[genome_id]['strain'] = ''
        else:
            importer.inputgenomes[genome_id]['strain'] = genome.strain.strain_id
        if genome.sample is None:
            importer.inputgenomes[genome_id]['sample'] = ''
        else:
            importer.inputgenomes[genome_id]['sample'] = genome.sample.sample_id
        
        script_rows = ['#!/bin/bash', 
                       'source ' + importer.config['cgcms.conda_path'],
                       'conda activate cgcms-jbrowse'
                      ]
        # Create working dir
        temp_dir = os.path.join(importer.config['cgcms.temp_dir'], genome_id + '_jbrowse')
        if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir)
        os.mkdir(temp_dir)
        os.mkdir(os.path.join(temp_dir, 'seq'))
        ref_seqs = []
        
        # Export genome sequence
        genome_fasta = os.path.join(temp_dir, genome_id + '.fna')
        if genome.gbk_filepath.endswith('.gz'):
            gbk_handle = gzip.open(genome.gbk_filepath, 'rt')
        else:
            gbk_handle = open(genome.gbk_filepath, 'r')
        parser = GenBank.parse(gbk_handle)
        with open(genome_fasta, 'w') as outfile:
            for gbk_record in parser:
                contig_sequence = gbk_record.sequence
                contig_id = gbk_record.locus
                contig_name = gbk_record.accession
                nucl_seq_uid = '>' + contig_id + '|' + genome_id
                outfile.write('>' + contig_id + '\n' + ''.join(contig_sequence) + '\n')
                ref_seqs.append({'name': contig_id, 'start': 0, 'end': gbk_record.size, 'length': gbk_record.size})
        gbk_handle.close()
        
        with open(os.path.join(temp_dir, 'seq', 'refSeqs.json'), 'w') as outfile:
            json.dump(ref_seqs, outfile)
            
        jbrowse_config = {'formatVersion' : 1,
                  'tracks': [],
                  'refseqs': genome_fasta + '.fai'
                 }
        jbrowse_track = {
                         'category' : 'Reference sequence',
                         'key' : 'Reference sequence',
                         'label' : 'DNA',
                         'seqType' : 'dna',
                         'storeClass' : 'JBrowse/Store/SeqFeature/IndexedFasta',
                         'type' : 'SequenceTrack',
                         'urlTemplate' : os.path.basename(genome_fasta)
                        }
        jbrowse_config['tracks'].append(jbrowse_track)
        # Index FASTA file
        script_rows.append('samtools faidx ' + genome_fasta)
            
        # Create GFF3 files
        feature_types = {'CDS': ('blue', 'CDSs', 'normal'),
                         'RNA': ('coral', 'RNAs', 'compact'),
                         'pseudogene': ('grey', 'Pseudogenes', 'compact'),
                         'operon': ('goldenrod', 'Operons', 'compact')
                         }
        for feature_type in feature_types:
            gff_file = importer.export_gff(genome_id, temp_dir, feature_type)
            script_rows.append('(grep ^"#" ' + gff_file + '; grep -v ^"#" ' + gff_file + ' | grep -v "^$" | grep "\t" | sort -k1,1 -k4,4n)|bgzip > ' + gff_file + '.gz')
            script_rows.append('tabix -p gff ' + gff_file + '.gz')
            script_rows.append('rm ' + gff_file)
            jbrowse_track = {'compress' : 0,
                             'key' : feature_types[feature_type][1],
                             'label' : feature_types[feature_type][1],
                             'storeClass' : 'JBrowse/Store/SeqFeature/GFF3Tabix',
                             'style' : {
                                'color': feature_types[feature_type][0],
                                },
                             'trackType' : 'CanvasFeatures',
                             'type' : 'JBrowse/View/Track/CanvasFeatures',
                             'displayMode': feature_types[feature_type][2],
                             'urlTemplate' : os.path.basename(gff_file) + '.gz',
                             
                             'onClick': {
                                 'label': '{Id}: {start}..{end}\n{product}',
                                 'title': '{name} {type}',
                                 'action': 'defaultDialog'
                                 }
                             }
            if feature_type == 'operon':
                jbrowse_track['fmtDetailValue_Name'] = "function(name, feature) { return '<a href=\"" + BASE_URL + "/operon/'+feature.get('internal_id')+'\" target=\"_blank\">' + name + ' ( '+ feature.get('internal_id') + ')</a>'; }"
            else:
                jbrowse_track['fmtDetailValue_Name'] = "function(name, feature) { return '<a href=\"" + BASE_URL + "/gene/'+feature.get('internal_id')+'\" target=\"_blank\">' + name + ' ( '+ feature.get('internal_id') + ')</a>'; }",
            jbrowse_config['tracks'].append(jbrowse_track)
            
        # Write config file
        with open(os.path.join(temp_dir, 'trackList.json'), 'w') as outfile:
            json.dump(jbrowse_config, outfile, indent = 2)
        
        script_path = os.path.join(importer.config['cgcms.temp_dir'], 'make_jbrowse_files.sh')
        with open(script_path, 'w') as outfile:
            outfile.write('\n'.join(script_rows))
            
        cmd = ['bash', script_path]
        with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
            for line in proc.stdout:
                print(line, end='')
        if proc.returncode != 0:
            # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
            # pylint: disable=no-member
            raise CalledProcessError(proc.returncode, proc.args)
        
        dest_dir = os.path.join(importer.config['cgcms.json_dir'], genome_id)
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)
        shutil.copytree(temp_dir, dest_dir)
        shutil.rmtree(temp_dir)
        
        # OLD: Generate static files
        # importer.export_jbrowse_genome_data(genome_id)
