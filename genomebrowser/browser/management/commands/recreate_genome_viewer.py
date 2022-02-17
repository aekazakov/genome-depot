import os
import gzip
from Bio import GenBank
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer

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
        
        # Export genome sequences
        genome_fasta = os.path.join(importer.config['cgcms.temp_dir'], genome_id + '.contigs.fasta')
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
        gbk_handle.close()
        
        # Generat static files
        importer.export_jbrowse_genome_data(genome_id)
