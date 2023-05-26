import os
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.annotator import Annotator


class Command(BaseCommand):
    help = 'For protein sequences from selected genomes uploaded to Django database, this program runs hmmsearch with PFAM and TIGRFAM HMMs, and writes domain annotations to database. Existing domain annotations are deleted.'

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            annotator = Annotator()
            genomes = {}
            with open(options['i'], 'r') as infile:
                for line in infile:
                    if line.startswith('#'):
                        continue
                    filepath, genome_name, _, _, _, _ = line.rstrip('\n\r').split('\t')
                    genomes[genome_name] = Genome.objects.get(name=genome_name).gbk_filepath
                    #genomes[genome_name] = filepath
            annotator.run_external_tools(genomes)
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
