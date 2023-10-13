import os
from django.core.management.base import BaseCommand
from browser.pipeline.genome_import import Importer

class Command(BaseCommand):
    help = '''Generates genome viewer static files for genomes
        from the input file and re-creates search databases'''

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            importer = Importer()
            print('Genomes file', options['i'])
            importer.generate_static_files(options['i'])
        else:
            raise ValueError('Genomes file not found.')
