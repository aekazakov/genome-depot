import os
from django.core.management.base import BaseCommand, CommandError
from browser.pipeline.genome_import import Importer

class Command(BaseCommand):
    help = '''Imports genomes from GBK files.

    Input is a tab-separated file with six columns:
    - path to Genbank file
    - genome ID (no spaces)
    - strain name (no spaces)
    - sample ID (no spaces)
    - URL (link to NCBI genome assembly etc.)
    - external ID (for example, "NCBI:GCF_000006945.2")
    '''

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            importer = Importer()
            print('Using genomes file', options['i'])
            lines = []
            with open(options['i'], 'r') as infile:
                print('Reading file of genomes', options['i'])
                for line in infile:
                    lines.append(line)
            result = importer.import_genomes(lines)
            print(result)
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
