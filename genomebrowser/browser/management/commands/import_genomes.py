import os
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer

class Command(BaseCommand):
    help = 'Imports genomes from GBK files'

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            importer = Importer()
            print('Using genomes file', options['i'])
            lines = []
            with open(options['i'], 'r') as infile:
                print('Reading file of genomes', in_file)
                for line in infile:
                    lines.append(line)
            result = importer.import_genomes(lines)
            print(result)
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
