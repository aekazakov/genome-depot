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
            importer.import_genomes(options['i'])
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
