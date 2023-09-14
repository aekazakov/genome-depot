import os
from django.core.management.base import BaseCommand
from browser.dataimport.importer import Importer

class Command(BaseCommand):
    help = 'Export genomes with gene annotations in genbank format'
    def add_arguments(self, parser):
        parser.add_argument('-g',
                            default='',
                            help='Comma-separated list of genome names'
                            )
        parser.add_argument('-o',
                            default='.',
                            help='Output directory'
                            )
    
    def handle(self, *args, **options):
        if not os.path.exists(options['o']):
            raise ValueError('Directory not exists:', options['o'])
        if options['g'] == '':
            genomes = []
        else:
            genomes = options['g'].split(',')
        importer = Importer()
        importer.export_genomes(options['o'], genomes)

