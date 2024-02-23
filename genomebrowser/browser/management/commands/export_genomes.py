from django.core.management.base import BaseCommand
from browser.util import export_genomes

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
        if options['g'] == '':
            genomes = []
        else:
            genomes = options['g'].split(',')
        export_genomes(options['o'], genomes)

