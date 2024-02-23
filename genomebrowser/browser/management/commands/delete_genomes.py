from django.core.management.base import BaseCommand
from browser.util import delete_genomes

class Command(BaseCommand):
    help = '''
    This command deletes all Annotation, Gene, Protein, Genome, Strain, Sample
    entries from the database for the genomes listed in the input file
    '''
    
    def add_arguments(self, parser):
        parser.add_argument('-i', default='', help='File with list of genomes')

    def handle(self, *args, **options):
        delete_genomes(options['i'])
