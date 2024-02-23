from django.core.management.base import BaseCommand
from browser.util import update_tags

class Command(BaseCommand):
    help = '''
    Assigns one or more tags to genomes listed
    in a text or tab-separated file.
    '''

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')
        parser.add_argument('-t', help='Tags, comma-separated')

    def handle(self, *args, **options):
        update_tags(options['i'], options['t'])
