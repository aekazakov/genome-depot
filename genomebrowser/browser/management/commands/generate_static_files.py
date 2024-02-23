from django.core.management.base import BaseCommand
from browser.util import generate_static_files

class Command(BaseCommand):
    help = '''Generates genome viewer static files for genomes
        from the input file and re-creates search databases'''

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        generate_static_files(options['i'])
