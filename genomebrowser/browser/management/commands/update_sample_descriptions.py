import os
from django.core.management.base import BaseCommand
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = '''
    For metaghenomic samples uploaded to Django database,
    this program changes sample descriptions.
    Input file must contain the following fields (tab-separated):
    1. Sample name  (as in the database).
    2. Full name.
    3. Description text
    '''

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to sample descriptions file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            annotator = Annotator()
            lines = []
            print('Reading descriptions from file')
            with open(options['i'], 'r') as in_file:
                for line in in_file:
                    lines.append(line)
            annotator.update_sample_descriptions(lines)
        else:
            raise FileNotFoundError('Input file not found.')
