from django.core.management.base import BaseCommand
from browser.pipeline.annotate import Annotator

class Command(BaseCommand):
    help = """Imports sample metadata from tab-separated file
    into the Django database.
    Sample records must be created before running this program.
    
    The input file must contain the following fields:
    1. Sample name (as in the database).
    2. Metadata source.
    3. Metadata URL (external).
    4. Metadata key.
    5. Metadata value.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to sample metadata file')

    def handle(self, *args, **options):
        lines = []
        annotator = Annotator()
        print('Reading metadata from file')
        with open(options['i'], 'r') as in_file:
            for line in in_file:
#               line = line.decode()
                lines.append(line)
        annotator.add_sample_metadata(lines)
