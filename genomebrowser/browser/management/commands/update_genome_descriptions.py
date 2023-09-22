import os
from django.core.management.base import BaseCommand
from browser.pipeline.annotate import Annotator

class Command(BaseCommand):
    help = """For genomes uploaded to Django database, 
    this program changes genome descriptions.
    Input file must contain the following fields:
    1. Genome name  (as in the database).
    2. New genome name (same as [1] if no change needed).
    3. Description text.
    4. External URL.
    5. External ID.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to genome descriptions file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            annotator = Annotator()
            annotator.update_genome_descriptions(options['i'])
        else:
            raise FileNotFoundError('Genomes file not found.')
