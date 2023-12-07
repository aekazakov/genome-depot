from django.core.management.base import BaseCommand
from browser.pipeline.annotate import Annotator

class Command(BaseCommand):
    help = """Imports annotations for protein-coding genes from 
    a tab-separated file and writes the annotations to the CGCMS database.
    The genes must be in the database.

    Annotations file must contain the following fields:
    1. Gene locus tag (as in the database).
    2. Genome name (as in the database).
    3. Annotation source.
    4. Annotation URL (external).
    5. Annotation key.
    6. Annotation value.
    7. Annotation note.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to tab-separated annotations file')

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.add_custom_annotations(options['i'])
