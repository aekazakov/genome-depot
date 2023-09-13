from django.core.management.base import BaseCommand
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = """For protein sequences uploaded to Django database, this program
    adds annotations from tab-separated file and writes the annotations to the database.

    Annotations file must contain the following fields:
    1. Gene locus tag (as in the database).
    2. Genome name (as in the database).
    3. Annotation source.
    4. Annotation URL (external).
    5. Annotataion key.
    6. Annotation value.
    7. Annotation note.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to tab-separated annotations file')

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.add_custom_annotations(options['i'])
