import os
from django.core.management.base import BaseCommand
from browser.models import *
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = """For samples uploaded into Django database, this program adds metadata from tab-separated file.
    Metadata file must contain the following fields:
    1. Sample name  (as in the database).
    2. Metadata source.
    3. Metadata URL (external).
    4. Metadata key.
    5. Metadata value.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to sample metadata file')

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.add_sample_metadata(options['i'])
