import os
from django.core.management.base import BaseCommand
from browser.models import *
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = """For strain uploaded to Django database, this program adds metadata from tab-separated file.
    Metadata file must contain the following fields:
    1. Strain name  (as in the database).
    2. Metadata source.
    3. Metadata URL (external).
    4. Metadata key.
    5. Metadata value.
    """

    def add_arguments(self, parser):
        parser.add_argument('-c', default='config.ini', help='Path to config file')
        parser.add_argument('-i', help='Path to strain metadata file')

    def handle(self, *args, **options):
        if os.path.exists(options['c']):
            print('Config file sent', options['c'])
            annotator = Annotator(options['c'])
            annotator.add_strain_metadata(options['i'])
        else:
            raise FileNotFoundError('Config file not found.')
