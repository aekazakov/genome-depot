import os
from django.core.management.base import BaseCommand
from browser.models import *
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = """For strain uploaded to Django database, this program updates metadata entries.
    The metadata is taken from two sources. The first is a spreadsheet provided with -i parameter. The second is isolates.genomics.lbl.gov API.
    The spreadsheet must contain strain identifier in the first column and names of metadata categories in the first row. All another non-empty cells will be considered values.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to strain metadata file')

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.update_strain_metadata(xlsx_path=options['i'], xlsx_file=None)
