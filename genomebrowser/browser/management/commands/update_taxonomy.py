import os
from django.core.management.base import BaseCommand
from browser.dataimport.taxonomy import update_taxonomy

class Command(BaseCommand):
    help = """Update taxonomy with new data from NCBI.
    """


    def handle(self, *args, **options):
        update_taxonomy()
