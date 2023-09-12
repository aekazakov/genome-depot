import os
from django.core.management.base import BaseCommand
from browser.dataimport.annotator import Annotator

class Command(BaseCommand):
    help = 'For protein sequences uploaded to Django database, this program finds FitnessBrowser genes from the same genome, and writes annotations to the database'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.add_fitbrowser_links()
