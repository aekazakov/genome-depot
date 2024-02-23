from django.core.management.base import BaseCommand
from browser.util import recreate_search_databases

class Command(BaseCommand):
    help = '''
    Deletes and re-cretes nucleotide and protein search
    databases. Use the function if the files are missing or corrupted,
    or if the genome import pipeline crashed before creating the 
    search database files.
    '''
    
    def handle(self, *args, **options):
        recreate_search_databases()
