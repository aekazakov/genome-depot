from django.core.management.base import BaseCommand
from browser.util import delete_all_data

class Command(BaseCommand):
    help = '''Deletes all data except configuration settings 
    from the GenomeDepot database
    '''
    def handle(self, *args, **options):
        delete_all_data()
