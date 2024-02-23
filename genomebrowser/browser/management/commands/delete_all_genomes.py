from django.core.management.base import BaseCommand
from browser.util import delete_all_genomes

class Command(BaseCommand):
    help = '''
    This command deletes all Annotation, Gene, Protein, Genome,
    Strain, Sample entries from the database
    '''
    
    def handle(self, *args, **options):
        delete_all_genomes()
