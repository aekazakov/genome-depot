import sys
from django.core.management.base import BaseCommand
from browser.util import delete_all_genomes

class Command(BaseCommand):
    help = '''
    This command deletes all Annotation, Gene, Protein, Genome,
    Strain, Sample entries from the database
    '''
    
    def handle(self, *args, **options):
        while True:
            answer = input('All genomes, strains and samples in the CGCMS database '
                'will be deleted. This action cannot be reversed. Continue? (y/n)'
            )
            if answer.lower() in ["y","yes"]:
                break
            elif answer.lower() in ["n","no"]:
                print('Exiting.')
                sys.exit()
            else:
                print('Please answer yes/y or no/n.')
