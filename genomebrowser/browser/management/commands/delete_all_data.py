import sys
from django.core.management.base import BaseCommand
from browser.util import delete_all_data

class Command(BaseCommand):
    help = '''Deletes all data except configuration settings 
    from the CGCMS database
    '''
    def handle(self, *args, **options):
        print('')
        while True:
            answer = input('All data in the CGCMS database will be deleted.'
                'This action cannot be reversed. Continue? (y/n)'
            )
            if answer.lower() in ["y","yes"]:
                break
            elif answer.lower() in ["n","no"]:
                print('Exiting.')
                sys.exit()
            else:
                print('Please answer yes/y or no/n.')
        delete_all_data()