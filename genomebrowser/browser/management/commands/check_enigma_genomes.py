import os
from django.core.management.base import BaseCommand
from django.core.mail import mail_admins
from genomebrowser.settings import BASE_DIR

class Command(BaseCommand):
    help = '''
    Checks ENIGMA data repository for new genomes.
    If you have no access to the repository, you don't have to run this command.
    '''

    def add_arguments(self, parser):
        parser.add_argument('username')
        parser.add_argument('password')

    def handle(self, *args, **options):
        if os.path.exists(os.path.join(BASE_DIR, 'enigma', 'enigma.py')):
            from enigma.enigma import check_enigma_repository
            message = check_enigma_repository(options['username'], options['password'])
            mail_admins('Check ENIGMA genomes repository', message)
        else:
            print('ENIGMA module not installed')

