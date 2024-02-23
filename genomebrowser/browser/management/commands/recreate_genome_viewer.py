from django.core.management.base import BaseCommand
from browser.util import regenerate_jbrowse_files

class Command(BaseCommand):
    help = '''Deletes and re-creates static files for one genome'''

    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        regenerate_jbrowse_files(options['g'])
