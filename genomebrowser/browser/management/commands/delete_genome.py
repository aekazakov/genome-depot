from django.core.management.base import BaseCommand
from browser.util import delete_genome

class Command(BaseCommand):
    help = 'Deletes one genome with all genes and annotations from the database'
    
    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        delete_genome(options['g'])
