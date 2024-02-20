import os
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.util import delete_genome

class Command(BaseCommand):
    help = 'Deletes one genome with all genes and annotations from the database'
    
    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        genome_name = options['g']
        delete_genome(genome_name)
