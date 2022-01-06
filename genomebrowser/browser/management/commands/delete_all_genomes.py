from django.core.management.base import BaseCommand
from browser.models import *

class Command(BaseCommand):
    help = 'Deletes all gene, protein, genome, strain, sample entries from the database'
    
    def handle(self, *args, **options):
        # Delete all mappings before deleting genes
        Annotation.objects.all().delete()
        # Delete genes before deleting proteins
        Gene.objects.all().delete()
        # Delete proteins
        Protein.objects.all().delete()
        # Delete genomes
        Genome.objects.all().delete()
        # Delete strains
        Strain_metadata.objects.all().delete()
        Strain.objects.all().delete()
        # Delete samples
        Sample_metadata.objects.all().delete()
        Sample.objects.all().delete()
