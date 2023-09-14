from django.core.management.base import BaseCommand
from browser.models import Annotation
from browser.models import Gene
from browser.models import Protein
from browser.models import Genome
from browser.models import Strain_metadata
from browser.models import Strain
from browser.models import Sample_metadata
from browser.models import Sample

class Command(BaseCommand):
    help = '''
    Deletes all Annotation, Gene, Protein, Genome,
    Strain, Sample entries from the database
    '''
    
    def handle(self, *args, **options):
        # Delete all mappings before deleting genes
        print("Deleting annotations...")
        Annotation.objects.all().delete()
        # Delete genes before deleting proteins
        print("Deleting genes...")
        Gene.objects.all().delete()
        # Delete proteins
        print("Deleting proteins...")
        Protein.objects.all().delete()
        # Delete genomes
        print("Deleting genomes...")
        Genome.objects.all().delete()
        # Delete strains
        print("Deleting strains...")
        Strain_metadata.objects.all().delete()
        Strain.objects.all().delete()
        # Delete samples
        print("Deleting samples...")
        Sample_metadata.objects.all().delete()
        Sample.objects.all().delete()
        print("...done.")
