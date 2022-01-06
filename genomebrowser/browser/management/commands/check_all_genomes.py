from django.core.management.base import BaseCommand
from browser.models import *

class Command(BaseCommand):
    help = 'Checks if any genomes have many genes without egg-nog mappings'
    
    def handle(self, *args, **options):
        genomes = Genome.objects.all()
        #annotated_proteins = Protein.objects.exclude(ortholog_groups = None).distinct()
        for genome in genomes:
            gene_count = Gene.objects.filter(genome = genome).exclude(protein__ortholog_groups=None).count()
            if gene_count < 1000:
                print(genome.name + '\t' + str(gene_count))

