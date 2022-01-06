from django.core.management.base import BaseCommand
from browser.models import *

class Command(BaseCommand):
    help = 'Checks if any genomes have many genes without egg-nog mappings'
    
    def handle(self, *args, **options):
        genomes = Genome.objects.select_related('strain')
        #annotated_proteins = Protein.objects.exclude(ortholog_groups = None).distinct()
        with open('genomes.txt', 'w') as outfile:
            for genome in genomes:
                gene_count = Gene.objects.filter(genome = genome).count()
                outfile.write('\t'.join([genome.name, genome.strain.strain_id, genome.strain.full_name, genome.external_id, genome.external_url, str(genome.size), str(genome.contigs), str(genome.genes), str(gene_count), str(genome.pub_date)]) + '\n')

