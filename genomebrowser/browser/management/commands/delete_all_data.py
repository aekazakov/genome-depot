from django.core.management.base import BaseCommand
from browser.models import *

class Command(BaseCommand):
    help = 'Deletes all genome browser entries from the database'
    def handle(self, *args, **options):
        # Delete all mappings before deleting genes
        Annotation.objects.all().delete()
        Eggnog_description.objects.all().delete()
        Cog_class.objects.all().delete()
        Kegg_reaction.objects.all().delete()
        Kegg_pathway.objects.all().delete()
        Kegg_ortholog.objects.all().delete()
        Go_term.objects.all().delete()
        Cazy_family.objects.all().delete()
        Ec_number.objects.all().delete()
        Tc_family.objects.all().delete()
        Ortholog_group.objects.all().delete()
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
        # Delete taxonomy entries
        Taxon.objects.all().delete()
