from django.core.management.base import BaseCommand
from browser.models import Annotation
from browser.models import Eggnog_description
from browser.models import Cog_class
from browser.models import Kegg_reaction
from browser.models import Kegg_pathway
from browser.models import Kegg_ortholog
from browser.models import Go_term
from browser.models import Cazy_family
from browser.models import Ec_number
from browser.models import Tc_family
from browser.models import Ortholog_group
from browser.models import Gene
from browser.models import Protein
from browser.models import Genome
from browser.models import Strain_metadata
from browser.models import Strain
from browser.models import Sample_metadata
from browser.models import Sample
from browser.models import Taxon

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
