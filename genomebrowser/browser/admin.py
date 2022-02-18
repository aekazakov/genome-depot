from django.contrib import admin
# Import your models here
from browser.models import Strain, Sample, Genome, Contig, Gene, Taxon, Cog_class, Kegg_reaction, Kegg_pathway, Kegg_ortholog, Go_term, Cazy_family, Ec_number, Ortholog_group, Eggnog_description, Tc_family, Strain_metadata, Sample_metadata, Protein, Annotation, Operon, Site, Regulon, Config, ChangeLog


# Register your models here.
admin.site.register(Strain)
admin.site.register(Sample)
admin.site.register(Genome)
admin.site.register(Gene)
admin.site.register(Taxon)
admin.site.register(Cog_class)
admin.site.register(Kegg_reaction)
admin.site.register(Kegg_pathway)
admin.site.register(Kegg_ortholog)
admin.site.register(Go_term)
admin.site.register(Cazy_family)
admin.site.register(Ec_number)
admin.site.register(Ortholog_group)
admin.site.register(Eggnog_description)
admin.site.register(Tc_family)
admin.site.register(Strain_metadata)
admin.site.register(Sample_metadata)
admin.site.register(Protein)
admin.site.register(Annotation)
admin.site.register(Operon)
admin.site.register(Site)
admin.site.register(Regulon)
admin.site.register(Config)
admin.site.register(Contig)
admin.site.register(ChangeLog)
