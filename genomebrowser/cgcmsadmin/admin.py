from django.contrib import admin
# Import your models here
from browser.models import Strain, Sample, Genome, Contig, Gene, Taxon, Cog_class, Kegg_reaction, Kegg_pathway, Kegg_ortholog, Go_term, Cazy_family, Ec_number, Ortholog_group, Eggnog_description, Tc_family, Strain_metadata, Sample_metadata, Protein, Annotation, Operon, Site, Regulon, Config, ChangeLog


class CgcmsAdminSite(admin.AdminSite):
    site_header = "CGCMS admin"
    site_title = "CGCMS Admin Portal"
    index_title = "Welcome to CGCMS administration interface"


cgcms_admin_site = CgcmsAdminSite(name='cgcms_admin')

# Register your models here.
class GenomeAdmin(admin.ModelAdmin):
    list_display = ['name', 'taxon', 'strain', 'size', 'contigs']
    ordering = ['name']
    
cgcms_admin_site.register(Genome, GenomeAdmin)

class GeneAdmin(admin.ModelAdmin):
    list_display = ['locus_tag', 'genome', 'contig', 'start', 'end']
    ordering = ['locus_tag', 'genome']

cgcms_admin_site.register(Gene, GeneAdmin)

class TaxonAdmin(admin.ModelAdmin):
    list_display = ['admin_name', 'parent_id', 'eggnog_taxid']
    ordering = ['name']

cgcms_admin_site.register(Taxon, TaxonAdmin)


class StrainAdmin(admin.ModelAdmin):
    list_display = ['strain_id', 'full_name', 'order', 'taxon']
    ordering = ['strain_id']

cgcms_admin_site.register(Strain, StrainAdmin)


class SampleAdmin(admin.ModelAdmin):
    list_display = ['sample_id', 'full_name']
    ordering = ['sample_id']

cgcms_admin_site.register(Sample, SampleAdmin)


class CogClassAdmin(admin.ModelAdmin):
    list_display = ['cog_id', 'description']
    ordering = ['cog_id']

cgcms_admin_site.register(Cog_class, CogClassAdmin)


class KeggReactionAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']

cgcms_admin_site.register(Kegg_reaction, KeggReactionAdmin)


class KeggPathwayAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']

cgcms_admin_site.register(Kegg_pathway, KeggPathwayAdmin)


class KeggOrthologAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']

cgcms_admin_site.register(Kegg_ortholog, KeggOrthologAdmin)


class GoTermAdmin(admin.ModelAdmin):
    list_display = ['go_id', 'go_namespace', 'description']
    ordering = ['go_id']

cgcms_admin_site.register(Go_term, GoTermAdmin)


class CazyFamilyAdmin(admin.ModelAdmin):
    list_display = ['cazy_id', 'description']
    ordering = ['cazy_id']

cgcms_admin_site.register(Cazy_family, CazyFamilyAdmin)


class EcNumberAdmin(admin.ModelAdmin):
    list_display = ['ec_number', 'description']
    ordering = ['ec_number']

cgcms_admin_site.register(Ec_number, EcNumberAdmin)


class TcFamilyAdmin(admin.ModelAdmin):
    list_display = ['tc_id', 'description']
    ordering = ['tc_id']

cgcms_admin_site.register(Tc_family, TcFamilyAdmin)


class OrthologGroupAdmin(admin.ModelAdmin):
    list_display = ['eggnog_id', 'taxon']
    ordering = ['eggnog_id']

cgcms_admin_site.register(Ortholog_group, OrthologGroupAdmin)


class EggnogDescriptionAdmin(admin.ModelAdmin):
    list_display = ['fingerprint', 'description']
    ordering = ['fingerprint']

cgcms_admin_site.register(Eggnog_description, EggnogDescriptionAdmin)


class StrainMetadataAdmin(admin.ModelAdmin):
    list_display = ['strain', 'key', 'source']
    ordering = ['strain', 'key']

cgcms_admin_site.register(Strain_metadata, StrainMetadataAdmin)


class SampleMetadataAdmin(admin.ModelAdmin):
    list_display = ['sample', 'key', 'source']
    ordering = ['sample', 'key']

cgcms_admin_site.register(Sample_metadata, SampleMetadataAdmin)


class ProteinAdmin(admin.ModelAdmin):
    list_display = ['name', 'protein_hash', 'length']
    ordering = ['protein_hash']

cgcms_admin_site.register(Protein, ProteinAdmin)


class AnnotationAdmin(admin.ModelAdmin):
    list_display = ['gene_id', 'key', 'value', 'source']
    ordering = ['gene_id', 'key']

cgcms_admin_site.register(Annotation, AnnotationAdmin)


class OperonAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome', 'contig', 'start', 'end']
    ordering = ['genome', 'name']

cgcms_admin_site.register(Operon, OperonAdmin)


class SiteAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome', 'contig', 'start', 'end', 'type']
    ordering = ['genome', 'name']

cgcms_admin_site.register(Site, SiteAdmin)


class RegulonAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome']
    ordering = ['genome', 'name']

cgcms_admin_site.register(Regulon, RegulonAdmin)


class ConfigAdmin(admin.ModelAdmin):
    list_display = ['param', 'value']
    ordering = ['param']

cgcms_admin_site.register(Config, ConfigAdmin)


class ContigAdmin(admin.ModelAdmin):
    list_display = ['contig_id', 'name', 'genome', 'size']
    ordering = ['genome', 'contig_id']

cgcms_admin_site.register(Contig, ContigAdmin)


class ChangeLogAdmin(admin.ModelAdmin):
    list_display = ['action', 'timestamp']
    ordering = ['timestamp']

cgcms_admin_site.register(ChangeLog, ChangeLogAdmin)
