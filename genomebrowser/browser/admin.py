import os
import csv
import shutil
from django import forms
from django.contrib import admin
from django.urls import path
from django.shortcuts import render
from django.shortcuts import redirect
from django.contrib import messages
# Import your models here
from browser.models import Strain, Sample, Genome, Contig, Gene, Taxon, Cog_class, Kegg_reaction, Kegg_pathway, Kegg_ortholog, Go_term, Cazy_family, Ec_number, Ortholog_group, Eggnog_description, Tc_family, Strain_metadata, Sample_metadata, Protein, Annotation, Operon, Site, Regulon, Config, ChangeLog
from browser.dataimport.importer import Importer
from cgcmsadmin.async_tasks import test_async_task, async_import_genomes, async_delete_genomes, async_import_sample_metadata, async_import_sample_descriptions, async_update_strain_metadata, async_import_annotations, async_import_regulon
from django_q.models import OrmQ


admin.site.site_header = "CGCMS admin"
admin.site.site_title = "CGCMS Admin Portal"
admin.site.index_title = "Welcome to CGCMS administration interface"


class CsvImportForm(forms.Form):
    csv_file = forms.FileField()

class ExcelImportForm(forms.Form):
    xlsx_file = forms.FileField()

@admin.action(description = 'Test sync action')
def test_task(self, request, queryset):
    test_async_task(request, queryset)
    
@admin.action(description = 'Delete genomes')
def delete_genomes(self, request, queryset):
    task_id = async_delete_genomes(request, queryset)
    messages.info(request, "Selected genomes are being deleted. Check queued task list for progress.")

def count_tasks(request):
    return str(OrmQ.objects.all().count())

# Register your models here.
class GenomeAdmin(admin.ModelAdmin):
    change_list_template = 'admin/genome_change_list.html'
    #actions = [delete_genome]
    actions = [test_task,
               delete_genomes]
    list_display = ['name', 'taxon', 'strain', 'sample', 'size', 'contigs']
    list_filter = (('strain', admin.EmptyFieldListFilter), ('sample', admin.EmptyFieldListFilter))
    ordering = ['name']
    search_fields = ['name']
    fields = (
        ('name', 'taxon'),
        ('strain', 'sample'),
        ('contigs', 'size', 'genes'),
        'description',
        ('json_url', 'pub_date'),
        ('external_id', 'external_url'),
        'gbk_filepath'
    )
    
    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('add/', self.import_genomes),
            path('delete-genomes/', self.delete_genomes),
        ]
        return my_urls + urls

    def import_genomes(self, request):
        if request.method == 'POST':
            print(request.FILES)
            csv_file = request.FILES["csv_file"]
            lines = []
            for line in csv_file:
                line = line.decode()
                print(line)
                lines.append(line)
            task_name = async_import_genomes(lines)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = CsvImportForm()
        payload = {'form': form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/csv_form.html", payload
        )
        
    def delete_genomes(self, request, queryset):
        self.message_user(request, "You asked for removal of " + str(len(queryset)) + " genomes")
        return redirect("..")

    def get_actions(self, request):
        actions = super().get_actions(request)
        if 'delete_selected' in actions:
            del actions['delete_selected']
        return actions
   
    
admin.site.register(Genome, GenomeAdmin)


class GeneAdmin(admin.ModelAdmin):
    list_display = ['locus_tag', 'genome', 'contig', 'start', 'end']
    list_filter = ['type']
    ordering = ['locus_tag', 'genome']
    search_fields = ['locus_tag', 'genome__name']

admin.site.register(Gene, GeneAdmin)


class TaxonAdmin(admin.ModelAdmin):
    list_display = ['admin_name', 'parent_id', 'eggnog_taxid']
    ordering = ['name']
    search_fields = ['name']

admin.site.register(Taxon, TaxonAdmin)


class StrainAdmin(admin.ModelAdmin):
    list_display = ['strain_id', 'full_name', 'order', 'taxon']
    ordering = ['strain_id']
    search_fields = ['strain_id', 'full_name']

admin.site.register(Strain, StrainAdmin)


class SampleAdmin(admin.ModelAdmin):
    list_display = ['sample_id', 'full_name']
    ordering = ['sample_id']
    search_fields = ['sample_id', 'full_name']

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import-descriptions/', self.import_sample_descriptions),
        ]
        return my_urls + urls

    def import_sample_descriptions(self, request):
        if request.method == 'POST':
            print(request.FILES)
            csv_file = request.FILES["csv_file"]
            lines = []
            for line in csv_file:
                line = line.decode()
                print(line)
                lines.append(line)
            task_name = async_import_sample_descriptions(lines)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = CsvImportForm()
        payload = {"form": form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_sample_descriptions_form.html", payload
        )

admin.site.register(Sample, SampleAdmin)


class CogClassAdmin(admin.ModelAdmin):
    list_display = ['cog_id', 'description']
    ordering = ['cog_id']
    search_fields = ['cog_id', 'description']

admin.site.register(Cog_class, CogClassAdmin)


class KeggReactionAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']
    search_fields = ['kegg_id', 'description']

admin.site.register(Kegg_reaction, KeggReactionAdmin)


class KeggPathwayAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']
    search_fields = ['kegg_id', 'description']

admin.site.register(Kegg_pathway, KeggPathwayAdmin)


class KeggOrthologAdmin(admin.ModelAdmin):
    list_display = ['kegg_id', 'description']
    ordering = ['kegg_id']
    search_fields = ['kegg_id', 'description']

admin.site.register(Kegg_ortholog, KeggOrthologAdmin)


class GoTermAdmin(admin.ModelAdmin):
    list_display = ['go_id', 'go_namespace', 'description']
    ordering = ['go_id']
    search_fields = ['go_id', 'go_namespace', 'description']

admin.site.register(Go_term, GoTermAdmin)


class CazyFamilyAdmin(admin.ModelAdmin):
    list_display = ['cazy_id', 'description']
    ordering = ['cazy_id']
    search_fields = ['cazy_id', 'description']

admin.site.register(Cazy_family, CazyFamilyAdmin)


class EcNumberAdmin(admin.ModelAdmin):
    list_display = ['ec_number', 'description']
    ordering = ['ec_number']
    search_fields = ['ec_number', 'description']

admin.site.register(Ec_number, EcNumberAdmin)


class TcFamilyAdmin(admin.ModelAdmin):
    list_display = ['tc_id', 'description']
    ordering = ['tc_id']
    search_fields = ['tc_id', 'description']

admin.site.register(Tc_family, TcFamilyAdmin)


class OrthologGroupAdmin(admin.ModelAdmin):
    list_display = ['eggnog_id', 'taxon']
    ordering = ['eggnog_id']
    search_fields = ['eggnog_id', 'taxon']

admin.site.register(Ortholog_group, OrthologGroupAdmin)


class EggnogDescriptionAdmin(admin.ModelAdmin):
    list_display = ['fingerprint', 'description']
    ordering = ['fingerprint']
    search_fields = ['fingerprint', 'description']

admin.site.register(Eggnog_description, EggnogDescriptionAdmin)


class StrainMetadataAdmin(admin.ModelAdmin):
    list_display = ['strain', 'key', 'source']
    ordering = ['strain', 'key']
    search_fields = ['strain', 'key', 'source']

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('add/', self.update_strain_metadata),
        ]
        return my_urls + urls

    def update_strain_metadata(self, request):
        if request.method == 'POST':
            print(request.FILES)
            xlsx_file = request.FILES["xlsx_file"]
            task_name = async_update_strain_metadata(xlsx_file)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = ExcelImportForm()
        payload = {"form": form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_strain_metadata_form.html", payload
        )

admin.site.register(Strain_metadata, StrainMetadataAdmin)


class SampleMetadataAdmin(admin.ModelAdmin):
    list_display = ['sample', 'key', 'source']
    ordering = ['sample', 'key']
    search_fields = ['sample', 'key', 'source']

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('add/', self.import_sample_metadata),
        ]
        return my_urls + urls

    def import_sample_metadata(self, request):
        if request.method == 'POST':
            print(request.FILES)
            csv_file = request.FILES["csv_file"]
            lines = []
            for line in csv_file:
                line = line.decode()
                print(line)
                lines.append(line)
            task_name = async_import_sample_metadata(lines)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = CsvImportForm()
        payload = {"form": form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_sample_metadata_form.html", payload
        )

admin.site.register(Sample_metadata, SampleMetadataAdmin)


class ProteinAdmin(admin.ModelAdmin):
    list_display = ['name', 'protein_hash', 'length']
    ordering = ['protein_hash']
    search_fields = ['name', 'protein_hash', 'length']

admin.site.register(Protein, ProteinAdmin)


class AnnotationAdmin(admin.ModelAdmin):
    list_display = ['gene_id', 'key', 'value', 'source']
    ordering = ['gene_id', 'key']
    search_fields = ['gene_id__locus_tag', 'key', 'value', 'source']

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import-annotations/', self.import_annotations),
        ]
        return my_urls + urls
    
    def import_annotations(self, request):
        if request.method == 'POST':
            print(request.FILES)
            csv_file = request.FILES["csv_file"]
            lines = []
            for line in csv_file:
                line = line.decode()
                print(line)
                lines.append(line)
            task_name = async_import_annotations(lines)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = CsvImportForm()
        payload = {"form": form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_annotations_form.html", payload
        )

admin.site.register(Annotation, AnnotationAdmin)


class OperonAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome', 'contig', 'start', 'end']
    ordering = ['genome', 'name']
    search_fields = ['name', 'genome', 'contig']

admin.site.register(Operon, OperonAdmin)


class SiteAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome', 'contig', 'start', 'end', 'type']
    ordering = ['genome', 'name']
    search_fields = ['name', 'genome', 'contig', 'type']

admin.site.register(Site, SiteAdmin)


class RegulonAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome']
    ordering = ['genome', 'name']
    search_fields = ['name', 'genome']

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import-regulons/', self.import_regulons),
        ]
        return my_urls + urls

    def import_regulons(self, request):
        if request.method == 'POST':
            print(request.FILES)
            csv_file = request.FILES["csv_file"]
            lines = []
            for line in csv_file:
                line = line.decode()
                print(line)
                lines.append(line)
            task_name = async_import_regulon(lines)
            # Do some staff
            self.message_user(request, "Your file was submitted for the processing with ID " + task_name)
            return redirect("..")
        form = CsvImportForm()
        payload = {"form": form}
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_regulons_form.html", payload
        )

admin.site.register(Regulon, RegulonAdmin)


class ConfigAdmin(admin.ModelAdmin):
    list_display = ['param', 'value']
    ordering = ['param']
    search_fields = ['param', 'value']

admin.site.register(Config, ConfigAdmin)


class ContigAdmin(admin.ModelAdmin):
    list_display = ['contig_id', 'name', 'genome', 'size']
    ordering = ['genome', 'contig_id']
    search_fields = ['contig_id', 'name', 'genome']

admin.site.register(Contig, ContigAdmin)


class ChangeLogAdmin(admin.ModelAdmin):
    list_display = ['action', 'timestamp']
    list_filter = ['timestamp']
    ordering = ['timestamp']

admin.site.register(ChangeLog, ChangeLogAdmin)

