import os
import uuid
import datetime
import zipfile
import logging
from pathlib import Path
from django.contrib import admin
from django.urls import path
from django.shortcuts import render
from django.shortcuts import redirect
from django.contrib import messages
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.urls import reverse_lazy
from django.db.models import Q
from django.forms import CheckboxSelectMultiple
from django_q.models import OrmQ
from django_q.monitor import Stat
# Import your models here
from browser.models import Strain
from browser.models import Sample
from browser.models import Genome
from browser.models import Contig
from browser.models import Gene
from browser.models import Taxon
from browser.models import Cog_class
from browser.models import Kegg_reaction
from browser.models import Kegg_pathway
from browser.models import Kegg_ortholog
from browser.models import Go_term
from browser.models import Cazy_family
from browser.models import Ec_number
from browser.models import Ortholog_group
from browser.models import Eggnog_description
from browser.models import Tc_family
from browser.models import Strain_metadata
from browser.models import Sample_metadata
from browser.models import Protein
from browser.models import Annotation
from browser.models import Operon
from browser.models import Site
from browser.models import Regulon
from browser.models import Config
from browser.models import ChangeLog
from browser.models import Tag
from browser.forms import TsvImportForm
from browser.forms import ExcelImportForm
from browser.forms import GenomeImportForm
from browser.forms import GenomeUploadForm
from browser.forms import GenomeDownloadForm
from browser.forms import TagModelForm
from browser.forms import AddTagForm
from browser.forms import ChooseAnnotationToolForm
from browser.async_tasks import async_import_genomes
from browser.async_tasks import async_delete_genomes
from browser.async_tasks import async_update_static_files
from browser.async_tasks import async_import_sample_metadata
from browser.async_tasks import async_import_sample_descriptions
from browser.async_tasks import async_update_strain_metadata
from browser.async_tasks import async_import_annotations
from browser.async_tasks import async_import_regulon
from browser.async_tasks import async_run_annotation_tools
from genomebrowser.settings import TITLE

logger = logging.getLogger("GenomeDepot")
admin.site.site_title = ""
admin.site.site_header = "GenomeDepot administration"
admin.site.index_title = TITLE + " administration"

@admin.action(description = 'Run annotation tools')
def run_annotation_tools(self, request, queryset):
    if 'do_action' in request.POST:
        form = ChooseAnnotationToolForm(request.POST)
        # this hack is necessary because the tools field is empty on form initialization to prevent calling db server on startup
        form.fields.get('tools').choices = [(item,item) for item in request.POST.getlist('tools')]
        if form.is_valid():
            tools = ['.'.join(item.split('.')[1:-1]) for item in form.cleaned_data['tools']]
            task_id = async_run_annotation_tools(request, queryset, tools)
            messages.info(request,
                          "The annotation pipeline is running for selected genomes. " +
                          "Check the status of the task " + task_id + 
                          " in the list of queued tasks for progress."
                          )
        else:
            logger.error(form.errors)
            
        return HttpResponseRedirect(request.get_full_path())
    else:
        plugins_enabled = [item.replace('.enabled', '.display_name') for item in Config.objects.filter(Q(param__startswith='plugins.')&Q(param__endswith='.enabled')&Q(value__in=('1','yes','Yes','y','Y'))).values_list('param', flat=True)]
        form = ChooseAnnotationToolForm()
        form.fields.get('tools').choices = Config.objects.filter(param__in=plugins_enabled).values_list('param','value')
        context = {'title': u'Choose tools', 'form': form}
        active_tasks = OrmQ.objects.all().count()
        context['genomes'] = queryset
        context['active_tasks'] = str(active_tasks)
        return render(request,
            'admin/choose_annotation_tool_form.html',
            context
        )


@admin.action(description = 'Delete genomes')
def delete_genomes(self, request, queryset):
    task_id = async_delete_genomes(request, queryset)
    messages.info(request,
                  "Selected genomes are being deleted. " +
                  "Check the status of the task " + task_id + 
                  " in the list of queued tasks for progress."
                  )


@admin.action(description = 'Update static files and re-build search databases')
def update_static_files(self, request, queryset):
    task_id = async_update_static_files(request, queryset)
    messages.info(request,
                  "Static files for selected genomes are being re-created. " +
                  "Check the status of the task " + task_id + 
                  " in the list of queued tasks for progress."
                  )


@admin.action(description = 'Add a tag')
def add_genome_tag(self, request, queryset):
    if 'do_action' in request.POST:
        form = AddTagForm(request.POST)
        if form.is_valid():
            tag = form.cleaned_data['tag']
            updated = 0
            for genome in queryset:
                genome.tags.add(tag)
                updated += 1
            messages.success(request, '{0} genomes were updated'.format(updated))
        else:
            logger.error(form.errors)
        redirect("..")
    else:
        form = AddTagForm()
        return render(request, 'admin/add_tag.html',
            {'title': u'Add tag',
             'objects': queryset,
             'form': form})


@admin.action(description = 'Remove a tag')
def remove_genome_tag(self, request, queryset):
    if 'do_action' in request.POST:
        form = AddTagForm(request.POST)
        if form.is_valid():
            tag = form.cleaned_data['tag']
            updated = 0
            for genome in queryset:
                genome.tags.remove(tag)
                updated += 1
            messages.success(request, '{0} genomes were updated'.format(updated))
        else:
            logger.error(form.errors)
        redirect("..")
    else:
        form = AddTagForm()
        return render(request, 'admin/remove_tag.html',
            {'title': u'Remove tag',
             'objects': queryset,
             'form': form})


def count_tasks(request):
    return str(OrmQ.objects.all().count())


def count_clusters(request):
    result = 0
    for stat in Stat.get_all():
        result += 1
    return str(result)


def pipeline_status(request):
    result = 'Off'
    workers = 0
    tasks = 0
    for stat in Stat.get_all():
        workers += 1
    if workers > 0:
        tasks = OrmQ.objects.all().count()
        if tasks > 0:
            result = 'Busy'
        else:
            result = 'Idle'
    return result


# Register your models here.
class GenomeAdmin(admin.ModelAdmin):
    change_list_template = 'admin/genome_change_list.html'
    #actions = [delete_genome]
    actions = [add_genome_tag,
               remove_genome_tag,
               run_annotation_tools,
               delete_genomes,
               update_static_files
               ]
    list_display = ['name',
                    'taxon',
                    'strain',
                    'sample',
                    'size',
                    'contigs',
                    'get_tags'
                    ]
    list_filter = (('strain', admin.EmptyFieldListFilter),
                   ('sample', admin.EmptyFieldListFilter),
                   ('tags')
                   )
    ordering = ['name']
    search_fields = ['name']
    fields = (
        ('name', 'taxon'),
        ('strain', 'sample'),
        ('contigs', 'size', 'genes'),
        'description',
        ('json_url', 'pub_date'),
        ('external_id', 'external_url'),
        'gbk_filepath',
        'tags'
    )
    readonly_fields=('name', )
    autocomplete_fields = ('strain', 'sample', 'taxon', )
    
    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('add/', self.import_genomes),
            path('delete-genomes/', self.delete_genomes),
            path('update-static/', self.update_static_files),
        ]
        return my_urls + urls

    def import_genomes(self, request):
        import_form = None
        upload_form = None
        download_form = None
        if request.method == 'POST':
            logger.debug(request.FILES)
            lines = []
            if request.POST['choice_field'] == 'import':
                import_form = GenomeImportForm(prefix='import', data=request.POST)
                if import_form.is_valid():
                    tsv_file = request.FILES.get("import_tsv_file", False)
                    if not tsv_file:
                        self.message_user(request,
                                          "Genome list file was not submitted",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                    for line in tsv_file:
                        line = line.decode()
                        line = line.rstrip('\n\r')
                        if line == '':
                            continue
                        row = line.split('\t')
                        if row[0].startswith('#'):
                            pass
                        elif not os.path.exists(row[0]):
                            logger.warning('%s not found in the filesystem', row[0])
                        lines.append('\t'.join(row))
                    if not lines:
                        self.message_user(request,
                                          "No valid entries in the list of genomes",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                else:
                    self.message_user(request,
                                      "Form is not valid",
                                      level=messages.ERROR
                                      )
                    return redirect("..")
            elif request.POST['choice_field'] == 'upload':
                upload_form = GenomeUploadForm(prefix='upload', data=request.POST)
                if upload_form.is_valid():
                    tsv_file = request.FILES.get("upload_tsv_file", False)
                    if not tsv_file:
                        self.message_user(request,
                                          "Genome list file was not submitted",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                    zip_file = request.FILES.get("zip_file",False)
                    if not zip_file:
                        self.message_user(request,
                                          "Zip file with genomes was not submitted",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                    temp_dir = Config.objects.get(param='core.temp_dir').value
                    upload_dir = os.path.join(temp_dir, str(uuid.uuid4()))
                    logger.debug(upload_dir)
                    Path(upload_dir).mkdir(parents=True, exist_ok=True)
                    os.chmod(upload_dir, 0o777)
                    zip_content = handle_zip_upload(zip_file, upload_dir)
                    for line in tsv_file:
                        line = line.decode()
                        line = line.rstrip('\n\r')
                        row = line.split('\t')
                        if row[0].startswith('#'):
                            pass
                        elif row[0] in zip_content:
                            row[0] = zip_content[row[0]]
                        elif row[0] == '':
                            logger.warning('File name is empty in the line: ' + line)
                        else:
                            logger.warning('%s not found in the zip archive', row[0])
                        lines.append('\t'.join(row))
                    if not lines:
                        self.message_user(request,
                                          "No valid entries in the list of genomes",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                else:
                    self.message_user(request,
                                      "Form is not valid",
                                      level=messages.ERROR
                                      )
                    return redirect("..")
            elif request.POST['choice_field'] == 'download':
                download_form = GenomeDownloadForm(prefix='download', data=request.POST)
                if download_form.is_valid():
                    tsv_file = request.FILES.get("download_tsv_file", False)
                    if not tsv_file:
                        self.message_user(request,
                                          "Genome list file was not submitted",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                    download_email = request.POST.get("download_email", False)
                    if not download_email:
                        self.message_user(request,
                                          "Email address required",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                    for line in tsv_file:
                        line = line.decode()
                        line = line.rstrip('\n\r')
                        if line == '':
                            continue
                        row = line.split('\t')
                        if row[0].startswith('#'):
                            pass
                        elif row[0] == '' and row[-1].startswith('NCBI:'):
                            lines.append('\t'.join(row))
                            # Assembly will be downloaded from NCBI later
                        elif not os.path.exists(row[0]):
                            logger.warning('%s not found', row[0])
                            #lines.append('\t'.join(row))
                    if not lines:
                        self.message_user(request,
                                          "No valid entries in the list of genomes",
                                          level=messages.ERROR
                                          )
                        return redirect("..")
                else:
                    self.message_user(request,
                                      "Form is not valid",
                                      level=messages.ERROR
                                      )
                    return redirect("..")
            else:
                self.message_user(request,
                                  "You must choose one of genome import options",
                                  level=messages.ERROR
                                  )
                return redirect("..")
            task_name = async_import_genomes(lines,
                                             request.POST['download_email']
                                             )
            # Notify the user
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
            return redirect("..")
        import_form = GenomeImportForm()
        upload_form = GenomeUploadForm()
        download_form = GenomeDownloadForm()
        payload = {
            'import_form': import_form,
            'upload_form': upload_form,
            'download_form': download_form
        }
        active_tasks = OrmQ.objects.all().count()
        payload['active_tasks'] = str(active_tasks)
        return render(
            request, "admin/import_genomes.html", payload
        )

    def delete_genomes(self, request, queryset):
        self.message_user(request,
                          "You asked for removal of " +
                          str(len(queryset)) +
                          " genomes")
        return redirect("..")

    def update_static_files(self, request, queryset):
        self.message_user(request,
                          "You asked for re-creating static files of " +
                          str(len(queryset)) +
                          " genomes")
        return redirect("..")
        
    def get_actions(self, request):
        actions = super().get_actions(request)
        if 'delete_selected' in actions:
            del actions['delete_selected']
        return actions
   
    def change_view(self, request, object_id=None, form_url='', extra_context=None):
        return super().change_view(request, object_id, form_url,
                               extra_context=dict(show_delete=False))    


admin.site.register(Genome, GenomeAdmin)


class GeneAdmin(admin.ModelAdmin):
    list_display = ['locus_tag', 'genome', 'contig', 'start', 'end']
    list_filter = ['type']
    ordering = ['locus_tag', 'genome']
    search_fields = ['locus_tag', 'genome__name']
    readonly_fields=('locus_tag', )
    autocomplete_fields = ('genome', 'contig', 'protein', 'operon', )

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
    autocomplete_fields = ('taxon',)

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
            logger.debug(request.FILES)
            tsv_file = request.FILES["tsv_file"]
            lines = []
            for line in tsv_file:
                line = line.decode()
                logger.debug(line)
                lines.append(line)
            task_name = async_import_sample_descriptions(lines)
            # Do some staff
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
            return redirect("..")
        form = TsvImportForm()
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
    search_fields = ['eggnog_id', 'taxon__name']
    autocomplete_fields = ('taxon', )
    
admin.site.register(Ortholog_group, OrthologGroupAdmin)


class EggnogDescriptionAdmin(admin.ModelAdmin):
    list_display = ['fingerprint', 'description']
    ordering = ['fingerprint']
    search_fields = ['fingerprint', 'description']

admin.site.register(Eggnog_description, EggnogDescriptionAdmin)


class StrainMetadataAdmin(admin.ModelAdmin):
    list_display = ['strain', 'key', 'source']
    ordering = ['strain', 'key']
    search_fields = ['strain__strain_id', 'key', 'source']
    autocomplete_fields = ('strain',)

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import/', self.admin_site.admin_view(self.update_strain_metadata)),
        ]
        return my_urls + urls

    def update_strain_metadata(self, request):
        if request.method == 'POST':
            logger.debug(request.FILES)
            xlsx_file = request.FILES["xlsx_file"]
            task_name = async_update_strain_metadata(xlsx_file)
            # Do some staff
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
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
    search_fields = ['sample__sample_id', 'key', 'source']
    autocomplete_fields = ('sample',)

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import/', self.admin_site.admin_view(self.import_sample_metadata)),
        ]
        return my_urls + urls

    def import_sample_metadata(self, request):
        if request.method == 'POST':
            logger.debug(request.FILES)
            tsv_file = request.FILES["tsv_file"]
            lines = []
            for line in tsv_file:
                line = line.decode()
                logger.debug(line)
                lines.append(line)
            task_name = async_import_sample_metadata(lines)
            # Do some staff
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
            return redirect("..")
        form = TsvImportForm()
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
    autocomplete_fields = ('taxonomy_id', )
    readonly_fields=('protein_hash', )
    raw_id_fields = ('ortholog_groups',
                     'cog_classes',
                     'kegg_reactions',
                     'kegg_pathways',
                     'kegg_orthologs',
                     'go_terms',
                     'cazy_families',
                     'ec_numbers',
                     'tc_families',
                     )
    

admin.site.register(Protein, ProteinAdmin)


class AnnotationAdmin(admin.ModelAdmin):
    list_display = ['gene_id', 'key', 'value', 'source']
    ordering = ['gene_id', 'key']
    search_fields = ['gene_id__locus_tag', 'key', 'value', 'source']
    autocomplete_fields = ('gene_id',)

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import-annotations/', self.import_annotations),
        ]
        return my_urls + urls
    
    def import_annotations(self, request):
        if request.method == 'POST':
            logger.debug(request.FILES)
            tsv_file = request.FILES["tsv_file"]
            lines = []
            for line in tsv_file:
                line = line.decode()
                logger.debug(line)
                lines.append(line)
            task_name = async_import_annotations(lines)
            # Do some staff
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
            return redirect("..")
        form = TsvImportForm()
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
    search_fields = ['name', 'genome__name', 'contig__contig_id']
    autocomplete_fields = ('genome', 'contig', )

admin.site.register(Operon, OperonAdmin)


class SiteAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome', 'contig', 'start', 'end', 'type']
    ordering = ['genome', 'name']
    search_fields = ['name', 'genome__name', 'contig__contig_id', 'type']
    autocomplete_fields = ('genome', 'contig', 'regulon',)
    raw_id_fields = ('operons', 'genes',)

admin.site.register(Site, SiteAdmin)


class RegulonAdmin(admin.ModelAdmin):
    list_display = ['name', 'genome']
    ordering = ['genome', 'name']
    search_fields = ['name', 'genome']
    autocomplete_fields = ('genome', )
    raw_id_fields = ('regulators', )

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('import-regulons/', self.import_regulons),
        ]
        return my_urls + urls

    def import_regulons(self, request):
        if request.method == 'POST':
            logger.debug(request.FILES)
            tsv_file = request.FILES["tsv_file"]
            lines = []
            for line in tsv_file:
                line = line.decode()
                logger.debug(line)
                lines.append(line)
            task_name = async_import_regulon(lines)
            # Do some staff
            self.message_user(request,
                              "Your file was submitted for the processing with ID " +
                              task_name
                              )
            return redirect("..")
        form = TsvImportForm()
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
    search_fields = ['contig_id', 'name', 'genome__name']
    autocomplete_fields = ('genome', )
    
admin.site.register(Contig, ContigAdmin)


class ChangeLogAdmin(admin.ModelAdmin):
    list_display = ['action', 'timestamp']
    list_filter = ['timestamp']
    ordering = ['timestamp']

admin.site.register(ChangeLog, ChangeLogAdmin)


class TagAdmin(admin.ModelAdmin):
    list_display = ['name', 'description', 'color', 'textcolor']
    list_filter = ['name']
    ordering = ['name']
    form = TagModelForm

admin.site.register(Tag, TagAdmin)

@login_required(login_url=reverse_lazy("admin:login"))
def clusters_view(request):
    cluster_count = 0
    clusters = []
    for stat in Stat.get_all():
        cluster_count += 1
        clusters.append({'id':str(stat.cluster_id),
                         'tob': str(stat.tob),
                         'uptime': str(datetime.timedelta(seconds=stat.uptime())),
                         'workers': ','.join([str(x) for x in stat.workers])
                         })
    if cluster_count > 0:
        context = {'cluster_count':str(cluster_count), 'clusters':clusters}
    else:
        context = {'cluster_count':'Off', 'clusters':clusters}
    return render(
        request, "admin/clusters.html", context
    )

@login_required(login_url=reverse_lazy("admin:login"))
def tools_view(request):
    context = {}
    context['tasks'] = count_tasks(request)
    context['clusters'] = pipeline_status(request)
    return render(
        request, "admin/tools.html", context
    )

def handle_zip_upload(zipf, upload_dir):
    result = {}
    with zipfile.ZipFile(zipf, "r", zipfile.ZIP_STORED) as openzip:
        filelist = openzip.infolist()
        for member in filelist:
            if member.is_dir():
                continue
            openzip.extract(member.filename, path=upload_dir)
            result[str(member.filename)] = os.path.join(upload_dir, member.filename)
            os.chmod(os.path.join(upload_dir, member.filename), 0o777)
    return result

