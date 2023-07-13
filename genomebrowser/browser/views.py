import csv
import time
import json
from django.shortcuts import render
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.template import loader
from django.http import Http404
from .models import *
from django.views import View, generic
from django.db.models import Q
from django.core.files.storage import default_storage
from django.core.exceptions import SuspiciousOperation
from django.urls import reverse
from browser.seqsearch import run_protein_search, run_nucleotide_search
from browser.comparative_analysis import get_scribl
from browser.conserved_regulon import build_conserved_regulon
from browser.taxonomy import generate_sunburst, get_taxon_children
# Create your views here.

class AnnotationSearchResultsSubView(generic.ListView):
    '''
        Sub-page for AJAX-based view of annotation search result
    '''
    context_object_name = 'annotationlist'
    template_name = 'browser/annotation_list_subpage.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        '''
            Generates context string
        '''
        context = super(AnnotationSearchResultsSubView,self).get_context_data(**kwargs)
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('annotation_query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        return context

    def get_queryset(self):
        '''
            Generates Annotation queryset
        '''
        annotation_query = self.request.GET.get('annotation_query')
        genome = self.request.GET.get('genome')
        if not annotation_query:
            object_list = Annotation.objects.none()
        elif genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon').prefetch_related('gene_id__genome__tags')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon').prefetch_related('gene_id__genome__tags')
        return object_list


class AnnotationSearchResultsAjaxView(View):
    '''
        AJAX-based view for annotation search results
    '''
    def get(self,request):
        '''
            Takes a GET request and returns a webpage that will send AJAX request
        '''
        context = {} #super(AnnotationSearchResultsAjaxView,self).get_context_data(**kwargs)
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('annotation_query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        for key, val in request.GET.items():
            context[key] = val
        return render(request,'browser/annotation_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, creates sub-page (AnnotationSearchResultsSubView) and sends it back as JSON
            
            For the testing of long-running tasks, uncomment sleep_timer=0 and set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()

        # print('REQUEST2')
        # for key, val in request.GET.items():
        #     print(key, val)
            
        # Sleep timer for testing to imitate long-running task
        # sleep_timer = 0
        # print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        # time.sleep(sleep_timer)

        context = {}
        sub_view = AnnotationSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()
        context['searchcontext'] = 'Search results for "' + request.GET.get('annotation_query') + '"'
        context['searchresult'] = sub_response.content.decode('utf-8')
        context['time'] = time.time()-start_time
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


class TaxonListView(generic.ListView):
    '''
        Class-based list view of Taxon. Displays a table of taxa.
    '''
    model = Taxon
    context_object_name = 'taxalist'
    paginate_by = 50

    def get_queryset(self):
        return Taxon.objects.order_by('name')


class TaxonSearchResultsView(generic.ListView):
    '''
        Returns results of search in strain names full names or taxonomic orders.
    '''
    model = Taxon
    context_object_name = 'taxalist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(TaxonSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Taxon.objects.filter(
                Q(name__icontains=query) | Q(taxonomy_id__icontains=query)
            ).order_by('name')
        else:
            object_list = Taxon.objects.none()
        return object_list


class StrainListView(generic.ListView):
    '''
        Class-based list view of Strain. Displays a table of strains.
    '''
    model = Strain
    context_object_name = 'strainlist'
    paginate_by = 50

    def get_queryset(self):
        return Strain.objects.order_by('strain_id')


class SampleListView(generic.ListView):
    '''
        Class-based list view of Sample. Displays a table of samples.
    '''
    model = Sample
    context_object_name = 'samplelist'
    paginate_by = 50

    def get_queryset(self):
        return Sample.objects.order_by('sample_id')


class GenomeListView(generic.ListView):
    '''
        Class-based list view of Genome. Displays a table of genomes.
    '''
    model = Genome
    context_object_name = 'genomelist'
    paginate_by = 50

    def get_queryset(self):
        return Genome.objects.order_by('name').select_related('strain', 'taxon').prefetch_related('tags')

    def get_context_data(self,**kwargs):
        context = super(GenomeListView,self).get_context_data(**kwargs)
        if not 'page_obj' in context or context['page_obj'].number == 1:
            sunburst = generate_sunburst()
            if sunburst:
                context['sunburst'] = sunburst
        return context
        
        
class OperonListView(generic.ListView):
    '''
        Class-based list view of Operon. Displays a table of operons.
    '''
    model = Operon
    template_name = 'operon_list.html'
    context_object_name = 'operonlist'
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        try:
            genome = Genome.objects.get(name = self.kwargs['genome'])
        except Genome.DoesNotExist:
            raise Http404('Genome ' + self.kwargs['genome'] + ' does not exist')
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '" in ' + genome.name + ' operons'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Operon.objects.filter(genome__name=genome).filter(
                Q(genes__name__icontains=query) | Q(genes__locus_tag__icontains=query)
            ).order_by('name').select_related('genome', 'contig').prefetch_related('genes', 'genome__tags')
        else:
            object_list = Operon.objects.filter(genome__name=genome).order_by('name').select_related('genome', 'contig').prefetch_related('genes', 'genome__tags')
        return object_list

        
class SiteListView(generic.ListView):
    '''
        Class-based list view of Site. Displays a table of sites.
    '''
    model = Site
    template_name = 'site_list.html'
    context_object_name = 'sitelist'
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        genome_name = self.kwargs['genome']
        genome = Genome.objects.get(name = genome_name)
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '" in ' + genome.name + ' sites'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Site.objects.filter(genome__name=genome).filter(
                Q(name__icontains=query) | Q(regulon__name__icontains=query)
            ).order_by('name').select_related('genome', 'contig').prefetch_related('genome__tags')
        else:
            object_list = Site.objects.filter(genome__name=genome).order_by('name').select_related('genome', 'contig').prefetch_related('genome__tags')
        return object_list


class RegulonListView(generic.ListView):
    '''
        Class-based list view of Regulon. Displays a table of regulons.
    '''
    model = Regulon
    template_name = 'regulon_list.html'
    context_object_name = 'regulonlist'
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        genome_name = self.kwargs['genome']
        genome = Genome.objects.get(name = genome_name)
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '" in ' + genome.name + ' regulons'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Regulon.objects.filter(genome__name=genome).filter(
                Q(name__icontains=query) | Q(description__icontains=query) | Q(regulators__locus_tag__icontains=query) | Q(regulators__name__icontains=query)
            ).order_by('name').select_related('genome').prefetch_related('genome__tags')
        else:
            object_list = Regulon.objects.filter(genome__name=genome).order_by('name').select_related('genome').prefetch_related('genome__tags')
        return object_list


class GeneListView(generic.ListView):
    '''
        Class-based list view of Gene. Displays a table of genes.
    '''
    model = Gene
    template_name = 'gene_list.html'
    context_object_name = 'genelist'
    paginate_by = 50

    def get_queryset(self):
        #return Gene.objects.order_by('locus_tag').select_related('genome__strain')
        #return Gene.objects.order_by('locus_tag').values('id', 'locus_tag', 'name', 'function', 'genome__strain__full_name')
        return Gene.objects.none()


class GeneSearchResultsSubView(generic.ListView):
    '''
        Sub-page for AJAX-based view of gene search result
    '''
    context_object_name = 'genelist'
    template_name = 'browser/gene_list_subpage.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        '''
            Generates context string
        '''
        context = super(GeneSearchResultsSubView,self).get_context_data(**kwargs)
        if self.request.GET.get('genome'):
            searchcontext = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'), self.request.GET.get('genome'))
        else:
            searchcontext = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'))
        context['searchcontext'] = searchcontext
        return context

    def get_queryset(self):
        '''
            Generates Gene queryset
        '''
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        query_type = self.request.GET.get('type')
        print(query, genome, query_type)
        if query_type == 'gene':
            if genome:
                object_list = Gene.objects.filter(locus_tag__exact=query).filter(genome__name=genome).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
                if not object_list:
                    object_list = Gene.objects.filter(genome__name=genome).filter(
                        Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                    ).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(locus_tag__exact=query).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
                if not object_list:
                    object_list = Gene.objects.filter(
                        Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                    ).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'og':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'ko_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'kp_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'kr_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'ec_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'tc_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'cazy_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'cog_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif query_type == 'go_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        elif genome:
            if query_type == 'ko':
                ko_ids = Kegg_ortholog.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'kp':
                kp_ids = Kegg_pathway.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon', 'protein').prefetch_related('protein__kegg_orthologs', 'genome__tags')
            elif query_type == 'kr':
                kr_ids = Kegg_reaction.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'ec':
                ec_ids = Ec_number.objects.filter(
                    Q(ec_number__icontains=query) | Q(description__icontains=query)
                ).values('ec_number')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'tc':
                tc_ids = Tc_family.objects.filter(
                    Q(tc_id__icontains=query) | Q(description__icontains=query)
                ).values('tc_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'cazy':
                cazy_ids = Cazy_family.objects.filter(
                    Q(cazy_id__icontains=query) | Q(description__icontains=query)
                ).values('cazy_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'cog':
                cog_ids = Cog_class.objects.filter(
                    Q(cog_id__icontains=query) | Q(description__icontains=query)
                ).values('cog_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            elif query_type == 'go':
                go_ids = Go_term.objects.filter(
                    Q(go_id__icontains=query) | Q(description__icontains=query)
                ).values('go_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
            else:
                object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag').select_related('genome', 'genome__taxon').prefetch_related('genome__tags')
        else:
            object_list = Gene.objects.none()
        return object_list


class GeneSearchResultsAjaxView(View):
    '''
        AJAX-based view for gene search results
    '''
    def get(self,request):
        '''
            Takes a GET request and returns a webpage that will send AJAX request
        '''
        context = {}
        if self.request.GET.get('genome'):
            searchcontext = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'), self.request.GET.get('genome'))
        else:
            searchcontext = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'))
        context['searchcontext'] = searchcontext

        # copy request paramters into context
        for key, val in request.GET.items():
            context[key] = val
        return render(request,'browser/gene_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, creates sub-page (GeneSearchResultsSubView) and sends it back as JSON
            
            For the testing of long-running tasks, uncomment sleep_timer=0 and set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()

        # Sleep timer for testing to imitate long-running task
        # sleep_timer = 0
        # print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        # time.sleep(sleep_timer)

        context = {}
        sub_view = GeneSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()

        if request.GET.get('genome'):
            context['searchcontext'] = generate_gene_search_context(request.GET.get('query'), request.GET.get('type'), request.GET.get('genome'))
            external = generate_external_link(request.GET.get('query'), request.GET.get('type'), request.GET.get('genome'))
        else:
            context['searchcontext'] = generate_gene_search_context(request.GET.get('query'), request.GET.get('type'))
            external = generate_external_link(request.GET.get('query'), request.GET.get('type'))

        if external != '':
            context['external'] = external

        context['searchresult'] = sub_response.content.decode('utf-8')
        context['time'] = time.time()-start_time
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


class GenomeSearchResultsView(generic.ListView):
    '''
        Returns results of search in genome names.
    '''
    # TODO: add search in taxon names
    model = Genome
    context_object_name = 'genomelist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(GenomeSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Genome.objects.filter(name__icontains=query).order_by('name').select_related('strain', 'sample', 'taxon').prefetch_related('tags')
        else:
            object_list = Genome.objects.none()
        return object_list


class StrainSearchResultsView(generic.ListView):
    '''
        Returns results of search in strain names full names or taxonomic orders.
    '''
    model = Strain
    context_object_name = 'strainlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(StrainSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Strain.objects.filter(
                Q(strain_id__icontains=query) | Q(full_name__icontains=query) | Q(order__icontains=query)
            ).order_by('strain_id')
        else:
            object_list = Strain.objects.none()
        return object_list


class KoSearchResultsView(generic.ListView):
    '''
        Returns results of search in KEGG Orthologs.
    '''
    model = Kegg_ortholog
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(KoSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if query:
            object_list = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=query) | Q(description__icontains=query)
            ).order_by('kegg_id')
        else:
            object_list = Kegg_ortholog.objects.order_by('kegg_id')
        return object_list


class KpSearchResultsView(generic.ListView):
    '''
        Returns results of search in KEGG pathways.
    '''
    model = Kegg_pathway
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(KpSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=query) | Q(description__icontains=query)
            ).order_by('kegg_id')
        else:
            object_list = Kegg_pathway.objects.order_by('kegg_id')
        return object_list


class KrSearchResultsView(generic.ListView):
    '''
        Returns results of search in KEGG reactions.
    '''
    model = Kegg_reaction
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(KrSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=query) | Q(description__icontains=query)
            ).order_by('kegg_id')
        else:
            object_list = Kegg_reaction.objects.order_by('kegg_id')
        return object_list


class EcSearchResultsView(generic.ListView):
    '''
        Returns results of search in EC numbers and descriptions.
    '''
    model = Ec_number
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(EcSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Ec_number.objects.filter(
                Q(ec_number__icontains=query) | Q(description__icontains=query)
            ).order_by('ec_number')
        else:
            object_list = Ec_number.objects.order_by('ec_number')
        return object_list


class TcSearchResultsView(generic.ListView):
    '''
        Returns results of search in TCDB numbers and descriptions.
    '''
    model = Tc_family
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(TcSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Tc_family.objects.filter(
                Q(tc_id__icontains=query) | Q(description__icontains=query)
            ).order_by('tc_id')
        else:
            object_list = Tc_family.objects.order_by('tc_id')
        return object_list


class CazySearchResultsView(generic.ListView):
    '''
        Returns results of search in CAZy identifiers and descriptions.
    '''
    model = Cazy_family
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(CazySearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Cazy_family.objects.filter(
                Q(cazy_id__icontains=query) | Q(description__icontains=query)
            ).order_by('cazy_id')
        else:
            object_list = Cazy_family.objects.order_by('cazy_id')
        return object_list


class GoSearchResultsView(generic.ListView):
    '''
        Returns results of search in GO terms.
    '''
    model = Go_term
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(GoSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if query:
            object_list = Go_term.objects.filter(
                Q(go_id__icontains=query) | Q(description__icontains=query)
            ).order_by('go_id')
        else:
            object_list = Go_term.objects.order_by('go_id')
        return object_list


class CogSearchResultsView(generic.ListView):
    '''
        Returns results of search in COG classes.
    '''
    model = Cog_class
    context_object_name = 'annotationlist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(CogSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Cog_class.objects.filter(
                Q(cog_id__icontains=query) | Q(description__icontains=query)
            ).order_by('cog_id')
        else:
            object_list = Cog_class.objects.order_by('cog_id')
        return object_list


def textsearch(request):
    '''
        Displays web-page with text search forms
    '''
    template = loader.get_template('browser/textsearch.html')
    context={}
    return HttpResponse(template.render(context, request))

def startpage(request):
    '''
        Displays home page
    '''
    template = loader.get_template('browser/index.html')  #('browser/landing.html')
    num_genomes = Genome.objects.all().count()
    context = {'num_genomes':num_genomes}
    return HttpResponse(template.render(context, request))

def show_help(request):
    '''
        Displays help page
    '''
    template = loader.get_template('browser/help.html')
    context = {}
    return HttpResponse(template.render(context, request))

'''
def tag_detail(request, name):
'''
#        Displays genome page.
'''
    try:
        genometag = Tag.objects.get(name = name)
    except Tag.DoesNotExist:
        raise Http404('Tag ' + name + ' does not exist')
    context = {'genometag': genometag}
    genomes = Genome.objects.filter(tags__id=genometag.id).order_by('name').select_related('strain', 'taxon').prefetch_related('tags')
    if genomes:
        context['genomes'] = genomes
    return render(request, 'browser/genometag.html', context)
'''

class TagView(generic.ListView):
    '''
        Returns results of search by tag name.
    '''
    model = Tag
    context_object_name = 'genomes'
    paginate_by = 50
    template_name = 'browser/genometag.html'

    def get_context_data(self,**kwargs):
        context = super().get_context_data(**kwargs)
        context['genometag'] = self.genometag
        return context

    def get_queryset(self): # new
        self.genometag = get_object_or_404(Tag, name=self.kwargs['name'])
        object_list = Genome.objects.filter(tags__id=self.genometag.id).order_by('name').select_related('strain', 'taxon').prefetch_related('tags')
        return object_list
    
def taxon_detail(request, taxonomy_id):
    '''
        Displays genome page.
    '''
    try:
        taxon = Taxon.objects.get(taxonomy_id = taxonomy_id)
        print(taxon.name, 'found')
    except Taxon.DoesNotExist:
        print('Taxonomy ID', taxonomy_id, 'not found')
        raise Http404('Taxon ' + taxonomy_id + ' does not exist')
    context = {'taxon': taxon}
    children = get_taxon_children(taxonomy_id)
    context['genomes'] = Genome.objects.filter(taxon__taxonomy_id__in=children)
    context['strains'] = Strain.objects.filter(taxon__taxonomy_id__in=children)
    context['sunburst'] = generate_sunburst(taxonomy_id)
    lineage = [taxon,]
    parent_id = taxon.parent_id
    iteration_count = 0
    while True:
        parent_taxon = Taxon.objects.get(taxonomy_id = parent_id)
        iteration_count += 1
        if parent_taxon.taxonomy_id == parent_taxon.parent_id or parent_id == '1':
            break
        lineage.append(parent_taxon)
        parent_id = parent_taxon.parent_id
    context['lineage'] = reversed(lineage)
    return render(request, 'browser/taxon.html', context)
    
def genome_detail(request, name):
    '''
        Displays genome page.
    '''
    try:
        genome = Genome.objects.get(name = name)
    except Genome.DoesNotExist:
        raise Http404('Genome ' + name + ' does not exist')
    context = {'genome': genome}
    context['operons'] = Operon.objects.filter(genome=genome).count()
    context['sites'] = Site.objects.filter(genome=genome).count()
    context['regulons'] = Regulon.objects.filter(genome=genome).count()
    if request.GET.get('contig'):
        context['highlight_start'] = request.GET.get('start')
        context['highlight_end'] = request.GET.get('end')
        display_offset = 5000
        contig = Contig.objects.get(contig_id = request.GET.get('contig'), genome=genome)
        context['contig'] = contig.contig_id
        display_start = int(request.GET.get('start')) - display_offset
        if display_start < 0:
            display_start = 0
        context['viewer_start'] = str(display_start)
        display_end = int(request.GET.get('end')) + display_offset
        if display_end > contig.size:
            display_end = contig.size
        context['viewer_end'] = str(display_end)
    return render(request, 'browser/genome.html', context)

def strain_detail(request, strain_id):
    '''
        Displays strain page.
    '''
    try:
        strain = Strain.objects.get(strain_id = strain_id)
        genomes = Genome.objects.filter(strain=strain).order_by('name').prefetch_related('tags')
        strain_metadata = Strain_metadata.objects.filter(strain=strain)
        metadata_entries = {}
        for metadata_entry in strain_metadata:
            if metadata_entry.source not in metadata_entries:
                metadata_entries[metadata_entry.source] = {}
                metadata_entries[metadata_entry.source]['entries'] = []
            metadata_entries[metadata_entry.source]['source'] = metadata_entry.source
            metadata_entries[metadata_entry.source]['url'] = metadata_entry.url
            metadata_entries[metadata_entry.source]['entries'].append({'key':metadata_entry.key,'value':metadata_entry.value})
        metadata = []
        for source in sorted(metadata_entries.keys()):
            metadata.append(metadata_entries[source])
        print(metadata)
    except Strain.DoesNotExist:
        return render(request, '404.html', {'searchcontext': 'Sample not found: ' + strain_id})
    return render(request, 'browser/strain.html', {'strain': strain, 'genomes':genomes, 'metadata':metadata})

def sample_detail(request, sample_id):
    '''
        Displays sample page.
    '''
    try:
        sample = Sample.objects.get(sample_id = sample_id)
        genomes = Genome.objects.filter(sample=sample).order_by('name').prefetch_related('tags')
        sample_metadata = Sample_metadata.objects.filter(sample=sample)
        metadata_entries = {}
        for metadata_entry in sample_metadata:
            if metadata_entry.source not in metadata_entries:
                metadata_entries[metadata_entry.source] = {}
                metadata_entries[metadata_entry.source]['entries'] = []
            metadata_entries[metadata_entry.source]['source'] = metadata_entry.source
            metadata_entries[metadata_entry.source]['url'] = metadata_entry.url
            metadata_entries[metadata_entry.source]['entries'].append({'key':metadata_entry.key,'value':metadata_entry.value})
        metadata = []
        for source in sorted(metadata_entries.keys()):
            metadata.append(metadata_entries[source])
        print(metadata)
    except Sample.DoesNotExist:
        return render(request, '404.html', {'searchcontext': 'Sample not found: ' + sample_id})
    return render(request, 'browser/sample.html', {'sample': sample, 'genomes':genomes, 'metadata':metadata})

def gene_detail(request, genome, locus_tag):
    '''
        Displays gene page.
    '''
    try:
        gene = Gene.objects.select_related(
            'genome', 'genome__taxon', 'protein', 'operon'
            ).prefetch_related(
            'protein__ortholog_groups__taxon', 'protein__kegg_orthologs', 'protein__kegg_pathways', 'protein__kegg_reactions', 'protein__ec_numbers', 'protein__go_terms', 'protein__tc_families', 'protein__cog_classes', 'genome__tags'
            ).get(genome__name=genome, locus_tag = locus_tag)
        annotations = Annotation.objects.filter(gene_id = gene).order_by('source', 'value')
        context = {'gene': gene}
        if annotations:
            context['annotations'] = annotations
        if gene.protein:
            context['protein'] = gene.protein
        viewer_start = gene.start - 5000
        if viewer_start < 0:
            viewer_start = 1
        viewer_end = gene.end + 5000
        if viewer_end > gene.contig.size:
            viewer_end = gene.contig.size
        context['viewer_start'] = str(viewer_start)
        context['viewer_end'] = str(viewer_end)
    except Gene.DoesNotExist:
        raise Http404('Gene not found')
    return render(request, 'browser/gene.html', context)

def operon_detail(request, name, genome=None):
    '''
        Displays operon page.
    '''
    try:
        operon = Operon.objects.select_related('contig', 'genome').prefetch_related('genome__tags').get(name=name)
        genes = Gene.objects.filter(operon = operon.id).select_related('contig', 'genome', 'genome__strain')
        contig = operon.contig
        context = {'operon': operon}
        context['genes'] = genes
        viewer_start = operon.start - 5000
        if viewer_start < 0:
            viewer_start = 1
        viewer_end = operon.end + 5000
        if viewer_end > contig.size:
            viewer_end = contig.size
        context['viewer_start'] = str(viewer_start)
        context['viewer_end'] = str(viewer_end)
        context['highlight_start'] = operon.start
        context['highlight_end'] = operon.end
        context['sites'] = operon.site_set.all()
    except Operon.DoesNotExist:
        raise Http404('Operon not found')
    return render(request, 'browser/operon.html', context)

def site_detail(request, genome, name):
    '''
        Displays site page.
    '''
    try:
        site = Site.objects.select_related(
            'genome', 'genome__strain', 'contig'
            ).prefetch_related('genome__tags').get(genome__name=genome, name = name)
        context = {'site': site}
        viewer_start = site.start - 5000
        if viewer_start < 0:
            viewer_start = 1
        viewer_end = site.end + 5000
        if viewer_end > site.contig.size:
            viewer_end = site.contig.size
        context['viewer_start'] = str(viewer_start)
        context['viewer_end'] = str(viewer_end)
        context['highlight_start'] = site.start
        context['highlight_end'] = site.end
        print(site.genes)
    except Site.DoesNotExist:
        raise Http404('Site not found')
    return render(request, 'browser/site.html', context)

def regulon_detail(request, genome, name):
    '''
        Displays regulon page.
    '''
    try:
        regulon = Regulon.objects.select_related(
            'genome', 'genome__strain'
            ).prefetch_related('regulators', 'genome__tags'
            ).get(genome__name=genome, name = name)
        context = {'regulon': regulon}
    except Regulon.DoesNotExist:
        raise Http404('Regulon not found')
    sites = Site.objects.filter(regulon = regulon).select_related(
        'contig').prefetch_related('genes', 'operons')
    context['sites'] = sites
    ortholog_groups = set()
    for regulator in regulon.regulators.all():
        for og in regulator.protein.ortholog_groups.all():
            ortholog_groups.add(og.id)
    context['ortholog_groups'] = Ortholog_group.objects.filter(id__in=list(ortholog_groups))
    print([x.name for x in sites])
    return render(request, 'browser/regulon.html', context)

def gene_byname(request):
    '''
        Displays gene page.
    '''
    if request.GET.get('genome') and request.GET.get('locus'):
        genome_name = request.GET.get('genome')
        locus_tag = request.GET.get('locus')
    else:
        raise Http404('Genome name or locus ID not provided')
    try:
        gene = Gene.objects.get(locus_tag = locus_tag, genome__name = genome_name).prefetch_related('genome__tags')
        annotations = Annotation.objects.filter(gene_id = gene).order_by('source', 'value')
        context = {'gene': gene}
        if annotations:
            context['annotations'] = annotations
        if gene.protein:
            context['protein'] = gene.protein
        reverse_strand = False
        if gene.strand == -1:
            reverse_strand = True
        viewer_start = gene.start - 5000
        if viewer_start < 0:
            viewer_start = 1
        viewer_end = gene.end + 5000
        if viewer_end > gene.contig.size:
            viewer_end = gene.contig.size
        context['reverse_strand'] = reverse_strand
        context['viewer_start'] = str(viewer_start)
        context['viewer_end'] = str(viewer_end)
    except Gene.DoesNotExist:
        raise Http404('Gene not found for genome ' + genome_name + ' and locus tag ' + locus_tag)
    except Gene.MultipleObjectsReturned:
        raise Http404('More than one gene was found')
    return render(request, 'browser/gene.html', context)

def protein_search_external(request):
    '''
        Returns protein search results for external query (from Fitness browser etc.).
    '''
    context = {}
    if request.GET.get("sequence"):
        result = {}
        sequence = request.GET.get("sequence")
        hits, searchcontext, query_len = run_protein_search(sequence)
        if searchcontext != '':
            context['searchcontext'] = searchcontext
        for row in hits:
            row=row.split('\t')
            print('Search for gene', row[1])
            if row[0] not in result:
                result[row[0]] = []
            unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
            query_cov = (query_len - unaligned_part) * 100.0 / query_len
            genes = Gene.objects.select_related('protein', 'genome', 'genome__taxon').prefetch_related('genome__tags').filter(protein__protein_hash = row[1])
            for target in genes:
                hit = [target.genome.name, target.locus_tag, target.genome.taxon.name,
                       target.function, '{:.1f}'.format(float(row[2])), row[3], '{:.1f}'.format(query_cov), row[10], row[11]]
                result[row[0]].append(hit)
        context['searchresult'] = result
    else:
        raise Http404('Provide protein sequence in FASTA format')
    return render(request, 'browser/proteinsearch.html', context)


class NsearchResultView(View):
    '''
        AJAX-based view for nucleotide similarity search
    '''
    def post(self,request):
        '''
            Takes a POST request and returns a webpage that will send AJAX request
        '''
        sequence = request.POST.get("sequence")
        #print('Sequence sent:', sequence)
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        for key, val in request.POST.items():
            context[key] = val
            #print(key, val)
        return render(request,'browser/nucleotidesearchajax.html', context)

    def get(self,request):
        '''
            Displays nucleotide sequence search form
        '''
        return render(request,'browser/nucleotidesearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, calls run_nucleotide_search and sends formatted results back as JSON
            
            For the testing of long-running tasks, uncomment sleep_timer=0 and set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()
        result = []
        #print('REQUEST2')
        #for key, val in request.POST.items():
        #    print(key, val)
        params = {'sequence': request.POST.get("sequence"), 'evalue': request.POST.get("evalue"), 'hitstoshow': request.POST.get("hitstoshow")}
        # Test sequence
        #params['sequence'] = '>test\nATGACGATGCACCCGGCAGCACCTCGAACTCCGCACTGGCGCTGCTTGCACCGGCTGGCATGGAGCCTGTCCTGCGCTGCCCTGTTGCCGCTCGCTGCACTGGCGCAGGACGTACCGTCCCGCGCCGTCACGCCGGTGTCCGCAGCGTCGCCGCCGCGCCAGTCGCAGGACGACGCCTGGTGGACCGGGCCGATGCTGGCGAACTCCGCCGCCACCCTGCCGCGCGGCCACGTCCTGATCGAGCCTTACGTCTACGACGTGTCCTCGCCGCACGCCGACGGCTACGGTTCGCTCACCTACATGCTCTACGGCCTCACCGACCGGCTGACGGTCGGCCTGATGCCGGTGCTGGGCTACAACCGCATGGATGGCCCGGGCGACAGCAGCGGGATCGGGCTGGGCGACGTCAGCGTGCAGGCGCAGTACCGGCTGACCGATGTGCCGGCGGGCAGTTCGCGGCCCACGGTCTCGCTGCAACTGCAGGAAACCCTGCCGACCGGCAAGTACGACCGGCTGGGCCGGCGACCCGGCAACGGCCTGGGCAGCGGCGCCACCACCACTACGCTGCAGGTCAACACGCAGACGTATTTCTGGTTGTCCAACGGCCGCATCCTGCGCATGCGCTTCAACGTGGCGCAATCATTCTCGACGCGGGCACGGGTCGAGGACATCAGCGTCTACGGCACCCCGGACGGCTTTCGCGGGCACGCCCGGCCGGGGCGTTCGTTCTTCGTCAATGCGGCCTGGGAGTACAGCCTCAGCCAGCGCTGGGTGCTGGCGCTCGACCTCACCTACCGGCGCAGCCACGGTGCCCGCGTGCGCGACGACGACCTCAATGGCGTGCCTGCCTTGCGTCTGGACGGCCGCTCCAGCGAGGCGTTCGGCTTTGCCCCGGCCATCGAGTACAGCTGGAGTCCGCGGCTCGGCGTGCTGTTCGGCACCCGCGTGATCACCGGCGGGCACAACACCGCGACCACGATCACGCCGGCGGTGGCCTTCAACTACGTGCACTGA'
        #print('Sequence received:', params['sequence'])
        #sleep_timer = 0
        #print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)
        hits, searchcontext, query_len = run_nucleotide_search(params)
        if hits:
            result.append('<table><thead><tr><th>Target contig (click to see hit)</th><th>Genome</th><th>%identity</th><th>Alignment length</th><th>%Query coverage</th><th>E-value</th><th>Bit-score</th></tr></thead><tbody>')
            for row in hits:
                row=row.split('\t')
                contig_name, genome_name = row[1].split('|')
                genome_name = genome_name.split('[')[0]
                print('Search for contig', contig_name, 'in genome', genome_name)
                target = Contig.objects.select_related('genome', 'genome__taxon').prefetch_related('genome__tags').get(contig_id = contig_name, genome__name = genome_name)
                ani = float(row[2])
                strand = '+'
                start = row[8]
                end = row[9]
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                if int(row[7]) < int(row[6]):
                    strand = '-'
                genome_tags = ''
                for genome_tag in target.genome.tags.all():
                    genome_tags += '<span class="genometag" style="background-color:' + genome_tag.color + '"><a href="' + reverse('tagdetails', args=(genome_tag.name,)) + '" style="color:' + genome_tag.textcolor + '" title="' + genome_tag.description + '">' + genome_tag.name + '</a></span>&nbsp;'
                hit = '<tr><td align=\"left\"><a href=\"' + reverse('genomedetails', args=(target.genome.name,)) + '?contig=' + target.contig_id + '&start=' + start + '&end=' + end + '\">' + contig_name + ': (' + strand + 'strand) ' + start + '..' + end + '</a></td><td align="left">' + target.genome.name + ' [' + target.genome.taxon.name + ']' + genome_tags + '</td><td>' + '{:.2f}'.format(float(row[2])) + '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) + '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                result.append(hit)
            result.append('</tbody></table>')
            context = {"searchresult":'\n'.join(result),"searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        elif searchcontext == '':
            context = {"searchresult":'',"searchcontext":'No hits found',"query_len":query_len,"time":time.time()-start_time}
        else:
            context = {"searchresult":"","searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        print(context)
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


def nucleotidesearchform(request):
    '''
        Displays nucleotide sequence search form
    '''
    return render(request,'browser/nucleotidesearchform.html')


class PsearchResultView(View):
    '''
        AJAX-based view for protein similarity search
    '''
    def post(self,request):
        '''
            Takes a POST request and returns a webpage that will send AJAX request
        '''
        #sequence = request.POST.get("sequence")
        #print('Sequence sent:', sequence)
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        #print('REQUEST1')
        for key, val in request.POST.items():
            context[key] = val
            #print(key, val)
        return render(request,'browser/proteinsearchajax.html', context)

    def get(self,request):
        '''
            Displays protein sequence search form
        '''
        return render(request,'browser/proteinsearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, calls run_protein_search and sends formatted results back as JSON
            
            For the testing of long-running tasks, uncomment sleep_timer=0 and set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()
        result = []
        #print('REQUEST2')
        #for key, val in request.POST.items():
        #    print(key, val)
        params = {'sequence': request.POST.get("sequence"), 'evalue': request.POST.get("evalue"), 'hitstoshow': request.POST.get("hitstoshow")}
        # Test sequence
        #params['sequence'] = 'MTKQVQEAYIVAATRTPVGKAPRGVFRNTRPDDMLAHVIRAVMAQAPGIDPHQIGDVIIGCAMPEAEQGMNVARIGLLLAGLPDTVPGVTVNRFCSSGLQSVAMAADRIRLGLDDLMLAGGTESMSMVPMMGHKIAMNPAIFNDENIGIAYGMGITAENVAKQWKVSREQQDAFSVESHRRALAAQAAGEFNDEISPFALDDHYPNLATRGIVTDSRRIDSDEGPRAGTTMEVLAKLKTVFRNGQFGGTVTAGNSSQMSDGAGAVLLASERAVKEYNLQPLARFVGFSVAGVPPEVMGIGPKEAIPKALKQAGLNRDQLDWIELNEAFAAQALAVMGDLGLDPDKVNPLGGAIALGHPLGATGAVRIATLVHGMRRRKQKYGMVTMCIGTGMGAAGIFEAL'
        #print('Sequence received:', sequence)
        #sleep_timer = 0
        #print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)
        hits, searchcontext, query_len = run_protein_search(params)

        if hits:
            result.append('<table><thead><tr><th>Target gene</th><th>Genome</th><th>%identity</th><th>Alignment length</th><th>%Query coverage</th><th>E-value</th><th>Bit-score</th></tr></thead><tbody>')
            for row in hits:
                row=row.split('\t')
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                genes = Gene.objects.select_related('protein', 'genome', 'genome__taxon').prefetch_related('genome__tags').filter(protein__protein_hash = row[1])
                for target in genes:
                    genome_tags = ''
                    for genome_tag in target.genome.tags.all():
                        genome_tags += '<span class="genometag" style="background-color:' + genome_tag.color + '"><a href="' + reverse('tagdetails', args=(genome_tag.name,)) + '" style="color:' + genome_tag.textcolor + '" title="' + genome_tag.description + '">' + genome_tag.name + '</a></span>&nbsp;'
                    hit = '<tr><td align=\"left\"><a href=\"' + reverse('genedetails', args=(target.genome.name, target.locus_tag)) + '\">' + target.locus_tag + '</a></td><td align="left">' + target.genome.name + ' [' + target.genome.taxon.name + ']' + genome_tags + '</td><td>' + '{:.1f}'.format(float(row[2])) + '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) + '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                    result.append(hit)
            result.append('</tbody></table>')
            context = {"searchresult":'\n'.join(result),"searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        elif searchcontext == '':
            context = {"searchresult":'',"searchcontext":'No hits found',"query_len":query_len,"time":time.time()-start_time}
        else:
            context = {"searchresult":"","searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


def proteinsearchform(request):
    '''
        Displays protein sequence search form
    '''
    return render(request,'browser/proteinsearchform.html')

def cregulon_view(request):
    """
    Conserved regulon view
    """
    og_id = request.GET.get('og')
    #locus_tag = request.GET.get('locus_tag')
    context = build_conserved_regulon(og_id)
    return render(request, 'browser/cregulon.html', context)

    
class ComparativeView(View):
    '''
        AJAX-based view for comparative gene plot page
    '''
    @staticmethod
    def verify_parameters(size, lines):
        '''
            Parameters validation
        '''
        try:
            if size is None:
                raise SuspiciousOperation('"size" parameter is missing')

            try:
                size = int(size)
            except TypeError:
                raise SuspiciousOperation("Unacceptable value '%s' for size parameter." % size)
                
            if int(size) not in [5, 10, 20, 40, 60, 80, 100]:
                raise SuspiciousOperation("Unacceptable value for locus size: '%s'" % size)
                
            if lines is None:
                raise SuspiciousOperation('"lines" parameter is missing')
                
            try:
                lines = int(lines)
            except TypeError:
                raise SuspiciousOperation("Unacceptable value '%s' for lines number." % lines)
                
            if int(lines) not in [10, 25, 50, 75, 100, 200]:
                raise SuspiciousOperation("Unacceptable value for lines number: '%s'" % lines)
                
        except Exception as e:
            raise SuspiciousOperation(str(e))
        
    def get(self,request):
        '''
            Displays a page that will send AJAX request
        '''
        locus_tag = request.GET.get('locus_tag')
        genome = request.GET.get('genome')
        og_id = request.GET.get('og')

        try:
            ComparativeView.verify_parameters(request.GET.get('size'), request.GET.get('lines'))
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
        
        #print('REQUEST1')
        #print('Request parameters:', og_id, genome, locus_tag, request.GET.get('size'), request.GET.get('lines'))
        context = {}
        for key, val in request.GET.items():
            context[key] = val
        try:
            gene = Gene.objects.select_related('protein', 'genome__strain', 'genome__sample', 'contig').get(locus_tag=locus_tag, genome__name=genome)
            og = Ortholog_group.objects.get(id=og_id)
        except Gene.DoesNotExist:
            raise Http404('Gene not found')
        except Ortholog_group.DoesNotExist:
            raise Http404('Ortholog group not found')
        context['gene'] = gene
        context['ortholog_group'] = og
        return render(request,'browser/scriblajax.html', context)

    @staticmethod
    def ajax_view(request):
        start_time = time.time()

        #print('REQUEST2')
        #for key, val in request.GET.items():
        #    print(key, val)
        try:
            ComparativeView.verify_parameters(request.GET.get('size'), request.GET.get('lines'))
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
            
        # Sleep timer for testing to imitate long-running task
        #sleep_timer = 0
        #print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)

        context = {}
        locus_tag = request.GET.get('locus_tag')
        genome = request.GET.get('genome')
        og_id = request.GET.get('og')

        try:
            gene = Gene.objects.select_related('protein', 'genome__strain', 'genome__sample', 'contig').get(locus_tag=locus_tag, genome__name=genome)
            og = Ortholog_group.objects.get(id=og_id)
        except Gene.DoesNotExist:
            raise Http404('Gene not found')
        except Ortholog_group.DoesNotExist:
            raise Http404('Ortholog group not found')
        scribl, tree_canvas, tree_newick, og_gene_count, plot_gene_count = get_scribl(gene, og, request)
        
        if og_gene_count == 1:
            scribl='<script type="text/javascript">\nalert("Comparative plot cannot be created for only one genome");</script>'
            context = {'scribl':scribl, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count, "time":time.time()-start_time}
        else:
            scribl='<script type="text/javascript">\nfunction draw(canvasName) {\nvar canvas = document.getElementById(canvasName);\nvar parent = document.getElementById("scribl-container");\nvar ctx = canvas.getContext("2d");\n' + scribl + '\nchart.draw();\nvar img = canvas.toDataURL("image/png");\ndocument.getElementById("pngexport").href = img;\n}\n</script>'
            plot_gene_count -= 1
            context = {'scribl':scribl, 'tree_canvas':tree_canvas, 'tree_newick':tree_newick, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count,"time":time.time()-start_time}
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")

    
def export_csv(request):
    '''
        Returns list of genes in CSV format
    '''
    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="export.tab"'

    writer = csv.writer(response, delimiter='\t')
    search_context = ('','')
    annotation_query = request.GET.get('annotation_query')
    query_type = request.GET.get('type')
    query = request.GET.get('query')
    genome = request.GET.get('genome')

    if annotation_query:
        search_context = ('Gene annotation query', annotation_query)
        if genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
    elif query_type == 'gene':
        if query:
            search_context = ('Query', query)
            if genome:
                object_list = Gene.objects.filter(genome__name=genome).filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query)
                ).order_by('locus_tag')
        elif genome:
            search_context = ('Genome', genome)
            object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
        else:
            object_list = Gene.objects.none()
    else:
        if query_type == 'og':
            search_context = ('Ortholog group', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=query).values('protein_hash')]
        elif query_type == 'ko_id':
            search_context = ('KEGG ortholog', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=query).values('protein_hash')]
        elif query_type == 'kp_id':
            search_context = ('KEGG pathway', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=query).values('protein_hash')]
        elif query_type == 'kr_id':
            search_context = ('KEGG reaction', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=query).values('protein_hash')]
        elif query_type == 'ec_id':
            search_context = ('EC number', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=query).values('protein_hash')]
        elif query_type == 'tc_id':
            search_context = ('TCDB family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=query).values('protein_hash')]
        elif query_type == 'cazy_id':
            search_context = ('CAZy family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=query).values('protein_hash')]
        elif query_type == 'cog_id':
            search_context = ('COG class', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=query).values('protein_hash')]
        elif query_type == 'go_id':
            search_context = ('GO term', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=query).values('protein_hash')]
        elif query_type == 'ko':
            search_context = ('KEGG ortholog query', query)
            ko_ids = Kegg_ortholog.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
        elif query_type == 'kp':
            search_context = ('KEGG pathway query', query)
            kp_ids = Kegg_pathway.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
        elif query_type == 'kr':
            search_context = ('KEGG reaction query', query)
            kr_ids = Kegg_reaction.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
        elif query_type == 'ec':
            search_context = ('EC number query', query)
            ec_ids = Ec_number.objects.filter(Q(ec_number__icontains=query) | Q(description__icontains=query)).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
        elif query_type == 'tc':
            search_context = ('TCDB family query', query)
            tc_ids = Tc_family.objects.filter(Q(tc_id__icontains=query) | Q(description__icontains=query)).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
        elif query_type == 'cazy':
            search_context = ('CAZy family query', query)
            cazy_ids = Cazy_family.objects.filter(Q(cazy_id__icontains=query) | Q(description__icontains=query)).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
        elif query_type == 'cog':
            search_context = ('COG class query', query)
            cog_ids = Cog_class.objects.filter(Q(cog_id__icontains=query) | Q(description__icontains=query)).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
        elif query_type == 'go':
            search_context = ('GO term query', query)
            go_ids = Go_term.objects.filter(Q(go_id__icontains=query) | Q(description__icontains=query)).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
        else:
            search_context = ('Unknown', query)
            proteins = []

        if genome:
            if not query_type:
                object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            if not query_type:
                object_list = Gene.objects.none()
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')

    if annotation_query:
        writer.writerow(['Locus tag', 'Name', 'Organism', 'Genome', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Function', 'Annotation_source', 'Annotation_type', 'Annotation_note'])
        for obj in object_list:
            writer.writerow([obj.gene_id.locus_tag, obj.gene_id.name, obj.gene_id.genome.taxon.name, obj.gene_id.genome.name, obj.gene_id.contig.contig_id, str(obj.gene_id.start), str(obj.gene_id.end), str(obj.gene_id.strand), obj.gene_id.type, obj.gene_id.function, obj.source, obj.key, obj.value, obj.note])
    elif search_context[0] == 'Genome':
        writer.writerow(['Locus tag', 'Name', 'Organism', 'Genome', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Function'])
        for gene in object_list:
            writer.writerow([gene.locus_tag, gene.name, gene.genome.taxon.name, gene.genome.name, gene.contig.contig_id, str(gene.start), str(gene.end), str(gene.strand), gene.type, gene.function])
    else:
        writer.writerow(['Locus tag', 'Name', 'Organism', 'Genome', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Function', search_context[0]])
        for gene in object_list:
            writer.writerow([gene.locus_tag, gene.name, gene.genome.taxon.name, gene.genome.name, gene.contig.contig_id, str(gene.start), str(gene.end), str(gene.strand), gene.type, gene.function, search_context[1]])

    return response

def export_fasta(request):
    '''
        Returns protein sequences in FASTA format
    '''
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="export.faa"'

    search_context = ('','')
    query_type = request.GET.get('type')
    query = request.GET.get('query')
    genome = request.GET.get('genome')
    annotation_query = request.GET.get('annotation_query')

    if annotation_query:
        search_context = ('Gene annotation query', annotation_query)
        if genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon', 'gene_id__protein')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon', 'gene_id__protein')
    elif query_type == 'gene':
        if query:
            search_context = ('Query', query)
            if genome:
                object_list = Gene.objects.filter(genome__name=genome).filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query)
                ).order_by('locus_tag')
        elif genome:
            object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
        else:
            object_list = Gene.objects.none()
    else:
        if query_type == 'og':
            search_context = ('Ortholog group', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=query).values('protein_hash')]
        elif query_type == 'ko_id':
            search_context = ('KEGG ortholog', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=query).values('protein_hash')]
        elif query_type == 'kp_id':
            search_context = ('KEGG pathway', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=query).values('protein_hash')]
        elif query_type == 'kr_id':
            search_context = ('KEGG reaction', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=query).values('protein_hash')]
        elif query_type == 'ec_id':
            search_context = ('EC number', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=query).values('protein_hash')]
        elif query_type == 'tc_id':
            search_context = ('TCDB family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=query).values('protein_hash')]
        elif query_type == 'cazy_id':
            search_context = ('CAZy family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=query).values('protein_hash')]
        elif query_type == 'cog_id':
            search_context = ('COG class', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=query).values('protein_hash')]
        elif query_type == 'go_id':
            search_context = ('GO term', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=query).values('protein_hash')]
        elif query_type == 'ko':
            search_context = ('KEGG ortholog query', query)
            ko_ids = Kegg_ortholog.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
        elif query_type == 'kp':
            search_context = ('KEGG pathway query', query)
            kp_ids = Kegg_pathway.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
        elif query_type == 'kr':
            search_context = ('KEGG reaction query', query)
            kr_ids = Kegg_reaction.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
        elif query_type == 'ec':
            search_context = ('EC number query', query)
            ec_ids = Ec_number.objects.filter(Q(ec_number__icontains=query) | Q(description__icontains=query)).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
        elif query_type == 'tc':
            search_context = ('TCDB family query', query)
            tc_ids = Tc_family.objects.filter(Q(tc_id__icontains=query) | Q(description__icontains=query)).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
        elif query_type == 'cazy':
            search_context = ('CAZy family query', query)
            cazy_ids = Cazy_family.objects.filter(Q(cazy_id__icontains=query) | Q(description__icontains=query)).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
        elif query_type == 'cog':
            search_context = ('COG class query', query)
            cog_ids = Cog_class.objects.filter(Q(cog_id__icontains=query) | Q(description__icontains=query)).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
        elif query_type == 'go':
            search_context = ('GO term query', query)
            go_ids = Go_term.objects.filter(Q(go_id__icontains=query) | Q(description__icontains=query)).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
        else:
            proteins = []

        if genome:
            if not query_type:
                object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            if not query_type:
                object_list = Gene.objects.none()
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')

    if annotation_query:
        for obj in object_list:
            if obj.gene_id.protein:
                response.write('>' + obj.gene_id.locus_tag + '|' + obj.gene_id.genome.name + ' [' + obj.gene_id.genome.taxon.name + ']\n' + obj.gene_id.protein.sequence + '\n')
    else:
        for gene in object_list:
            if gene.protein:
                response.write('>' + gene.locus_tag + '|' + gene.genome.name + ' [' + gene.genome.taxon.name + ']\n' + gene.protein.sequence + '\n')

    return response
    
def handler404(request, exception):
    '''
        Returns 404 page
    '''
    return render(request, '404.html', status=404)
    
def generate_gene_search_context(query, query_type, genome=None):
    '''
        Generates search context string and external link for various query types
    '''
    searchcontext = ''
    if query is None:
        query = ''
    if genome is None:
        if query_type=='gene':
            searchcontext = 'Search results for "' + query + '"'
        elif query_type=='og':
            eggnog_og = Ortholog_group.objects.get(id=query)
            searchcontext = 'Genes from Ortholog Group ' + eggnog_og.eggnog_id + '[' + eggnog_og.taxon.name + ']'
        elif query_type=='ko_id':
            searchcontext = 'Genes from KEGG Ortholog Group ' + query
        elif query_type=='kp_id':
            searchcontext = 'Genes from KEGG pathway ' + query
        elif query_type=='kr_id':
            searchcontext = 'Genes from KEGG reaction ' + query
        elif query_type=='ec_id':
            searchcontext = 'Genes with EC number ' + query
        elif query_type=='tc_id':
            searchcontext = 'Genes from TCDB family ' + query
        elif query_type=='cazy_id':
            searchcontext = 'Genes from CAZy family ' + query
        elif query_type=='cog_id':
            searchcontext = 'Genes from COG class ' + query
        elif query_type=='go_id':
            searchcontext = 'Genes assigned to GO term ' + query
        elif query_type=='ko':
            searchcontext = 'Genes assigned to KEGG Ortholog groups containing "' + query + '"'
        elif query_type=='kp':
            searchcontext = 'Genes assigned to KEGG pathways containing "' + query + '"'
        elif query_type=='kr':
            searchcontext = 'Genes assigned to KEGG reactions containing "' + query + '"'
        elif query_type=='ec':
            searchcontext = 'Genes assigned to EC numbers containing "' + query + '"'
        elif query_type=='tc':
            searchcontext = 'Genes assigned to TCDB families containing "' + query + '"'
        elif query_type=='cazy':
            searchcontext = 'Genes assigned to CAZy families containing "' + query + '"'
        elif query_type=='cog':
            searchcontext = 'Genes assigned to COG class containing "' + query + '"'
        elif query_type=='go':
            searchcontext = 'Genes assigned to GO terms containing "' + query + '"'
    else:
        if query_type=='ko':
            searchcontext = 'Genes from genome ' + genome + ' assigned to KEGG Ortholog groups containing "' + query + '"'
        elif query_type=='kp':
            searchcontext = 'Genes from genome ' + genome + ' assigned to KEGG pathways containing "' + query + '"'
        elif query_type=='kr':
            searchcontext = 'Genes from genome ' + genome + ' assigned to KEGG reactions containing "' + query + '"'
        elif query_type=='ec':
            searchcontext = 'Genes from genome ' + genome + ' assigned to EC numbers containing "' + query + '"'
        elif query_type=='tc':
            searchcontext = 'Genes from genome ' + genome + ' assigned to TCDB families containing "' + query + '"'
        elif query_type=='cazy':
            searchcontext = 'Genes from genome ' + genome + ' assigned to CAZy families containing "' + query + '"'
        elif query_type=='cog':
            searchcontext = 'Genes from genome ' + genome + ' assigned to COG class containing "' + query + '"'
        elif query_type=='go':
            searchcontext = 'Genes from genome ' + genome + ' assigned to GO terms containing "' + query + '"'
        elif query_type=='gene':
            searchcontext = 'Genes from genome ' + genome
    return searchcontext

def generate_external_link(query, query_type, genome=None):
    '''
        Generates search context string and external link for various query types
    '''
    external = ''
    if query is None:
        query = ''
    if genome is None:
        if query_type=='ko_id':
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='kp_id':
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='kr_id':
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='ec_id':
            external = 'https://www.kegg.jp/dbget-bin/www_bget?ec:' + query
        elif query_type=='tc_id':
            external = 'http://www.tcdb.org/search/result.php?tc=' + query + '#' + query
        elif query_type=='cazy_id':
            external = 'http://www.cazy.org/' + query + '.html'
        elif query_type=='cog_id':
            external = 'https://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/list' + query + '.html'
        elif query_type=='go_id':
            external = 'https://www.ebi.ac.uk/QuickGO/search/' + query
    else:
        if query_type=='kp':
            kp_ids = Kegg_pathway.objects.filter(Q(kegg_id__icontains=query) | Q(description__icontains=query)).values('kegg_id')
            ext_kegg_map_url = ''
            if kp_ids.count() == 1:
                gene_list = Gene.objects.filter(genome__name=genome, protein__kegg_pathways__kegg_id__in=kp_ids).select_related('protein').prefetch_related('protein__kegg_orthologs')
                print(kp_ids[0]['kegg_id'])
                kegg_map_url = 'https://www.kegg.jp/pathway/' + kp_ids[0]['kegg_id'] + '+'
                print(kegg_map_url)
                ko_ids = set()
                for gene in gene_list:
                    for ko in gene.protein.kegg_orthologs.all():
                        ko_ids.add(ko.kegg_id)
                if ko_ids:
                    ext_kegg_map_url = kegg_map_url + '+'.join(sorted(list(ko_ids)))
            if ext_kegg_map_url != '':
                external = ext_kegg_map_url
    if external != '':
        external = '<a href="' + external + '">External link</a>'
    return external

