import csv
import time
import json
import logging
from django.shortcuts import render
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.template import loader
from django.views import View, generic
from django.db.models import Q
from django.core.exceptions import SuspiciousOperation
from django.urls import reverse
from genomebrowser.settings import TITLE
from browser.models import Annotation
from browser.models import Cazy_family
from browser.models import Contig
from browser.models import Cog_class
from browser.models import Ec_number
from browser.models import Gene
from browser.models import Genome
from browser.models import Go_term
from browser.models import Kegg_ortholog
from browser.models import Kegg_pathway
from browser.models import Kegg_reaction
from browser.models import Operon
from browser.models import Ortholog_group
from browser.models import Protein
from browser.models import Regulon
from browser.models import Tag
from browser.models import Taxon
from browser.models import Tc_family
from browser.models import Sample
from browser.models import Sample_metadata
from browser.models import Site
from browser.models import Strain
from browser.models import Strain_metadata
from browser.colors import COLORS
from browser.treemap import generate_og_treemap, generate_genes_treemap
from browser.seqsearch import run_protein_search, run_nucleotide_search
from browser.comparative_analysis import get_scribl
from browser.conserved_regulon import build_conserved_regulon
from browser.conserved_operon import build_conserved_operon 
from browser.taxonomy import generate_genome_sunburst, get_taxon_children, generate_genes_sunburst
# Create your views here.
logger = logging.getLogger("GenomeDepot")


class AnnotationSearchResultsSubView(generic.ListView):
    '''
        Sub-page for AJAX-based view of annotation search result
    '''
    context_object_name = 'annotationlist'
    template_name = 'browser/annotation_list_subpage.html'
    # For testing in synchronous mode, uncomment the template below
    # template_name = 'browser/annotation_list.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        '''
            Generates context string
        '''
        context = super(AnnotationSearchResultsSubView,self)\
        .get_context_data(**kwargs)
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('annotation_query') + '"'
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
                (Q(source__icontains=annotation_query) | 
                 Q(value__icontains=annotation_query) | 
                 Q(note__icontains=annotation_query)) & 
                 Q(gene_id__genome__name=genome)
            ).order_by(
                'gene_id__locus_tag'
            ).select_related(
                'gene_id', 'gene_id__genome', 'gene_id__genome__taxon'
            ).prefetch_related('gene_id__genome__tags').distinct()
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) |
                Q(value__icontains=annotation_query) |
                Q(note__icontains=annotation_query)
            ).order_by(
                'gene_id__locus_tag'
            ).select_related(
                'gene_id', 'gene_id__genome', 'gene_id__genome__taxon'
            ).prefetch_related('gene_id__genome__tags').distinct()
        return object_list


class AnnotationSearchResultsAjaxView(View):
    '''
        AJAX-based view for annotation search results
    '''
    def get(self,request):
        '''
            Takes a GET request and returns a webpage that will send AJAX request
        '''
        context = {}
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('annotation_query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        for key, val in request.GET.items():
            context[key] = val
        return render(request,'browser/annotation_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, creates sub-page (AnnotationSearchResultsSubView)
            and sends it back as JSON
            
        '''
        start_time = time.time()

        context = {}
        sub_view = AnnotationSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()
        context['searchcontext'] = 'Search results for "' + \
        request.GET.get('annotation_query') + '"'
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
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Taxon.objects.filter(
                Q(name__icontains=query) | Q(taxonomy_id__exact=query)
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
        return Genome.objects.order_by('name').select_related(
            'strain', 'sample', 'taxon'
            ).prefetch_related('tags')

    def get_context_data(self,**kwargs):
        context = super(GenomeListView,self).get_context_data(**kwargs)
        if 'page_obj' not in context or context['page_obj'].number == 1:
            # Call generate_sunburst without parameters to let it
            # choose the root node
            sunburst = generate_genome_sunburst()
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
        genome_name = self.kwargs['genome']
        try:
            genome = Genome.objects.prefetch_related('tags').get(
                name = genome_name
                )
        except Genome.DoesNotExist:
            return {'searchcontext': 'Genome ' + genome_name + ' does not exist'}
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('query') + '" in ' + genome.name + ' operons'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Operon.objects.filter(genome__name=genome).filter(
                Q(genes__name__icontains=query) | 
                Q(genes__locus_tag__icontains=query) |
                Q(name__icontains=query)
            ).order_by(
                'name'
            ).select_related(
                'genome', 'contig'
            ).prefetch_related(
                'genes', 'genome__tags'
            ).distinct()
        else:
            object_list = Operon.objects.filter(genome__name=genome).order_by(
                'name'
            ).select_related(
                'genome', 'contig'
            ).prefetch_related(
                'genes', 'genome__tags'
            ).distinct()
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
        try:
            genome = Genome.objects.prefetch_related('tags').get(name = genome_name)
        except Genome.DoesNotExist:
            return {'searchcontext': 'Genome ' + genome_name + ' does not exist'}
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('query') + '" in ' + genome.name + ' sites'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Site.objects.filter(genome__name=genome).filter(
                Q(name__icontains=query) |
                Q(regulon__name__icontains=query)
            ).order_by(
                'name'
            ).select_related(
                'genome', 'contig'
            ).prefetch_related(
                'genome__tags'
            ).distinct()
        else:
            object_list = Site.objects.filter(genome__name=genome).order_by(
                'name'
            ).select_related(
                'genome', 'contig'
            ).prefetch_related(
                'genome__tags'
            ).distinct()
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
        try:
            genome = Genome.objects.prefetch_related('tags').get(name = genome_name)
        except Genome.DoesNotExist:
            return {'searchcontext': 'Genome ' + genome_name + ' does not exist'}
        context['genome'] = genome
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + \
            self.request.GET.get('query') + '" in ' + genome.name + ' regulons'
        return context

    def get_queryset(self):
        genome = self.kwargs['genome']
        if self.request.GET.get('query'):
            query = self.request.GET.get('query')
            object_list = Regulon.objects.filter(
                genome__name=genome
            ).filter(
                Q(name__icontains=query) |
                Q(description__icontains=query) |
                Q(regulators__locus_tag__icontains=query) |
                Q(regulators__name__icontains=query)
            ).order_by(
                'name'
            ).select_related(
                'genome'
            ).prefetch_related(
                'genome__tags'
            ).distinct()
        else:
            object_list = Regulon.objects.filter(
                genome__name=genome
            ).order_by(
                'name'
            ).select_related(
                'genome'
            ).prefetch_related(
                'genome__tags'
            ).distinct()
        return object_list


class GeneSearchResultsSubView(generic.ListView):
    '''
        Sub-page for AJAX-based view of gene search result
    '''
    context_object_name = 'genelist'
    template_name = 'browser/gene_list_subpage.html'
    # For testing in synchronous mode, uncomment the template below
    # template_name = 'browser/gene_list.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        '''
            Generates context string
        '''
        context = super(GeneSearchResultsSubView,self).get_context_data(**kwargs)
        genome = self.request.GET.get('genome')
        if genome:
            searchcontext = generate_gene_search_context(
                self.request.GET.get('query'),
                self.request.GET.get('type'),
                genome
            )
        else:
            searchcontext = generate_gene_search_context(
                self.request.GET.get('query'),
                self.request.GET.get('type')
            )
        context['searchcontext'] = searchcontext
        return context

    def get_queryset(self):
        '''
            Generates Gene queryset
        '''
        query = self.request.GET.get('query')
        query = query.strip()
        genome = self.request.GET.get('genome')
        query_type = self.request.GET.get('type')
        #logger.debug(', '.join([query, genome, query_type]))
        if query_type == 'gene':
            if genome:
                object_list = Gene.objects.filter(
                    locus_tag__exact=query
                ).filter(
                    genome__name=genome
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
                if not object_list:
                    object_list = Gene.objects.filter(
                        genome__name=genome
                    ).filter(
                        Q(name__icontains=query) |
                        Q(locus_tag__icontains=query) |
                        Q(function__icontains=query) |
                        Q(contig__contig_id__icontains=query)
                    ).order_by(
                        'locus_tag'
                    ).select_related(
                        'genome', 'genome__taxon'
                    ).prefetch_related(
                        'genome__tags'
                    ).distinct()
            elif query == '':
                object_list = Gene.objects.none()
            else:
                object_list = Gene.objects.filter(
                    locus_tag__exact=query
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
                if not object_list:
                    object_list = Gene.objects.filter(
                        Q(name__icontains=query) |
                        Q(locus_tag__icontains=query) |
                        Q(function__icontains=query) |
                        Q(genome__name__icontains=query) |
                        Q(contig__contig_id__icontains=query)
                    ).order_by(
                        'locus_tag'
                    ).select_related(
                        'genome', 'genome__taxon'
                    ).prefetch_related(
                        'genome__tags'
                    ).distinct()
        elif query_type == 'og':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__in=Protein.objects.filter(
                            ortholog_groups__id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                            ortholog_groups__id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'ko_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome,protein__in=Protein.objects.filter(
                        kegg_orthologs__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                        kegg_orthologs__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'kp_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome,protein__in=Protein.objects.filter(
                    kegg_pathways__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                            kegg_pathways__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'kr_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__in=Protein.objects.filter(
                        kegg_reactions__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                        kegg_reactions__kegg_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'ec_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__in=Protein.objects.filter(
                        ec_numbers__ec_number=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                        ec_numbers__ec_number=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'tc_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__in=Protein.objects.filter(
                        tc_families__tc_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                        tc_families__tc_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'cazy_id':
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__in=Protein.objects.filter(
                        cazy_families__cazy_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    protein__in=Protein.objects.filter(
                        cazy_families__cazy_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'cog_id':
            if genome:
                object_list = Gene.objects.filter(
                genome__name=genome, protein__in=Protein.objects.filter(
                        cog_classes__cog_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                protein__in=Protein.objects.filter(
                        cog_classes__cog_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
        elif query_type == 'go_id':
            if genome:
                object_list = Gene.objects.filter(
                genome__name=genome, protein__in=Protein.objects.filter(
                        go_terms__go_id=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                protein__in=Protein.objects.filter(
                        go_terms__go_id=query)
                ).order_by(
                'locus_tag'
                ).select_related(
                'genome', 'genome__taxon'
                ).prefetch_related(
                'genome__tags'
                )
        elif genome:
            if query_type == 'ko':
                ko_ids = Kegg_ortholog.objects.filter(
                    Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'kegg_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'kp':
                kp_ids = Kegg_pathway.objects.filter(
                    Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'kegg_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                ).prefetch_related(
                    'protein__kegg_orthologs', 'genome__tags'
                )
            elif query_type == 'kr':
                kr_ids = Kegg_reaction.objects.filter(
                    Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'kegg_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'ec':
                ec_ids = Ec_number.objects.filter(
                    Q(ec_number__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'ec_number'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            ec_numbers__ec_number__in=ec_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'tc':
                tc_ids = Tc_family.objects.filter(
                    Q(tc_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'tc_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            tc_families__tc_id__in=tc_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'cazy':
                cazy_ids = Cazy_family.objects.filter(
                    Q(cazy_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'cazy_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            cazy_families__cazy_id__in=cazy_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'cog':
                cog_ids = Cog_class.objects.filter(
                    Q(cog_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'cog_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            cog_classes__cog_id__in=cog_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            elif query_type == 'go':
                go_ids = Go_term.objects.filter(
                    Q(go_id__icontains=query) |
                    Q(description__icontains=query)
                ).values(
                    'go_id'
                )
                proteins = [item['protein_hash'] for item in Protein.objects.filter(
                            go_terms__go_id__in=go_ids).values('protein_hash')
                            ]
                object_list = Gene.objects.filter(
                genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
            else:
                object_list = Gene.objects.filter(
                    genome__name=genome
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                )
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
            searchcontext = generate_gene_search_context(self.request.GET.get('query'),
                                                         self.request.GET.get('type'),
                                                         self.request.GET.get('genome')
                                                         )
        else:
            searchcontext = generate_gene_search_context(self.request.GET.get('query'),
                                                         self.request.GET.get('type')
                                                         )
        context['searchcontext'] = searchcontext

        # copy request paramters into context
        for key, val in request.GET.items():
            context[key] = val
        return render(request,'browser/gene_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, creates sub-page (GeneSearchResultsSubView)
            and sends it back as JSON
            
            N.B.: For the testing of long-running tasks, uncomment sleep_timer=0
            and set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()

        # Sleep timer for testing only 
        # Uncomment next three lines to imitate long-running task with 30 sec waiting
        # sleep_timer = 30
        # logger.debug('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        # time.sleep(sleep_timer)

        context = {}
        sub_view = GeneSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()

        if request.GET.get('genome'):
            context['searchcontext'] = generate_gene_search_context(
                request.GET.get('query'),
                request.GET.get('type'),
                request.GET.get('genome')
                )
            external = generate_external_link(request.GET.get('query'),
                                              request.GET.get('type'),
                                              request.GET.get('genome')
                                              )
        else:
            context['searchcontext'] = generate_gene_search_context(
                request.GET.get('query'),
                request.GET.get('type')
                )
            external = generate_external_link(request.GET.get('query'),
                                              request.GET.get('type')
                                              )

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
            context['searchcontext'] = 'Search results for "' + \
                                       self.request.GET.get('query') + '"'
        elif self.request.GET.get('taxon'):
            taxon = Taxon.objects.get(taxonomy_id = self.request.GET.get('taxon'))
            context['searchcontext'] = 'Search results for ' + \
                                       taxon.name + ' [' + taxon.rank + ']'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        taxon = self.request.GET.get('taxon')
        if query:
            object_list = Genome.objects.filter(
                name__icontains=query
            ).distinct().order_by(
                'name'
            ).select_related(
                'strain', 'sample', 'taxon'
            ).prefetch_related(
                'tags'
            )
        elif taxon:
            children = get_taxon_children(taxon)
            object_list = Genome.objects.filter(
                taxon__taxonomy_id__in=children
            ).distinct().order_by(
                'name'
            ).select_related(
                'strain', 'sample', 'taxon'
            ).prefetch_related(
                'tags'
            )
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
            context['searchcontext'] = 'Search results for "' + \
                                       self.request.GET.get('query') + '"'
        elif self.request.GET.get('taxon'):
            taxon = Taxon.objects.get(taxonomy_id = self.request.GET.get('taxon'))
            context['searchcontext'] = 'Search results for ' + \
                                       taxon.name + ' [' + taxon.rank + ']'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        taxon = self.request.GET.get('taxon')
        if query:
            object_list = Strain.objects.filter(
                Q(strain_id__icontains=query) |
                Q(full_name__icontains=query) |
                Q(order__icontains=query)|
                Q(strain_metadata__value__icontains=query)
            ).distinct().order_by('strain_id')
        elif taxon:
            children = get_taxon_children(taxon)
            object_list = Strain.objects.filter(
                taxon__taxonomy_id__in=children
            ).distinct().order_by('strain_id')
        else:
            object_list = Strain.objects.none()
        return object_list


class SampleSearchResultsView(generic.ListView):
    '''
        Returns results of search in sample names, descriptions or metadata.
    '''
    model = Sample
    context_object_name = 'samplelist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(SampleSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + \
                                       self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        if query:
            object_list = Sample.objects.filter(
                Q(sample_id__icontains=query) |
                Q(full_name__icontains=query) |
                Q(description__icontains=query)|
                Q(sample_metadata__value__icontains=query)
            ).distinct().order_by('sample_id')
        else:
            object_list = Sample.objects.none()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self):
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Kegg_ortholog.objects.filter(
                    (Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('kegg_id').distinct()
            else:
                object_list = Kegg_ortholog.objects.filter(
                    protein__in = proteins
                ).order_by('kegg_id').distinct()
        elif query:
            object_list = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).order_by('kegg_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Kegg_pathway.objects.filter(
                    (Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('kegg_id').distinct()
            else:
                object_list = Kegg_pathway.objects.filter(
                    protein__in = proteins
                ).order_by('kegg_id').distinct()
        elif query:
            object_list = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=query) | Q(description__icontains=query)
            ).order_by('kegg_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Kegg_reaction.objects.filter(
                    (Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('kegg_id').distinct()
            else:
                object_list = Kegg_reaction.objects.filter(
                    protein__in = proteins
                ).order_by('kegg_id').distinct()
        elif query:
            object_list = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).order_by('kegg_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Ec_number.objects.filter(
                    (Q(ec_number__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('ec_number').distinct()
            else:
                object_list = Ec_number.objects.filter(
                    protein__in = proteins
                ).order_by('ec_number').distinct()
        elif query:
            object_list = Ec_number.objects.filter(
                Q(ec_number__icontains=query) |
                Q(description__icontains=query)
                ).order_by('ec_number').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Tc_family.objects.filter(
                    (Q(tc_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('tc_id').distinct()
            else:
                object_list = Tc_family.objects.filter(
                    protein__in = proteins
                ).order_by('tc_id').distinct()
        elif query:
            object_list = Tc_family.objects.filter(
                Q(tc_id__icontains=query) | Q(description__icontains=query)
            ).order_by('tc_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Cazy_family.objects.filter(
                    (Q(cazy_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('cazy_id').distinct()
            else:
                object_list = Cazy_family.objects.filter(
                    protein__in = proteins
                ).order_by('cazy_id').distinct()
        elif query:
            object_list = Cazy_family.objects.filter(
                Q(cazy_id__icontains=query) |
                Q(description__icontains=query)
            ).order_by('cazy_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Go_term.objects.filter(
                    (Q(go_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('go_id').distinct()
            else:
                object_list = Go_term.objects.filter(
                    protein__in = proteins
                ).order_by('go_id').distinct()
        elif query:
            object_list = Go_term.objects.filter(
                Q(go_id__icontains=query) |
                Q(description__icontains=query)
            ).order_by('go_id').distinct()
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
        if self.request.GET.get('genome'):
            context['genome'] = self.request.GET.get('genome')
            context['searchcontext'] = 'Search results for ' \
                                       + self.request.GET.get('genome') + ' genome'
        if self.request.GET.get('query'):
            if self.request.GET.get('genome'):
                context['searchcontext'] = 'Search results for "' \
                                           + self.request.GET.get('query') + '" in ' \
                                           + self.request.GET.get('genome') + ' genome'
            else:
                context['searchcontext'] = 'Search results for "' + \
                                           self.request.GET.get('query') + '"'
        return context

    def get_queryset(self): # new
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        if genome:
            proteins = Gene.objects.filter(
                genome__name = genome
            ).values_list("protein__id", flat=True).distinct()
            if query:
                object_list = Cog_class.objects.filter(
                    (Q(cog_id__icontains=query) |
                    Q(description__icontains=query)) &
                    Q(protein__in = proteins)
                ).order_by('cog_id').distinct()
            else:
                object_list = Cog_class.objects.filter(
                    protein__in = proteins
                ).order_by('cog_id').distinct()
        elif query:
            object_list = Cog_class.objects.filter(
                Q(cog_id__icontains=query) |
                Q(description__icontains=query)
            ).order_by('cog_id').distinct()
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
    template = loader.get_template('browser/index.html')
    #num_genomes = Genome.objects.all().count()
    context = {'site_title':TITLE}
    return HttpResponse(template.render(context, request))

def show_help(request):
    '''
        Displays help page
    '''
    template = loader.get_template('browser/help.html')
    context = {'site_title':TITLE}
    return HttpResponse(template.render(context, request))


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
        object_list = Genome.objects.filter(
            tags__id=self.genometag.id
        ).order_by(
            'name'
        ).select_related(
            'strain', 'taxon'
        ).prefetch_related(
            'tags'
        )
        return object_list
    
def taxon_detail(request, taxonomy_id):
    '''
        Displays taxon page.
    '''
    try:
        taxon = Taxon.objects.get(taxonomy_id = taxonomy_id)
    except Taxon.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Taxon ' + taxonomy_id + ' does not exist'}
                      )
    context = {'taxon': taxon}
    children = get_taxon_children(taxonomy_id)
    context['sunburst'] = generate_genome_sunburst(taxonomy_id, children)
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
        genome = Genome.objects.select_related(
            'taxon', 'strain', 'sample'
        ).prefetch_related(
            'tags'
        ).get(
            name = name
        )
    except Genome.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Genome ' + name + ' does not exist'}
                      )
    context = {'genome': genome}
    context['operons'] = Operon.objects.filter(genome=genome).count()
    context['sites'] = Site.objects.filter(genome=genome).count()
    context['regulons'] = Regulon.objects.filter(genome=genome).count()
    lineage = [genome.taxon,]
    parent_id = genome.taxon.parent_id
    iteration_count = 0
    while True:
        try:
            parent_taxon = Taxon.objects.get(taxonomy_id = parent_id)
        except Taxon.DoesNotExist:
            logger.error('Taxonomy ID ' + parent_id + ' does not exist. Update taxonomy records.')
            break
        iteration_count += 1
        if parent_taxon.taxonomy_id == parent_taxon.parent_id or parent_id == '1':
            break
        lineage.append(parent_taxon)
        parent_id = parent_taxon.parent_id
    context['lineage'] = reversed(lineage)
    
    if request.GET.get('contig'):
        context['highlight_start'] = request.GET.get('start')
        context['highlight_end'] = request.GET.get('end')
        display_offset = 5000
        contig = Contig.objects.get(contig_id = request.GET.get('contig'),
                                    genome=genome
                                    )
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
    context = {}
    try:
        strain = Strain.objects.get(id = strain_id)
        context['strain'] = strain
        genomes = Genome.objects.filter(
            strain=strain
        ).order_by(
            'name'
        ).prefetch_related(
            'tags'
        )
        context['genomes'] = genomes
        strain_metadata = Strain_metadata.objects.filter(strain=strain)
        metadata_entries = {}
        for metadata_entry in strain_metadata:
            if metadata_entry.source not in metadata_entries:
                metadata_entries[metadata_entry.source] = {}
                metadata_entries[metadata_entry.source]['entries'] = []
            metadata_entries[metadata_entry.source]['source'] = metadata_entry.source
            metadata_entries[metadata_entry.source]['url'] = metadata_entry.url
            metadata_entries[metadata_entry.source]['entries']\
            .append({'key':metadata_entry.key,'value':metadata_entry.value})
        metadata = []
        for source in sorted(metadata_entries.keys()):
            metadata.append(metadata_entries[source])
        context['metadata'] = metadata
        lineage = [strain.taxon,]
        parent_id = strain.taxon.parent_id
        iteration_count = 0
        while True:
            parent_taxon = Taxon.objects.get(taxonomy_id = parent_id)
            iteration_count += 1
            if parent_taxon.taxonomy_id == parent_taxon.parent_id or parent_id == '1':
                break
            lineage.append(parent_taxon)
            parent_id = parent_taxon.parent_id
        context['lineage'] = reversed(lineage)
    except Strain.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Sample not found: ' + str(strain_id)}
                      )
    return render(request,
                  'browser/strain.html',
                  context
                  )

def sample_detail(request, sample_id):
    '''
        Displays sample page.
    '''
    try:
        sample = Sample.objects.get(id = sample_id)
        genomes = Genome.objects.filter(
            sample=sample
        ).order_by(
            'name'
        ).prefetch_related(
            'tags'
        )
        sample_metadata = Sample_metadata.objects.filter(sample=sample)
        metadata_entries = {}
        for metadata_entry in sample_metadata:
            if metadata_entry.source not in metadata_entries:
                metadata_entries[metadata_entry.source] = {}
                metadata_entries[metadata_entry.source]['entries'] = []
            metadata_entries[metadata_entry.source]['source'] = metadata_entry.source
            metadata_entries[metadata_entry.source]['url'] = metadata_entry.url
            metadata_entries[metadata_entry.source]['entries']\
            .append({'key':metadata_entry.key,'value':metadata_entry.value})
        metadata = []
        for source in sorted(metadata_entries.keys()):
            metadata.append(metadata_entries[source])
    except Sample.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Sample not found: ' + str(sample_id)}
                      )
    return render(request,
                  'browser/sample.html',
                  {'sample': sample, 'genomes':genomes, 'metadata':metadata}
                  )

def gene_detail(request, genome, locus_tag):
    '''
        Displays gene page.
    '''
    try:
        gene = Gene.objects.select_related(
            'genome', 'genome__taxon', 'protein', 'operon'
        ).prefetch_related(
            'protein__ortholog_groups__taxon',
            'protein__kegg_orthologs',
            'protein__kegg_pathways',
            'protein__kegg_reactions',
            'protein__ec_numbers',
            'protein__go_terms',
            'protein__tc_families',
            'protein__cog_classes',
            'genome__tags'
        ).get(
            genome__name=genome, locus_tag = locus_tag
        )
        annotations = Annotation.objects.filter(
            gene_id = gene
        ).order_by(
            'source', 'value'
        )
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
        return render(request,
                      '404.html',
                      {'searchcontext': 'Gene ' + locus_tag + ' does not exist'}
                      )
    return render(request, 'browser/gene.html', context)

def operon_detail(request, name, genome=None):
    '''
        Displays operon page.
    '''
    try:
        operon = Operon.objects.select_related(
            'contig', 'genome'
        ).prefetch_related(
            'genome__tags'
        ).get(
            name=name
        )
        genes = Gene.objects.filter(
            operon = operon.id
        ).select_related(
            'contig', 'genome', 'genome__strain'
        )
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
        return render(request,
                      '404.html',
                      {'searchcontext': 'Operon ' + name + ' does not exist'}
                      )
    return render(request, 'browser/operon.html', context)

def site_detail(request, genome, name):
    '''
        Displays site page.
    '''
    try:
        site = Site.objects.select_related(
            'genome', 'genome__strain', 'contig'
        ).prefetch_related(
            'genome__tags'
        ).get(
            genome__name=genome, name = name
        )
        context = {'site': site}
        viewer_start = site.start - 500
        if viewer_start < 0:
            viewer_start = 1
        viewer_end = site.end + 500
        if viewer_end > site.contig.size:
            viewer_end = site.contig.size
        context['viewer_start'] = str(viewer_start)
        context['viewer_end'] = str(viewer_end)
        context['highlight_start'] = site.start
        context['highlight_end'] = site.end
    except Site.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Site ' + name + ' does not exist'}
                      )
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
        return render(request,
                      '404.html',
                      {'searchcontext': 'Regulon ' + name + ' does not exist'}
                      )
    sites = Site.objects.filter(regulon = regulon).select_related(
        'contig').prefetch_related('genes', 'operons')
    context['sites'] = sites
    ortholog_groups = set()
    for regulator in regulon.regulators.all():
        for og in regulator.protein.ortholog_groups.all():
            ortholog_groups.add(og.id)
    context['ortholog_groups'] = Ortholog_group.objects.filter(
        id__in=list(ortholog_groups)
        )
    return render(request, 'browser/regulon.html', context)


def og_detail(request, og_id):
    context = {}
    try:
        ortholog_group = Ortholog_group.objects.get(id=og_id)
    except Ortholog_group.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Ortholog group does not exist'}
                      )
    context['ortholog_group'] = ortholog_group
    return render(request, 'browser/family.html', context)


def gene_byname(request):
    '''
        Displays gene page.
    '''
    if request.GET.get('genome') and request.GET.get('locus'):
        genome_name = request.GET.get('genome')
        locus_tag = request.GET.get('locus')
    else:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Genome name or locus ID not provided'}
                      )
    try:
        gene = Gene.objects.get(
            locus_tag = locus_tag, genome__name = genome_name
        )
        annotations = Annotation.objects.filter(
            gene_id = gene
        ).order_by(
            'source', 'value'
        )
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
        return render(request,
                      '404.html',
                      {'searchcontext': 'Gene not found for genome ' + genome_name +
                      ' and locus tag ' + locus_tag}
                      )
    except Gene.MultipleObjectsReturned:
        return render(request,
                      '404.html',
                      {'searchcontext': 'More than one gene was found'}
                      )
    return render(request, 'browser/gene.html', context)

def protein_search_external(request):
    '''
        Returns protein search results for external query (from Fitness browser etc.).
    '''
    context = {}
    if request.GET.get("sequence"):
        result = {}
        sequence = request.GET.get("sequence")
        hits, searchcontext, query_len, _ = run_protein_search(sequence)
        if searchcontext != '':
            context['searchcontext'] = searchcontext
        for row in hits:
            row=row.split('\t')
            logger.info('Search for gene %s', row[1])
            if row[0] not in result:
                result[row[0]] = []
            unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
            query_cov = (query_len - unaligned_part) * 100.0 / query_len
            genes = Gene.objects.select_related(
                'protein', 'genome', 'genome__taxon'
            ).prefetch_related(
                'genome__tags'
            ).filter(
                protein__protein_hash = row[1]
            )
            for target in genes:
                hit = [target.genome.name,
                       target.locus_tag,
                       target.genome.taxon.name,
                       target.function,
                       '{:.1f}'.format(float(row[2])),
                       row[3],
                       '{:.1f}'.format(query_cov),
                       row[10],
                       row[11]
                       ]
                result[row[0]].append(hit)
        context['searchresult'] = result
    else:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Provide protein sequence in FASTA format'}
                      )
    return render(request, 'browser/proteinsearch.html', context)


class NsearchResultView(View):
    '''
        AJAX-based view for nucleotide similarity search
    '''
    def post(self,request):
        '''
            Takes a POST request and returns a webpage that will send AJAX request
        '''
        #logger.debug(('Sequence sent: %s', request.POST.get("sequence"))
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        for key, val in request.POST.items():
            context[key] = val
            #logger.debug('%s:%s', key, val)
        return render(request,'browser/nucleotidesearchajax.html', context)

    def get(self,request):
        '''
            Displays nucleotide sequence search form
        '''
        return render(request,'browser/nucleotidesearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, calls run_nucleotide_search and 
            sends formatted results back as JSON
            
            N.B.: For the testing of long-running tasks, uncomment sleep_timer=0 and 
            set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()
        result = []
        #logger.debug('REQUEST2')
        #for key, val in request.POST.items():
        #    logger.debug('%s:%s',key, val)
        params = {'sequence': request.POST.get("sequence"),
                  'evalue': request.POST.get("evalue"),
                  'hitstoshow': request.POST.get("hitstoshow")
                  }
        #sleep_timer = 0
        #logger.debug('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)
        hits, searchcontext, query_len, query_name = run_nucleotide_search(params)
        if hits:
            result.append('<table><thead><tr><th>Target contig (click to see hit)' +
                          '</th><th>Genome</th><th>%identity</th><th>Alignment length'+
                          '</th><th>%Query coverage</th><th>E-value</th><th>Bit-score' +
                          '</th></tr></thead><tbody>'
                          )
            for row in hits:
                row=row.split('\t')
                contig_name, genome_name = row[1].split('|')
                genome_name = genome_name.split('[')[0]
                target = Contig.objects.select_related(
                    'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                ).get(
                    contig_id = contig_name, genome__name = genome_name
                )
                strand = '+'
                start = row[8]
                end = row[9]
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                if int(row[7]) < int(row[6]):
                    strand = '-'
                genome_tags = ''
                for genome_tag in target.genome.tags.all():
                    genome_tags += '<span class="genometag" style="background-color:'\
                    + genome_tag.color + '"><a href="' + \
                    reverse('tagdetails', args=(genome_tag.name,)) + \
                    '" style="color:' + genome_tag.textcolor + \
                    '" title="' + genome_tag.description + '">' + \
                    genome_tag.name + '</a></span>&nbsp;'
                hit = '<tr><td align=\"left\"><a href=\"' + \
                      reverse('genomedetails', args=(target.genome.name,)) + \
                      '?contig=' + target.contig_id + '&start=' + start + '&end=' + \
                      end + '\">' + contig_name + ': (' + strand + 'strand) ' + \
                      start + '..' + end + '</a></td><td align="left">' + \
                      target.genome.name + ' [' + target.genome.taxon.name + ']' + \
                      genome_tags + '</td><td>' + '{:.2f}'.format(float(row[2])) + \
                      '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) \
                      + '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                result.append(hit)
            result.append('</tbody></table>')
            context = {"searchresult":'\n'.join(result),
                       "searchcontext":searchcontext,
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' bp',
                       "time":time.time()-start_time
                       }
        elif searchcontext == '':
            context = {"searchresult":'',"searchcontext":'No hits found',
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' bp',
                       "time":time.time()-start_time
                       }
        else:
            context = {"searchresult":"",
                       "searchcontext":searchcontext,
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' bp',
                       "time":time.time()-start_time
                       }
        #logger.debug(context)
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
        sequence = request.POST.get("sequence")
        #logger.debug('Sequence sent: %s', sequence)
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        #logger.debug('REQUEST1')
        for key, val in request.POST.items():
            context[key] = val
            #logger.debug('%s:%s', key, val)
        return render(request,'browser/proteinsearchajax.html', context)

    def get(self,request):
        '''
            Displays protein sequence search form
        '''
        return render(request,'browser/proteinsearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        '''
            Takes AJAX request, calls run_protein_search and 
            sends formatted results back as JSON
            
            N.B.: For the testing of long-running tasks, uncomment sleep_timer=0 and
            set the number of seconds to delay the AJAX response
        '''
        start_time = time.time()
        result = []
        #logger.debug('REQUEST2')
        #for key, val in request.POST.items():
        #    logger.debug('%s,%s', key, val)
        params = {'sequence': request.POST.get("sequence"),
                  'evalue': request.POST.get("evalue"),
                  'hitstoshow': request.POST.get("hitstoshow")
                  }
        #sleep_timer = 0
        #logger.debug('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)
        hits, searchcontext, query_len, query_name = run_protein_search(params)

        if hits:
            result.append('<table><thead><tr><th>Target gene</th><th>Genome</th>' +\
                          '<th>Function</th><th>%identity</th><th>Alignment length</th><th>' +\
                          '%Query coverage</th><th>E-value</th><th>Bit-score</th>' +\
                          '</tr></thead><tbody>'
                          )
            for row in hits:
                row=row.split('\t')
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                genes = Gene.objects.select_related(
                    'protein', 'genome', 'genome__taxon'
                ).prefetch_related(
                    'genome__tags'
                ).filter(
                    protein__protein_hash = row[1]
                )
                for target in genes:
                    genome_tags = ''
                    for genome_tag in target.genome.tags.all():
                        genome_tags += '<span class="genometag" ' +\
                        'style="background-color:' + genome_tag.color + \
                        '"><a href="' + \
                        reverse('tagdetails', args=(genome_tag.name,)) + \
                        '" style="color:' + genome_tag.textcolor + \
                        '" title="' + genome_tag.description + '">' + \
                        genome_tag.name + '</a></span>&nbsp;'
                    hit = '<tr><td align=\"left\"><a href=\"' + \
                    reverse('genedetails',
                            args=(target.genome.name, target.locus_tag)) + \
                    '\">' + target.locus_tag + '</a></td><td align="left">' + \
                     '<a href="' + reverse('genomedetails', args=(target.genome.name,)) + \
                     '" title="' + target.genome.taxon.name + '">' + target.genome.name +  '</a> ' +\
                     genome_tags + '</td><td>' + target.function + '</td><td>' + '{:.1f}'.format(float(row[2])) + \
                    '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) + \
                    '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                    result.append(hit)
            result.append('</tbody></table>')
            context = {"searchresult":'\n'.join(result),
                       "searchcontext":searchcontext,
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' aa',
                       "time":time.time()-start_time
                       }
        elif searchcontext == '':
            context = {"searchresult":'',
                       "searchcontext":'No hits found',
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' aa',
                       "time":time.time()-start_time
                       }
        else:
            context = {"searchresult":"",
                       "searchcontext":searchcontext,
                       "query_len":query_len,
                       "query_name":'Query: ' + query_name + ', ' +  str(query_len) + ' aa',
                       "time":time.time()-start_time
                       }
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


def proteinsearchform(request):
    '''
        Displays protein sequence search form
    '''
    return render(request,'browser/proteinsearchform.html')

def cregulon_view(request):
    '''
        Conserved regulon view
    '''
    og_id = request.GET.get('og')
    #locus_tag = request.GET.get('locus_tag')
    context = build_conserved_regulon(og_id)
    return render(request, 'browser/cregulon.html', context)

    
def conserved_operon_view(request, operon_id):
    '''
        Conserved operon view
    '''
    context = {}
    try:
        context['start_operon'] = Operon.objects.get(id=operon_id)
    except Operon.DoesNotExist:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Operon ' + str(operon_id) + ' does not exist'}
                      )
    return render(request, 'browser/coperon.html', context)


def conserved_operon_data(request, operon_id):
    '''
        Returns conserved operon data for Ajax request
    '''
    context = {}
    treemap, sunburst, operon_ids, functional_profile_tsv = build_conserved_operon(operon_id)
    context['treemap'] = treemap
    context['sunburst'] = sunburst
    context['tsv_profile'] = functional_profile_tsv
    object_list = Operon.objects.filter(id__in=operon_ids).order_by(
        'name'
    ).select_related(
        'genome', 'contig'
    ).prefetch_related(
        'genes', 'genome__tags'
    )
    context['operonlist'] = loader.render_to_string('browser/operon_list_subview.html', {'operonlist':object_list})
    data = json.dumps(context)
    return HttpResponse(data,content_type="application/json")

    
def pathway_view(request):
    context = {}
    items = []
    
    if request.GET.get('genome') and request.GET.get('pathway'):
        genome = request.GET.get('genome')
        query = request.GET.get('pathway')
        kp = Kegg_pathway.objects.get(kegg_id=query)
        #kp_id = kp.kegg_id
        ext_kegg_map_url = ''
        context['pathway'] = kp
        context['genome'] = genome
        gene_list = Gene.objects.filter(
                    genome__name=genome, protein__kegg_pathways__kegg_id=kp.kegg_id
                    ).select_related(
                        'protein', 'genome'
                    ).prefetch_related(
                        'protein__kegg_orthologs'
                    )
        #logger.debug(kp.kegg_id)
        kegg_map_url = 'https://www.kegg.jp/kegg-bin/show_pathway?' + \
                       kp.kegg_id + '/'
        #logger.debug(kegg_map_url)
        ko_ids = {}
        for gene in gene_list:
            for ko in gene.protein.kegg_orthologs.all():
                if ko.kegg_id not in ko_ids:
                    ko_ids[ko.kegg_id] = []
                ko_ids[ko.kegg_id].append(gene)
        if ko_ids:
            ko_count = 0
            for ko_id in sorted(ko_ids.keys()):
                item = {}
                item['ko'] = ko_id
                item['genes'] = ko_ids[ko_id]
                items.append(item)
                color = COLORS[ko_count % len(COLORS)]
                kegg_map_url += ko_id + '%09%23' + color[0] + ',%23' + color[1] + '/'
                item['bg'] = color[0]
                item['fg'] = color[1]
                ko_count += 1
        context['iframeheight'] = 850
        if ko_count > 18:
            context['iframeheight'] = 50*ko_count
        kegg_map_url += 'default%3dpink'
        
        if len(kegg_map_url) > 2559:
            context['error'] = 'The pathway URL cannot be displayed by KEGG mapper tool because it is too long (' + str(len(kegg_map_url)) + ' symbols). Try another pathway with smaller number of genes.'
            return render(request, 'browser/pathway.html', context)
        context['external'] = kegg_map_url
        context['items'] = items

    return render(request, 'browser/pathway.html', context)
    
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
                raise SuspiciousOperation(
                    "Unacceptable value '%s' for size parameter." % size
                    )
                
            if int(size) not in [5, 10, 20, 40, 60, 80, 100]:
                raise SuspiciousOperation(
                    "Unacceptable value for locus size: '%s'" % size
                    )
                
            if lines is None:
                raise SuspiciousOperation('"lines" parameter is missing')
                
            try:
                lines = int(lines)
            except TypeError:
                raise SuspiciousOperation(
                    "Unacceptable value '%s' for lines number." % lines
                    )
                
            if int(lines) not in [10, 25, 50, 75, 100, 200]:
                raise SuspiciousOperation(
                    "Unacceptable value for lines number: '%s'" % lines
                    )
                
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
            ComparativeView.verify_parameters(request.GET.get('size'),
                                              request.GET.get('lines')
                                              )
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
        
        #logger.debug('REQUEST1')
        #logger.debug('Request parameters: %s %s %s %s %s', og_id, genome,
        #             locus_tag, request.GET.get('size'), request.GET.get('lines')
        #             )
        context = {}
        for key, val in request.GET.items():
            context[key] = val
        try:
            gene = Gene.objects.select_related(
                'protein', 'genome__strain', 'genome__sample', 'contig'
            ).get(
                locus_tag=locus_tag, genome__name=genome
            )
            og = Ortholog_group.objects.get(id=og_id)
        except Gene.DoesNotExist:
            return render(request,
                          '404.html',
                          {'searchcontext': 'Gene not found'}
                          )
        except Ortholog_group.DoesNotExist:
            return render(request,
                          '404.html',
                          {'searchcontext': 'Ortholog group not found'}
                          )
        context['gene'] = gene
        context['size' + request.GET.get('size')] = '1'
        context['lines' + request.GET.get('lines')] = '1'
        context['ortholog_group'] = og
        return render(request,'browser/scriblajax.html', context)

    @staticmethod
    def ajax_view(request):
        start_time = time.time()

        #logger.debug('REQUEST2')
        #for key, val in request.GET.items():
        #    logger.debug('%s:%s', key, val)
        try:
            ComparativeView.verify_parameters(request.GET.get('size'),
                                              request.GET.get('lines')
                                              )
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
            
        # Sleep timer for testing to imitate long-running task
        #sleep_timer = 0
        #logger.debug('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        #time.sleep(sleep_timer)

        context = {}
        locus_tag = request.GET.get('locus_tag')
        genome = request.GET.get('genome')
        og_id = request.GET.get('og')

        try:
            gene = Gene.objects.select_related(
                'protein', 'genome__strain', 'genome__sample', 'contig'
            ).get(
                locus_tag=locus_tag, genome__name=genome
            )
            og = Ortholog_group.objects.get(id=og_id)
        except Gene.DoesNotExist:
            return render(request,
                          '404.html',
                          {'searchcontext': 'Gene not found'}
                          )
        except Ortholog_group.DoesNotExist:
            return render(request,
                          '404.html',
                          {'searchcontext': 'Ortholog group not found'}
                          )
        scribl, tree_canvas, tree_newick, og_gene_count, plot_gene_count, treemap_gene_ids = \
            get_scribl(gene, og, request)
            
        treemap = ''
        if treemap_gene_ids:
            #print(len(treemap_gene_ids), 'genes for treemap generation')
            treemap, functional_profile = generate_genes_treemap(treemap_gene_ids)

        if og_gene_count == 1:
            scribl='<script type="text/javascript">\nalert("Comparative plot cannot' +\
            ' be created for only one genome");</script>'
            context = {'scribl':scribl,
                       'og_gene_count':og_gene_count,
                       'plot_gene_count':plot_gene_count,
                       "time":time.time()-start_time
                       }
        else:
            scribl='<script type="text/javascript">\nfunction draw(canvasName) ' +\
            '{\nvar canvas = document.getElementById(canvasName);\nvar parent = ' +\
            'document.getElementById("scribl-container");\nvar ctx = canvas.' +\
            'getContext("2d");\n' + scribl + '\nchart.draw();\nvar img = canvas.' +\
            'toDataURL("image/png");\ndocument.getElementById("pngexport").href = ' +\
            'img;\n}\n</script>'
            plot_gene_count -= 1
            context = {'scribl':scribl,
                       'tree_canvas':tree_canvas,
                       'tree_newick':tree_newick,
                       'og_gene_count':og_gene_count,
                       'plot_gene_count':plot_gene_count,
                       'time':time.time()-start_time,
                       'treemap':treemap,
                       'tsv_profile': functional_profile
                       }
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")

    
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
            if query == '':
                searchcontext = 'Query string is empty'
            else:
                searchcontext = 'Search results for "' + query + '"'
        elif query_type=='og':
            eggnog_og = Ortholog_group.objects.get(id=query)
            searchcontext = 'Genes from Ortholog Group ' + eggnog_og.eggnog_id + \
                            '[' + eggnog_og.taxon.name + ']'
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
            searchcontext = 'Genes assigned to KEGG Ortholog groups containing "' +\
                            query + '"'
        elif query_type=='kp':
            searchcontext = 'Genes assigned to KEGG pathways containing "' + \
                            query + '"'
        elif query_type=='kr':
            searchcontext = 'Genes assigned to KEGG reactions containing "' + \
                            query + '"'
        elif query_type=='ec':
            searchcontext = 'Genes assigned to EC numbers containing "' + \
                            query + '"'
        elif query_type=='tc':
            searchcontext = 'Genes assigned to TCDB families containing "' + \
                            query + '"'
        elif query_type=='cazy':
            searchcontext = 'Genes assigned to CAZy families containing "' + \
                            query + '"'
        elif query_type=='cog':
            searchcontext = 'Genes assigned to COG class containing "' + query + '"'
        elif query_type=='go':
            searchcontext = 'Genes assigned to GO terms containing "' + query + '"'
        else:
            searchcontext = 'This query type is not supported'
    else:
        if query_type=='gene':
            searchcontext = 'Genes from genome ' + genome
        elif query_type=='ko_id':
            searchcontext = 'Genes from KEGG Ortholog Group ' + query + \
                            ' in genome ' +  genome
        elif query_type=='kp_id':
            searchcontext = 'Genes from KEGG pathway ' + query + \
                            ' in genome ' +  genome
        elif query_type=='kr_id':
            searchcontext = 'Genes from KEGG reaction ' + query + \
                            ' in genome ' +  genome
        elif query_type=='ec_id':
            searchcontext = 'Genes with EC number ' + query + \
                            ' in genome ' +  genome
        elif query_type=='tc_id':
            searchcontext = 'Genes from TCDB family ' + query + \
                            ' in genome ' +  genome
        elif query_type=='cazy_id':
            searchcontext = 'Genes from CAZy family ' + query + \
                            ' in genome ' +  genome
        elif query_type=='cog_id':
            searchcontext = 'Genes from COG class ' + query + \
                            ' in genome ' +  genome
        elif query_type=='go_id':
            searchcontext = 'Genes assigned to GO term ' + query + \
                            ' in genome ' +  genome
        else:
            searchcontext = 'This query type is not supported'
    return searchcontext

    
def get_og_data(request):
    og_id = request.GET.get('og')
    #print('AJAX request for', og_id)
    ortholog_group = Ortholog_group.objects.get(id=og_id)
    treemap, functional_profile, gene_ids = generate_og_treemap(ortholog_group)
    context = {'treemap':treemap, 'og_gene_count':str(len(gene_ids)), 'tsv_profile':functional_profile}
    sunburst = generate_genes_sunburst(gene_ids)
    context['sunburst'] = sunburst
    data = json.dumps(context)
    return HttpResponse(data,content_type="application/json")

    
def generate_external_link(query, query_type, genome=None):
    '''
        Generates search context string and external link for various query types
    '''
    external = ''
    link_text = ''
    if query is None:
        query = ''
    else:
        if query_type=='kp' or query_type=='kp_id':
            if genome is None:
                if query_type=='kp_id':
                    external = '<a href="https://www.kegg.jp/dbget-bin/www_bget?' + query  + '" target="blank_">Search for "' + query + '" in KEGG</a>'
            else:
                kp_ids = Kegg_pathway.objects.filter(
                    Q(kegg_id__icontains=query) |
                    Q(description__icontains=query)
                ).values('kegg_id')
                ext_kegg_map_url = ''
                if kp_ids.count() == 1:
                    gene_list = Gene.objects.filter(
                        genome__name=genome, protein__kegg_pathways__kegg_id__in=kp_ids
                    ).select_related(
                        'protein'
                    ).prefetch_related(
                        'protein__kegg_orthologs'
                    )
                    #logger.debug(kp_ids[0]['kegg_id'])
                    kegg_map_url = 'https://www.kegg.jp/pathway/' + \
                                   kp_ids[0]['kegg_id'] + '+'
                    #logger.debug(kegg_map_url)
                    ko_ids = set()
                    for gene in gene_list:
                        for ko in gene.protein.kegg_orthologs.all():
                            ko_ids.add(ko.kegg_id)
                    if ko_ids:
                        ext_kegg_map_url = kegg_map_url + \
                                           '+'.join(sorted(list(ko_ids)))
                    link_text = 'View KEGG map ' + kp_ids[0]['kegg_id'] + ' for these functions'
                if ext_kegg_map_url != '':
                    external = ext_kegg_map_url
                if external != '':
                    external = '<a href="'+ reverse('pathway') + '?genome=' + genome + '&pathway=' + kp_ids[0]['kegg_id'] + '">' + link_text + '</a>'  #'<a href="' + external + '" target="blank_">' + link_text + '</a>'
        elif query_type=='ko_id':
            external = '<a href="https://www.kegg.jp/dbget-bin/www_bget?' + query + '" target="blank_">Search for "' + query + '" in KEGG</a>'
        elif query_type=='kp_id':
            external = '<a href="https://www.kegg.jp/dbget-bin/www_bget?' + query + '" target="blank_">Search for "' + query + '" in KEGG</a>'
        elif query_type=='kr_id':
            external = '<a href="https://www.kegg.jp/dbget-bin/www_bget?' + query + '" target="blank_">Search for "' + query + '" in KEGG</a>'
        elif query_type=='ec_id':
            external = '<a href="https://www.kegg.jp/dbget-bin/www_bget?ec:' + query + '" target="blank_">Search for "' + query + '" in KEGG</a>'
        elif query_type=='tc_id':
            external = '<a href="http://www.tcdb.org/search/result.php?tc=' + query + '#' + query + '" target="blank_">Search for "' + query + '" in TCDB</a>'
        elif query_type=='cazy_id':
            external = '<a href="http://www.cazy.org/' + query + '.html" target="blank_">View "' + query + '" page in CAZy</a>'
        elif query_type=='go_id':
            external = '<a href="https://www.ebi.ac.uk/QuickGO/search/' + query + '" target="blank_">Search for "' + query + '" in Gene Onthology</a>'
        else:
            external = ''
    return external

