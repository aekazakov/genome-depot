import csv
import time
import json
from django.shortcuts import render
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

# Create your views here.

class AnnotationSearchResultsSubView(generic.ListView):
    context_object_name = 'annotationlist'
    template_name = 'browser/annotation_list_subpage.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(AnnotationSearchResultsSubView,self).get_context_data(**kwargs)
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('annotation_query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        return context

    def get_queryset(self):
        annotation_query = self.request.GET.get('annotation_query')
        genome = self.request.GET.get('genome')
        if not annotation_query:
            object_list = Annotation.objects.none()
        elif genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
        return object_list


class AnnotationSearchResultsAjaxView(View):
    def get(self,request):
        context = {} #super(AnnotationSearchResultsAjaxView,self).get_context_data(**kwargs)
        #annotation_query = self.request.GET.get('annotation_query')
        if self.request.GET.get('annotation_query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('annotation_query') + '"'
        else:
            context['searchcontext'] = 'Query string is empty'
        # print('REQUEST1')
        # print('Request parameters:', request.GET.get('annotation_query'), request.GET.get('genome'))
        for key, val in request.GET.items():
            # print(key, val)
            context[key] = val
        return render(request,'browser/annotation_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        start_time = time.time()

        # print('REQUEST2')
        # for key, val in request.GET.items():
        #     print(key, val)
            
        # Sleep timer for testing to imitate long-running task
        # sleep_timer = 0
        # print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        # time.sleep(sleep_timer)

        context = {}
        #genome = request.GET.get('genome')
        #annotation_query = request.GET.get('query')
        sub_view = AnnotationSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()
        context['searchcontext'] = 'Search results for "' + request.GET.get('annotation_query') + '"'
        context['searchresult'] = sub_response.content.decode('utf-8')
        context['time'] = time.time()-start_time
        #print(context['searchresult'])
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


class StrainListView(generic.ListView):
    model = Strain
    context_object_name = 'strainlist'
    paginate_by = 50

    def get_queryset(self):
        return Strain.objects.order_by('strain_id')


class SampleListView(generic.ListView):
    model = Sample
    context_object_name = 'samplelist'
    paginate_by = 50

    def get_queryset(self):
        return Sample.objects.order_by('sample_id')


class GenomeListView(generic.ListView):
    model = Genome
    context_object_name = 'genomelist'
    paginate_by = 50

    def get_queryset(self):
        return Genome.objects.order_by('name').select_related('strain', 'taxon')


class OperonListView(generic.ListView):
    model = Operon
    template_name = 'operon_list.html'
    context_object_name = 'operonlist'
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        genome = Genome.objects.get(name = self.kwargs['genome'])
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
            ).order_by('name').select_related('genome', 'contig').prefetch_related('genes')
        else:
            object_list = Operon.objects.filter(genome__name=genome).order_by('name').select_related('genome', 'contig').prefetch_related('genes')
        return object_list

        
class SiteListView(generic.ListView):
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
            ).order_by('name').select_related('genome', 'contig')
        else:
            object_list = Site.objects.filter(genome__name=genome).order_by('name').select_related('genome', 'contig')
        return object_list


class RegulonListView(generic.ListView):
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
            ).order_by('name').select_related('genome')
        else:
            object_list = Regulon.objects.filter(genome__name=genome).order_by('name').select_related('genome')
        return object_list


class GeneListView(generic.ListView):
    model = Gene
    template_name = 'gene_list.html'
    context_object_name = 'genelist'
    paginate_by = 50

    def get_queryset(self):
        #return Gene.objects.order_by('locus_tag').select_related('genome__strain')
        #return Gene.objects.order_by('locus_tag').values('id', 'locus_tag', 'name', 'function', 'genome__strain__full_name')
        return Gene.objects.none()


class GeneSearchResultsSubView(generic.ListView):
    context_object_name = 'genelist'
    template_name = 'browser/gene_list_subpage.html'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(GeneSearchResultsSubView,self).get_context_data(**kwargs)
        if self.request.GET.get('genome'):
            searchcontext, external = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'), self.request.GET.get('genome'))
        else:
            searchcontext, external = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'))
        context['searchcontext'] = searchcontext
        if external != '':
            context['external'] = external
        return context

    def get_queryset(self):
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        query_type = self.request.GET.get('type')
        if query_type == 'gene':
            if genome:
                object_list = Gene.objects.filter(genome__name=genome).filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'og':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')

        elif query_type == 'ko_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'kp_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'kr_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'ec_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'tc_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'cazy_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'cog_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif query_type == 'go_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=query).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif genome:
            if query_type == 'ko':
                ko_ids = Kegg_ortholog.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'kp':
                kp_ids = Kegg_pathway.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon', 'protein').prefetch_related('protein__kegg_orthologs')
                self.kegg_map_url = None
                if kp_ids.count() == 1:
                    print(kp_ids[0]['kegg_id'])
                    kegg_map_url = 'https://www.kegg.jp/pathway/' + kp_ids[0]['kegg_id'] + '+'
                    print(kegg_map_url)
                    ko_ids = set()
                    for gene in object_list:
                        for ko in gene.protein.kegg_orthologs.all():
                            ko_ids.add(ko.kegg_id)
                    if ko_ids:
                        self.kegg_map_url = kegg_map_url + '+'.join(sorted(list(ko_ids)))
            elif query_type == 'kr':
                kr_ids = Kegg_reaction.objects.filter(
                    Q(kegg_id__icontains=query) | Q(description__icontains=query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'ec':
                ec_ids = Ec_number.objects.filter(
                    Q(ec_number__icontains=query) | Q(description__icontains=query)
                ).values('ec_number')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'tc':
                tc_ids = Tc_family.objects.filter(
                    Q(tc_id__icontains=query) | Q(description__icontains=query)
                ).values('tc_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'cazy':
                cazy_ids = Cazy_family.objects.filter(
                    Q(cazy_id__icontains=query) | Q(description__icontains=query)
                ).values('cazy_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'cog':
                cog_ids = Cog_class.objects.filter(
                    Q(cog_id__icontains=query) | Q(description__icontains=query)
                ).values('cog_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif query_type == 'go':
                go_ids = Go_term.objects.filter(
                    Q(go_id__icontains=query) | Q(description__icontains=query)
                ).values('go_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag').select_related('genome', 'genome__taxon')
        else:
            object_list = Gene.objects.none()
        return object_list


class GeneSearchResultsView(generic.ListView):
    model = Gene
    context_object_name = 'genelist'
    paginate_by = 50

    def get_context_data(self,**kwargs):
        context = super(GeneSearchResultsView,self).get_context_data(**kwargs)
        if self.request.GET.get('query'):
            context['searchcontext'] = 'Search results for "' + self.request.GET.get('query') + '"'
        elif self.request.GET.get('og'):
            eggnog_og = Ortholog_group.objects.get(id=self.request.GET.get('og'))
            context['searchcontext'] = 'Genes from Ortholog Group ' + eggnog_og.eggnog_id + '[' + eggnog_og.taxon.name + ']'
        elif self.request.GET.get('ko'):
            context['searchcontext'] = 'Genes from KEGG Ortholog Group ' + self.request.GET.get('ko')
            context['external'] = 'https://www.kegg.jp/dbget-bin/www_bget?' + self.request.GET.get('ko')
        elif self.request.GET.get('kp'):
            context['searchcontext'] = 'Genes from KEGG pathway ' + self.request.GET.get('kp')
            context['external'] = 'https://www.kegg.jp/dbget-bin/www_bget?' + self.request.GET.get('kp')
        elif self.request.GET.get('kr'):
            context['searchcontext'] = 'Genes from KEGG reaction ' + self.request.GET.get('kr')
            context['external'] = 'https://www.kegg.jp/dbget-bin/www_bget?' + self.request.GET.get('kr')
        elif self.request.GET.get('ec'):
            context['searchcontext'] = 'Genes with EC number ' + self.request.GET.get('ec')
            context['external'] = 'https://www.kegg.jp/dbget-bin/www_bget?ec:' + self.request.GET.get('ec')
        elif self.request.GET.get('tc'):
            #tc_id = Tc_family.objects.get(id=self.request.GET.get('tc')).tc_id
            context['searchcontext'] = 'Genes from TCDB family ' + self.request.GET.get('tc')
            context['external'] = 'http://www.tcdb.org/search/result.php?tc=' + self.request.GET.get('tc') + '#' + self.request.GET.get('tc')
        elif self.request.GET.get('cazy'):
            context['searchcontext'] = 'Genes from CAZy family ' + self.request.GET.get('cazy')
            context['external'] = 'http://www.cazy.org/' + self.request.GET.get('cazy') + '.html'
        elif self.request.GET.get('cog'):
            context['searchcontext'] = 'Genes from COG class ' + self.request.GET.get('cog')
            context['external'] = 'https://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/list' + self.request.GET.get('cog') + '.html'
        elif self.request.GET.get('go'):
            context['searchcontext'] = 'Genes assigned to GO term ' + self.request.GET.get('go')
            context['external'] = 'https://www.ebi.ac.uk/QuickGO/search/' + self.request.GET.get('go')
        elif self.request.GET.get('genome'):
            print(self.request.GET)
            genome = self.request.GET.get('genome')
            if self.request.GET.get('ko_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to KEGG Ortholog groups containing "' + self.request.GET.get('ko_query') + '"'
            elif self.request.GET.get('kp_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to KEGG pathways containing "' + self.request.GET.get('kp_query') + '"'
                if self.kegg_map_url is not None:
                    context['external'] = self.kegg_map_url
            elif self.request.GET.get('kr_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to KEGG reactions containing "' + self.request.GET.get('kr_query') + '"'
            elif self.request.GET.get('ec_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to EC numbers containing "' + self.request.GET.get('ec_query') + '"'
            elif self.request.GET.get('tc_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to TCDB families containing "' + self.request.GET.get('tc_query') + '"'
            elif self.request.GET.get('cazy_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to CAZy families containing "' + self.request.GET.get('cazy_query') + '"'
            elif self.request.GET.get('cog_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to COG class containing "' + str(self.request.GET.get('cog_query')) + '"'
            elif self.request.GET.get('go_query'):
                context['searchcontext'] = 'Genes from genome ' + genome + ' assigned to GO term containing "' + str(self.request.GET.get('go_query')) + '"'
            else:
                context['searchcontext'] = 'Genes from genome ' + genome
        return context

    def get_queryset(self):
        query = self.request.GET.get('query')
        genome = self.request.GET.get('genome')
        og = self.request.GET.get('og')
        ko = self.request.GET.get('ko')
        kp = self.request.GET.get('kp')
        kr = self.request.GET.get('kr')
        ec = self.request.GET.get('ec')
        tc = self.request.GET.get('tc')
        cazy = self.request.GET.get('cazy')
        cog = self.request.GET.get('cog')
        go = self.request.GET.get('go')
        og_query = self.request.GET.get('og_query')
        ko_query = self.request.GET.get('ko_query')
        kp_query = self.request.GET.get('kp_query')
        kr_query = self.request.GET.get('kr_query')
        ec_query = self.request.GET.get('ec_query')
        tc_query = self.request.GET.get('tc_query')
        cazy_query = self.request.GET.get('cazy_query')
        cog_query = self.request.GET.get('cog_query')
        go_query = self.request.GET.get('go_query')
        if query:
            if genome:
                object_list = Gene.objects.filter(genome__name=genome).filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) | Q(locus_tag__icontains=query) | Q(function__icontains=query)
                ).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif og:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=og).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif ko:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=ko).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif kp:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=kp).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif kr:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=kr).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif ec:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=ec).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif tc:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=tc).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif cazy:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=cazy).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif cog:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=cog).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif go:
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=go).values('protein_hash')]
            if genome:
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
        elif genome:
            if ko_query:
                ko_ids = Kegg_ortholog.objects.filter(
                    Q(kegg_id__icontains=ko_query) | Q(description__icontains=ko_query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif kp_query:
                kp_ids = Kegg_pathway.objects.filter(
                    Q(kegg_id__icontains=kp_query) | Q(description__icontains=kp_query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon', 'protein').prefetch_related('protein__kegg_orthologs')
                self.kegg_map_url = None
                if kp_ids.count() == 1:
                    print(kp_ids[0]['kegg_id'])
                    kegg_map_url = 'https://www.kegg.jp/pathway/' + kp_ids[0]['kegg_id'] + '+'
                    print(kegg_map_url)
                    ko_ids = set()
                    for gene in object_list:
                        for ko in gene.protein.kegg_orthologs.all():
                            ko_ids.add(ko.kegg_id)
                    if ko_ids:
                        self.kegg_map_url = kegg_map_url + '+'.join(sorted(list(ko_ids)))
            elif kr_query:
                kr_ids = Kegg_reaction.objects.filter(
                    Q(kegg_id__icontains=kr_query) | Q(description__icontains=kr_query)
                ).values('kegg_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif ec_query:
                ec_ids = Ec_number.objects.filter(
                    Q(ec_number__icontains=ec_query) | Q(description__icontains=ec_query)
                ).values('ec_number')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif tc_query:
                tc_ids = Tc_family.objects.filter(
                    Q(tc_id__icontains=tc_query) | Q(description__icontains=tc_query)
                ).values('tc_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif cazy_query:
                cazy_ids = Cazy_family.objects.filter(
                    Q(cazy_id__icontains=cazy_query) | Q(description__icontains=cazy_query)
                ).values('cazy_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif cog_query:
                cog_ids = Cog_class.objects.filter(
                    Q(cog_id__icontains=cog_query) | Q(description__icontains=cog_query)
                ).values('cog_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            elif go_query:
                go_ids = Go_term.objects.filter(
                    Q(go_id__icontains=go_query) | Q(description__icontains=go_query)
                ).values('go_id')
                proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
                object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag').select_related('genome', 'genome__taxon')
            else:
                object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag').select_related('genome', 'genome__taxon')
        else:
            object_list = Gene.objects.none()
        return object_list


class GeneSearchResultsAjaxView(View):
    def get(self,request):
        context = {}
        if self.request.GET.get('genome'):
            searchcontext, external = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'), self.request.GET.get('genome'))
        else:
            searchcontext, external = generate_gene_search_context(self.request.GET.get('query'), self.request.GET.get('type'))
        context['searchcontext'] = searchcontext
        if external != '':
            context['external'] = external

        print('REQUEST1')
        print('Request parameters:', request.GET.get('annotation_query'), request.GET.get('genome'))
        for key, val in request.GET.items():
            print(key, val)
            context[key] = val
        return render(request,'browser/gene_list_ajax.html', context)

    @staticmethod
    def ajax_view(request):
        start_time = time.time()

        print('REQUEST2')
        for key, val in request.GET.items():
            print(key, val)
            
        # Sleep timer for testing to imitate long-running task
        # sleep_timer = 0
        # print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        # time.sleep(sleep_timer)

        context = {}
        #genome = request.GET.get('genome')
        #annotation_query = request.GET.get('query')
        sub_view = GeneSearchResultsSubView()
        sub_view.setup(request)
        sub_response = sub_view.dispatch(request)
        sub_response.render()
        if request.GET.get('genome'):
            context['searchcontext'] = generate_gene_search_context(request.GET.get('query'), request.GET.get('type'), request.GET.get('genome'))
        else:
            context['searchcontext'] = generate_gene_search_context(request.GET.get('query'), request.GET.get('type'))
        context['searchresult'] = sub_response.content.decode('utf-8')
        context['time'] = time.time()-start_time
        #print(context['searchresult'])
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


class GenomeSearchResultsView(generic.ListView):
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
            object_list = Genome.objects.filter(name__icontains=query).order_by('name').select_related('strain', 'sample', 'taxon')
        else:
            object_list = Genome.objects.none()
        return object_list


class StrainSearchResultsView(generic.ListView):
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
    template = loader.get_template('browser/textsearch.html')
    context={}
    return HttpResponse(template.render(context, request))


def startpage(request):
    template = loader.get_template('browser/landing.html')
    num_genomes = Genome.objects.all().count()
    context = {'num_genomes':num_genomes}
    return HttpResponse(template.render(context, request))


def show_help(request):
    template = loader.get_template('browser/help.html')
    context = {}
    return HttpResponse(template.render(context, request))


def genome_detail(request, name):
    try:
        genome = Genome.objects.get(name = name)
    except Genome.DoesNotExist:
        raise Http404('Genome not found: ' + name)
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
    try:
        strain = Strain.objects.get(strain_id = strain_id)
        genomes = Genome.objects.filter(strain=strain).order_by('name')
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
    try:
        sample = Sample.objects.get(sample_id = sample_id)
        genomes = Genome.objects.filter(sample=sample).order_by('name')
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
    try:
        gene = Gene.objects.select_related(
            'genome', 'genome__taxon', 'protein', 'operon'
            ).prefetch_related(
            'protein__ortholog_groups__taxon', 'protein__kegg_orthologs', 'protein__kegg_pathways', 'protein__kegg_reactions', 'protein__ec_numbers', 'protein__go_terms', 'protein__tc_families', 'protein__cog_classes'
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
    try:
        operon = Operon.objects.select_related('contig', 'genome').get(name=name)
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
    try:
        site = Site.objects.select_related(
            'genome', 'genome__strain', 'contig'
            ).get(genome__name=genome, name = name)
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
    try:
        regulon = Regulon.objects.select_related(
            'genome', 'genome__strain'
            ).prefetch_related('regulators'
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
    if request.GET.get('genome') and request.GET.get('locus'):
        genome_name = request.GET.get('genome')
        locus_tag = request.GET.get('locus')
    else:
        raise Http404('Genome name or locus ID not provided')
    try:
        gene = Gene.objects.get(locus_tag = locus_tag, genome__name = genome_name)
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


def genomes_index(request):
    template = loader.get_template('browser/genomes.html')
    strains = Strain.objects.order_by('strain_id')
    genomes = Genome.objects.order_by('name')
    context = {
        'strains': strains,
        'genomes': genomes,
    }
    return HttpResponse(template.render(context, request))


def protein_search_external(request):
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
            genes = Gene.objects.select_related('protein', 'genome', 'genome__taxon').filter(protein__protein_hash = row[1])
            for target in genes:
                hit = [target.genome.name, target.locus_tag, target.genome.taxon.name,
                       target.function, '{:.1f}'.format(float(row[2])), row[3], '{:.1f}'.format(query_cov), row[10], row[11]]
                result[row[0]].append(hit)
        context['searchresult'] = result
    else:
        raise Http404('Provide protein sequence in FASTA format')
    return render(request, 'browser/proteinsearch.html', context)


class NsearchResultView(View):
    def post(self,request):
        sequence = request.POST.get("sequence")
        print('Sequence sent:', sequence)
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        print('REQUEST1')
        for key, val in request.POST.items():
            context[key] = val
            print(key, val)
        return render(request,'browser/nucleotidesearchajax.html', context)

    def get(self,request):
        return render(request,'browser/nucleotidesearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        start_time = time.time()
        result = []
        print('REQUEST2')
        for key, val in request.POST.items():
            print(key, val)
        #sequence = request.POST.get("sequence")
        params = {'sequence': request.POST.get("sequence"), 'evalue': request.POST.get("evalue"), 'hitstoshow': request.POST.get("hitstoshow")}
        #params['sequence'] = '>test\nATGACGATGCACCCGGCAGCACCTCGAACTCCGCACTGGCGCTGCTTGCACCGGCTGGCATGGAGCCTGTCCTGCGCTGCCCTGTTGCCGCTCGCTGCACTGGCGCAGGACGTACCGTCCCGCGCCGTCACGCCGGTGTCCGCAGCGTCGCCGCCGCGCCAGTCGCAGGACGACGCCTGGTGGACCGGGCCGATGCTGGCGAACTCCGCCGCCACCCTGCCGCGCGGCCACGTCCTGATCGAGCCTTACGTCTACGACGTGTCCTCGCCGCACGCCGACGGCTACGGTTCGCTCACCTACATGCTCTACGGCCTCACCGACCGGCTGACGGTCGGCCTGATGCCGGTGCTGGGCTACAACCGCATGGATGGCCCGGGCGACAGCAGCGGGATCGGGCTGGGCGACGTCAGCGTGCAGGCGCAGTACCGGCTGACCGATGTGCCGGCGGGCAGTTCGCGGCCCACGGTCTCGCTGCAACTGCAGGAAACCCTGCCGACCGGCAAGTACGACCGGCTGGGCCGGCGACCCGGCAACGGCCTGGGCAGCGGCGCCACCACCACTACGCTGCAGGTCAACACGCAGACGTATTTCTGGTTGTCCAACGGCCGCATCCTGCGCATGCGCTTCAACGTGGCGCAATCATTCTCGACGCGGGCACGGGTCGAGGACATCAGCGTCTACGGCACCCCGGACGGCTTTCGCGGGCACGCCCGGCCGGGGCGTTCGTTCTTCGTCAATGCGGCCTGGGAGTACAGCCTCAGCCAGCGCTGGGTGCTGGCGCTCGACCTCACCTACCGGCGCAGCCACGGTGCCCGCGTGCGCGACGACGACCTCAATGGCGTGCCTGCCTTGCGTCTGGACGGCCGCTCCAGCGAGGCGTTCGGCTTTGCCCCGGCCATCGAGTACAGCTGGAGTCCGCGGCTCGGCGTGCTGTTCGGCACCCGCGTGATCACCGGCGGGCACAACACCGCGACCACGATCACGCCGGCGGTGGCCTTCAACTACGTGCACTGA'
        print('Sequence received:', params['sequence'])
        sleep_timer = 0
        print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        time.sleep(sleep_timer)
        hits, searchcontext, query_len = run_nucleotide_search(params)
        if hits:
            result.append('<table><thead><tr><th>Target contig (click to see hit)</th><th>Genome</th><th>%identity</th><th>Alignment length</th><th>%Query coverage</th><th>E-value</th><th>Bit-score</th></tr></thead><tbody>')
            for row in hits:
                row=row.split('\t')
#                if row[0] not in result:
#                    result[row[0]] = []
                contig_name, genome_name = row[1].split('|')
                genome_name = genome_name.split('[')[0]
                print('Search for contig', contig_name, 'in genome', genome_name)
                target = Contig.objects.select_related('genome', 'genome__taxon').get(contig_id = contig_name, genome__name = genome_name)
                ani = float(row[2])
                strand = '+'
                start = row[8]
                end = row[9]
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                if int(row[7]) < int(row[6]):
                    strand = '-'
#                hit = [contig_name + ': (' + strand + 'strand) ' + start + '..' + end,
#                        target, '{:.2f}'.format(float(row[2])), row[3], '{:.1f}'.format(query_cov), row[10], row[11], start, end]
                hit = '<tr><td align=\"left\"><a href=\"' + reverse('genomedetails', args=(target.genome.name,)) + '?contig=' + target.contig_id + '&start=' + start + '&end=' + end + '\">' + contig_name + ': (' + strand + 'strand) ' + start + '..' + end + '</a></td><td align="left">' + target.genome.name + ' [' + target.genome.taxon.name + ']</td><td>' + '{:.2f}'.format(float(row[2])) + '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) + '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                result.append(hit)
            result.append('</tbody></table>')
        elif searchcontext == '':
            result = ['<h5>No hits found.</h5>']
        else:
            context = {"searchcontext":searchcontext}
            return render(request,'404.html', context)
        context = {"searchresult":'\n'.join(result),"searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        print(context)
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


def nucleotidesearchform(request):
    return render(request,'browser/nucleotidesearchform.html')


class PsearchResultView(View):

    def post(self,request):
        sequence = request.POST.get("sequence")
        print('Sequence sent:', sequence)
        context = {'csrfmiddlewaretoken': request.POST.get('csrfmiddlewaretoken')}
        print('REQUEST1')
        for key, val in request.POST.items():
            context[key] = val
            print(key, val)
        return render(request,'browser/proteinsearchajax.html', context)

    def get(self,request):
        return render(request,'browser/proteinsearchajax.html')
    
    @staticmethod
    def ajax_view(request):
        start_time = time.time()
        result = []
        print('REQUEST2')
        for key, val in request.POST.items():
            print(key, val)
        sequence = request.POST.get("sequence")
        params = {'sequence': request.POST.get("sequence"), 'evalue': request.POST.get("evalue"), 'hitstoshow': request.POST.get("hitstoshow")}
        #sequence = 'MTKQVQEAYIVAATRTPVGKAPRGVFRNTRPDDMLAHVIRAVMAQAPGIDPHQIGDVIIGCAMPEAEQGMNVARIGLLLAGLPDTVPGVTVNRFCSSGLQSVAMAADRIRLGLDDLMLAGGTESMSMVPMMGHKIAMNPAIFNDENIGIAYGMGITAENVAKQWKVSREQQDAFSVESHRRALAAQAAGEFNDEISPFALDDHYPNLATRGIVTDSRRIDSDEGPRAGTTMEVLAKLKTVFRNGQFGGTVTAGNSSQMSDGAGAVLLASERAVKEYNLQPLARFVGFSVAGVPPEVMGIGPKEAIPKALKQAGLNRDQLDWIELNEAFAAQALAVMGDLGLDPDKVNPLGGAIALGHPLGATGAVRIATLVHGMRRRKQKYGMVTMCIGTGMGAAGIFEAL'
        print('Sequence received:', sequence)
        sleep_timer = 0
        print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        time.sleep(sleep_timer)
        hits, searchcontext, query_len = run_protein_search(params)

        if hits:
            result.append('<table><thead><tr><th>Target gene</th><th>Genome</th><th>%identity</th><th>Alignment length</th><th>%Query coverage</th><th>E-value</th><th>Bit-score</th></tr></thead><tbody>')
            for row in hits:
                row=row.split('\t')
                unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
                query_cov = (query_len - unaligned_part) * 100.0 / query_len
                genes = Gene.objects.select_related('protein', 'genome', 'genome__taxon').filter(protein__protein_hash = row[1])
                for target in genes:
                    hit = '<tr><td align=\"left\"><a href=\"' + reverse('genedetails', args=(target.genome.name, target.locus_tag)) + '\">' + target.locus_tag + '</a></td><td align="left">' + target.genome.name + ' [' + target.genome.taxon.name + ']</td><td>' + '{:.1f}'.format(float(row[2])) + '</td><td>' + row[3] + '</td><td>' + '{:.1f}'.format(query_cov) + '</td><td>' + row[10] + '</td><td>' + row[11] + '</td></tr>'
                    result.append(hit)
            result.append('</tbody></table>')
        elif searchcontext == '':
            result = ['<h5>No hits found.</h5>']
        else:
            context = {"searchcontext":searchcontext}
            return render(request,'404.html', context)
        context = {"searchresult":'\n'.join(result),"searchcontext":searchcontext,"query_len":query_len,"time":time.time()-start_time}
        print(context)
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")


def proteinsearchform(request):
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

    @staticmethod
    def verify_parameters(size, lines):
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
        locus_tag = request.GET.get('locus_tag')
        genome = request.GET.get('genome')
        og_id = request.GET.get('og')

        try:
            ComparativeView.verify_parameters(request.GET.get('size'), request.GET.get('lines'))
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
        
        print('REQUEST1')
        print('Request parameters:', og_id, genome, locus_tag, request.GET.get('size'), request.GET.get('lines'))
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

        print('REQUEST2')
        for key, val in request.GET.items():
            print(key, val)
        try:
            ComparativeView.verify_parameters(request.GET.get('size'), request.GET.get('lines'))
        except SuspiciousOperation as e:
            return render(request, '404.html', {'searchcontext': str(e)})
            
        # Sleep timer for testing to imitate long-running task
        sleep_timer = 0
        print('DELAY FOR ' + str(sleep_timer) + ' SECONDS')
        time.sleep(sleep_timer)

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
            #context = {'gene': gene, 'ortholog_group':og, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count, 'scribl_message':''}
            scribl='<script type="text/javascript">\nalert("Comparative plot cannot be created for only one genome");</script>'
            context = {'scribl':scribl, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count, "time":time.time()-start_time}
        else:
            #context = {'gene': gene, 'ortholog_group':og, 'scribl':scribl, 'tree_canvas':tree_canvas, 'tree_newick':tree_newick, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count}
            scribl='<script type="text/javascript">\nfunction draw(canvasName) {\nvar canvas = document.getElementById(canvasName);\nvar parent = document.getElementById("scribl-container");\nvar ctx = canvas.getContext("2d");\n' + scribl + '\nchart.draw();\nvar img = canvas.toDataURL("image/png");\ndocument.getElementById("pngexport").href = img;\n}\n</script>'
            plot_gene_count -= 1
            context = {'scribl':scribl, 'tree_canvas':tree_canvas, 'tree_newick':tree_newick, 'og_gene_count':og_gene_count, 'plot_gene_count':plot_gene_count,"time":time.time()-start_time}
        # print(context)
        data = json.dumps(context)
        return HttpResponse(data,content_type="application/json")

    
def export_csv(request):
    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="export.tab"'

    writer = csv.writer(response, delimiter='\t')
    search_context = ('','')
    query = request.GET.get('query')
    genome = request.GET.get('genome')
    og = request.GET.get('og')
    ko = request.GET.get('ko')
    kp = request.GET.get('kp')
    kr = request.GET.get('kr')
    ec = request.GET.get('ec')
    tc = request.GET.get('tc')
    cazy = request.GET.get('cazy')
    cog = request.GET.get('cog')
    go = request.GET.get('go')
    og_query = request.GET.get('og_query')
    ko_query = request.GET.get('ko_query')
    kp_query = request.GET.get('kp_query')
    kr_query = request.GET.get('kr_query')
    ec_query = request.GET.get('ec_query')
    tc_query = request.GET.get('tc_query')
    cazy_query = request.GET.get('cazy_query')
    cog_query = request.GET.get('cog_query')
    go_query = request.GET.get('go_query')
    annotation_query = request.GET.get('annotation_query')
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
    elif og:
        search_context = ('Ortholog group', og)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=og).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif ko:
        search_context = ('KEGG ortholog', ko)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=ko).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif kp:
        search_context = ('KEGG pathway', kp)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=kp).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif kr:
        search_context = ('KEGG reaction', kr)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=kr).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif ec:
        search_context = ('EC number', ec)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=ec).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif tc:
        search_context = ('TCDB family', tc)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=tc).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif cazy:
        search_context = ('CAZy family', cazy)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=cazy).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif cog:
        search_context = ('COG class', cog)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=cog).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif go:
        search_context = ('GO term', go)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=go).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif annotation_query:
        search_context = ('Gene annotation query', annotation_query)
        if genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon')
    elif genome:
        if ko_query:
            search_context = ('KEGG ortholog query', ko_query)
            ko_ids = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=ko_query) | Q(description__icontains=ko_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif kp_query:
            search_context = ('KEGG pathway query', kp_query)
            kp_ids = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=kp_query) | Q(description__icontains=kp_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif kr_query:
            search_context = ('KEGG reaction query', kr_query)
            kr_ids = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=kr_query) | Q(description__icontains=kr_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif ec_query:
            search_context = ('EC number query', ec_query)
            ec_ids = Ec_number.objects.filter(
                Q(ec_number__icontains=ec_query) | Q(description__icontains=ec_query)
            ).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif tc_query:
            search_context = ('TCDB family query', tc_query)
            tc_ids = Tc_family.objects.filter(
                Q(tc_id__icontains=tc_query) | Q(description__icontains=tc_query)
            ).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif cazy_query:
            search_context = ('CAZy family query', cazy_query)
            cazy_ids = Cazy_family.objects.filter(
                Q(cazy_id__icontains=cazy_query) | Q(description__icontains=cazy_query)
            ).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif cog_query:
            search_context = ('COG class query', cog_query)
            cog_ids = Cog_class.objects.filter(
                Q(cog_id__icontains=cog_query) | Q(description__icontains=cog_query)
            ).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif go_query:
            search_context = ('GO term query', go_query)
            go_ids = Go_term.objects.filter(
                Q(go_id__icontains=go_query) | Q(description__icontains=go_query)
            ).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
    else:
        object_list = Gene.objects.none()

    if annotation_query:
        writer.writerow(['Locus tag', 'Name', 'Organism', 'Genome', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Function', 'Annotation_source', 'Annotation_type', 'Annotation_note'])
        for obj in object_list:
            writer.writerow([obj.gene_id.locus_tag, obj.gene_id.name, obj.gene_id.genome.taxon.name, obj.gene_id.genome.name, obj.gene_id.contig.contig_id, str(obj.gene_id.start), str(obj.gene_id.end), str(obj.gene_id.strand), obj.gene_id.type, obj.gene_id.function, obj.source, obj.key, obj.value, obj.note])
    else:
        writer.writerow(['Locus tag', 'Name', 'Organism', 'Genome', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Function', search_context[0]])
        for gene in object_list:
            writer.writerow([gene.locus_tag, gene.name, gene.genome.taxon.name, gene.genome.name, gene.contig.contig_id, str(gene.start), str(gene.end), str(gene.strand), gene.type, gene.function, search_context[1]])

    return response

def export_fasta(request):
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="export.faa"'

    search_context = ('','')
    query = request.GET.get('query')
    genome = request.GET.get('genome')
    og = request.GET.get('og')
    ko = request.GET.get('ko')
    kp = request.GET.get('kp')
    kr = request.GET.get('kr')
    ec = request.GET.get('ec')
    tc = request.GET.get('tc')
    cazy = request.GET.get('cazy')
    cog = request.GET.get('cog')
    go = request.GET.get('go')
    og_query = request.GET.get('og_query')
    ko_query = request.GET.get('ko_query')
    kp_query = request.GET.get('kp_query')
    kr_query = request.GET.get('kr_query')
    ec_query = request.GET.get('ec_query')
    tc_query = request.GET.get('tc_query')
    cazy_query = request.GET.get('cazy_query')
    cog_query = request.GET.get('cog_query')
    go_query = request.GET.get('go_query')
    annotation_query = request.GET.get('annotation_query')
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
    elif og:
        search_context = ('Ortholog group', og)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(ortholog_groups__id=og).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif ko:
        search_context = ('KEGG ortholog', ko)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id=ko).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif kp:
        search_context = ('KEGG pathway', kp)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id=kp).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif kr:
        search_context = ('KEGG reaction', kr)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id=kr).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif ec:
        search_context = ('EC number', ec)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number=ec).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif tc:
        search_context = ('TCDB family', tc)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id=tc).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif cazy:
        search_context = ('CAZy family', cazy)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id=cazy).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif cog:
        search_context = ('COG class', cog)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id=cog).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif go:
        search_context = ('GO term', go)
        proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id=go).values('protein_hash')]
        if genome:
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(protein__protein_hash__in=proteins).order_by('locus_tag')
    elif annotation_query:
        search_context = ('Gene annotation query', annotation_query)
        if genome:
            object_list = Annotation.objects.filter(
                (Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)) & Q(gene_id__genome__name=genome)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon', 'gene_id__protein')
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) | Q(value__icontains=annotation_query) | Q(note__icontains=annotation_query)
            ).order_by('gene_id__locus_tag').select_related('gene_id', 'gene_id__genome', 'gene_id__genome__taxon', 'gene_id__protein')
    elif genome:
        if ko_query:
            search_context = ('KEGG ortholog query', ko_query)
            ko_ids = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=ko_query) | Q(description__icontains=ko_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif kp_query:
            search_context = ('KEGG pathway query', kp_query)
            kp_ids = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=kp_query) | Q(description__icontains=kp_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif kr_query:
            search_context = ('KEGG reaction query', kr_query)
            kr_ids = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=kr_query) | Q(description__icontains=kr_query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif ec_query:
            search_context = ('EC number query', ec_query)
            ec_ids = Ec_number.objects.filter(
                Q(ec_number__icontains=ec_query) | Q(description__icontains=ec_query)
            ).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif tc_query:
            search_context = ('TCDB family query', tc_query)
            tc_ids = Tc_family.objects.filter(
                Q(tc_id__icontains=tc_query) | Q(description__icontains=tc_query)
            ).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(tc_families__tc_id__in=tc_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif cazy_query:
            search_context = ('CAZy family query', cazy_query)
            cazy_ids = Cazy_family.objects.filter(
                Q(cazy_id__icontains=cazy_query) | Q(description__icontains=cazy_query)
            ).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cazy_families__cazy_id__in=cazy_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif cog_query:
            search_context = ('COG class query', cog_query)
            cog_ids = Cog_class.objects.filter(
                Q(cog_id__icontains=cog_query) | Q(description__icontains=cog_query)
            ).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(cog_classes__cog_id__in=cog_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        elif go_query:
            search_context = ('GO term query', go_query)
            go_ids = Go_term.objects.filter(
                Q(go_id__icontains=go_query) | Q(description__icontains=go_query)
            ).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(go_terms__go_id__in=go_ids).values('protein_hash')]
            object_list = Gene.objects.filter(genome__name=genome, protein__protein_hash__in=proteins).order_by('locus_tag')
        else:
            object_list = Gene.objects.filter(genome__name=genome).order_by('locus_tag')
    else:
        object_list = Gene.objects.none()

    if annotation_query:
        for obj in object_list:
            if obj.gene_id.protein:
                response.write('>' + obj.gene_id.locus_tag + '|' + obj.gene_id.genome.name + ' [' + obj.gene_id.genome.taxon.name + ']\n' + obj.gene_id.protein.sequence + '\n')
    else:
        for gene in object_list:
            response.write('>' + gene.locus_tag + '|' + gene.genome.name + ' [' + gene.genome.taxon.name + ']\n' + gene.protein.sequence + '\n')

    return response
    
def handler404(request, exception):
    return render(request, '404.html', status=404)
    
def generate_gene_search_context(query, query_type, genome=None):
    searchcontext = ''
    external = ''
    if genome is None:
        if query_type=='gene':
            searchcontext = 'Search results for "' + query + '"'
        elif query_type=='og':
            eggnog_og = Ortholog_group.objects.get(id=query)
            searchcontext = 'Genes from Ortholog Group ' + eggnog_og.eggnog_id + '[' + eggnog_og.taxon.name + ']'
        elif query_type=='ko_id':
            searchcontext = 'Genes from KEGG Ortholog Group ' + query
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='kp_id':
            searchcontext = 'Genes from KEGG pathway ' + query
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='kr_id':
            searchcontext = 'Genes from KEGG reaction ' + query
            external = 'https://www.kegg.jp/dbget-bin/www_bget?' + query
        elif query_type=='ec_id':
            searchcontext = 'Genes with EC number ' + query
            external = 'https://www.kegg.jp/dbget-bin/www_bget?ec:' + query
        elif query_type=='tc_id':
            searchcontext = 'Genes from TCDB family ' + query
            external = 'http://www.tcdb.org/search/result.php?tc=' + query + '#' + query
        elif query_type=='cazy_id':
            searchcontext = 'Genes from CAZy family ' + query
            external = 'http://www.cazy.org/' + query + '.html'
        elif query_type=='cog_id':
            searchcontext = 'Genes from COG class ' + query
            external = 'https://ftp.ncbi.nih.gov/pub/COG/COG2014/static/lists/list' + query + '.html'
        elif query_type=='go_id':
            searchcontext = 'Genes assigned to GO term ' + query
            external = 'https://www.ebi.ac.uk/QuickGO/search/' + query
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
            searchcontext = 'Genes from genome ' + genome + ' assigned to GO term containing "' + query + '"'
        elif query_type=='gene':
            searchcontext = 'Genes from genome ' + genome
    return searchcontext, external
