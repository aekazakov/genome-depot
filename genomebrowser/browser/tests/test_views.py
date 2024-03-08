import os
import hashlib
import gzip
from unittest import skip
from django.test import TestCase
from django.test import Client
from django.test import RequestFactory
from django.utils import timezone
from browser.models import Operon
from browser.models import Site
from browser.views import StrainListView
from browser.views import SampleListView
from browser.views import GenomeListView
from browser.views import AnnotationSearchResultsSubView
from browser.views import RegulonListView
from browser.views import OperonListView
from browser.views import SiteListView
from genomebrowser.settings import BASE_DIR
from browser.comparative_analysis import _get_color, make_muscle_alignment
from browser.seqsearch import _sanitize_sequence

from browser.util import recreate_search_databases
# Create your tests here.

class BrowserViewsTest(TestCase):
    '''
        Views testing
    '''
    fixtures = ['minigenomes.testdata.json']
    
    @classmethod
    def setUpTestData(cls):
        cls.client = Client()
        cls.csrf_client = Client(enforce_csrf_checks=True)
        
    def test_start_page(self):
        '''
            Testing Start page view
            path('', views.startpage, name='index')
        '''
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<div class="logo">GenomeDepot</div>')

    def test_genomes_page(self):
        '''
            Testing Genomes list page view
            path('genomes/', views.GenomeListView.as_view(), name='genome_list')
        '''
        #response = self.client.get('/genomes/')
        request = RequestFactory().get('/genomes/')
        view = GenomeListView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Genomes</h2>')
        self.assertContains(response, 'E_coli_BW2952')
        self.assertContains(response, 'window.PLOTLYENV=window.PLOTLYENV')

    def test_strains_page(self):
        '''
            Testing Strains list page view
            path('strains/', views.StrainListView.as_view(), name='strain_list')
        '''
        request = RequestFactory().get('/strains/')
        view = StrainListView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Strains</h2>')
        self.assertContains(response, 'BW2952')

    def test_samples_page(self):
        '''
            Testing Samples list page view
            path('samples/', views.SampleListView.as_view(), name='sample_list')
        '''
        #response = self.client.get('/samples/')
        request = RequestFactory().get('/samples/')
        view = SampleListView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Samples</h2>')
        self.assertContains(response, 'test_sample')

    def test_taxa_page(self):
        '''
            Testing Taxon list page view
            path('taxa/', views.TaxonListView.as_view(), name='taxa_list')
        '''
        response = self.client.get('/taxa/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Taxonomy')
        self.assertContains(response, 'Escherichia coli')
        
    def test_taxon_page(self):
        '''
            Testing Taxon page view
            path('taxonomy/<str:taxonomy_id>', views.taxon_detail, name='taxondetails'),
        '''
        response = self.client.get('/taxonomy/562/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'NCBI Taxonomy ID')
        self.assertContains(response, 'Escherichia coli')
        # Taxon does not exist
        response = self.client.get('/taxonomy/9999999999/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')
        
    def test_sample_page(self):
        '''
            Testing Sample page view
            path('sample/<int:sample_id>', views.sample_detail, name='sampledetails')
        '''
        response = self.client.get('/sample/1/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Sample</h2>')
        self.assertContains(response, 'Test sample')
        self.assertContains(response, 'This metadata record is for testing only')

        response = self.client.get('/sample/000000000000/')
        self.assertEqual(response.status_code, 200)

    def test_strain_page(self):
        '''
            Testing Strain page view
            path('strain/<int:strain_id>/', views.strain_detail, name='straindetails')
        '''
        response = self.client.get('/strain/157/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Strain</h2>')
        self.assertContains(response, 'E_coli')

        response = self.client.get('/strain/000000000000/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sample not found')

    def test_genome_page(self):
        '''
            Testing Genome page view
            path('genome/<str:name>', views.genome_detail, name='genomedetails')
        '''
        response = self.client.get('/genome/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'E_coli_BW2952')
        
        response = self.client.get('/genome/XXXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')
            
        response = self.client.get('/genome/E_coli_CFT073/', {'contig':'NC_004431', 'start':'1048', 'end':'3510'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'E_coli_CFT073')

    def test_tag_page(self):
        '''
            Testing Tag page view
            path('tag/<str:name>/', views.TagView.as_view(), name='tagdetails')
        '''
        response = self.client.get('/tag/test/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Tag')
        self.assertContains(response, 'Test genome')

    def test_gene_page(self):
        '''
            Testing Gene page view
            path('gene/<str:genome>/<str:locus_tag>/', views.gene_detail, name='genedetails')
        '''
        response = self.client.get('/gene/E_coli_BW2952/BWG_RS00020/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'BWG_RS00020')

        response = self.client.get('/gene/E_coli_BW2952/XXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_operonlist_page(self):
        '''
            Testing Operons list page view
            path('operons/<str:genome>/', views.OperonListView.as_view(), name='operonlist')
        '''
        #response = self.client.get('/operons/E_coli_BW2952/')
        # with genome only
        request = RequestFactory().get('/operons/E_coli_BW2952/')
        view = OperonListView.as_view()
        #view.setup(request)
        response = view(request, **{'genome':'E_coli_BW2952'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operons in')
        self.assertContains(response, 'NC_012759: 190..5020')
        # with genome and text query
        response = view(request, **{'genome':'E_coli_BW2952', 'query':'BWG_RS00005'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operons in')
        self.assertContains(response, 'NC_012759: 190..5020')
        # if genome does not exist
        response = view(request, **{'genome':'XXXXXXXXXXXX', 'query':'BWG_RS00005'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_sitelist_page(self):
        '''
            Testing Sites list page view
            path('sites/<str:genome>/', views.SiteListView.as_view(), name='sitelist')
        '''
        #response = self.client.get('/sites/E_coli_BW2952/')
        request = RequestFactory().get('/sites/E_coli_BW2952/')
        view = SiteListView.as_view()
        #view.setup(request)
        response = view(request, **{'genome':'E_coli_BW2952'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sites in')
        self.assertContains(response, 'NC_012759: complement(70130..70146)')
        # if genome does not exist
        request = RequestFactory().get('/sites/XXXXXXXXXXXX/')
        view = SiteListView.as_view()
        response = view(request, **{'genome':'XXXXXXXXXXXX'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_regulonlist_page(self):
        '''
            Testing Regulons list page view
            path('regulons/<str:genome>/', views.RegulonListView.as_view(), name='regulonlist')
        '''
        #response = self.client.get('/regulons/E_coli_BW2952/')
        request = RequestFactory().get('/regulons/E_coli_BW2952/')
        view = RegulonListView.as_view()
        #view.setup(request)
        response = view(request, **{'genome':'E_coli_BW2952'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulons in')
        self.assertContains(response, 'AraC')
        # if genome does not exist
        response = view(request, **{'genome':'XXXXXXXXXXXX'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_gene_search_page(self):
        '''
            Testing Genes list for a genome view with "genome" parameter
            path('searchgene/', views.GeneSearchResultsAjaxView.as_view(), name='searchgene')
        '''
        # Gene list for a genome
        response = self.client.get('/searchgene/',
                                   {'genome':'E_coli_BW2952', 'type':'gene'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')
        self.assertContains(response,
                            'data: {\'query\': "", \'type\': "gene", ' +\
                            '\'genome\': "E_coli_BW2952", \'page\': "" }'
                            )
        # search in one genome 
        response = self.client.get('/searchgene/',
                                   {'genome':'E_coli_BW2952', 'type':'gene', 'query':'BWG_RS00010'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')
        self.assertContains(response,
                            'data: {\'query\': "BWG_RS00010", \'type\': "gene", ' +\
                            '\'genome\': "E_coli_BW2952", \'page\': "" }'
                            )
        response = self.client.get('/searchgene/',
                                   {'genome':'E_coli_BW2952', 'type':'gene', 'query':'thrl'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')
        self.assertContains(response,
                            'data: {\'query\': "thrl", \'type\': "gene", ' +\
                            '\'genome\': "E_coli_BW2952", \'page\': "" }'
                            )
        response = self.client.get('/searchgene/',
                                   {'genome':'E_coli_BW2952', 'type':'gene', 'query':''}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')

    def test_gene_search_byregulator_page(self):
        '''
            Testing Genes list for a regulator view without "genome" parameter
            path('searchgene/', views.GeneSearchResultsAjaxView.as_view(), name='searchgene')
        '''
        # Gene list for "regulator" text query
        response = self.client.get('/searchgene/',
                                   {'query':'regulator', 'type':'gene'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response,
                            'data: {\'query\': "regulator", \'type\': ' +\
                            '"gene", \'genome\': "", \'page\': "" }'
                            )

    def test_loading_gene_search_bygenome_page(self):
        '''
            Testing Genes search results with "genome" parameter
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with empty query
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'gene', 'query':''}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')
        self.assertContains(response, 'bifunctional aspartate kinase')
        # with text query in one genome
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'gene', 'query':'thr'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genes from genome E_coli_BW2952')
        self.assertContains(response, 'thr operon leader peptide')
        # with text query in all genomes
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'gene', 'query':'thrl'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response, 'thr operon leader peptide')
        # with empty query string
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'gene', 'query':''}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Query string is empty')

    def test_loading_gene_search_by_og_page(self):
        '''
            Testing Genes search results for og type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'og', 'query':102482}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertNotContains(response, 'C_RS00015')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'og', 'query':102482}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

    def test_loading_gene_search_by_ko_id_page(self):
        '''
            Testing Genes search results for ko_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'ko_id', 'query':'K12524'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertNotContains(response, 'C_RS00015')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'ko_id', 'query':'K12524'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

    def test_loading_gene_search_by_kp_id_page(self):
        '''
            Testing Genes search results for kp_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'kp_id', 'query':'map00300'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertNotContains(response, 'C_RS00015')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'kp_id', 'query':'map00300'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

    def test_loading_gene_search_by_kr_id_page(self):
        '''
            Testing Genes search results for kp_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'kr_id', 'query':'R00480'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertNotContains(response, 'C_RS00015')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'kr_id', 'query':'R00480'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

    def test_loading_gene_search_by_ec_id_page(self):
        '''
            Testing Genes search results for ec_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'ec_id', 'query':'1.1.1.3'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertNotContains(response, 'C_RS00015')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'ec_id', 'query':'1.1.1.3'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

    def test_loading_gene_search_by_tc_id_page(self):
        '''
            Testing Genes search results for ec_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'tc_id', 'query':'1.A.33.1'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00070')
        self.assertNotContains(response, 'C_RS00080')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'tc_id', 'query':'1.A.33.1'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00070')
        self.assertContains(response, 'C_RS00080')

    def test_loading_gene_search_by_cazy_id_page(self):
        '''
            Testing Genes search results for cazy_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_CFT073', 'type':'cazy_id', 'query':'GH94'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'C_RS00340')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'cazy_id', 'query':'GH94'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'C_RS00340')


    def test_loading_gene_search_by_cog_id_page(self):
        '''
            Testing Genes search results for cog_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'cog_id', 'query':'C'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00200')
        self.assertNotContains(response, 'C_RS00220')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'cog_id', 'query':'C'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00200')
        self.assertContains(response, 'C_RS00220')

    def test_loading_gene_search_by_go_id_page(self):
        '''
            Testing Genes search results for go_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        # with genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'go_id', 'query':'GO:0000028'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00110')
        self.assertNotContains(response, 'C_RS00120')
        # without genome keyword
        response = self.client.get('/loadinggenesearch/',
                                   {'type':'go_id', 'query':'GO:0000028'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00110')
        self.assertContains(response, 'C_RS00120')

    def test_loading_gene_search_by_ko_page(self):
        '''
            Testing Genes search results for ko type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'ko', 'query':'thrA'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')

    def test_loading_gene_search_by_kp_page(self):
        '''
            Testing Genes search results for kp type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'kp', 'query':'Lysine'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_kp_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')

    def test_loading_gene_search_by_kr_page(self):
        '''
            Testing Genes search results for kr type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'kr', 'query':'ATP:L-aspartate'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_kr_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')

    def test_loading_gene_search_by_tc_page(self):
        '''
            Testing Genes search results for tc type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'tc', 'query':'Hsp70'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_tc_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00070')

    def test_loading_gene_search_by_ec_page(self):
        '''
            Testing Genes search results for ec type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'ec', 'query':'HSDH'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_ec_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')

    def test_loading_gene_search_by_cazy_page(self):
        '''
            Testing Genes search results for cazy type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_CFT073', 'type':'cazy', 'query':'cellobiose'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_cazy_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'C_RS00340')

    def test_loading_gene_search_by_cog_page(self):
        '''
            Testing Genes search results for cog type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'cog', 'query':'Energy'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_cog_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00200')

    def test_loading_gene_search_by_go_page(self):
        '''
            Testing Genes search results for go_id type
            path('loadinggenesearch/',
                 views.GeneSearchResultsAjaxView.ajax_view,
                 name='loadinggenesearch'
                 )
        '''
        response = self.client.get('/loadinggenesearch/',
                                   {'genome':'E_coli_BW2952', 'type':'go', 'query':'homoserine dehydrogenase activity'}
                                   )
        self.assertEqual(response.status_code, 200)
        #print('test_loading_gene_search_by_go_page\n', response.content)
        self.assertContains(response, 'Gene ID')
        self.assertContains(response, 'BWG_RS00010')

    def test_genome_search_page(self):
        '''
            Testing Genome search results view
            path('searchgenome/', views.GenomeSearchResultsView.as_view(), name='searchgenome')
        '''
        # with query keyword
        response = self.client.get('/searchgenome/',{'query':'E_coli_BW2952'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genomes')
        self.assertContains(response, 'Escherichia coli BW2952')
        # with taxon keyword
        response = self.client.get('/searchgenome/',{'taxon':'562'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genomes')
        self.assertContains(response, 'Escherichia coli BW2952')
        # with empty query
        response = self.client.get('/searchgenome/',{'query':''})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genomes')
        self.assertContains(response, 'No genomes found.')

    def test_strain_search_page(self):
        '''
            Testing Strain search results view
            path('searchstrain/', views.StrainSearchResultsView.as_view(), name='searchstrain')
        '''
        # with query keyword
        response = self.client.get('/searchstrain/',{'query':'BW2952'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Strain')
        self.assertContains(response, 'Escherichia coli BW2952')
        # with taxon keyword
        response = self.client.get('/searchstrain/',{'taxon':'562'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Strain')
        self.assertContains(response, 'Escherichia coli BW2952')
        # with empty query
        response = self.client.get('/searchstrain/',{'query':''})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Strain')
        self.assertContains(response, 'No strains found.')

    def test_sample_search_page(self):
        '''
            Testing Sample search results view
            path('searchsample/', views.SampleSearchResultsView.as_view(), name='searchsample')
        '''
        response = self.client.get('/searchsample/',{'query':'test_sample'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sample')
        self.assertContains(response, 'Test sample')
        # with empty query
        response = self.client.get('/searchsample/',{'query':''})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sample')
        self.assertContains(response, 'No samples found.')

    def test_taxon_search_page(self):
        '''
            Testing Taxon search results view
            path('searchtaxon/', views.TaxonSearchResultsView.as_view(), name='searchtaxon')
        '''
        response = self.client.get('/searchtaxon/',{'query':'Escherichia'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Taxon')
        self.assertContains(response, 'Escherichia coli')
        response = self.client.get('/searchtaxon/',{})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Query string is empty')

    def test_loading_og_treemap_page(self):
        '''
            Testing Genes list for a genome view with "genome" parameter
            path('loadingogtreemap/',  # No test!
                 views.get_og_data,
                 name='loadingogtreemap'
                 )
        '''
        response = self.client.get('/loadingogtreemap/', {'og':102483})
        self.assertEqual(response.status_code, 200)
        #print('test_loading_og_treemap_page response\n', response.content)
        self.assertContains(response, 'tsv_profile')
        self.assertContains(response, 'GHMP_kinases_C')

    def test_annotationlist_page(self):
        '''
            Testing Annotation list view
            path('searchannotation/',
                views.AnnotationSearchResultsAjaxView.as_view(),
                name='searchannotation'
                ),
        '''
        # Annotation list for "Pfam" text query
        response = self.client.get('/searchannotation/',
                                   {'annotation_query':'Pfam'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response,
            'data: {\'annotation_query\': "Pfam", \'genome\': "", \'page\': "", \'type\': "annotation" }')
        # Annotation list without query
        response = self.client.get('/searchannotation/',
                                   {}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Query string is empty')

    def test_annotationlist_genome_page(self):
        '''
            Testing Annotation list view for a  genome
            path('searchannotation/',
                views.AnnotationSearchResultsAjaxView.as_view(),
                name='searchannotation'
                ),
        '''
        # Annotation list for "Pfam" text query in one genome
        response = self.client.get('/searchannotation/',
                                   {'annotation_query':'Pfam', 'genome':'E_coli_BW2952'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response,
            'data: {\'annotation_query\': "Pfam", \'genome\': "E_coli_BW2952", \'page\': "", \'type\': "annotation" }')

    def test_annotationlist_subpage(self):
        '''
            Testing AnnotationSearchResultsSubView class
        '''
        # with text query only
        request = RequestFactory().get('/loadingtextsearch/', {'annotation_query':'Pfam'})
        view = AnnotationSearchResultsSubView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Pfam database')
        # with text query and genome
        request = RequestFactory().get('/loadingtextsearch/', {'annotation_query':'Pfam', 'genome':'E_coli_BW2952'})
        view = AnnotationSearchResultsSubView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Pfam database')
        # without text query
        request = RequestFactory().get('/loadingtextsearch/')
        view = AnnotationSearchResultsSubView()
        view.setup(request)
        response = view.get(request)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'No annotations found.')

    def test_operon_page(self):
        '''
            Testing Operon page view
            path('operon/<str:genome>/<str:name>/', views.operon_detail, name='operondetails')
        '''
        genome_id = 'E_coli_BW2952'
        operon_id = Operon.objects.values_list('name', flat=True)[0]
        response = self.client.get('/operon/' + genome_id + '/' + operon_id + '/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operon information')
        self.assertContains(response, operon_id)
        response = self.client.get('/operon/' + genome_id + '/XXXXXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_conserved_operon_page(self):
        '''
            Testing Conserved operon page view
            path('coperon/<int:operon_id>/', views.conserved_operon_view, name='coperon')
        '''
        response = self.client.get('/coperon/65903/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Conserved operon</h2>')
        self.assertContains(response, 'conservoperondata/65903/')

    def test_conserved_operon_data_page(self):
        '''
            Testing Conserved operon data page view
            path('conservoperondata/<int:operon_id>/', views.conserved_operon_data, name='conservoperondata')
        '''
        response = self.client.get('/conservoperondata/65903/')
        #print('test_conserved_operon_data_page\n', response.content[:1000])
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'tsv_profile')
        self.assertContains(response, 'E_coli_BW2952_operon_49')

    def test_pathway_map_page(self):
        '''
            Testing KEGG pathway map view
            path('pathway/', views.pathway_view, name='pathway')
        '''
        response = self.client.get('/pathway/', {'genome':'E_coli_BW2952', 'pathway':'map00300'})
        #print('test_pathway_map_page\n', response.content[:1000])
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Pathway map</h2>')
        self.assertContains(response, 'BWG_RS00010')

    def test_site_page(self):
        '''
            Testing Site page view
            path('site/<str:genome>/<str:name>/', views.site_detail, name='sitedetails')
        '''
        genome_id = 'E_coli_BW2952'
        site_id = Site.objects.filter(
            genome__name=genome_id
        ).values_list(
            'name', flat=True
        )[0]
        response = self.client.get('/site/' + genome_id + '/' + site_id + '/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Site information')
        self.assertContains(response, site_id)

        response = self.client.get('/site/' + genome_id + '/XXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_regulon_page(self):
        '''
            Testing Regulon page view
            path('regulon/<str:genome>/<str:name>/', views.regulon_detail, name='regulondetails')
        '''
        genome_id = 'E_coli_BW2952'
        response = self.client.get('/regulon/' + genome_id + '/AraC/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulon information')
        self.assertContains(response, 'AraC')

        response = self.client.get('/regulon/' + genome_id + '/XXXXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_conserved_regulon_page(self):
        '''
            Testing Conserved regulon page view
            path('cregulon/', views.cregulon_view, name='cregulon')
        '''
        genome_id = 'E_coli_BW2952'
        response = self.client.get('/cregulon/', {'og':102420})
        #print('test_conserved_regulon_page\n',response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Conserved regulon</h2>')
        self.assertContains(response, 'AraC')

    def test_ogroup_page(self):
        '''
            Testing Orthogroup page view
            path('ogroup/<int:og_id>/', views.og_detail, name='ogdetails'),
        '''
        response = self.client.get('/ogroup/102405/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Ortholog group')
        self.assertContains(response, '1RIFA')

        response = self.client.get('/ogroup/0000000000/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'does not exist')

    def test_textsearch_page(self):
        '''
            Testing search page view'
            path('textsearch/', views.textsearch, name='textsearch'),
        '''
        response = self.client.get('/textsearch/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Find genes by gene, contig, genome identifier or product')

    def test_loading_textsearch_ajax_page(self):
        '''
            Testing search page view'
            path('loadingtextsearch/',
                 views.AnnotationSearchResultsAjaxView.ajax_view,
                 name="loadingtextsearch"
                 )
        '''
        response = self.client.get('/loadingtextsearch/?annotation_query=Pfam&genome=E_coli_BW2952')
        #print('test_loading_textsearch_ajax_page response\n',response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Annotation')
        self.assertContains(response, 'Pfam database')

    def test_kolist_page(self):
        '''
            Testing KEGG orthologs list page view
            path('kos/', views.KoSearchResultsView.as_view(), name='kos')
        '''
        response = self.client.get('/kos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_kolist_query_page(self):
        '''
            Testing KEGG orthologs list page view with a text query
            path('kos/', views.KoSearchResultsView.as_view(), name='kos')
        '''
        response = self.client.get('/kos/', {'query':'isopropylmalate'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

        response = self.client.get('/kos/', {'query':'isopropylmalate', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')


    def test_kolist_genome_page(self):
        '''
            Testing KEGG orthologs list page view for a genome
            path('kos/', views.KoSearchResultsView.as_view(), name='kos')
        '''
        response = self.client.get('/kos/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_pathwayslist_page(self):
        '''
            Testing KEGG pathways list page view
            path('pathways/', views.KpSearchResultsView.as_view(), name='kps')
        '''
        response = self.client.get('/pathways/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_pathwayslist_query_page(self):
        '''
            Testing KEGG pathways list page view with a text query
            path('pathways/', views.KpSearchResultsView.as_view(), name='kps')
        '''
        response = self.client.get('/pathways/', {'query':'Pentose'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

        response = self.client.get('/pathways/', {'query':'Pentose', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_pathwayslist_genome_page(self):
        '''
            Testing KEGG pathways list page view for a genome
            path('pathways/', views.KpSearchResultsView.as_view(), name='kps')
        '''
        response = self.client.get('/pathways/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_reactionslist_page(self):
        '''
            Testing KEGG reactions list page view
            path('reactions/', views.KrSearchResultsView.as_view(), name='krs')
        '''
        response = self.client.get('/reactions/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R00006')

    def test_reactionslist_query_page(self):
        '''
            Testing KEGG reactions list page view with a text query
            path('reactions/', views.KrSearchResultsView.as_view(), name='krs')
        '''
        response = self.client.get('/reactions/', {'query':'tetrahydrobiopterin'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R11765')

        response = self.client.get('/reactions/', {'query':'tetrahydrobiopterin', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R11765')

    def test_reactionslist_genome_page(self):
        '''
            Testing KEGG reactions list page view for a genome
            path('reactions/', views.KrSearchResultsView.as_view(), name='krs')
        '''
        response = self.client.get('/reactions/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R00006')

    def test_enzymeslist_page(self):
        '''
            Testing EC numbers list page view
            path('enzymes/', views.EcSearchResultsView.as_view(), name='enzymes')
        '''
        response = self.client.get('/enzymes/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.262')

    def test_enzymeslist_query_page(self):
        '''
            Testing EC numbers list page view with a text query
            path('enzymes/', views.EcSearchResultsView.as_view(), name='enzymes')
        '''
        response = self.client.get('/enzymes/', {'query':'homoserine'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.3')

        response = self.client.get('/enzymes/', {'query':'homoserine', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.3')

    def test_enzymeslist_genome_page(self):
        '''
            Testing EC numbers list page view for a genome
            path('enzymes/', views.EcSearchResultsView.as_view(), name='enzymes')
        '''
        response = self.client.get('/enzymes/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.3')

    def test_transporterslist_page(self):
        '''
            Testing TC families list page view
            path('transporters/', views.TcSearchResultsView.as_view(), name='transporters')
        '''
        response = self.client.get('/transporters/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.A.33.1')

    def test_transporterslist_query_page(self):
        '''
            Testing TC families list page view with a text query
            path('transporters/', views.TcSearchResultsView.as_view(), name='transporters')
        '''
        response = self.client.get('/transporters/', {'query':'ABC'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '3.A.1')

        response = self.client.get('/transporters/', {'query':'ABC', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '3.A.1')

    def test_transporterslist_genome_page(self):
        '''
            Testing TC families list page view for a genome
            path('transporters/', views.TcSearchResultsView.as_view(), name='transporters')
        '''
        response = self.client.get('/transporters/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '3.A.1')

    def test_cazylist_page(self):
        '''
            Testing CAZy families list page view
            path('cazy/', views.CazySearchResultsView.as_view(), name='cazy')
        '''
        response = self.client.get('/cazy/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')
        self.assertContains(response, 'GH94')

    def test_cazylist_query_page(self):
        '''
            Testing CAZy families list page view with a text query
            path('cazy/', views.CazySearchResultsView.as_view(), name='cazy')
        '''
        response = self.client.get('/cazy/', {'query':'GH94'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')
        self.assertContains(response, 'cellobiose')

        response = self.client.get('/cazy/', {'query':'GH94', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')
        self.assertContains(response, 'cellobiose')

    def test_cazylist_genome_page(self):
        '''
            Testing CAZy families list page view for a genome
            path('cazy/', views.CazySearchResultsView.as_view(), name='cazy')
        '''
        response = self.client.get('/cazy/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')
        self.assertContains(response, 'cellobiose')

    def test_cogslist_page(self):
        '''
            Testing COG classes list page view
            path('cogs/', views.CogSearchResultsView.as_view(), name='cogs'),
        '''
        response = self.client.get('/cogs/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Energy production and conversion')

    def test_cogslist_query_page(self):
        '''
            Testing COG classes list page view with a text query
            path('cogs/', views.CogSearchResultsView.as_view(), name='cogs'),
        '''
        response = self.client.get('/cogs/', {'query':'I'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Lipid transport and metabolism')

        response = self.client.get('/cogs/', {'query':'I', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Lipid transport and metabolism')

    def test_cogslist_genome_page(self):
        '''
            Testing COG classes list page view for a genome
            path('cogs/', views.CogSearchResultsView.as_view(), name='cogs'),
        '''
        response = self.client.get('/cogs/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Energy production and conversion')

    def test_golist_page(self):
        '''
            Testing GO terms list page view
            path('gos/', views.GoSearchResultsView.as_view(), name='gos')
        '''
        response = self.client.get('/gos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'GO:0000003')

    def test_golist_query_page(self):
        '''
            Testing GO terms list page view with a text query
            path('gos/', views.GoSearchResultsView.as_view(), name='gos')
        '''
        response = self.client.get('/gos/', {'query':'GO:0003824'})
        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, 'GO:0000015')
        self.assertContains(response, 'catalytic activity')

        response = self.client.get('/gos/', {'query':'GO:0003824', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, 'GO:0000015')
        self.assertContains(response, 'catalytic activity')

    def test_golist_genome_page(self):
        '''
            Testing GO terms list page view with a text query
            path('gos/', views.GoSearchResultsView.as_view(), name='gos')
        '''
        response = self.client.get('/gos/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'GO:0003824')

    def test_proteinsearch_page(self):
        '''
            Testing protein search form page view
            path('protsearchform/',views.proteinsearchform,name="proteinsearchform")
        '''
        response = self.client.get('/protsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sequence search')

    def test_proteinsearchresults_page(self):
        '''
            Testing protein search results page view
            path('protsearch/',views.PsearchResultView.as_view(),name="proteinsearch")
        '''
        recreate_search_databases()
        sequence = '''>C_RS00020
MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERF
CQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVA
PCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYS
RQPELAANLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPDTAQRVADWLGKNYL
QNQEGFVHICRLDTAGARVLEN'''
        page = self.csrf_client.get('/protsearchform/')
        token = page.context.get("csrf_token")
        response = self.csrf_client.post('/protsearch/', {'sequence':sequence, 'evalue':'0.000001', 'hitstoshow': '10', 'csrfmiddlewaretoken':token})
        #print('Protein search\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Sequence search</h2>')
        self.assertContains(response, 'url: "/loadingprotsearch/"')
        self.assertContains(response, 'EC_RS00020')

    def test_loading_proteinsearch_results_page(self):
        '''
            Testing actual protein search results page view
            path('loadingprotsearch/',views.PsearchResultView.ajax_view,name="loadingprotsearch")
        '''
        recreate_search_databases()
        sequence = '''>C_RS00020
MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERF
CQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVA
PCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYS
RQPELAANLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPDTAQRVADWLGKNYL
QNQEGFVHICRLDTAGARVLEN'''
        page = self.csrf_client.get('/protsearchform/')
        token = page.context.get("csrf_token")
        response = self.csrf_client.post('/loadingprotsearch/', {'sequence':sequence, 'evalue':'0.000001', 'hitstoshow': '10', 'csrfmiddlewaretoken':token})
        #print('test_loading_proteinsearch_results_page\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<th>Target gene</th>')
        self.assertContains(response, 'C_RS00020')


    def test_protein_external_search_page(self):
        '''
            Testing external protein search results page view
            path('seqsearch', views.protein_search_external, name='proteinsearchexternal')
        '''
        recreate_search_databases()
        sequence = '''>C_RS00020%0AMVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERF\nCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVA\nPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYS\nRQPELAANLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPDTAQRVADWLGKNYL\nQNQEGFVHICRLDTAGARVLEN'''
        response = self.client.get('/seqsearch/', {'sequence':sequence})
        #print('External protein search\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Sequence search</h2>')
        self.assertContains(response, 'Target gene')
        self.assertContains(response, 'E_coli_CFT073')

    def test_nucleotidesearch_page(self):
        '''
            Testing nucleotide sequence search page view
            path('nuclsearchform/',views.nucleotidesearchform,name="nucleotidesearchform")
        '''
        response = self.client.get('/nuclsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response,
                            'Enter one nucleotide sequence in FASTA format'
                            )

    def test_nucleotidesearchresults_page(self):
        '''
            Testing nucleotide sequence search page view
            path('nuclsearch/',views.NsearchResultView.as_view(),name="nucleotidesearch")
        '''
        sequence = '''>NC_004431:3512..4444
ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCT
GTTGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCGGCAGAGACATTCAGTCTCAACAACCTCGGACGCTTT
GCCGATAAGCTGCCGTCAGAACCTCGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCAGGAGCTGGGTAAG
CAAATTCCAGTGGCGATGACTCTGGAAAAGAATATGCCAATCGGCTCGGGCTTAGGCTCCAGTGCCTGTTCGGTGGTC
GCGGCGCTGATGGCGATGAATGAACACTGTGGCAAGCCGCTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTG
GAAGGACGAATCTCCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTTGGTGGTATGCAGTTGATGATCGAA
GAAAACGACATCATCAGCCAGCAAGTGCCTGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGGATTAAAGTC
TCGACGGCAGAAGCCAGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATTTGGCA
GGCTTCATTCACGCCTGCTATTCCCGTCAGCCTGAGCTTGCCGCGAATCTGATGAAAGATGTTATCGCTGAACCCTAC
CGTGAACGGTTACTGCCTGGTTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCC
GGCTCCGGCCCGACTTTGTTCGCTCTGTGTGACAAGCCGGATACCGCCCAGCGCGTTGCCGACTGGTTGGGTAAGAAC
TACCTGCAAAATCAGGAAGGTTTTGTTCATATTTGCCGGCTGGATACGGCGGGCGCACGAGTACTGGAAAACTAA'''
        page = self.csrf_client.get('/nuclsearchform/')
        token = page.context.get("csrf_token")
        response = self.csrf_client.post('/nuclsearch/', {'sequence':sequence, 'evalue':'0.000001', 'hitstoshow': '10', 'csrfmiddlewaretoken':token})
        #print('Nucleotide search\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Sequence search</h2>')
        self.assertContains(response, 'url: "/loadingnuclsearch/"')
        self.assertContains(response, 'NC_004431:3512..4444')

    def test_loading_nuclsearch_results_page(self):
        '''
            Testing actual nucleotide sequence search results view
            path('loadingnuclsearch/',views.NsearchResultView.ajax_view,name="loadingnuclsearch")
        '''
        recreate_search_databases()
        sequence = '''>NC_004431:3512..4444
ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCT
GTTGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCGGCAGAGACATTCAGTCTCAACAACCTCGGACGCTTT
GCCGATAAGCTGCCGTCAGAACCTCGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCAGGAGCTGGGTAAG
CAAATTCCAGTGGCGATGACTCTGGAAAAGAATATGCCAATCGGCTCGGGCTTAGGCTCCAGTGCCTGTTCGGTGGTC
GCGGCGCTGATGGCGATGAATGAACACTGTGGCAAGCCGCTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTG
GAAGGACGAATCTCCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTTGGTGGTATGCAGTTGATGATCGAA
GAAAACGACATCATCAGCCAGCAAGTGCCTGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGGATTAAAGTC
TCGACGGCAGAAGCCAGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATTTGGCA
GGCTTCATTCACGCCTGCTATTCCCGTCAGCCTGAGCTTGCCGCGAATCTGATGAAAGATGTTATCGCTGAACCCTAC
CGTGAACGGTTACTGCCTGGTTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCC
GGCTCCGGCCCGACTTTGTTCGCTCTGTGTGACAAGCCGGATACCGCCCAGCGCGTTGCCGACTGGTTGGGTAAGAAC
TACCTGCAAAATCAGGAAGGTTTTGTTCATATTTGCCGGCTGGATACGGCGGGCGCACGAGTACTGGAAAACTAA'''
        page = self.csrf_client.get('/nuclsearchform/')
        token = page.context.get("csrf_token")
        response = self.csrf_client.post('/loadingnuclsearch/', {'sequence':sequence, 'evalue':'0.000001', 'hitstoshow': '10', 'csrfmiddlewaretoken':token})
        #print('test_loading_nuclsearch_results_page\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Target contig')
        self.assertContains(response, 'NC_004431: (+strand) 3512..4444')

    def test_help_page(self):
        '''
            Testing help page view
            path('help/', views.show_help, name='help')
        '''
        response = self.client.get('/help/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'About this site')

    def test_comparative_view(self):
        '''
            Testing comparative plot view
            path('comparative/', views.ComparativeView.as_view(), name='comparative')
        '''
        response = self.client.get('/comparative/',
                                   {'genome':'E_coli_BW2952',
                                    'locus_tag':'BWG_RS00020',
                                    'og':'102492',
                                    'size':'10',
                                    'lines':'50'
                                    }
                                   )
        #print('test_comparative_view\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Comparative analysis</h2>')

    def test_loading_comparative_view(self):
        '''
            Testing loading comparative plot view
            path('loadingscribl/',views.ComparativeView.ajax_view,name="loadingscribl")
        '''
        response = self.client.get('/loadingscribl/',
                                   {'genome':'E_coli_BW2952',
                                    'locus_tag':'BWG_RS00020',
                                    'og':'102492',
                                    'size':'10',
                                    'lines':'50'
                                    }
                                   )
        #print('test_loading_comparative_view\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<script type=\\"text/javascript\\">')
        self.assertContains(response, 'E_coli_CFT073')

    def test_pagenotfound(self):
        '''
            Test for non-existent URL view
            Invokes views.handler404
        '''
        response = self.client.get('/nonexistingpage')
        self.assertEqual(response.status_code, 404)

    def test_getgene_view(self):
        '''
            Testing get gene view
            path('getgene/', views.gene_byname, name='genebyname')
        '''
        response = self.client.get('/getgene/', {'genome':'E_coli_CFT073','locus':'C_RS00015'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'C_RS00015')

    def test_export_genome_csv(self):
        '''
            Testing CSV export of genomes table
            path('export/', export_text.export_csv, name='export')
        '''
        response = self.client.get('/export/', {'query':'E_coli_CFT073','type':'genome'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Name\tTags\t')
        self.assertContains(response, 'E_coli_CFT073')

        response = self.client.get('/export/', {'type':'genome'})
        #print('Genome list export\n', response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Name\tTags\t')
        self.assertContains(response, 'E_coli_CFT073')

    def test_export_annotation_csv(self):
        '''
            Testing CSV export of annotations table
            path('export/', export_text.export_csv, name='export')
        '''
        response = self.client.get('/export/', {'annotation_query':'Pfam','type':'annotation'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '\tAnnotation_source\t')
        self.assertContains(response, 'E_coli_CFT073')

        response = self.client.get('/export/', {'annotation_query':'Pfam','type':'annotation', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '\tAnnotation_source\t')
        self.assertContains(response, 'E_coli_CFT073')

    def test_export_gene_csv(self):
        '''
            Testing CSV export of genes table
            path('export/', export_text.export_csv, name='export')
        '''
        response = self.client.get('/export/', {'type':'gene', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'E_coli_CFT073')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'query':'C_RS00015','type':'gene'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'E_coli_CFT073')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'query':'C_RS00015','type':'gene', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'E_coli_CFT073')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'og', 'query':102482})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'ko_id', 'query':'K12524'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'kp_id', 'query':'map00300'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'kr_id', 'query':'R00480'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'ec_id', 'query':'1.1.1.3'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'tc_id', 'query':'1.A.33.1'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00070')
        self.assertContains(response, 'C_RS00080')

        response = self.client.get('/export/', {'type':'cazy_id', 'query':'GH94'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'C_RS00340')
        
        response = self.client.get('/export/', {'type':'cog_id', 'query':'C'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00200')
        self.assertContains(response, 'C_RS00220')

        response = self.client.get('/export/', {'type':'go_id', 'query':'GO:0000028'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00110')
        self.assertContains(response, 'C_RS00120')

        response = self.client.get('/export/', {'type':'ko', 'query':'thrA'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'kp', 'query':'Lysine'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'kr', 'query':'ATP:L-aspartate'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')
        self.assertContains(response, 'C_RS00015')

        response = self.client.get('/export/', {'type':'tc', 'query':'Hsp70'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00070')

        response = self.client.get('/export/', {'type':'ec', 'query':'HSDH'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')

        response = self.client.get('/export/', {'type':'cazy', 'query':'cellobiose'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'C_RS00340')

        response = self.client.get('/export/', {'type':'cog', 'query':'Energy'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00200')

        response = self.client.get('/export/', {'type':'go', 'query':'homoserine dehydrogenase activity'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome')
        self.assertContains(response, 'BWG_RS00010')


    def test_export_fasta(self):
        '''
            Testing export of protein fasta
            path('exportfasta/', export_text.export_fasta, name='exportfasta')
        '''
        response = self.client.get('/exportfasta/', {'query':'C_RS00015','type':'gene', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'og', 'query':102482})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'ko_id', 'query':'K12524'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'kp_id', 'query':'map00300'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'kr_id', 'query':'R00480'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'ec_id', 'query':'1.1.1.3'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'tc_id', 'query':'1.A.33.1'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00070')
        self.assertContains(response, 'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEG')

        response = self.client.get('/exportfasta/', {'type':'cazy_id', 'query':'GH94'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00340')
        
        response = self.client.get('/exportfasta/', {'type':'cog_id', 'query':'C'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00200')
        self.assertContains(response, '>C_RS00220')
        self.assertContains(response, 'MDHLPMPKFGPLAGLRVVFSGIEIAGPFAGQ')

        response = self.client.get('/exportfasta/', {'type':'go_id', 'query':'GO:0000028'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00110')
        self.assertContains(response, '>C_RS00120')

        response = self.client.get('/exportfasta/', {'type':'ko', 'query':'thrA'})
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'kp', 'query':'Lysine'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'kr', 'query':'ATP:L-aspartate'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'tc', 'query':'Hsp70'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00070')
        self.assertContains(response, 'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEG')

        response = self.client.get('/exportfasta/', {'type':'ec', 'query':'HSDH'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')

        response = self.client.get('/exportfasta/', {'type':'cazy', 'query':'cellobiose'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00340')

        response = self.client.get('/exportfasta/', {'type':'cog', 'query':'Energy'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00200')

        response = self.client.get('/exportfasta/', {'type':'go', 'query':'homoserine dehydrogenase activity'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>BWG_RS00010')
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'query':'K12524','type':'ko_id', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'annotation_query':'AA_kinase', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'annotation_query':'AA_kinase'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'query':'C_RS00015','type':'gene'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'type':'gene', 'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

        response = self.client.get('/exportfasta/', {'genome':'E_coli_CFT073'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '>C_RS00015')
        self.assertContains(response, 'MRVLKFGGTSVANAERFLRVADILESNARQ')

    def test_export_gbk(self):
        '''
            Testing GBK file export
            path('exportgbk/<str:name>/', export_text.export_gbk, name='exportgbk')
        '''
        response = self.client.get('/exportgbk/E_coli_CFT073/')
        self.assertEqual(response.status_code, 200)

        response = self.client.get('/exportgbk/XXXXXXXXXXX/')
        self.assertEqual(response.status_code, 200)

    '''
        Accessory functions testing
    '''
    def test_comparative_get_color(self):
        print('Testing get_color function for comparative plot')
        color = _get_color(6)
        self.assertEqual(color,
            '.setColorGradient(\'rgb(255, 182, 199)\', \'rgb(204, 102, 119)\')'
            )

    def test_make_muscle_alignment(self):
        print('Testing multiple alignment with muscle')
        proteins = '>1\nMKRISTTITTTTITTGNGAG\n' +\
        '>2\nMKRISTTITTTITITTGNGAG\n' +\
        '>3\nMKRISTTITTTITITTGNGAG'
        outfasta = make_muscle_alignment(proteins)
        outfasta_lines = outfasta.split('\n')
        self.assertEqual(outfasta_lines[-6], 'MKRISTTITTT-TITTGNGAG')

    def test_sanitize_sequence(self):
        print('Testing _sanitize_sequence for sequence similarity search')
        protein = '>1\nMKRISTTITTTTITTGNGAG\n'
        seq_id, seq = _sanitize_sequence(protein)
        self.assertEqual(seq_id, '1')
        self.assertEqual(seq, 'MKRISTTITTTTITTGNGAG')
        
        protein = '>qqq\n MKRISTTITTTTI\nTTGNGAG*\n'
        seq_id, seq = _sanitize_sequence(protein)
        self.assertEqual(seq_id, 'qqq')
        self.assertEqual(seq, 'MKRISTTITTTTITTGNGAG')

        nucleotide = '''>qqq
ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCT'''
        seq_id, seq = _sanitize_sequence(nucleotide)
        self.assertEqual(seq_id, 'qqq')
        self.assertEqual(seq, 'ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCT')

        nucleotide = '''>qqq
        ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGG
        GGCGGCGGTGACACCT* '''
        seq_id, seq = _sanitize_sequence(nucleotide)
        self.assertEqual(seq_id, 'qqq')
        self.assertEqual(seq, 'ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCT')


