import os
import hashlib
from unittest import skip
#from contextlib import contextmanager

from django.test import TestCase
from django.test import Client
from django.utils import timezone
#from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.db.utils import IntegrityError

#from browser.models import Taxon
#from browser.models import Config
#from browser.models import Strain
#from browser.models import Sample
#from browser.models import Strain_metadata
#from browser.models import Genome
#from browser.models import Annotation
#from browser.models import Regulon
#from browser.models import Contig
#from browser.models import Sample_metadata
#from browser.models import Kegg_ortholog
#from browser.models import Kegg_pathway
#from browser.models import Kegg_reaction
#from browser.models import Go_term
#from browser.models import Cog_class
#from browser.models import Ec_number
#from browser.models import Cazy_family
from browser.models import Operon
#from browser.models import Ortholog_group
#from browser.models import Eggnog_description
#from browser.models import Tc_family
#from browser.models import Protein
#from browser.models import Gene
from browser.models import Site
from genomebrowser.settings import BASE_DIR
from browser.pipeline.genome_import import Importer
from browser.comparative_analysis import _get_color, make_muscle_alignment

# Create your tests here.

class BrowserViewsTest(TestCase):
    '''
        Views testing
    '''
    fixtures = ['minigenomes.testdata.json']
    
    @classmethod
    def setUpTestData(cls):
        cls.client = Client()
        
    def test_start_page(self):
        print('Testing Start page view')
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Welcome')

    def test_genomes_page(self):
        print('Testing Genomes list page view')
        response = self.client.get('/genomes/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Genomes</h2>')
        self.assertContains(response, 'E_coli_BW2952')

    def test_strains_page(self):
        print('Testing Strains list page view')
        response = self.client.get('/strains/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Strains</h2>')
        self.assertContains(response, 'BW2952')

    def test_samples_page(self):
        print('Testing Samples list page view')
        response = self.client.get('/samples/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Samples</h2>')
        self.assertContains(response, 'No samples found.')

    def test_genome_page(self):
        print('Testing Genome page view')
        response = self.client.get('/genome/E_coli_BW2952')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'E_coli_BW2952')

    def test_gene_page(self):
        print('Testing Gene page view')
        response = self.client.get('/gene/E_coli_BW2952/BWG_RS00020/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'BWG_RS00020')

    def test_operonlist_page(self):
        print('Testing Operons list page view')
        response = self.client.get('/operons/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operons in')
        self.assertContains(response, 'NC_012759: 190..5020')

    def test_sitelist_page(self):
        print('Testing Sites list page view')
        response = self.client.get('/sites/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sites in')
        self.assertContains(response, 'NC_012759: complement(70130..70146)')

    def test_regulonlist_page(self):
        print('Testing Regulons list page view')
        response = self.client.get('/regulons/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulons in')
        self.assertContains(response, 'AraC')

    def test_genelist1_page(self):
        print('Testing Genes list for a genome view')
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

    def test_gene_search_byregulator_page(self):
        print('Testing Genes list for a regulator view')
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

    def test_annotationlist_page(self):
        print('Testing Annotation list view')
        # Annotation list for "Pfam" text query
        response = self.client.get('/searchannotation/',
                                   {'annotation_query':'Pfam'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response,
            'data: {\'annotation_query\': "Pfam", \'genome\': "", \'page\': "", \'type\': "annotation" }')

    def test_operon_page(self):
        print('Testing Operon page view')
        genome_id = 'E_coli_BW2952'
        operon_id = Operon.objects.values_list('name', flat=True)[0]
        response = self.client.get('/operon/' + genome_id + '/' + operon_id + '/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operon information')
        self.assertContains(response, operon_id)

    def test_site_page(self):
        print('Testing Site page view')
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

    def test_regulon_page(self):
        print('Testing Regulon page view')
        genome_id = 'E_coli_BW2952'
        response = self.client.get('/regulon/' + genome_id + '/AraC/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulon information')
        self.assertContains(response, 'AraC')

    def test_textsearch_page(self):
        print('Testing search page view')
        response = self.client.get('/textsearch/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Find genes by gene, contig, genome identifier or product')

    def test_kolist_page1(self):
        print('Testing KEGG orthologs list page view')
        response = self.client.get('/kos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_kolist_page2(self):
        print('Testing KEGG orthologs list page view with a text query')
        response = self.client.get('/kos/', {'query':'isopropylmalate'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_pathwayslist_page1(self):
        print('Testing KEGG pathways list page view')
        response = self.client.get('/pathways/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_pathwayslist_page2(self):
        print('Testing KEGG pathways list page view with a text query')
        response = self.client.get('/pathways/', {'query':'Pentose'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_reactionslist_page1(self):
        print('Testing KEGG reactions list page view')
        response = self.client.get('/reactions/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R00006')

    def test_reactionslist_page2(self):
        print('Testing KEGG reactions list page view with a text query')
        response = self.client.get('/reactions/', {'query':'tetrahydrobiopterin'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R11765')

    def test_enzymeslist_page1(self):
        print('Testing EC numbers list page view')
        response = self.client.get('/enzymes/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.262')

    def test_enzymeslist_page2(self):
        print('Testing EC numbers list page view with a text query')
        response = self.client.get('/enzymes/', {'query':'homoserine'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.3')

    def test_transporterslist_page1(self):
        print('Testing TC families list page view')
        response = self.client.get('/transporters/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.A.33.1')

    def test_cazylist_page(self):
        print('Testing CAZy families list page view')
        response = self.client.get('/cazy/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')

    def test_cogslist_page1(self):
        print('Testing COG classes list page view')
        response = self.client.get('/cogs/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Energy production and conversion')

    def test_golist_page1(self):
        print('Testing GO terms list page view')
        response = self.client.get('/gos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'GO:0000003')

    def test_proteinsearch_page(self):
        print('Testing protein sequence search page view')
        response = self.client.get('/protsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sequence search')

    def test_nucleotidesearch_page(self):
        print('Testing nucleotide sequence search page view')
        response = self.client.get('/nuclsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response,
                            'Enter one nucleotide sequence in FASTA format'
                            )

    def test_help_page(self):
        print('Testing help page view')
        response = self.client.get('/help/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'About this site')

    def test_comparative_view(self):
        print('Testing non-existant page view')
        response = self.client.get('/comparative/',
                                   {'genome':'E_coli_BW2952',
                                    'locus_tag':'BWG_RS00020',
                                    'og':'102492',
                                    'size':'10',
                                    'lines':'50'
                                    }
                                   )
        #print(response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Comparative analysis</h2>')

    def test_pagenotfound(self):
        response = self.client.get('/nonexistingpage')
        self.assertEqual(response.status_code, 404)

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


