import time
from unittest import skip
from contextlib import contextmanager
from django.test import TransactionTestCase
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium.webdriver.firefox.webdriver import WebDriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

from browser.dataimport.importer import Importer
# Create your tests here.

class ImporterTestCase(TransactionTestCase):
    fixtures = ['testdata.json']
    
    def setUp(self):
        self.importer = Importer()

    @skip("passed")
    def test_find_taxonomic_order(self):
        """Taxonomic order correctly identified"""
        pseudomonadales = self.importer.get_taxonomic_order('286')
        self.assertEqual(pseudomonadales, 'Pseudomonadales')
        pseudomonadales = self.importer.get_taxonomic_order('72274')
        self.assertEqual(pseudomonadales, 'Pseudomonadales')
        unknown = self.importer.get_taxonomic_order('0')
        self.assertEqual(unknown, 'Unknown')
        unknown = self.importer.get_taxonomic_order('randomnonsense')
        self.assertEqual(unknown, 'Unknown')
        unknown = self.importer.get_taxonomic_order('2')
        self.assertEqual(unknown, 'Unknown')
        unknown = self.importer.get_taxonomic_order('0')
        self.assertEqual(unknown, 'Unknown')

    @skip("passed")
    def test_generate_strain_data(self):
        """Read GBK file and return correct strain data"""
        p_aeruginosa = self.importer.generate_strain_data('/mnt/data2/Bacteriocins/genomes/Pseudomonas_aeruginosa_PAO1.gb', 'PAO1', '')
        self.assertEqual(p_aeruginosa.strain_id, 'PAO1')
        self.assertEqual(p_aeruginosa.full_name, 'Pseudomonas aeruginosa PAO1')
        self.assertEqual(p_aeruginosa.taxon.taxonomy_id, '208964')
        self.assertEqual(p_aeruginosa.order, 'Pseudomonadales')
        isolate = self.importer.generate_strain_data('/mnt/data2/ENIGMA/genome_files/genbank/DP16D-E2.genome.gbff.gz', 'DP16D-E2', '')
        self.assertEqual(isolate.strain_id, 'DP16D-E2')
        self.assertEqual(isolate.full_name, 'Environmental isolate DP16D-E2')
        self.assertEqual(isolate.taxon.taxonomy_id, '48479')
        self.assertEqual(isolate.order, 'Unknown')

    @skip("fix perl version incompatibility later")
    def test_importer(self):
        lines = []
        with open('/mnt/data2/CGCMS/test_data/test_genome_import.txt', 'r') as infile:
            for line in infile:
                lines.append(line.rstrip('\n\r'))
        result = self.importer.import_genomes(lines)
        self.assertEqual(len(self.importer.inputgenomes), 1)
        self.assertEqual(result, 'Done!')


class BrowserTestCase(StaticLiveServerTestCase):
    fixtures = ['testdata.json']
    
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.selenium = WebDriver()
        cls.selenium.implicitly_wait(10)

    @classmethod
    def tearDownClass(cls):
        cls.selenium.quit()
        super().tearDownClass()
        
    #@skip("test passed")
    def test_search_page(self):
        print('Search in annotations')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('name', 'annotation_query')
        text_input.send_keys('tetracycline')
        self.selenium.find_element('name', 'annotation_query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert 'I6K11_00030' in self.selenium.page_source
        assert 'LC317985_1_410_-1' not in self.selenium.page_source
        
        print('Search in genes')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('id', 'gene-query')
        text_input.send_keys('tetracycline')
        self.selenium.find_element('id', 'gene-query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert 'LC317985_1_410_-1' in self.selenium.page_source
        assert 'LC317985_3092_3715_-1' not in self.selenium.page_source

        print('Search in KEGG orthologs')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('id', 'ko-query')
        text_input.send_keys('gltS')
        self.selenium.find_element('id', 'ko-query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert 'K03312' in self.selenium.page_source
        assert 'K00348' not in self.selenium.page_source
        
        print('Search in enzymes')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('id', 'ec-query')
        text_input.send_keys('1.6.5.8')
        self.selenium.find_element('id', 'ec-query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert '7.2.1.1' in self.selenium.page_source
        
        print('Search in transporters')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('id', 'tc-query')
        text_input.send_keys('iron')
        self.selenium.find_element('id', 'tc-query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert '2.A.108.1' in self.selenium.page_source
        assert '2.A.1.2.38' not in self.selenium.page_source
        
        print('Search in COG classes')
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('id', 'cog-query')
        text_input.send_keys('K')
        self.selenium.find_element('id', 'cog-query').send_keys(Keys.RETURN)
        element = WebDriverWait(self.selenium, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element.text)
        assert 'Transcription' in self.selenium.page_source
        assert 'Energy' not in self.selenium.page_source

    def test_search_genome_page(self):
        print('Search for unnamed in genomes list')
        self.selenium.get('%s%s' % (self.live_server_url, '/genomes/'))
        text_input = self.selenium.find_element('name', 'query')
        text_input.send_keys('unnamed')
        self.selenium.find_element('name', 'query').send_keys(Keys.RETURN)
        with self.wait_for_page_load(timeout=10):
            element = self.selenium.find_element(By.CLASS_NAME, 'table-wrapper')
            print(element.text)
            assert 'plasmid_unnamed1' in self.selenium.page_source
            assert 'Ecoli_plasmid_p15' not in self.selenium.page_source

    @contextmanager
    def wait_for_page_load(self, timeout=30):
        old_page = self.selenium.find_element(By.TAG_NAME, "html")
        yield WebDriverWait(self.selenium, timeout).until(
            EC.staleness_of(old_page)
        )
