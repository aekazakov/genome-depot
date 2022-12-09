import time
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
        
    def test_search_page(self):
        self.selenium.get('%s%s' % (self.live_server_url, '/textsearch/'))
        text_input = self.selenium.find_element('name', 'annotation_query')
        text_input.send_keys('tetracycline')
        self.selenium.find_element('name', 'annotation_query').send_keys(Keys.RETURN)
        # time.sleep(10)
        element = WebDriverWait(self.selenium, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "table-wrapper")))
        print(element)
        assert 'I6K11_00030' in self.selenium.page_source
