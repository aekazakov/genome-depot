import time
from unittest import skip
from contextlib import contextmanager
from django.test import TransactionTestCase
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from browser.models import *
from browser.dataimport.importer import Importer
# Create your tests here.

class BrowserTestCase(TransactionTestCase):
    def setUp(self):
        pass
        
    def test_config(self):
        config = Config(param='parameter name', value='parameter value')
        self.assertEqual(config.param, 'parameter name')
        config.save()
        config_saved = Config.objects.get(param='parameter name')
        self.assertEqual(config_saved.value, 'parameter value')

    def test_taxon(self):
        taxon = Taxon(taxonomy_id='10', eggnog_taxid='10', rank='genus', parent_id='1706371', name='Cellvibrio')
        self.assertEqual(taxon.name, 'Cellvibrio')
        self.assertEqual(taxon.taxonomy_id, '10')
        self.assertNotEqual(taxon.taxonomy_id, 10)
        taxon.save()
        taxon_saved = Taxon.objects.get(taxonomy_id='10')
        self.assertEqual(taxon_saved.name, 'Cellvibrio')
        self.assertEqual(taxon_saved.taxonomy_id, '10')
        
    def test_strain(self):
        taxon = Taxon(taxonomy_id='666685', eggnog_taxid='666685', rank='species', parent_id='75309', name='Rhodanobacter denitrificans')
        taxon.save()
        strain = Strain(strain_id='FW104-10B01', full_name='Rhodanobacter denitrificans str. FW104-10B01', order='Xanthomonadales', taxon=taxon)
        strain.save()
        self.assertEqual(strain.strain_id, 'FW104-10B01')
        self.assertEqual(strain.taxon.name, 'Rhodanobacter denitrificans')
        strain_saved = Strain.objects.get(strain_id='FW104-10B01')
        self.assertEqual(strain_saved.strain_id, 'FW104-10B01')
        self.assertEqual(strain_saved.taxon.name, 'Rhodanobacter denitrificans')
        
    def test_sample(self):
        sample = Sample(sample_id='FW106-02', full_name='FW106 groundwater metagenome')
        sample.description = 'Metagenomic sample from groundwater of FW106 well, site Y-12 West, collection date 2014-06-09, 0.2 micron filter'
        sample.save()
        self.assertEqual(sample.sample_id, 'FW106-02')
        sample_saved = Sample.objects.get(sample_id='FW106-02')
        self.assertEqual(sample_saved.sample_id, 'FW106-02')

    def test_strain_metadata(self):
        taxon = Taxon(taxonomy_id='666685', eggnog_taxid='666685', rank='species', parent_id='75309', name='Rhodanobacter denitrificans')
        taxon.save()
        strain = Strain(strain_id='FW104-10B01', full_name='Rhodanobacter denitrificans str. FW104-10B01', order='Xanthomonadales', taxon=taxon)
        strain.save()
        strain_metadata = Strain_metadata(strain=strain, source='abcdef', url='https://nih.gov/', key='ghijkl', value='mnopqr')
        self.assertEqual(strain_metadata.source, 'abcdef')
        self.assertEqual(strain_metadata.strain.strain_id, 'FW104-10B01')
        strain_metadata.save()
        strain_metadata_saved = Strain_metadata.objects.filter(source='abcdef')
        self.assertEqual(len(list(strain_metadata_saved.all())), 1)
        self.assertEqual(list(strain_metadata_saved.all())[0].strain.strain_id, 'FW104-10B01')


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


