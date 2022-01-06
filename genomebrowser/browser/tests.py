from django.test import TestCase

from browser.dataimport.importer import Importer
# Create your tests here.

class ImporterTestCase(TestCase):
    def setUp(self):
        self.importer = Importer('/mnt/data2/ENIGMA/genomes/dev/config.ini')

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
        p_aeruginosa = self.importer.generate_strain_data('/mnt/data2/Bacteriocins/genomes/Pseudomonas_aeruginosa_PAO1.gb', 'PAO1')
        self.assertEqual(p_aeruginosa.strain_id, 'PAO1')
        self.assertEqual(p_aeruginosa.full_name, 'Pseudomonas aeruginosa PAO1')
        self.assertEqual(p_aeruginosa.taxon.taxonomy_id, '208964')
        self.assertEqual(p_aeruginosa.order, 'Pseudomonadales')
        isolate = self.importer.generate_strain_data('/mnt/data2/ENIGMA/genome_files/genbank/DP16D-E2.genome.gbff.gz', 'DP16D-E2')
        self.assertEqual(isolate.strain_id, 'DP16D-E2')
        self.assertEqual(isolate.full_name, 'Environmental isolate DP16D-E2')
        self.assertEqual(isolate.taxon.taxonomy_id, '48479')
        self.assertEqual(isolate.order, 'Unknown')

    def test_importer(self):
        self.importer.import_genomes('/mnt/data2/ENIGMA/genomes/dev/genomes.txt')
        self.assertEqual(len(self.importer.inputgenomes), 6)
