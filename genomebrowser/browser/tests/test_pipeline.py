import os
import hashlib
from Bio import GenBank
from unittest import skip

from django.test import TestCase
#from django.db.utils import IntegrityError

from genomebrowser.settings import BASE_DIR
from browser.pipeline.genome_import import Importer
from browser.models import Taxon
from browser.models import Strain
from browser.models import Sample

# Create your tests here.


class ImporterTestCase(TestCase):
    '''
        Testing genome import pipeline
    '''
    
    #fixtures = ['testdata.json']
    fixtures = ['minigenomes.testdata.json']
    
    @classmethod
    def setUpTestData(cls):
        cls.importer = Importer()

    def test_sanitize_genome_id(self):
        '''
            Tests genome id conversion
        '''
        print('Testing sanitize_genome_id function')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli_CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli.CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli..CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli_α_CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli___CFT073'), 'E_coli_CFT073')
 
    def test_load_genome_list(self):
        '''
            Read tab-separated file with a list of genomes to import
        '''
        print('Testing load_genome_list function')
        in_file = '../testdata/import_minigenomes.txt'
        self.importer.load_genome_list(in_file)
        self.assertEqual(self.importer.inputgenomes['E_coli_CFT073']['strain'], 'CFT073') # first line imported
        self.assertEqual(self.importer.inputgenomes['E_coli_BW2952']['strain'], 'BW2952') # last line imported
        
    def test_check_genomes(self):
        '''
            Test the function checking if genome names are unique
        '''
        in_file = '../testdata/import_minigenomes.txt'
        self.importer.load_genome_list(in_file)
        self.assertRaises(ValueError,self.importer.check_genomes)

    def create_tag(self):
        '''
            Test the function creating genome tags
        '''
        current_date = str(timezone.localdate(timezone.now()))
        tag_name = 'imported_' + str(timezone.localdate(timezone.now()))
        tag_description = 'Genomes imported on ' + current_date
        self.importer.create_tag()
        self.assertEqual(Tag.objects.get(name=tag_name).description, tag_description)

    def test_get_taxonomic_order(self):
        '''
            Test if taxonomic order is correctly identified
        '''
        print('Testing find_taxonomic_order function')
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

    def test_get_gbk_organism(self):
        '''
            Test if organism name is correctly identified
        '''
        in_file = '../testdata/E_coli_BW2952.100000.gbk'
        out = self.importer.get_gbk_organism(in_file)
        self.assertEqual(out['organism'], 'Escherichia coli BW2952')
        
    def test_create_parent_taxa(self):
        print('Testing create_parent_taxa function')
        taxonomy_id = '817'
        self.assertRaises(Taxon.DoesNotExist, Taxon.objects.get, taxonomy_id='816')
        self.importer.create_parent_taxa(taxonomy_id)
        self.assertEqual(Taxon.objects.get(taxonomy_id='816').taxonomy_id, '816')

    def test_generate_strain_data(self):
        '''
            Read GBK file and return correct strain data
        '''
        print('Testing generate_strain_data function')
        # Test NCBI genome
        out = self.importer.generate_strain_data(
            '../testdata/E_coli_BW2952.100000.gbk',
            'BW2952',
            ''
            )
        self.assertEqual(out.strain_id, 'BW2952')
        self.assertEqual(out.full_name, 'Escherichia coli BW2952')
        self.assertEqual(out.taxon.taxonomy_id, '595496')
        self.assertEqual(out.order, 'Enterobacterales')

    def test_prepare_strain_data(self):
        '''
            Test if a strain record is created
        '''
        self.importer.inputgenomes['Test_genome']['strain'] = 'Test'
        self.importer.inputgenomes['Test_genome']['gbk'] = '../testdata/E_coli_BW2952.100000.gbk'
        self.importer.prepare_strain_data()
        self.assertEqual(self.importer.strain_instances['Test'].strain_id, 'Test')
        
    def test_prepare_sample_data(self):
        '''
            Test if a sample record is created
        '''
        in_file = '../testdata/import_minigenomes.txt'
        # Prepare the data
        self.importer.inputgenomes['Test']['sample'] = 'Test'
        self.importer.prepare_sample_data()
        self.assertEqual(self.importer.sample_instances['Test'].sample_id, 'Test')
        
    def test_parse_location(self):
        '''
            Test the parse_location function
        '''
        in_file = '../testdata/import_minigenomes.txt'
        # Prepare the data
        location = self.importer.parse_location('190..255', 100000)
        self.assertEqual(location[0], 190)
        self.assertEqual(location[1], 255)
        self.assertEqual(location[2], 1)
        location = self.importer.parse_location('complement(5683..6459)', 100000)
        self.assertEqual(location[0], 5683)
        self.assertEqual(location[1], 6459)
        self.assertEqual(location[2], -1)
        location = self.importer.parse_location('complement(19658..>19795)', 100000)
        self.assertEqual(location[0], 19658)
        self.assertEqual(location[1], 19795)
        self.assertEqual(location[2], -1)
        
    def test_process_feature(self):
        '''
            Test the process_feature function
        '''
        in_file = '../testdata/E_coli_BW2952.100000.gbk'
        gbk_handle = open(in_file, 'r')
        gbk_record = GenBank.read(gbk_handle)
        feature=gbk_record.features[2]
        gbk_handle.close()
        print(feature)
        feature_uid, locus_tag = self.importer.process_feature(feature, 'E_coli_BW2952', 'NC_012759', [], 100000)
        self.assertEqual(locus_tag, 'BWG_RS00005')
        self.assertEqual(self.importer.gene_data[feature_uid]['locus_tag'], 'BWG_RS00005')
        
    def test_export_contigs(self):
        '''
            Test the export_contigs function
        '''
        self.importer.export_contigs()
        self.assertEqual(len(self.importer.staticfiles[self.importer.config['cgcms.search_db_dir']]), 1)


    @skip("this is a very long test")
    def test_genome_import_pipeline(self):
        '''
            This test runs the entire genome import pipeline for three minigenomes
        '''
        print('Testing the entire genome import pipeline')
        lines = []
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), 'r') as infile:
            for line in infile:
                lines.append(line.rstrip('\n\r'))
                print(line)
        result = self.importer.import_genomes(lines)
        self.assertEqual(len(self.importer.inputgenomes), 3)
        self.assertEqual(result, 'Done!')

