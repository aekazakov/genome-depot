import os
import gzip
import shutil
from Bio import GenBank
from unittest import skip

from django.test import TestCase, TransactionTestCase

from genomebrowser.settings import BASE_DIR
from browser.models import Config
from browser.models import Kegg_reaction
from browser.models import Genome
from browser.models import Tag

from browser.util import export_genome
from browser.util import export_genomes
from browser.util import download_ncbi_assembly
from browser.util import delete_all_data
from browser.util import delete_all_genomes
from browser.util import delete_genome
from browser.util import delete_genomes
from browser.util import generate_static_files
from browser.util import import_config
from browser.util import export_config
from browser.util import regenerate_jbrowse_files
from browser.util import recreate_search_databases
from browser.util import update_tags
# Create your tests here.


class UtilTestCase(TestCase):
    '''
        Testing utils
    '''
    
    fixtures = ['minigenomes.testdata.json']
    
    def setUp(self) -> None:
        self.temp_dir = Config.objects.get(param='core.temp_dir').value

    def test_export_genome(self):
        '''
            Tests export_genome function
        '''
        genome_id = 'E_coli_BW2952'
        genome = Genome.objects.get(name = genome_id)
        out_file = os.path.join(self.temp_dir, 'test.gbk') 
        with open(out_file, 'w') as output_buffer:
            export_genome(genome, output_buffer)

        gbk_handle = open(out_file, 'r')
        gbk_record = GenBank.read(gbk_handle)
        gbk_handle.close()
        self.assertEqual(gbk_record.accession[0], 'NC_012759')
 
    def test_export_genomes(self):
        '''
            Test genome export
        '''
        genome = 'E_coli_BW2952'
        export_genomes(self.temp_dir, [genome, ])
        gbk_handle = gzip.open(os.path.join(self.temp_dir, genome + '.gbff.gz'), 'rt')
        gbk_record = GenBank.read(gbk_handle)
        gbk_handle.close()
        self.assertEqual(gbk_record.accession[0], 'NC_012759')

    def test_download_ncbi_assembly(self):
        '''
            Test assembly download
        '''
        assembly_id = 'GCF_020892305.1'
        email = 'aekazakov@lbl.gov'
        outfile = download_ncbi_assembly(assembly_id, email, self.temp_dir)
        
        gbk_handle = gzip.open(outfile, 'rt')
        gbk_record = GenBank.read(gbk_handle)
        gbk_handle.close()
        self.assertEqual(gbk_record.accession[0], 'NC_073066')

    #@skip("skip for now")
    def test_delete_all_data(self):
        '''
            Test database wipe out
        '''
        delete_all_data(confirm=False)
        self.assertEqual(Kegg_reaction.objects.all().count(), 0)

    #@skip("skip for now")
    def test_delete_all_genomes(self):
        '''
            Test deleting all_genomes
        '''
        delete_all_genomes(confirm=False)
        self.assertEqual(Genome.objects.all().count(), 0)

    def test_delete_genome(self):
        '''
            Test deleting one genome
        '''
        genome = 'E_coli_BW2952'
        delete_genome(genome)
        self.assertEqual(Genome.objects.filter(name=genome).count(), 0)

    def test_delete_genomes(self):
        '''
            Test deleting genomes by list
        '''
        genome = 'E_coli_BW2952'
        test_file = os.path.join(self.temp_dir, 'test.txt') 
        with open(test_file, 'w') as outfile:
            outfile.write('E_coli_BW2952\n')
        delete_genomes(test_file)
        self.assertEqual(Genome.objects.filter(name=genome).count(), 0)

    def test_generate_static_files(self):
        '''
            Test re-creating jbrowse files for a list of genomes from file
        '''
        genome = 'E_coli_BW2952'
        json_dir = os.path.join(Config.objects.get(param='core.json_dir').value, genome)
        if os.path.exists(json_dir):
            shutil.rmtree(json_dir)
        test_file = os.path.join(self.temp_dir, 'test.txt') 
        with open(test_file, 'w') as outfile:
            outfile.write('../testdata/E_coli_BW2952.100000.gbk\tE_coli_BW2952\tBW2952\t\thttps://www.ncbi.nlm.nih.gov/assembly/GCF_000022345.1\tNCBI:GCF_000022345.1')
        generate_static_files(test_file)
        self.assertTrue(os.path.exists(json_dir))

    def test_import_config(self):
        '''
            Test config import
        '''
        test_file = os.path.join(self.temp_dir, 'test.txt') 
        with open(test_file, 'w') as outfile:
            outfile.write('test_param=test_value\n')
        import_config(test_file)
        self.assertEqual(Config.objects.get(param='test_param').value, 'test_value')

    def test_export_config(self):
        '''
            Test config export
        '''
        test_file = os.path.join(self.temp_dir, 'test.txt') 
        if os.path.exists(test_file):
            os.remove(test_file)
        export_config(test_file)
        with open(test_file, 'r') as infile:
            param, value = infile.readline().rstrip('\n').split('=')
            self.assertEqual(Config.objects.get(param=param).value, value)

    def test_regenerate_jbrowse_files(self):
        '''
            Test re-creating jbrowse files for one genome
        '''
        genome = 'E_coli_BW2952'
        json_dir = os.path.join(Config.objects.get(param='core.json_dir').value, genome)
        if os.path.exists(json_dir):
            shutil.rmtree(json_dir)
        regenerate_jbrowse_files(genome)
        self.assertTrue(os.path.exists(json_dir))

    def test_recreate_search_databases(self):
        '''
            Test re-creating search databases
        '''
        os.remove(os.path.join(Config.objects.get(param='core.search_db_dir').value, 'blast_nucl.ndb'))
        os.remove(os.path.join(Config.objects.get(param='core.search_db_dir').value, 'blast_prot.pdb'))
        recreate_search_databases()
        self.assertTrue(os.path.exists(os.path.join(Config.objects.get(param='core.search_db_dir').value, 'blast_nucl.ndb')))
        self.assertTrue(os.path.exists(os.path.join(Config.objects.get(param='core.search_db_dir').value, 'blast_prot.pdb')))

    def test_update_tags(self):
        '''
            Test updating genome tags
        '''
        genome = 'E_coli_BW2952'
        test_file = os.path.join(self.temp_dir, 'test.txt') 
        with open(test_file, 'w') as outfile:
            outfile.write('../testdata/E_coli_BW2952.100000.gbk\tE_coli_BW2952\tBW2952\t\thttps://www.ncbi.nlm.nih.gov/assembly/GCF_000022345.1\tNCBI:GCF_000022345.1')
        update_tags(test_file, 'test_tag_1,test_tag_2')
        self.assertTrue(Tag.objects.filter(name__icontains='test_tag').count(),2)
