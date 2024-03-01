import os
import hashlib
from Bio import GenBank
from Bio.SeqIO.FastaIO import SimpleFastaParser
from unittest import skip
from collections import defaultdict

from django.test import TestCase, TransactionTestCase
from django.utils import timezone
#from django.db.utils import IntegrityError

from genomebrowser.settings import BASE_DIR
from browser.models import Annotation
from browser.models import Taxon
from browser.models import Genome
from browser.models import Strain
from browser.models import Sample
from browser.models import Strain_metadata
from browser.models import Sample_metadata
from browser.models import Regulon
from browser.models import Tag

from browser.pipeline.genome_import import Importer
from browser.pipeline.annotate import Annotator
from browser.pipeline.util import export_proteins_bygenome,export_nucl_bygenome,export_proteins
from browser.pipeline.taxonomy import update_taxonomy

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
        cls.annotator= Annotator()

    def test_sanitize_genome_id(self):
        '''
            Tests genome id conversion
        '''
        print('Testing sanitize_genome_id function')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli_CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli!CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli.CFT073'), 'E_coli.CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli_α_CFT073'), 'E_coli_CFT073')
        self.assertEqual(self.importer.sanitize_genome_id('E_coli___CFT073'), 'E_coli_CFT073')
 
    def test_load_genome_list(self):
        '''
            Read tab-separated file with a list of genomes to import
        '''
        print('Testing load_genome_list function')
        in_file = '../testdata/import_minigenomes.txt'
        self.importer.inputgenomes = defaultdict(dict)
        self.importer.load_genome_list(in_file)
        self.assertEqual(self.importer.inputgenomes['E_coli_CFT073']['strain'], 'CFT073') # first line imported
        self.assertEqual(self.importer.inputgenomes['E_coli_BW2952']['strain'], 'BW2952') # last line imported
        
    def test_check_genomes(self):
        '''
            Test the function checking if genome names are unique
        '''
        print('Testing check_genomes function')
        in_file = '../testdata/import_minigenomes.txt'
        self.importer.load_genome_list(in_file)
        self.assertRaises(ValueError,self.importer.check_genomes)

    def test_create_tag(self):
        '''
            Test the function creating genome tags
        '''
        print('Testing create_tag function')
        current_date = str(timezone.localdate(timezone.now()))
        tag_name = 'imported_' + str(timezone.localdate(timezone.now()))
        tag_description = 'Genomes imported on ' + current_date
        self.importer.create_tag()
        self.assertEqual(Tag.objects.get(name=tag_name).description, tag_description)

    def test_get_taxonomic_order(self):
        '''
            Test if taxonomic order is correctly identified
        '''
        print('Testing get_taxonomic_order function')
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
        print('Testing get_gbk_organism function')
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
        print('Testing prepare_strain_data function')
        self.importer.inputgenomes['Test_genome']['strain'] = 'Test'
        self.importer.inputgenomes['Test_genome']['gbk'] = '../testdata/E_coli_BW2952.100000.gbk'
        self.importer.prepare_strain_data()
        self.assertEqual(self.importer.strain_instances['Test'].strain_id, 'Test')
        
    def test_prepare_sample_data(self):
        '''
            Test if a sample record is created
        '''
        print('Testing prepare_sample_data function')
        in_file = '../testdata/import_minigenomes.txt'
        # Prepare the data
        self.importer.inputgenomes['Test']['sample'] = 'Test'
        self.importer.prepare_sample_data()
        self.assertEqual(self.importer.sample_instances['Test'].sample_id, 'Test')
        
    def test_parse_location(self):
        '''
            Test the parse_location function
        '''
        print('Testing parse_location function')
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
        print('Testing process_feature function')
        in_file = '../testdata/E_coli_BW2952.100000.gbk'
        gbk_handle = open(in_file, 'r')
        gbk_record = GenBank.read(gbk_handle)
        feature=gbk_record.features[2]
        gbk_handle.close()
        print(feature)
        feature_uid, locus_tag = self.importer.process_feature(feature, 'E_coli_BW2952', 'NC_012759', [], 100000)
        self.assertEqual(locus_tag, 'BWG_RS00005')
        self.assertEqual(self.importer.gene_data[feature_uid]['locus_tag'], 'BWG_RS00005')
        
    #@skip("skip for now")
    def test_export_contigs(self):
        '''
            Test the export_contigs function
        '''
        print('Testing export_contigs function')
        self.importer.export_contigs()
        self.assertEqual(len(self.importer.staticfiles[self.importer.config['cgcms.search_db_dir']]), 1)

    #@skip("skip for now")
    def test_export_proteins(self):
        '''
            Test the export_proteins function in util
            It only checks the output file is not empty
        '''
        print('Testing export_proteins function')
        genome_id = 155
        test_file = os.path.join(self.importer.config['cgcms.temp_dir'], 'testfile.faa')

        export_proteins([genome_id,], test_file)

        with open(test_file, 'r') as infile:
            for title, seq in SimpleFastaParser(infile):
                self.assertNotEqual(title, '')
                self.assertNotEqual(seq, '')

    #@skip("skip for now")
    def test_export_proteins_bygenome(self):
        '''
            Test the export_proteins_bygenome function in util
        '''
        print('Testing export_proteins_bygenome function')
        genome_name = 'E_coli_BW2952'
        genome = Genome.objects.get(name=genome_name)
        out_dir = self.importer.config['cgcms.temp_dir']
        
        export_proteins_bygenome({genome_name:genome.gbk_filepath}, out_dir)

        out_file = os.path.join(out_dir, str(genome.id) + '.faa')
        with open(out_file, 'r') as infile:
            for title, seq in SimpleFastaParser(infile):
                self.assertEqual(title[:3], 'BWG')
                self.assertNotEqual(seq, '')
        #os.remove(out_file)

    #@skip("skip for now")
    def test_export_nucl_bygenome(self):
        '''
            Test the export_nucl_bygenome function in util
            It only checks the output file is not empty
        '''
        print('Testing export_nucl_bygenome function')
        genome_name = 'E_coli_BW2952'
        genome = Genome.objects.get(name=genome_name)
        test_file = os.path.join(self.importer.config['cgcms.temp_dir'], 'testfile.fna')

        out_dir = self.importer.config['cgcms.temp_dir']
        
        export_nucl_bygenome({genome_name:genome.gbk_filepath}, out_dir)

        out_file = os.path.join(out_dir, str(genome.id) + '.fna')
        with open(out_file, 'r') as infile:
            for title, seq in SimpleFastaParser(infile):
                self.assertNotEqual(title, '')
                self.assertNotEqual(seq, '')
        #os.remove(out_file)

    def test_add_regulons(self):
        '''
            Test the add_regulons function in Annotator
        '''
        print('Testing add_regulons function')
        genome_name = 'E_coli_BW2952'
        regulon_name = 'Test'
        reg_gene_ids = 'BWG_RS00020'
        target_gene_id = 'BWG_RS00005'
        contig = 'NC_012759'
        start = '187'
        end = '189'
        strand = '1'
        sequence = 'TCC'
        # Prepare the data
        test_data = []
        test_line = '\t'.join([regulon_name, genome_name, reg_gene_ids, target_gene_id, contig, \
            start, end, strand, sequence]) + '\n'
        test_data.append(test_line)
        #Existing regulon
        test_line = '\t'.join(['AraC', genome_name, reg_gene_ids, target_gene_id, contig, \
            start, end, strand, sequence]) + '\n'
        test_data.append(test_line)
        #Non-existant gene
        test_line = '\t'.join([regulon_name, genome_name, reg_gene_ids, 'fakegene', contig, \
            start, end, strand, sequence]) + '\n'
        test_data.append(test_line)
        
        self.annotator.add_regulons(test_data)
        saved_regulon = Regulon.objects.get(name='Test')
        self.assertEqual(saved_regulon.genome.name, genome_name)

    def test_add_custom_annotations(self):
        '''
            Test the add_custom_annotations function in Annotator
        '''
        print('Testing add_custom_annotations function')
        genome_id = 'E_coli_BW2952'
        # Prepare the data
        test_file = os.path.join(self.importer.config['cgcms.temp_dir'], 'testfile_annot.tsv')
        with open(test_file, 'w') as outfile:
            outfile.write('#comment line\n')
            outfile.write('\t'.join([
                'BWG_RS00005', 
                genome_id,
                'Test source',
                'Test_URL',
                'Test',
                'Test_value',
                'Test description'
                ]))
        
        self.annotator.add_custom_annotations(test_file)
        os.remove(test_file)
        saved_annotation = Annotation.objects.get(key='Test')
        self.assertEqual(saved_annotation.note, 'Test description')
        
    def test_add_strain_metadata(self):
        '''
            Test the add_custom_annotations function in Annotator
        '''
        print('Testing add_strain_metadata function')
        strain_id = 'CFT073'
        # Prepare the data
        test_file = os.path.join(self.importer.config['cgcms.temp_dir'], 'testfile_meta.tsv')
        with open(test_file, 'w') as outfile:
            outfile.write('#comment line\n')
            outfile.write('\t'.join([
                strain_id,
                'Test source',
                'Test_URL',
                'Test',
                'Test_value'
                ]))
        
        self.annotator.add_strain_metadata(test_file)
        os.remove(test_file)
        saved_metadata = Strain_metadata.objects.get(key='Test')
        self.assertEqual(saved_metadata.value, 'Test_value')

    def test_update_strain_metadata(self):
        '''
            Test the update_strain_metadata function in Annotator
        '''
        print('Testing update_strain_metadata function')
        strain_id = 'BW2952'
        test_file = 'test_strain_metadata.xlsx'
        self.annotator.update_strain_metadata(xlsx_path=test_file)
        saved_metadata = Strain_metadata.objects.get(key='Phylogenetic Order')
        self.assertEqual(saved_metadata.strain.strain_id, strain_id)
        self.assertEqual(saved_metadata.value, 'Enterobacterales')

    def test_add_sample_metadata(self):
        '''
            Test the add_sample_metadata function in Annotator
        '''
        print('Testing add_sample_metadata function')
        # Prepare the data
        test_line = '\t'.join(['test_sample', 'Test source', 'Test_URL', 'Test', 'Test_value']) + '\n'
        test_data = [test_line, ]
        
        self.annotator.add_sample_metadata(test_data)
        saved_metadata = Sample_metadata.objects.get(key='Test')
        self.assertEqual(saved_metadata.value, 'Test_value')

    def test_update_genome_descriptions(self):
        '''
            Test the update_genome_descriptions function in Annotator
        '''
        print('Testing update_genome_descriptions function')
        genome_id = 'E_coli_BW2952'
        test_description = 'Updated genome description'
        # Prepare the data
        test_file = os.path.join(self.importer.config['cgcms.temp_dir'], 'testfile_genome.tsv')
        with open(test_file, 'w') as outfile:
            outfile.write('#comment line\n')
            outfile.write('\t'.join([
                genome_id,
                test_description
                ]) + '\n')
        
        self.annotator.update_genome_descriptions(test_file)
        os.remove(test_file)
        genome = Genome.objects.get(name=genome_id)
        self.assertEqual(genome.description, test_description)

    def test_update_sample_descriptions(self):
        '''
            Test the update_sample_descriptions function in Annotator
        '''
        print('Testing update_sample_descriptions function')
        sample_id = 'test_sample'
        test_fullname = 'Updated sample fullname'
        test_description = 'Updated sample description'
        # Prepare the data
        test_line = '\t'.join([sample_id, test_fullname, test_description]) + '\n'
        test_data = [test_line, ]

        self.annotator.update_sample_descriptions(test_data)
        
        sample = Sample.objects.get(sample_id=sample_id)
        self.assertEqual(sample.full_name, test_fullname)
        self.assertEqual(sample.description, test_description)

    #@skip("skip for now")
    def test_update_taxonomy(self):
        '''
            Test the update_taxonomy function
        '''
        print('Testing update_taxonomy function')
        taxon = Taxon.objects.get(taxonomy_id="1224")
        self.assertEqual(taxon.name, "Proteobacteria")
        update_taxonomy()
        updated_taxon = Taxon.objects.get(taxonomy_id="1224")
        self.assertEqual(updated_taxon.name, "Pseudomonadota")


class PipelineTestCase(TransactionTestCase):
    '''
        Testing pipeline functions that require TransactionTestCase
    '''
    fixtures = ['minigenomes.testdata.json']
    
    @classmethod
    def setUp(self):
        self.importer = Importer()
        self.importer.config['cgcms.generate_names_command'] = '/mnt/data/work/CGCMS/external_tools/jbrowse/bin/generate-names.pl'
        self.importer.config['plugins.macsyfinder.conda_env'] = 'cgcms-macsyfinder'
        self.importer.config['plugins.macsyfinder.model'] = 'TXSScan'
        self.importer.config['plugins.macsyfinder.models_dir'] = '/mnt/data/work/CGCMS/external_refdata/macsyfinder/data'
        self.importer.config['plugins.defensefinder.conda_env'] = 'cgcms-defensefinder'
        self.importer.config['plugins.defensefinder.defensefinder_models_dir'] = '/mnt/data/work/CGCMS/external_refdata/defensefinder/data'
        self.importer.config['plugins.hmmsearch_tigrfam.conda_env'] = 'cgcms_hmmsearch'
        self.importer.config['plugins.hmmsearch_tigrfam.hmmsearch_command'] = 'hmmsearch'
        self.importer.config['plugins.hmmsearch_tigrfam.display_name'] = 'TIGRFAM database'
        self.importer.config['plugins.hmmsearch_tigrfam.hmm_lib'] = '/mnt/data/work/CGCMS/external_refdata/tigrfam/TIGRFAM.HMM'
        self.importer.config['plugins.hmmsearch_tigrfam.ref_data'] = '/mnt/data/work/CGCMS/external_refdata/tigrfam/ref_tigrfam.txt'
        self.importer.config['plugins.hmmsearch_pfam.conda_env'] = 'cgcms_hmmsearch'
        self.importer.config['plugins.hmmsearch_pfam.hmmsearch_command'] = 'hmmsearch'
        self.importer.config['plugins.hmmsearch_pfam.display_name'] = 'Pfam database'
        self.importer.config['plugins.hmmsearch_pfam.hmm_lib'] = '/mnt/data/work/CGCMS/external_refdata/pfam/Pfam-A.hmm'
        self.importer.config['plugins.hmmsearch_pfam.ref_data'] = '/mnt/data/work/CGCMS/external_refdata/pfam/ref_pfam.txt'
        self.importer.config['plugins.gapmind.enabled'] = '1'
        self.importer.config['plugins.gapmind.conda_env'] = 'cgcms-gapmind'
        self.importer.config['plugins.gapmind.gapmind_dir'] = '/mnt/data/work/CGCMS/external_tools/PaperBLAST'
        self.importer.config['plugins.gapmind.threads'] = '8'
        self.importer.config['plugins.phispy.conda_env'] = 'cgcms-phispy'
        self.importer.config['plugins.fama.fama_dir'] = '/mnt/data/work/CGCMS/external_tools/fama/py'
        self.importer.config['plugins.fama.fama_config'] = '/mnt/data/work/CGCMS/external_tools/fama/config.ini'
        self.importer.config['plugins.antismash.conda_env'] = 'cgcms-antismash'
        self.importer.config['plugins.antismash.antismash_cmd'] = 'antismash'
        self.importer.config['plugins.amrfinder.threads'] = '8'
        self.importer.config['plugins.amrfinder.conda_env'] = 'cgcms-amrfinder'
        self.importer.config['plugins.amrfinder.display_name'] = 'AMRFinderPlus'

        self.importer.config['cgcms.poem_dir'] = '/mnt/data/work/sandbox/poem_py3k/POEM_py3k'

        self.annotator= Annotator()

    #@skip("skip for now")
    def test_run_external_tools(self):
        '''
            Test the run_external_tools function in Annotator.
            Gapmind should find 5 hits in the E_coli_BW2952 minigenome
        '''
        print('Testing run_external_tools function')
        self.annotator.config['plugins.gapmind.enabled'] = '1'
        self.annotator.config['plugins.gapmind.conda_env'] = 'cgcms-gapmind'
        self.annotator.config['plugins.gapmind.gapmind_dir'] = '/mnt/data/work/CGCMS/external_tools/PaperBLAST'
        self.annotator.config['plugins.gapmind.threads'] = '8'
        self.annotator.config['plugins.gapmind.display_name'] = 'GapMind'
        
        genome_id = 'E_coli_BW2952'
        genome = Genome.objects.get(name=genome_id)
        test_plugin = 'gapmind'
        self.annotator.run_external_tools({genome_id:genome.gbk_filepath})
        saved_annotations = list(Annotation.objects.filter(gene_id__genome__name=genome_id, source='GapMind'))
        print('Annotations:', saved_annotations)
        self.assertEqual(len(saved_annotations), 5)

    #@skip("this is a very long test")
    def test_genome_import_pipeline(self):
        '''
            This test runs the entire genome import pipeline for three minigenomes
        '''
        print('Testing genome_import_pipeline function')
        Genome.objects.get(name = 'E_coli_CFT073').delete()
        Protein.objects.filter(gene=None).delete()
        print('Testing the entire genome import pipeline')
        print('Genomes before import', str(Genome.objects.values_list('name', flat=True)))
        lines = []
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                row[1] = row[1] + '.test'
                lines.append('\t'.join(row))
                print('Test input ' + '\t'.join(row))
        result = self.importer.import_genomes(lines)
        self.assertEqual(len(self.importer.inputgenomes), 3)
        self.assertEqual(result, 'Done!')
        print('Genomes after import', str(Genome.objects.values_list('name', flat=True)))
        test_genome = Genome.objects.get(name = 'E_coli_CFT073.test')
        self.assertEqual(test_genome.size, 100000)

    def test_predict_operons(self):
        '''
            Test the predict_operons function
        '''
        print('Testing predict_operons function')
        self.importer.inputgenomes['E_coli_BW2952']['gbk'] = '../testdata/E_coli_BW2952.100000.gbk'
        operons_data = self.importer.predict_operons()
        print(operons_data)
        self.assertTrue('E_coli_BW2952' in operons_data)
        self.assertEqual(len(operons_data['E_coli_BW2952']), 22)
