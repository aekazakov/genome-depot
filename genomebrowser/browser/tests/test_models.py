import os
import hashlib
from unittest import skip
#from contextlib import contextmanager

from django.test import TestCase
from django.utils import timezone
#from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.db import transaction
from django.db.utils import IntegrityError

from browser.models import Taxon
from browser.models import Config
from browser.models import Strain
from browser.models import Sample
from browser.models import Strain_metadata
from browser.models import Tag
from browser.models import Genome
from browser.models import Annotation
from browser.models import Regulon
from browser.models import Contig
from browser.models import Sample_metadata
from browser.models import Kegg_ortholog
from browser.models import Kegg_pathway
from browser.models import Kegg_reaction
from browser.models import Go_term
from browser.models import Cog_class
from browser.models import Ec_number
from browser.models import Cazy_family
from browser.models import Operon
from browser.models import Ortholog_group
from browser.models import Eggnog_description
from browser.models import Tc_family
from browser.models import Protein
from browser.models import Gene
from browser.models import Site

# Create your tests here.

class BrowserModelsTest(TestCase):
    '''
        Models testing
    '''
    @classmethod
    def setUpTestData(cls):
        cls.config = Config.objects.create(param='parameter name', value='parameter value')
        cls.taxon = Taxon.objects.create(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        cls.strain = Strain.objects.create(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=cls.taxon
                        )
        cls.sample = Sample.objects.create(sample_id='FW106-02',
                        full_name='FW106 groundwater metagenome'
                        )
        cls.strain_metadata = Strain_metadata.objects.create(strain=cls.strain,
                                          source='abcdef',
                                          url='https://nih.gov/',
                                          key='ghijkl',
                                          value='mnopqr'
                                          )
        cls.sample_metadata = Sample_metadata.objects.create(sample=cls.sample,
                                          source='abcdef',
                                          url='https://www.lbl.gov/',
                                          key='ghijkl',
                                          value='mnopqr'
                                          )
        cls.tag = Tag.objects.create(name='Test tag',
            description = 'Test tag description'
        )

        cls.genome = Genome.objects.create(name='FW104-10B01',
                        description='Test description',
                        strain=cls.strain,
                        sample=cls.sample,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=cls.taxon)
        cls.contig = Contig.objects.create(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=cls.genome)
        cls.genome.contigs = 1
        cls.genome.genes = 1
        cls.genome.size += cls.contig.size
        cls.genome.save()
        cls.operon = Operon.objects.create(name='operon_1',
                        start=10,
                        end=990,
                        strand=-1,
                        genome=cls.genome,
                        contig=cls.contig
                        )
        cls.cog_class = Cog_class.objects.create(cog_id='X', description='Test description')
        cls.kegg_reaction = Kegg_reaction.objects.create(kegg_id='R00001', description='Test description')        
        cls.kegg_pathway = Kegg_pathway.objects.create(kegg_id='map00001', description='Test description')
        cls.kegg_ortholog = Kegg_ortholog.objects.create(kegg_id='KO0001', description='Test description')
        cls.go_term = Go_term.objects.create(go_id='GO:0000001',
                          go_namespace='Test name',
                          description='Test description'
                          )
        cls.cazy_family = Cazy_family.objects.create(cazy_id='GT99', description='Test description')
        cls.ec_number = Ec_number.objects.create(ec_number='99.1.1.1', description='Test description')
        cls.tc_family = Tc_family.objects.create(tc_id='GT99', description='Test description')
        cls.ortholog_group = Ortholog_group.objects.create(eggnog_id='COG001', taxon=cls.taxon)
        description='Test description'
        cls.eggnog_description = Eggnog_description.objects.create(
            fingerprint=hashlib.md5(description.encode('utf-8')).hexdigest(),
            description=description
            )
        sequence = 'MAAAAPPPTTTTT'
        cls.protein = Protein.objects.create(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(
                                                   sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=cls.taxon
                          )
        cls.protein.cog_classes.add(cls.cog_class)
        cls.protein.kegg_reactions.add(cls.kegg_reaction)
        cls.protein.kegg_pathways.add(cls.kegg_pathway)
        cls.protein.kegg_orthologs.add(cls.kegg_ortholog)
        cls.protein.go_terms.add(cls.go_term)
        cls.protein.cazy_families.add(cls.cazy_family)
        cls.protein.ec_numbers.add(cls.ec_number)
        cls.protein.tc_families.add(cls.tc_family)
        cls.protein.ortholog_groups.add(cls.ortholog_group)
        cls.protein.eggnog_description = cls.eggnog_description
        cls.protein.save()
        cls.gene = Gene.objects.create(name='genE',
                    locus_tag= 'Aaa_0001',
                    contig=cls.contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=cls.genome,
                    function='hypothetical protein',
                    protein=cls.protein,
                    operon=cls.operon)
        cls.regulon = Regulon.objects.create(name='TesT',
            genome=cls.genome,
            description='Test regulon')
        cls.regulon.regulators.add(cls.gene)
        cls.regulon.save()
        cls.site = Site.objects.create(name='Test site 1',
                    type='TFBS',
                    start=991,
                    end=999,
                    strand=-1,
                    sequence='AAAGGGGTTT',
                    regulon=cls.regulon,
                    genome=cls.genome,
                    contig=cls.contig
                    )
        cls.site.operons.add(cls.operon)
        cls.site.genes.add(cls.gene)
        cls.site.save()
        cls.annotation = Annotation.objects.create(gene_id=cls.gene,
                                source='Personal communication',
                                url='https://nih.gov',
                                key='group',
                                value='group_name',
                                note='Test description'
                                )
        
    def test_config(self):
        print('Testing Config object creation')
        self.assertEqual(self.config.param, 'parameter name')
        self.assertEqual(self.config.value, 'parameter value')
        self.assertEqual(str(self.config), 'parameter name:parameter value')

    def test_taxon(self):
        print('Testing Taxon object creation')
        self.assertEqual(self.taxon.name, 'Rhodanobacter denitrificans')
        self.assertEqual(self.taxon.taxonomy_id, '666685')
        self.assertEqual(self.taxon.parent_id, '75309')
        self.assertEqual(self.taxon.rank, 'species')
        self.assertEqual(self.taxon.eggnog_taxid, '666685')
        self.assertEqual(str(self.taxon), 'Rhodanobacter denitrificans (666685)')

    def test_strain(self):
        print('Testing Strain object creation')
        self.assertEqual(self.strain.strain_id, 'FW104-10B01')
        self.assertEqual(self.strain.full_name, 'Rhodanobacter denitrificans str. FW104-10B01')
        self.assertEqual(self.strain.taxon.name, 'Rhodanobacter denitrificans')
        self.assertEqual(str(self.strain), 'FW104-10B01 (Xanthomonadales)')

    def test_sample(self):
        print('Testing Sample object creation')
        self.sample.description = 'Metagenomic sample from groundwater of FW106 well'
        self.sample.save()
        self.assertEqual(self.sample.sample_id, 'FW106-02')
        self.assertEqual(str(self.sample), 'FW106-02 (FW106 groundwater metagenome)')

    def test_strain_metadata(self):
        print('Testing Strain object creation')
        self.assertEqual(self.strain_metadata.source, 'abcdef')
        self.assertEqual(self.strain_metadata.strain.strain_id, 'FW104-10B01')
        self.assertEqual(str(self.strain_metadata), 'abcdef: mnopqr')

    def test_sample_metadata(self):
        print('Testing Sample_metadata object creation')
        self.assertEqual(self.sample_metadata.source, 'abcdef')
        self.assertEqual(self.sample_metadata.sample.sample_id, 'FW106-02')
        sample_metadata_saved = Sample_metadata.objects.get(source='abcdef')
        self.assertEqual(sample_metadata_saved.sample.sample_id, 'FW106-02')
        self.assertEqual(str(sample_metadata_saved), 'abcdef: mnopqr')

    def test_tag(self):
        print('Testing Tag model')
        self.assertEqual(self.tag.name, 'Test tag')
        tag_saved = Tag.objects.get(name='Test tag')
        self.assertEqual(tag_saved.description, 'Test tag description')
        self.assertEqual(str(tag_saved), 'Test tag')

    def test_genome(self):
        print('Testing Genome object creation')
        self.assertEqual(self.genome.name, 'FW104-10B01')
        genome_saved = Genome.objects.get(name='FW104-10B01')
        self.assertEqual(genome_saved.name, 'FW104-10B01')
        self.assertEqual(genome_saved.contigs, 1)
        self.assertEqual(genome_saved.size, 1000)
        self.assertEqual(genome_saved.genes, 1)
        self.assertEqual(str(genome_saved), 'FW104-10B01')
        tag_saved = Tag.objects.get(name='Test tag')
        genome_saved.tags.add(tag_saved)
        self.assertEqual(genome_saved.get_tags(), 'Test tag')
        
    def test_contig(self):
        print('Testing Contig object creation')
        self.assertEqual(self.contig.name, 'scaffold0001')
        contig_saved = Contig.objects.get(contig_id='scaffold1')
        self.assertEqual(contig_saved.name, 'scaffold0001')
        self.assertEqual(contig_saved.genome.name, 'FW104-10B01')
        self.assertEqual(contig_saved.genome.size, 1000)
        self.assertEqual(contig_saved.genome.contigs, 1)
        self.assertEqual(str(contig_saved), 'scaffold1')

    def test_operon(self):
        print('Testing Operon object creation')
        self.assertEqual(self.operon.name, 'operon_1')
        operon_saved = Operon.objects.get(name='operon_1')
        self.assertEqual(operon_saved.name, 'operon_1')
        self.assertEqual(operon_saved.genome.name, 'FW104-10B01')
        self.assertEqual(operon_saved.contig.name, 'scaffold0001')
        self.assertEqual(str(operon_saved), 'operon_1')

    def test_cog_class(self):
        print('Testing Cog_class object creation')
        self.assertEqual(self.cog_class.cog_id, 'X')
        cog_class_saved = Cog_class.objects.get(cog_id='X')
        self.assertEqual(cog_class_saved.description, 'Test description')
        self.assertEqual(str(cog_class_saved), 'X')

    def test_kegg_reaction(self):
        print('Testing Kegg_reaction object creation')
        self.assertEqual(self.kegg_reaction.kegg_id, 'R00001')
        kegg_reaction_saved = Kegg_reaction.objects.get(kegg_id='R00001')
        self.assertEqual(kegg_reaction_saved.description, 'Test description')
        self.assertEqual(str(kegg_reaction_saved), 'R00001')

    def test_kegg_pathway(self):
        print('Testing Kegg_pathway object creation')
        self.assertEqual(self.kegg_pathway.kegg_id, 'map00001')
        kegg_pathway_saved = Kegg_pathway.objects.get(kegg_id='map00001')
        self.assertEqual(kegg_pathway_saved.description, 'Test description')
        self.assertEqual(str(kegg_pathway_saved), 'map00001')

    def test_kegg_ortholog(self):
        print('Testing Kegg_ortholog object creation')
        self.assertEqual(self.kegg_ortholog.kegg_id, 'KO0001')
        kegg_ortholog_saved = Kegg_ortholog.objects.get(kegg_id='KO0001')
        self.assertEqual(kegg_ortholog_saved.description, 'Test description')
        self.assertEqual(str(kegg_ortholog_saved), 'KO0001')

    def test_go_term(self):
        print('Testing Go_term object creation')
        self.assertEqual(self.go_term.go_namespace, 'Test name')
        go_term_saved = Go_term.objects.get(go_id='GO:0000001')
        self.assertEqual(go_term_saved.description, 'Test description')
        self.assertEqual(str(go_term_saved), 'GO:0000001')

    def test_cazy_family(self):
        print('Testing Cazy_family object creation')
        self.assertEqual(self.cazy_family.cazy_id, 'GT99')
        cazy_family_saved = Cazy_family.objects.get(cazy_id='GT99')
        self.assertEqual(cazy_family_saved.description, 'Test description')
        self.assertEqual(str(cazy_family_saved), 'GT99')

    def test_ec_number(self):
        print('Testing Ec_number object creation')
        ec_number = Ec_number(ec_number='99.1.1.1', description='Test description')
        self.assertEqual(self.ec_number.ec_number, '99.1.1.1')
        ec_number_saved = Ec_number.objects.get(ec_number='99.1.1.1')
        self.assertEqual(ec_number_saved.description, 'Test description')
        self.assertEqual(str(ec_number_saved), '99.1.1.1')

    def test_tc_family(self):
        print('Testing Tc_family object creation')
        self.assertEqual(self.tc_family.tc_id, 'GT99')
        tc_family_saved = Tc_family.objects.get(tc_id='GT99')
        self.assertEqual(tc_family_saved.description, 'Test description')
        self.assertEqual(str(tc_family_saved), 'GT99')

    def test_ortholog_group(self):
        print('Testing Ortholog_group object creation')
        self.assertEqual(self.ortholog_group.eggnog_id, 'COG001')
        ortholog_group_saved = Ortholog_group.objects.get(eggnog_id='COG001')
        self.assertEqual(ortholog_group_saved.taxon.name,
                         'Rhodanobacter denitrificans'
                         )
        self.assertEqual(str(ortholog_group_saved), 'COG001')

    def test_eggnog_description(self):
        print('Testing Eggnog_description object creation')
        description='Test description'
        self.assertEqual(self.eggnog_description.description, 'Test description')
        eggnog_description_saved = Eggnog_description.objects.get(
            fingerprint=hashlib.md5(description.encode('utf-8')).hexdigest()
            )
        self.assertEqual(eggnog_description_saved.description, 'Test description')
        self.assertEqual(str(eggnog_description_saved), 'Test description')

    def test_protein(self):
        print('Testing Protein object creation')
        sequence = 'MAAAAPPPTTTTT'
        self.assertEqual(self.protein.name, 'ProT')
        protein_saved = Protein.objects.get(name='ProT')
        self.assertEqual(protein_saved.name, 'ProT')
        self.assertEqual(protein_saved.sequence, sequence)
        self.assertEqual(protein_saved.taxonomy_id.taxonomy_id, '666685')
        description='Test description'
        protein_saved = Protein.objects.get(name='ProT')
        self.assertEqual(protein_saved.cog_classes.all()[0].cog_id, 'X')
        self.assertEqual(protein_saved.kegg_reactions.all()[0].kegg_id, 'R00001')
        self.assertEqual(protein_saved.kegg_pathways.all()[0].kegg_id, 'map00001')
        self.assertEqual(protein_saved.kegg_orthologs.all()[0].kegg_id, 'KO0001')
        self.assertEqual(protein_saved.go_terms.all()[0].go_id, 'GO:0000001')
        self.assertEqual(protein_saved.cazy_families.all()[0].cazy_id, 'GT99')
        self.assertEqual(protein_saved.ec_numbers.all()[0].ec_number, '99.1.1.1')
        self.assertEqual(protein_saved.tc_families.all()[0].tc_id, 'GT99')
        self.assertEqual(protein_saved.ortholog_groups.all()[0].eggnog_id, 'COG001')
        self.assertEqual(protein_saved.eggnog_description.description,
             'Test description'
        )
        self.assertEqual(str(protein_saved), hashlib.md5(
            sequence.encode('utf-8')).hexdigest()
        )
        # Tests for error on incomplete init args
        with transaction.atomic():
            with self.assertRaises(IntegrityError):
                # no protein name provided
                protein = Protein(length=100,
                                  protein_hash=hashlib.md5(
                                                           sequence.encode('utf-8')
                                                           ).hexdigest(),
                                  sequence=sequence,
                                  taxonomy_id=self.taxon
                                  )
                protein.save()
        with transaction.atomic():
            with self.assertRaises(IntegrityError):
                # no protein length provided
                protein = Protein(name='ProT',
                    protein_hash=hashlib.md5(
                        'AAAAAAAAAAAAAAWWWWWWWWWWWWWWWW'.encode('utf-8')
                        ).hexdigest(),
                    sequence=sequence,
                    taxonomy_id=self.taxon
                    )
                protein.save()
        with transaction.atomic():
            with self.assertRaises(ValueError):
                # string taxonomy ID instead of Taxon
                protein = Protein(name='ProT',
                    length=100,
                    protein_hash=hashlib.md5(
                        'MMMMMMMMMMMMMGGGGGGGGGGGGGGGG'.encode('utf-8')
                        ).hexdigest(),
                    sequence=sequence,
                    taxonomy_id='0123'
                    )
                protein.save()

    def test_gene(self):
        print('Testing Gene object creation')
        self.assertEqual(self.gene.name, 'genE')
        gene_saved = Gene.objects.get(name='genE')
        self.assertEqual(gene_saved.locus_tag, 'Aaa_0001')
        self.assertEqual(gene_saved.genome.name, 'FW104-10B01')
        self.assertEqual(gene_saved.contig.name, 'scaffold0001')
        self.assertEqual(gene_saved.protein.name, 'ProT')
        self.assertEqual(gene_saved.operon.name, 'operon_1')
        self.assertEqual(str(gene_saved), 'Aaa_0001')

    def test_regulon(self):
        print('Testing Regulon object creation')
        self.assertEqual(self.regulon.name, 'TesT')
        regulon_saved = Regulon.objects.get(name='TesT')
        self.assertEqual(regulon_saved.description, 'Test regulon')
        self.assertEqual(regulon_saved.genome.name, 'FW104-10B01')
        self.assertEqual(str(regulon_saved), 'TesT(FW104-10B01)')

    def test_site(self):
        print('Testing Site object creation')
        self.assertEqual(self.site.name, 'Test site 1')
        site_saved = Site.objects.get(name='Test site 1')
        self.assertEqual(site_saved.genome.name, 'FW104-10B01')
        self.assertEqual(site_saved.regulon.name, 'TesT')
        self.assertEqual(str(site_saved), 'Test site 1')

    def test_annotation(self):
        print('Testing Annotation object creation')
        self.assertEqual(self.annotation.value, 'group_name')
        annotation_saved = Annotation.objects.get(value='group_name')
        self.assertEqual(annotation_saved.gene_id.locus_tag, 'Aaa_0001')
        self.assertEqual(annotation_saved.key, 'group')
        self.assertEqual(str(annotation_saved), 'Personal communication: group_name')
