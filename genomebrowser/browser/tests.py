import hashlib
from unittest import skip
#from contextlib import contextmanager

from django.test import TransactionTestCase
from django.test import Client
from django.utils import timezone
#from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.db.utils import IntegrityError

from browser.models import Taxon
from browser.models import Config
from browser.models import Strain
from browser.models import Sample
from browser.models import Strain_metadata
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
from browser.dataimport.importer import Importer
from browser.comparative_analysis import _get_color, make_muscle_alignment

# Create your tests here.

class BrowserTestCase(TransactionTestCase):
    '''
        Models testing
    '''
    fixtures = ['minigenomes.testdata.json']
    
    def setUp(self):
        self.client = Client()
        
    def test_config(self):
        config = Config(param='parameter name', value='parameter value')
        self.assertEqual(config.param, 'parameter name')
        config.save()
        config_saved = Config.objects.get(param='parameter name')
        self.assertEqual(config_saved.value, 'parameter value')

    def test_taxon(self):
        taxon = Taxon(taxonomy_id='10',
                      eggnog_taxid='10',
                      rank='genus',
                      parent_id='1706371',
                      name='Cellvibrio'
                      )
        self.assertEqual(taxon.name, 'Cellvibrio')
        self.assertEqual(taxon.taxonomy_id, '10')
        self.assertNotEqual(taxon.taxonomy_id, 10)
        taxon.save()
        taxon_saved = Taxon.objects.get(taxonomy_id='10')
        self.assertEqual(taxon_saved.name, 'Cellvibrio')
        self.assertEqual(taxon_saved.taxonomy_id, '10')
        
    def test_strain(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        self.assertEqual(strain.strain_id, 'FW104-10B01')
        self.assertEqual(strain.taxon.name, 'Rhodanobacter denitrificans')
        strain_saved = Strain.objects.get(strain_id='FW104-10B01')
        self.assertEqual(strain_saved.strain_id, 'FW104-10B01')
        self.assertEqual(strain_saved.taxon.name, 'Rhodanobacter denitrificans')
        
    def test_sample(self):
        sample = Sample(sample_id='FW106-02',
                        full_name='FW106 groundwater metagenome'
                        )
        sample.description = 'Metagenomic sample from groundwater of FW106 well, ' +\
        'site Y-12 West, collection date 2014-06-09, 0.2 micron filter'
        sample.save()
        self.assertEqual(sample.sample_id, 'FW106-02')
        sample_saved = Sample.objects.get(sample_id='FW106-02')
        self.assertEqual(sample_saved.sample_id, 'FW106-02')

    def test_strain_metadata(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        strain_metadata = Strain_metadata(strain=strain,
                                          source='abcdef',
                                          url='https://nih.gov/',
                                          key='ghijkl',
                                          value='mnopqr'
                                          )
        self.assertEqual(strain_metadata.source, 'abcdef')
        self.assertEqual(strain_metadata.strain.strain_id, 'FW104-10B01')
        strain_metadata.save()
        strain_metadata_saved = Strain_metadata.objects.filter(source='abcdef')
        self.assertEqual(len(list(strain_metadata_saved.all())), 1)
        self.assertEqual(list(strain_metadata_saved.all())[0].strain.strain_id,
                         'FW104-10B01'
                         )

    def test_sample_metadata(self):
        sample = Sample(sample_id='FW106-02', full_name='FW106 groundwater metagenome')
        sample.description = 'Metagenomic sample from groundwater of FW106 well, ' +\
        'site Y-12 West, collection date 2014-06-09, 0.2 micron filter'
        sample.save()
        sample_metadata = Sample_metadata(sample=sample,
                                          source='abcdef',
                                          url='https://www.lbl.gov/',
                                          key='ghijkl',
                                          value='mnopqr'
                                          )
        self.assertEqual(sample_metadata.source, 'abcdef')
        self.assertEqual(sample_metadata.sample.sample_id, 'FW106-02')
        sample_metadata.save()
        sample_metadata_saved = Sample_metadata.objects.filter(source='abcdef')
        self.assertEqual(len(list(sample_metadata_saved.all())), 1)
        self.assertEqual(list(sample_metadata_saved.all())[0].sample.sample_id,
                         'FW106-02'
                         )

    def test_genome(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        sample = Sample(sample_id='FW106-02',
                        full_name='FW106 groundwater metagenome'
                        )
        sample.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        sample=sample,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        self.assertEqual(genome.name, 'FW104-10B01')
        genome.save()
        genome_saved = Genome.objects.get(name='FW104-10B01')
        self.assertEqual(genome_saved.name, 'FW104-10B01')
        self.assertEqual(genome_saved.contigs, 0)
        self.assertEqual(genome_saved.size, 0)
        self.assertEqual(genome_saved.genes, 0)
        
    def test_contig(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome)
        genome.contigs = 1
        genome.size += contig.size
        self.assertEqual(contig.name, 'scaffold0001')
        genome.save()
        contig.save()
        contig_saved = Contig.objects.get(contig_id='scaffold1')
        self.assertEqual(contig_saved.name, 'scaffold0001')
        self.assertEqual(contig_saved.genome.name, 'FW104-10B01')
        self.assertEqual(contig_saved.genome.size, 1000)
        self.assertEqual(contig_saved.genome.contigs, 1)

    def test_operon(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome
                        )
        genome.contigs = 1
        genome.size += contig.size
        genome.save()
        contig.save()
        operon = Operon(name='operon_1',
                        start=10,
                        end=990,
                        strand=-1,
                        genome=genome,
                        contig=contig
                        )
        self.assertEqual(operon.name, 'operon_1')
        operon.save()
        operon_saved = Operon.objects.get(name='operon_1')
        self.assertEqual(operon_saved.name, 'operon_1')
        self.assertEqual(operon_saved.genome.name, 'FW104-10B01')
        self.assertEqual(operon_saved.contig.name, 'scaffold0001')

    def test_cog_class(self):
        cog_class = Cog_class(cog_id='X', description='Test description')
        self.assertEqual(cog_class.cog_id, 'X')
        cog_class.save()
        cog_class_saved = Cog_class.objects.get(cog_id='X')
        self.assertEqual(cog_class_saved.description, 'Test description')

    def test_kegg_reaction(self):
        kegg_reaction = Kegg_reaction(kegg_id='R00001', description='Test description')
        self.assertEqual(kegg_reaction.kegg_id, 'R00001')
        kegg_reaction.save()
        kegg_reaction_saved = Kegg_reaction.objects.get(kegg_id='R00001')
        self.assertEqual(kegg_reaction_saved.description, 'Test description')

    def test_kegg_pathway(self):
        kegg_pathway = Kegg_pathway(kegg_id='map00001', description='Test description')
        self.assertEqual(kegg_pathway.kegg_id, 'map00001')
        kegg_pathway.save()
        kegg_pathway_saved = Kegg_pathway.objects.get(kegg_id='map00001')
        self.assertEqual(kegg_pathway_saved.description, 'Test description')

    def test_kegg_ortholog(self):
        kegg_ortholog = Kegg_ortholog(kegg_id='KO0001', description='Test description')
        self.assertEqual(kegg_ortholog.kegg_id, 'KO0001')
        kegg_ortholog.save()
        kegg_ortholog_saved = Kegg_ortholog.objects.get(kegg_id='KO0001')
        self.assertEqual(kegg_ortholog_saved.description, 'Test description')

    def test_go_term(self):
        go_term = Go_term(go_id='GO:0000001',
                          go_namespace='Test name',
                          description='Test description'
                          )
        self.assertEqual(go_term.go_namespace, 'Test name')
        go_term.save()
        go_term_saved = Go_term.objects.get(go_id='GO:0000001')
        self.assertEqual(go_term_saved.description, 'Test description')

    def test_cazy_family(self):
        cazy_family = Cazy_family(cazy_id='GT99', description='Test description')
        self.assertEqual(cazy_family.cazy_id, 'GT99')
        cazy_family.save()
        cazy_family_saved = Cazy_family.objects.get(cazy_id='GT99')
        self.assertEqual(cazy_family_saved.description, 'Test description')

    def test_ec_number(self):
        ec_number = Ec_number(ec_number='99.1.1.1', description='Test description')
        self.assertEqual(ec_number.ec_number, '99.1.1.1')
        ec_number.save()
        ec_number_saved = Ec_number.objects.get(ec_number='99.1.1.1')
        self.assertEqual(ec_number_saved.description, 'Test description')

    def test_tc_family(self):
        tc_family = Tc_family(tc_id='GT99', description='Test description')
        self.assertEqual(tc_family.tc_id, 'GT99')
        tc_family.save()
        tc_family_saved = Tc_family.objects.get(tc_id='GT99')
        self.assertEqual(tc_family_saved.description, 'Test description')

    def test_ortholog_group(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        ortholog_group = Ortholog_group(eggnog_id='COG001', taxon=taxon)
        self.assertEqual(ortholog_group.eggnog_id, 'COG001')
        ortholog_group.save()
        ortholog_group_saved = Ortholog_group.objects.get(eggnog_id='COG001')
        self.assertEqual(ortholog_group_saved.taxon.name,
                         'Rhodanobacter denitrificans'
                         )

    def test_eggnog_description(self):
        description='Test description'
        eggnog_description = Eggnog_description(
            fingerprint=hashlib.md5(description.encode('utf-8')).hexdigest(),
            description=description
            )
        self.assertEqual(eggnog_description.description, 'Test description')
        eggnog_description.save()
        eggnog_description_saved = Eggnog_description.objects.get(
            fingerprint=hashlib.md5(description.encode('utf-8')).hexdigest()
            )
        self.assertEqual(eggnog_description_saved.description, 'Test description')

    def test_protein(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        sequence = 'MAAAAPPPTTTTT'
        protein = Protein(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(
                                                   sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=taxon
                          )
        self.assertEqual(protein.name, 'ProT')
        protein.save()
        protein_saved = Protein.objects.get(name='ProT')
        self.assertEqual(protein_saved.name, 'ProT')
        self.assertEqual(protein_saved.sequence, 'MAAAAPPPTTTTT')
        self.assertEqual(protein_saved.taxonomy_id.taxonomy_id, '666685')
        description='Test description'
        cog_class = Cog_class(cog_id='X', description=description)
        cog_class.save()
        kegg_reaction = Kegg_reaction(kegg_id='R00001', description='Test description')
        kegg_reaction.save()
        kegg_pathway = Kegg_pathway(kegg_id='map00001', description='Test description')
        kegg_pathway.save()
        kegg_ortholog = Kegg_ortholog(kegg_id='KO0001', description='Test description')
        kegg_ortholog.save()
        go_term = Go_term(go_id='GO:0000001',
                          go_namespace='Test name',
                          description='Test description'
                          )
        go_term.save()
        cazy_family = Cazy_family(cazy_id='GT99', description='Test description')
        cazy_family.save()
        ec_number = Ec_number(ec_number='99.1.1.1', description='Test description')
        ec_number.save()
        tc_family = Tc_family(tc_id='GT99', description='Test description')
        tc_family.save()
        ortholog_group = Ortholog_group(eggnog_id='COG001', taxon=taxon)
        ortholog_group.save()
        eggnog_description = Eggnog_description(
            fingerprint=hashlib.md5(description.encode('utf-8')).hexdigest(),
            description=description
            )
        eggnog_description.save()
        protein.cog_classes.add(cog_class)
        protein.kegg_reactions.add(kegg_reaction)
        protein.kegg_pathways.add(kegg_pathway)
        protein.kegg_orthologs.add(kegg_ortholog)
        protein.go_terms.add(go_term)
        protein.cazy_families.add(cazy_family)
        protein.ec_numbers.add(ec_number)
        protein.tc_families.add(tc_family)
        protein.ortholog_groups.add(ortholog_group)
        protein.eggnog_description = eggnog_description
        protein.save()
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
        # Test for error on incomplete init args
        protein = Protein(length=100,
                          taxonomy_id=taxon)
        protein.save()
        with self.assertRaises(IntegrityError):
            # no protein name provided
            protein = Protein(length=100,
                              protein_hash=hashlib.md5(
                                                       sequence.encode('utf-8')
                                                       ).hexdigest(),
                              sequence=sequence,
                              taxonomy_id=taxon
                              )
            protein.save()
        with self.assertRaises(IntegrityError):
            # no protein length provided
            protein = Protein(name='ProT',
                protein_hash=hashlib.md5(
                    'AAAAAAAAAAAAAAWWWWWWWWWWWWWWWW'.encode('utf-8')
                    ).hexdigest(),
                sequence=sequence,
                taxonomy_id=taxon
                )
            protein.save()
        with self.assertRaises(IntegrityError):
            # no protein hash provided
            protein = Protein(name='ProT',
                              length=100,
                              sequence=sequence,
                              taxonomy_id=taxon
                              )
            protein.save()
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
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome
                        )
        genome.contigs = 1
        genome.size += contig.size
        genome.save()
        contig.save()
        operon = Operon(name='operon_1',
                        start=10,
                        end=990,
                        strand=-1,
                        genome=genome,
                        contig=contig
                        )
        operon.save()
        sequence = 'MAAAAPPPTTTTT'
        protein = Protein(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(
                                                   sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=taxon
                          )
        protein.save()
        gene = Gene(name='genE',
                    locus_tag= 'Aaa_0001',
                    contig=contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=genome,
                    function='hypothetical protein',
                    protein=protein,
                    operon=operon)
        self.assertEqual(gene.name, 'genE')
        gene.save()
        gene_saved = Gene.objects.get(name='genE')
        self.assertEqual(gene_saved.locus_tag, 'Aaa_0001')
        self.assertEqual(gene_saved.genome.name, 'FW104-10B01')
        self.assertEqual(gene_saved.contig.name, 'scaffold0001')
        self.assertEqual(gene_saved.protein.name, 'ProT')
        self.assertEqual(gene_saved.operon.name, 'operon_1')
        gene = Gene(name='genE',
                    locus_tag= 'Aaa_0001',
                    #contig=contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=genome,
                    function='hypothetical protein',
                    protein=protein,
                    operon=operon)

    def test_regulon(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome
                        )
        genome.contigs = 1
        genome.size += contig.size
        genome.save()
        contig.save()
        operon = Operon(name='operon_1',
                        start=10,
                        end=990,
                        strand=-1,
                        genome=genome,
                        contig=contig
                        )
        operon.save()
        sequence = 'MAAAAPPPTTTTT'
        protein = Protein(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=taxon)
        protein.save()
        gene = Gene(name='genE',
                    locus_tag= 'Aaa_0001',
                    contig=contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=genome,
                    function='hypothetical protein',
                    protein=protein,
                    operon=operon)
        self.assertEqual(gene.name, 'genE')
        gene.save()
        regulon = Regulon(name='TesT', genome=genome, description='Test regulon')
        regulon.save()
        regulon.regulators.add(gene)
        self.assertEqual(regulon.name, 'TesT')
        regulon.save()
        regulon_saved = Regulon.objects.get(name='TesT')
        self.assertEqual(regulon_saved.description, 'Test regulon')
        self.assertEqual(regulon_saved.genome.name, 'FW104-10B01')

    def test_site(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans'
                      )
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome
                        )
        genome.contigs = 1
        genome.size += contig.size
        genome.save()
        contig.save()
        operon = Operon(name='operon_1',
                        start=10,
                        end=990,
                        strand=-1,
                        genome=genome,
                        contig=contig
                        )
        operon.save()
        sequence = 'MAAAAPPPTTTTT'
        protein = Protein(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(
                                                   sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=taxon)
        protein.save()
        gene = Gene(name='genE',
                    locus_tag= 'Aaa_0001',
                    contig=contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=genome,
                    function='hypothetical protein',
                    protein=protein,
                    operon=operon)
        self.assertEqual(gene.name, 'genE')
        gene.save()
        regulon = Regulon(name='TesT', genome=genome, description='Test regulon')
        regulon.save()
        regulon.regulators.add(gene)
        regulon.save()

        site = Site(name='Test site 1',
                    type='TFBS',
                    start=991,
                    end=999,
                    strand=-1,
                    sequence='AAAGGGGTTT',
                    regulon=regulon,
                    genome=genome,
                    contig=contig
                    )
        site.save()
        site.operons.add(operon)
        site.genes.add(gene)
        self.assertEqual(site.name, 'Test site 1')
        site.save()
        site_saved = Site.objects.get(name='Test site 1')
        self.assertEqual(site_saved.genome.name, 'FW104-10B01')
        self.assertEqual(site_saved.regulon.name, 'TesT')

    def test_annotation(self):
        taxon = Taxon(taxonomy_id='666685',
                      eggnog_taxid='666685',
                      rank='species',
                      parent_id='75309',
                      name='Rhodanobacter denitrificans')
        taxon.save()
        strain = Strain(strain_id='FW104-10B01',
                        full_name='Rhodanobacter denitrificans str. FW104-10B01',
                        order='Xanthomonadales',
                        taxon=taxon
                        )
        strain.save()
        genome = Genome(name='FW104-10B01',
                        description='Test description',
                        strain=strain,
                        json_url='https://example.com/genomes/FW104-10B01',
                        external_url="https://nih.gov",
                        external_id="123456",
                        gbk_filepath='',
                        pub_date=timezone.now(),
                        contigs=0,
                        genes=0,
                        size=0,
                        taxon=taxon)
        contig = Contig(contig_id='scaffold1',
                        name='scaffold0001',
                        size=1000,
                        genome=genome
                        )
        genome.contigs = 1
        genome.size += contig.size
        genome.save()
        contig.save()
        sequence = 'MAAAAPPPTTTTT'
        protein = Protein(name='ProT',
                          length=100,
                          protein_hash=hashlib.md5(sequence.encode('utf-8')
                                                   ).hexdigest(),
                          sequence=sequence,
                          taxonomy_id=taxon)
        protein.save()
        gene = Gene(name='genE',
                    locus_tag= 'Aaa_0001',
                    contig=contig,
                    type='CDS',
                    start=10,
                    end=310,
                    strand=-1,
                    genome=genome,
                    function='hypothetical protein',
                    protein=protein)
        gene.save()
        
        annotation = Annotation(gene_id=gene,
                                source='Personal communication',
                                url='https://nih.gov',
                                key='group',
                                value='group_name',
                                note='Test description'
                                )
        self.assertEqual(annotation.value, 'group_name')
        annotation.save()
        annotation_saved = Annotation.objects.get(value='group_name')
        self.assertEqual(annotation_saved.gene_id.locus_tag, 'Aaa_0001')
        self.assertEqual(annotation_saved.key, 'group')
        
    def test_homepage(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<a href="/genomes/">Click here for a list of')

    def test_genomes_page(self):
        response = self.client.get('/genomes/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Genomes</h2>')
        self.assertContains(response, 'E_coli_BW2952')

    def test_strains_page(self):
        response = self.client.get('/strains/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Strains</h2>')
        self.assertContains(response, 'BW2952')

    def test_samples_page(self):
        response = self.client.get('/samples/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Samples</h2>')
        self.assertContains(response, 'No samples found.')

    def test_genome_page(self):
        response = self.client.get('/genome/E_coli_BW2952')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'E_coli_BW2952')

    def test_gene_page(self):
        response = self.client.get('/gene/E_coli_BW2952/BWG_RS00020/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Genome viewer')
        self.assertContains(response, 'BWG_RS00020')

    def test_operonlist_page(self):
        response = self.client.get('/operons/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operons in')
        self.assertContains(response, 'NC_012759: 190..5020')

    def test_sitelist_page(self):
        response = self.client.get('/sites/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sites in')
        self.assertContains(response, 'NC_012759: complement(70130..70146)')

    def test_regulonlist_page(self):
        response = self.client.get('/regulons/E_coli_BW2952/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulons in')
        self.assertContains(response, 'AraC')

    def test_genelist1_page(self):
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

    def test_genelist2_page(self):
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
        # Annotation list for "Pfam" text query
        response = self.client.get('/searchannotation/',
                                   {'annotation_query':'Pfam'}
                                   )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search results for')
        self.assertContains(response,
                            'data: {\'annotation_query\': "Pfam", ' +\
                            '\'genome\': "", \'page\': "" }'
                            )

    def test_operon_page(self):
        operon_id = Operon.objects.values_list('name', flat=True)[0]
        response = self.client.get('/operon/' + operon_id + '/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Operon information')
        self.assertContains(response, operon_id)

    def test_site_page(self):
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
        genome_id = 'E_coli_BW2952'
        response = self.client.get('/regulon/' + genome_id + '/AraC/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Regulon information')
        self.assertContains(response, 'AraC')

    def test_textsearch_page(self):
        response = self.client.get('/textsearch/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Search genes by name or function')

    def test_kolist_page1(self):
        response = self.client.get('/kos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_kolist_page2(self):
        response = self.client.get('/kos/', {'query':'isopropylmalate'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'K00052')

    def test_pathwayslist_page1(self):
        response = self.client.get('/pathways/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_pathwayslist_page2(self):
        response = self.client.get('/pathways/', {'query':'Pentose'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'map00030')

    def test_reactionslist_page1(self):
        response = self.client.get('/reactions/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R00006')

    def test_reactionslist_page2(self):
        response = self.client.get('/reactions/', {'query':'tetrahydrobiopterin'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'R11765')

    def test_enzymeslist_page1(self):
        response = self.client.get('/enzymes/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.262')

    def test_enzymeslist_page2(self):
        response = self.client.get('/enzymes/', {'query':'homoserine'})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.1.1.3')

    def test_transporterslist_page1(self):
        response = self.client.get('/transporters/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '1.A.33.1')

    def test_cazylist_page(self):
        response = self.client.get('/cazy/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'CAZymes')

    def test_cogslist_page1(self):
        response = self.client.get('/cogs/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Energy production and conversion')

    def test_golist_page1(self):
        response = self.client.get('/gos/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '[GO:0000003]')

    def test_proteinsearch_page(self):
        response = self.client.get('/protsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Sequence search')

    def test_nucleotidesearch_page(self):
        response = self.client.get('/nuclsearchform/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response,
                            'Enter one nucleotide sequence in FASTA format'
                            )

    def test_help_page(self):
        response = self.client.get('/help/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'ENIGMA genome browser')

    def test_comparative_view(self):
        response = self.client.get('/comparative/',
                                   {'genome':'E_coli_BW2952',
                                    'locus_tag':'BWG_RS00020',
                                    'og':'102492',
                                    'size':'10',
                                    'lines':'50'
                                    }
                                   )
        print(response.content)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '<h2>Comparative analysis</h2>')

    def test_comparative_get_color(self):
        color = _get_color(6)
        self.assertEqual(color,
            '.setColorGradient(\'rgb(255, 182, 199)\', \'rgb(204, 102, 119)\')'
            )

    def test_make_muscle_alignment(self):
        proteins = '>1\nMKRISTTITTTTITTGNGAG\n' +\
        '>2\nMKRISTTITTTITITTGNGAG\n' +\
        '>3\nMKRISTTITTTITITTGNGAG'
        outfasta = make_muscle_alignment(proteins)
        outfasta_lines = outfasta.split('\n')
        self.assertEqual(outfasta_lines[-6], 'MKRISTTITTT-TITTGNGAG')


class ImporterTestCase(TransactionTestCase):
    fixtures = ['testdata.json']
    
    def setUp(self):
        self.importer = Importer()

    @skip("working")
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

    @skip("working")
    def test_generate_strain_data(self):
        """Read GBK file and return correct strain data"""
        p_aeruginosa = self.importer.generate_strain_data(
            '/mnt/data2/Bacteriocins/genomes/Pseudomonas_aeruginosa_PAO1.gb',
            'PAO1',
            ''
            )
        self.assertEqual(p_aeruginosa.strain_id, 'PAO1')
        self.assertEqual(p_aeruginosa.full_name, 'Pseudomonas aeruginosa PAO1')
        self.assertEqual(p_aeruginosa.taxon.taxonomy_id, '208964')
        self.assertEqual(p_aeruginosa.order, 'Pseudomonadales')
        isolate = self.importer.generate_strain_data(
            '/mnt/data2/ENIGMA/genome_files/genbank/DP16D-E2.genome.gbff.gz',
            'DP16D-E2',
            ''
            )
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


