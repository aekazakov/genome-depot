import os
from unittest import skip

from django.test import TransactionTestCase

from browser.pipeline.annotate import Annotator
from browser.models import Genome
from browser.pipeline.plugins.genomedepot_amrfinder import application as amrfinder_plugin
from browser.pipeline.plugins.genomedepot_antismash import application as antismash_plugin
from browser.pipeline.plugins.genomedepot_defensefinder import application as defensefinder_plugin
from browser.pipeline.plugins.genomedepot_ecis_screen import application as ecis_screen_plugin
from browser.pipeline.plugins.genomedepot_fama import application as fama_plugin
from browser.pipeline.plugins.genomedepot_gapmind import application as gapmind_plugin
from browser.pipeline.plugins.genomedepot_hmmsearch_pfam import application as hmmsearch_pfam_plugin
from browser.pipeline.plugins.genomedepot_hmmsearch_tigrfam import application as hmmsearch_tigrfam_plugin
from browser.pipeline.plugins.genomedepot_macsyfinder import application as macsyfinder_plugin
from browser.pipeline.plugins.genomedepot_phispy import application as phispy_plugin
from browser.pipeline.plugins.genomedepot_genomad import application as genomad_plugin


class PluginsTestCase(TransactionTestCase):

    fixtures = ['minigenomes.testdata.json']

    @classmethod
    def setUp(self) -> None:
        self.genome_id = 'E_coli_BW2952'
        genome = Genome.objects.get(name=self.genome_id)
        self.test_genomes = {self.genome_id:genome.gbk_filepath}
        self.annotator= Annotator()
        self.annotator.config['plugins.macsyfinder.conda_env'] = 'genomedepot-macsyfinder'
        self.annotator.config['plugins.macsyfinder.model'] = 'TXSScan'
        self.annotator.config['plugins.macsyfinder.models_dir'] = '/mnt/data/work/CGCMS/external_refdata/macsyfinder/data'
        self.annotator.config['plugins.defensefinder.conda_env'] = 'genomedepot-defensefinder'
        self.annotator.config['plugins.defensefinder.defensefinder_models_dir'] = '/mnt/data/work/CGCMS/external_refdata/defensefinder/data'
        self.annotator.config['plugins.hmmsearch_tigrfam.conda_env'] = 'genomedepot-hmmsearch'
        self.annotator.config['plugins.hmmsearch_tigrfam.hmmsearch_command'] = 'hmmsearch'
        self.annotator.config['plugins.hmmsearch_tigrfam.display_name'] = 'TIGRFAM database'
        self.annotator.config['plugins.hmmsearch_tigrfam.hmm_lib'] = '/mnt/data/work/CGCMS/external_refdata/tigrfam/TIGRFAM.HMM'
        self.annotator.config['plugins.hmmsearch_tigrfam.ref_data'] = '/mnt/data/work/CGCMS/external_refdata/tigrfam/ref_tigrfam.txt'
        self.annotator.config['plugins.hmmsearch_pfam.conda_env'] = 'genomedepot-hmmsearch'
        self.annotator.config['plugins.hmmsearch_pfam.hmmsearch_command'] = 'hmmsearch'
        self.annotator.config['plugins.hmmsearch_pfam.display_name'] = 'Pfam database'
        self.annotator.config['plugins.hmmsearch_pfam.hmm_lib'] = '/mnt/data/work/CGCMS/external_refdata/pfam/Pfam-A.hmm'
        self.annotator.config['plugins.hmmsearch_pfam.ref_data'] = '/mnt/data/work/CGCMS/external_refdata/pfam/ref_pfam.txt'
        self.annotator.config['plugins.gapmind.conda_env'] = 'genomedepot-gapmind'
        self.annotator.config['plugins.gapmind.gapmind_dir'] = '/mnt/data/work/CGCMS/external_tools/PaperBLAST'
        self.annotator.config['plugins.gapmind.threads'] = '8'
        self.annotator.config['plugins.phispy.conda_env'] = 'genomedepot-phispy'
        self.annotator.config['plugins.fama.fama_dir'] = '/mnt/data/work/CGCMS/external_tools/fama/py'
        self.annotator.config['plugins.fama.fama_config'] = '/mnt/data/work/CGCMS/external_tools/fama/config.ini'
        self.annotator.config['plugins.antismash.conda_env'] = 'genomedepot-antismash'
        self.annotator.config['plugins.antismash.antismash_cmd'] = 'antismash'
        self.annotator.config['plugins.amrfinder.enabled'] = '1'
        self.annotator.config['plugins.amrfinder.threads'] = '8'
        self.annotator.config['plugins.amrfinder.conda_env'] = 'genomedepot-amrfinder'
        self.annotator.config['plugins.amrfinder.display_name'] = 'AMRFinderPlus'
        self.annotator.config['plugins.ecis_screen.ecis-screen_cmd'] = '/mnt/data/work/CGCMS/external_tools/eCIS-screen/HMMsearch_genomesII_fast.pl'
        self.annotator.config['plugins.ecis_screen.ecis_hmm'] = '/mnt/data/work/CGCMS/external_tools/eCIS-screen/eCIS.hmm'
        self.annotator.config['plugins.ecis_screen.enabled'] = '1'
        self.annotator.config['plugins.genomad.enabled'] = '1'
        self.annotator.config['plugins.genomad.display_name'] = 'geNomad'
        self.annotator.config['plugins.genomad.conda_env'] = 'genomedepot-genomad'
        self.annotator.config['plugins.genomad.ref_db'] = '/mnt/data/work/CGCMS/external_refdata/geNomad/genomad_db'

        #@skip("skip test")
    def test_amrfinder_plugin(self):
        print('test_amrfinder_plugin')
        plugin_outfile = amrfinder_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'amrfinder-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'amrfinder-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_antismash_plugin(self):
        print('test_antismash_plugin')
        plugin_outfile = antismash_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'antismash-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'antismash-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_defensefinder_plugin(self):
        print('test_defensefinder_plugin')
        plugin_outfile = defensefinder_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'defensefinder-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'defensefinder-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_ecis_screen_plugin(self):
        print('test_ecis_screen_plugin')
        plugin_outfile = ecis_screen_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'ecis-screen-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'ecis-screen-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_fama_plugin(self):
        print('test_fama_plugin')
        plugin_outfile = fama_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'fama-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'fama-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_gapmind_plugin(self):
        print('test_gapmind_plugin')
        plugin_outfile = gapmind_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'gapmind-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'gapmind-plugin-output.txt'
                               )))

    #@skip("skip test")
    def test_genomad_plugin(self):
        print('test_genomad_plugin')
        plugin_outfile = genomad_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'genomad-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'genomad-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_hmmsearch_pfam_plugin(self):
        print('test_hmmsearch_pfam_plugin')
        plugin_outfile = hmmsearch_pfam_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'hmmsearch-pfam-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'hmmsearch-pfam-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_hmmsearch_tigrfam_plugin(self):
        print('test_hmmsearch_tigrfam_plugin')
        plugin_outfile = hmmsearch_tigrfam_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'hmmsearch-tigrfam-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'hmmsearch-tigrfam-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_macsyfinder_plugin(self):
        print('test_macsyfinder_plugin')
        plugin_outfile = macsyfinder_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'macsyfinder-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'macsyfinder-plugin-output.txt'
                               )))
        
    #@skip("skip test")
    def test_phispy_plugin(self):
        print('test_phispy_plugin')
        plugin_outfile = phispy_plugin(self.annotator, self.test_genomes)
        self.assertEqual(plugin_outfile, os.path.join(self.annotator.config['core.temp_dir'],
                               'phispy-plugin-output.txt'
                               ))
        self.assertTrue(os.path.exists(os.path.join(self.annotator.config['core.temp_dir'],
                               'phispy-plugin-output.txt'
                               )))
