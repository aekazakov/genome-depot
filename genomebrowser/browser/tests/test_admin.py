import os
from unittest import skip

from django.test import TestCase
from django.test import TransactionTestCase
from django.test import Client
from django.test import RequestFactory

from django.contrib.auth.models import User
from django.contrib.admin.sites import AdminSite
from django.shortcuts import reverse

from genomebrowser.settings import BASE_DIR
from browser.models import Config
from browser.models import Genome
from browser.models import Annotation
from browser.models import Strain
from browser.models import Sample
from browser.models import Tag
from browser.models import Sample_metadata
from browser.models import Strain_metadata
from browser.models import Regulon

from browser.admin import GenomeAdmin
from browser.tasks import test_task_impl
from browser.tasks import import_genomes_impl
from browser.tasks import delete_genomes_impl
from browser.tasks import update_static_files_impl
from browser.tasks import import_sample_metadata_impl
from browser.tasks import import_sample_descriptions_impl
from browser.tasks import update_strain_metadata_impl
from browser.tasks import import_annotations_impl
from browser.tasks import import_regulon_impl
from browser.tasks import run_annotation_pipeline_impl


class AdminTestCase(TestCase):

    fixtures = ['minigenomes.testdata.json']

    def setUp(self) -> None:
        self.site = AdminSite()
        self.username = 'new_user'
        self.password = 'password'
        email = 'aekazakov@lbl.gov'  #'test@example.com'
        self.superuser = User.objects.create_superuser(
            username=self.username, email=email, password=self.password,
        )

    def test_delete_genome_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        data = {'action': 'delete_genomes',
                '_selected_action': Genome.objects.filter(name='E_coli_CFT073').values_list('pk', flat=True),
                }
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_test_task_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        data = {'action': 'test_task',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                }
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_run_annotation_tools_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        change_url = reverse("admin:browser_genome_changelist")
        # send GET request
        data = {'action': 'run_annotation_tools',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
               }
        response = self.client.get(change_url, data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        data = {'action': 'run_annotation_tools',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                'do_action': 'yes',
                'tools': ['gapmind', ], }
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    #@skip("skip for now")
    def test_import_genomes_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'import_genomes',
            '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
            }
        response = self.client.get('/admin/browser/genome/add/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'import_genomes',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                'do_action': 'yes',
                'download_email': '',
                'tsv_file': fp, }
            response = self.client.post('/admin/browser/genome/add/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_update_static_files_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        data = {'action': 'update_static_files',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),}
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_add_genome_tag_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'add_genome_tag',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                }
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.get(change_url, data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        data = {'action': 'add_genome_tag',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                'do_action': 'yes',
                'tag': Tag.objects.filter(name='test').values_list('pk', flat=True),}
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_remove_genome_tag_action(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'remove_genome_tag',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                }
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.get(change_url, data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        data = {'action': 'remove_genome_tag',
                '_selected_action': Genome.objects.filter(name='E_coli_BW2952').values_list('pk', flat=True),
                'do_action': 'yes',
                'tag': Tag.objects.filter(name='test').values_list('pk', flat=True),}
        change_url = reverse("admin:browser_genome_changelist")
        response = self.client.post(change_url, data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_import_sample_descriptions(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'import_sample_descriptions',
                '_selected_action': Sample.objects.all().values_list('pk', flat=True),
                }
        response = self.client.get('/admin/browser/sample/import-descriptions/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'import_sample_descriptions',
                    '_selected_action': Sample.objects.all().values_list('pk', flat=True),
                    'tsv_file': fp, }
            response = self.client.post('/admin/browser/sample/import-descriptions/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_update_strain_metadata(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'update_strain_metadata',
                '_selected_action': Strain.objects.all().values_list('pk', flat=True),
                }
        response = self.client.get('/admin/browser/strain_metadata/import/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'update_strain_metadata',
                    '_selected_action': Strain.objects.all().values_list('pk', flat=True),
                    'xlsx_file': fp, }
            response = self.client.post('/admin/browser/strain_metadata/import/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_import_sample_metadata(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'import_sample_metadata',
                '_selected_action': Sample.objects.all().values_list('pk', flat=True),
                }
        response = self.client.get('/admin/browser/sample_metadata/import/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'import_sample_metadata',
                    '_selected_action': Sample.objects.all().values_list('pk', flat=True),
                    'tsv_file': fp, }
            #change_url = reverse("admin:browser_sample_import")
            response = self.client.post('/admin/browser/sample_metadata/import/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_import_annotations(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'import_annotations',
                '_selected_action': [],
                }
        response = self.client.get('/admin/browser/annotation/import-annotations/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'import_annotations',
                    '_selected_action': [],
                    'tsv_file': fp, }
            #change_url = reverse("admin:browser_annotation_import-annotations")
            response = self.client.post('/admin/browser/annotation/import-annotations/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    def test_import_regulons(self):
        # This test only submits data to django_Q, not running actual backend
        self.client.login(username=self.username, password=self.password)
        # send GET request
        data = {'action': 'import_regulons',
                '_selected_action': [],
                }
        response = self.client.get('/admin/browser/regulon/import-regulons/', data, follow=True)
        self.assertEqual(response.status_code, 200)
        # send POST request
        with open(os.path.join(BASE_DIR, '../testdata/import_minigenomes.txt'), "rb") as fp:
            data = {'action': 'import_regulons',
                    '_selected_action': [],
                    'tsv_file': fp, }
            #change_url = reverse("admin:browser_regulon_import-regulons")
            response = self.client.post('/admin/browser/regulon/import-regulons/', data, follow=True)
        self.client.logout()
        self.assertEqual(response.status_code, 200)

    # Test task implementations
    
    def test_test_task_impl(self):
        # Tests implementation of test_task_impl task
        request = RequestFactory().get('/admin/')
        genomes = 'genome1, genome2'
        out = test_task_impl(request, genomes)
        self.assertTrue(out.startswith('Genomes:'))
        
    def test_run_annotation_pipeline_impl(self):
        # Tests implementation of run_annotation_pipeline task
        genomes = ['E_coli_BW2952',]
        plugins = ['gapmind',]
        args = (genomes, plugins)
        out = run_annotation_pipeline_impl(args)
        self.assertEqual(out, 'Annotation pipeline finished.')
        
    def test_update_static_files_impl(self):
        # Tests implementation of update_static_files task
        genomes = ['E_coli_BW2952',]
        update_static_files_impl(genomes)
        json_dir = os.path.join(Config.objects.get(param='cgcms.json_dir').value, genomes[0])
        self.assertTrue(os.path.exists(json_dir))

    def test_delete_genomes_impl(self):
        # Tests implementation of delete_genomes task
        genome = 'E_coli_BW2952'
        genomes = [genome,]
        delete_genomes_impl(genomes)
        self.assertEqual(Genome.objects.filter(name=genome).count(), 0)

    def test_import_sample_metadata_impl(self):
        # Tests implementation of import_sample_metadata task
        lines = ['test_sample\ttest\ttest\ttest\ttest',]
        import_sample_metadata_impl(lines)
        self.assertEqual(Sample_metadata.objects.get(source='test').value, 'test')

    def test_import_sample_descriptions_impl(self):
        # Tests implementation of import_sample_metadata task
        lines = ['test_sample\ttest\ttest description',]
        import_sample_descriptions_impl(lines)
        self.assertEqual(Sample.objects.get(sample_id='test_sample').description, 'test description')

    def test_import_annotations_impl(self):
        # Tests implementation of import_annotations task
        lines = ['BWG_RS00070\tE_coli_BW2952\ttest\ttest\ttest\ttest\ttest',]
        import_annotations_impl(lines)
        self.assertEqual(Annotation.objects.get(gene_id__locus_tag='BWG_RS00070', source='test').note, 'test')

    def test_import_regulon_impl(self):
        # Tests implementation of import_regulon task
        lines = ['TEST\tE_coli_BW2952\tBWG_RS00070\tBWG_RS00070\tNC_012759\t1\t3\t1\tATT',]
        import_regulon_impl(lines)
        self.assertEqual(Regulon.objects.get(name='TEST').genome.name, 'E_coli_BW2952')

    def test_update_strain_metadata_impl(self):
        # Tests implementation of update_strain_metadata task
        strain_id = 'BW2952'
        test_file = 'test_strain_metadata.xlsx'
        with open(test_file, "rb") as fp:
            update_strain_metadata_impl(fp)
        saved_metadata = Strain_metadata.objects.get(key='Phylogenetic Order')
        self.assertEqual(saved_metadata.strain.strain_id, strain_id)
        self.assertEqual(saved_metadata.value, 'Enterobacterales')


class AdminTransactionTestCase(TransactionTestCase):
    '''
        Testing admin functions that require TransactionTestCase
    '''

    fixtures = ['minigenomes.testdata.json']
    
    def setUp(self):
        self.site = AdminSite()

    def test_import_genomes_impl(self):
        # Tests implementation of import_genomes task
        lines = ['../testdata/E_coli_BW2952.100000.gbk\tE_coli_BW2952.test\tBW2952\t\thttps://www.ncbi.nlm.nih.gov/assembly/GCF_000022345.1\tNCBI:GCF_000022345.1',]
        email = 'test@example.com'
        args = (lines, email)
        result = import_genomes_impl(args)
        self.assertEqual(result, 'Done!')
        
