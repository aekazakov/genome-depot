import os
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer

class Command(BaseCommand):
    help = 'Deletes a genome with all genes and annotations from the database'
    
    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        genome_name = options['g']
        if genome_name == '':
            raise CommandError('Genome name required')
        print('Looking for genome', genome_name)
        genome_set = Genome.objects.filter(name=genome_name)
        if genome_set.count() == 0:
            print('Genome ' + genome_name + ' not found')
            raise CommandError()
        elif genome_set.count() > 1:
            print('Not unique genome name: ' + genome_name)
            raise CommandError()
        importer = Importer()
        print('Note 1: Genome deletion removes contigs, genes and gene annotations for this genome. It also removes proteins not linked to any gene. But it does not remove a strain associated with this genome if the strain has other genomes.')
        print('Deleting genome...')
        genome_set.delete()
        print('Deleting proteins not linked to genes...')
        Protein.objects.filter(gene=None).delete()
        print('Deleting strains not linked to genomes...')
        Strain.objects.filter(genome=None).delete()
        print('Deleting samples not linked to genomes...')
        Sample.objects.filter(genome=None).delete()
        importer.export_proteins()
        importer.export_contigs()
        importer.delete_search_databases()
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_nucl'])), importer.config['cgcms.search_db_nucl'])
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_prot'])), importer.config['cgcms.search_db_prot'])
        importer.create_search_databases()
        shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'], genome_name))
        os.remove(importer.config['cgcms.search_db_nucl'])
        os.remove(importer.config['cgcms.search_db_prot'])
        # importer.cleanup()
        print('Done!')
        
