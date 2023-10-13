import os
import sys
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.models import Protein
from browser.models import Genome
from browser.models import Strain
from browser.models import Sample
from browser.pipeline.genome_import import Importer

class Command(BaseCommand):
    help = 'Deletes a genome with all genes and annotations from the database'
    
    def add_arguments(self, parser):
        parser.add_argument('-i', default='', help='File with list of genomes')

    def handle(self, *args, **options):
        while True:
            answer = input('All genomes listed in the input file '
                'will be deleted. This action cannot be reversed. Continue? (y/n)'
            )
            if answer.lower() in ["y","yes"]:
                break
            elif answer.lower() in ["n","no"]:
                print('Exiting.')
                sys.exit()
            else:
                print('Please answer yes or no.')
        genomes_file = options['i']
        if not os.path.exists(genomes_file):
            raise CommandError(genomes_file + ' not found')
        print('Deleting genomes...')
        importer = Importer()
        with open(genomes_file, 'r') as infile:
            for line in infile:
                genome_name = line.rstrip('\n\r')
                if genome_name == '':
                    continue
                print('Looking for genome', genome_name)
                genome_set = Genome.objects.filter(name=genome_name)
                if genome_set.count() == 0:
                    print('Genome ' + genome_name + ' not found')
                    continue
                elif genome_set.count() > 1:
                    print('Not unique genome name: ' + genome_name)
                    continue
                if os.path.exists(os.path.join(importer.config['cgcms.json_dir'],
                                               genome_name
                                               )):
                    shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'],
                                               genome_name
                                               ))
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
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                        os.path.basename(importer.config['cgcms.search_db_nucl'])
                        ),
                        importer.config['cgcms.search_db_nucl']
                        )
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                        os.path.basename(importer.config['cgcms.search_db_prot'])),
                        importer.config['cgcms.search_db_prot']
                        )
        importer.create_search_databases()
        os.remove(importer.config['cgcms.search_db_nucl'])
        os.remove(importer.config['cgcms.search_db_prot'])
        # importer.cleanup()
        print('Done!')
