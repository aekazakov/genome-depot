from django.core.management.base import BaseCommand, CommandError
from browser.models import Genome
from browser.pipeline.genome_import import Importer

class Command(BaseCommand):
    help = '''Deletes and re-creates static files for one genome'''

    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        genome_id = options['g']
        # Check genome ID
        if genome_id == '':
            raise CommandError('Genome name required')
        print('Looking for genome', genome_id)
        genome_set = Genome.objects.filter(name=genome_id)
        if genome_set.count() == 0:
            print('Genome ' + genome_id + ' not found')
            raise CommandError()
        elif genome_set.count() > 1:
            print('Non-unique genome name: ' + genome_id)
            raise CommandError()
        print('Genome found:', genome_id)
        genome = genome_set[0]
        
        # Configure importer
        importer = Importer()
        importer.inputgenomes[genome_id]['gbk'] = genome.gbk_filepath
        importer.inputgenomes[genome_id]['url'] = genome.external_url
        importer.inputgenomes[genome_id]['external_id'] = genome.external_id
        if genome.strain is None:
            importer.inputgenomes[genome_id]['strain'] = ''
        else:
            importer.inputgenomes[genome_id]['strain'] = genome.strain.strain_id
        if genome.sample is None:
            importer.inputgenomes[genome_id]['sample'] = ''
        else:
            importer.inputgenomes[genome_id]['sample'] = genome.sample.sample_id
        
        importer.export_jbrowse_genome_data(genome_id)
