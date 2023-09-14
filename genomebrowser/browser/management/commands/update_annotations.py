import os
from django.core.management.base import BaseCommand, CommandError
from browser.models import Genome
from browser.dataimport.annotator import Annotator


class Command(BaseCommand):
    help = '''
    For protein sequences from selected genomes uploaded to Django database,
    this program runs one of annotation tool plugins.
    '''

    def add_arguments(self, parser):
        parser.add_argument('-i',
                            default='genomes.txt',
                            help='Path to input file'
                            )
        parser.add_argument('-t',
                            help='Annotation tools'
                            )

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            annotator = Annotator()
            tool = options['t']
            if tool == 'hmmsearch':
                genome_names = []
                with open(options['i'], 'r') as infile:
                    for line in infile:
                        if line.startswith('#'):
                            continue
                        filepath, genome_name, _, _, _, _ = \
                        line.rstrip('\n\r').split('\t')
                        genome_names.append(genome_name)
                genomes=[item['id'] for item in Genome.objects.filter(
                         name__in = genome_names).values('id')
                         ]
                annotator.update_pfam_domains(genomes)
                annotator.update_tigrfam_domains(genomes)
            else:
                plugins_available = set()
                for param in annotator.config:
                    if param.startswith('plugins.') and annotator.config[param] != '':
                        plugin = param.split('.')[1]
                        plugins_available.add(plugin)
                if tool not in plugins_available:
                    plugins_available.add('hmmsearch')
                    print('Available tools:', ','.join(list(plugins_available)))
                    raise CommandError('Annotation tool', tool, 'not found.')
                    #genomes[genome_name] = filepath
                genomes = {}
                with open(options['i'], 'r') as infile:
                    for line in infile:
                        if line.startswith('#'):
                            continue
                        filepath, genome_name, _, _, _, _ = \
                        line.rstrip('\n\r').split('\t')
                        genomes[genome_name] = Genome.objects.get(
                                                                  name=genome_name
                                                                  ).gbk_filepath
                annotator.run_external_tools(genomes, plugin=tool)
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
