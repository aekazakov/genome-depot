import os
import logging
from django.core.management.base import BaseCommand, CommandError
from browser.models import Genome
from browser.pipeline.annotate import Annotator
logger = logging.getLogger("GenomeDepot")


class Command(BaseCommand):
    help = '''
    For protein-coding genes from input genomes, this command runs one,
    several or all annotation tools.
    Existing annotations for these tools in the input genomes will be deleted.

    Input is either a file with list of genome names or a tab-separated 
    file with six columns:
    - path to Genbank file
    - genome ID (no spaces)
    - strain name (no spaces)
    - sample ID (no spaces)
    - URL (link to NCBI genome assembly etc.)
    - external ID (for example, "NCBI:GCF_000006945.2")
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='genomes.txt',
            help='Path to input file'
        )
        parser.add_argument(
            '-t',
            help='Annotation tools'
        )
        parser.add_argument(
            '--all',
            action='store_true',
            help='Ignore -t option and run all annotation tools.'
        )

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            genomes = {}
            with open(options['i'], 'r') as infile:
                for line in infile:
                    if line.startswith('#'):
                        continue
                    row = line.rstrip('\n\r').split('\t')
                    if len(row) == 1:
                        genome_name = row[0]
                    else:
                        genome_name = row[1]
                    try:
                        genomes[genome_name] = Genome.objects.get(
                            name=genome_name
                        ).gbk_filepath
                    except Genome.DoesNotExist:
                        logger.error('Genome not found:' + str(genome_name))
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
        annotator = Annotator()
        if options['all']:
            annotator.run_annotation_pipeline(genomes)
        else:
            tool = options['t']
            plugins_available = set()
            for param in annotator.config:
                if param.startswith('plugins.') and annotator.config[param] != '':
                    plugin = param.split('.')[1]
                    plugins_available.add(plugin)
            if tool not in plugins_available:
                print('Available tools:', ','.join(list(plugins_available)))
                raise CommandError('Annotation tool', tool, 'not found.')
            annotator.run_external_tools(genomes, plugin_name=tool)
