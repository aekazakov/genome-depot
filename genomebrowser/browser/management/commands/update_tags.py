import os
from django.core.management.base import BaseCommand, CommandError
from browser.models import Tag, Genome

class Command(BaseCommand):
    help = '''
    For protein sequences from selected genomes uploaded to Django
    database, this program runs one of annotation tool plugins.
    '''

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')
        parser.add_argument('-t', help='Tags, comma-separated')

    def handle(self, *args, **options):
        if not os.path.exists(options['i']):
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
        
        # create tags if missing
        tag_names = options['t'].split(',')
        tags = {}
        for tag_name in tag_names:
            tag_name = tag_name.strip()
            if not tag_name.isalnum():
                raise CommandError('Tag must contain only alphabet letter (a-z)' +
                                    'and numbers (0-9). Correct the tag ' + tag_name)
            try:
                tag = Tag.objects.get(name=tag_name)
                tags[tag_name] = tag
            except Tag.DoesNotExist:
                tag = Tag(name = tag_name, description = '')
                tag.save()
                tags[tag_name] = tag
        # add tags to genome
        with open(options['i'], 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                row = line.rstrip('\n\r').split('\t')
                if len(row) == 1:
                    genome_name = row[0]
                else:
                    filepath, genome_name, _, _, _, _ = line.rstrip('\n\r').split('\t')
                try:
                    genome = Genome.objects.get(name = genome_name)
                    for tag_name, tag in tags.items():
                        genome.tags.add(tag)
                    genome.save()
                except Genome.DoesNotExist:
                    print('Genome not found:', genome_name)
                    continue
