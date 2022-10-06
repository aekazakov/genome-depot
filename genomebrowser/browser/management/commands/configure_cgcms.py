import os
from django.core.management.base import BaseCommand
from browser.models import Config
from browser.dataimport.importer import Importer
from genomebrowser.settings import BASE_URL

class Command(BaseCommand):
    help = 'Imports config settings from file'

    def add_arguments(self, parser):
        parser.add_argument('-i', default='configs.txt', help='Path to parameters file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            print('Importing parameters from', options['i'])
            configs = {}
            with open(options['i'], 'r') as infile:
                for line in infile:
                    key, val = line.rstrip('\n\r').split('=')
                    configs[key.strip()] = val.strip()
            configs_saved = set()
            for item in Config.objects.all():
                if item.param in configs:
                    configs_saved.add(item.param)
                    if item.value != configs[item.param]:
                        item.value = configs[item.param]
                        item.save()
            for param in configs:
                if param not in configs_saved:
                    config = Config.objects.create(param=param, value=configs[param])
                    
            tracklist_template = 'trackList.json.template'
            tracklist_file = os.path.join(Config.objects.get(param='cgcms.json_dir').value, 'trackList.json')
            print('Writing', tracklist_file)
            with open(tracklist_file, 'w') as outfile:
                with open(tracklist_template, 'r') as infile:
                    for line in infile:
                        outfile.write(line.replace('https://example.com/', BASE_URL))
        else:
            raise ValueError('Parameters file not found.')
