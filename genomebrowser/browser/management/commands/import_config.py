import os
from django.core.management.base import BaseCommand
from browser.models import Config

class Command(BaseCommand):
    help = ''''Imports settings from a config file
    (text file with key/value entries separated by "=" symbol)
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='configs.txt',
            help='Path to parameters file'
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Overwrite existing config value'
        )
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
                    if options['overwrite']:
                        if item.value != configs[item.param]:
                            item.value = configs[item.param]
                            item.save()
            for param in configs:
                if param not in configs_saved:
                    _ = Config.objects.create(param=param, value=configs[param])
        else:
            raise ValueError('Parameters file not found.')
