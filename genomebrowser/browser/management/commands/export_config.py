import os
import sys
from django.core.management.base import BaseCommand
from browser.models import Config

class Command(BaseCommand):
    help = '''This command exports configuration parameters into 
    a config file (text file with key/value entries separated by "=" symbol)'''

    def add_arguments(self, parser):
        parser.add_argument('-o',
                            default='exported_config.txt',
                            help='Output file name'
                            )

    def handle(self, *args, **options):
        if os.path.exists(options['o']):
            print('Output file already exists')
            sys.exit()
        with open(options['o'], 'w') as outfile:
            for item in Config.objects.values_list('param', 'value'):
                outfile.write('='.join((str(x) for x in item)) + '\n')
