from django.core.management.base import BaseCommand
from browser.util import export_config

class Command(BaseCommand):
    help = '''Exports configuration parameters into 
    a config file (text file with key/value entries separated by "=" symbol)'''

    def add_arguments(self, parser):
        parser.add_argument('-o',
                            default='exported_config.txt',
                            help='Output file name'
                            )

    def handle(self, *args, **options):
        export_config(options['o'])
