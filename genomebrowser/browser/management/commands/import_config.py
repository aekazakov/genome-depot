from django.core.management.base import BaseCommand
from browser.util import import_config

class Command(BaseCommand):
    help = '''Imports settings from a config file
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
        import_config(options['i'], options['overwrite'])
