from django.core.management.base import BaseCommand
from browser.pipeline.taxonomy import update_taxonomy

class Command(BaseCommand):
    help = """
    Updates taxonomy with new data from NCBI.
    """

    def handle(self, *args, **options):
        update_taxonomy()
