from django.core.management.base import BaseCommand
from browser.pipeline.annotate import Annotator

class Command(BaseCommand):
    help = '''
    For protein sequences uploaded to Django database, 
    this program runs hmmsearch with PFAM and TIGRFAM HMMs,
    and writes domain annotations to database.
    Existing domain annotations are deleted.'''

    def handle(self, *args, **options):
        annotator = Annotator()
        annotator.add_pfam_domains()
        annotator.add_tigrfam_domains()
