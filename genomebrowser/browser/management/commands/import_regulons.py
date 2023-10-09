from django.core.management.base import BaseCommand
from browser.pipeline.annotate import Annotator

class Command(BaseCommand):
    help = """This command imports regulons from a tab-separated file.
    Genomes must be uploaded into the CGCMS database in advance. 

    The input file must contain the following fields:
    1. Regulon name.
    2. Genome name  (as in the database).
    3. Regulatory gene locus_tag.
    4. Target gene locus_tag.
    5. Contig ID (as in the database).
    6. Site start.
    7. Site end.
    8. Site strand.
    9. Site sequence.
    """

    def add_arguments(self, parser):
        parser.add_argument('-i', help='Path to the input file')

    def handle(self, *args, **options):
        annotator = Annotator()
        lines = []
        with open(options['i'], 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                lines.append(line)
        annotator.add_regulons(lines)
