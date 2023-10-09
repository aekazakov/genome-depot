import os
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.models import Contig
from browser.pipeline.genome_import import Importer

class Command(BaseCommand):
    help = '''
    This command regenerates nucleotide and protein search
    databases if the database files are missing or corrupted,
    or if the genome import pipeline crashed before creating the 
    search database files.
    '''
    
    def handle(self, *args, **options):
        importer = Importer()
        gbff_dir = os.path.join(importer.config['cgcms.static_dir'], 'gbff')
        nucl_db_fasta = importer.config['cgcms.search_db_nucl']
        os.remove(nucl_db_fasta)
        with open(nucl_db_fasta, 'w') as outfile:
            for genome_data in Genome.objects.values_list('name', 'gbk_filepath'):
                if genome_data[1].endswith('.gz'):
                    gbk_handle = gzip.open(genome_data[1], 'rt')
                else:
                    gbk_handle = open(genome_data[1], 'r')
                parser = GenBank.parse(gbk_handle)
                for gbk_record in parser:
                    contig_sequence = gbk_record.sequence
                    contig_id = gbk_record.locus
                    outfile.write('>' + contig_id + '|' +
                                  genome_data[0] + '\n' +
                                  ''.join(contig_sequence) + '\n'
                                  )
        importer.export_proteins()
        importer.delete_search_databases()
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                        os.path.basename(importer.config['cgcms.search_db_nucl'])),
                        importer.config['cgcms.search_db_nucl']
                        )
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                        os.path.basename(importer.config['cgcms.search_db_prot'])),
                        importer.config['cgcms.search_db_prot']
                        )
        importer.create_search_databases()
        print('Done!')
