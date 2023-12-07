import gzip
import os
import shutil
from Bio import GenBank
from django.core.management.base import BaseCommand
from browser.pipeline.genome_import import Importer
from browser.models import Genome

class Command(BaseCommand):
    help = '''
    Removes and re-generates nucleotide and protein search
    databases. Use it if the files are missing or corrupted,
    or if the genome import pipeline crashed before creating the 
    search database files.
    '''
    
    def handle(self, *args, **options):
        importer = Importer()
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
