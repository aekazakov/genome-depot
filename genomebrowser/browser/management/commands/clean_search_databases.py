import os
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer

class Command(BaseCommand):
    help = 'This command cleans nucleotide and protein search databases after deletion of a genome from admin interface'
    
    def handle(self, *args, **options):
        importer = Importer()
        self._clean_nucl_database(importer)
        importer.export_proteins()
        importer.delete_search_databases()
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_nucl'])), importer.config['cgcms.search_db_nucl'])
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_prot'])), importer.config['cgcms.search_db_prot'])
        importer.create_search_databases()
        print('Done!')
        
        
    def _clean_nucl_database(self, importer):
        """Copies nucleotide database from static directory for appending. 
        Returns unique sequence IDs.
        """
        contigs = set()
        for item in Contig.objects.values('contig_id', 'genome__name'):
            contig_uid = item['contig_id'] + '|' + item['genome__name']
            contigs.add(contig_uid)

        nucl_db_fasta = importer.config['cgcms.search_db_nucl']
        backup_file = os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(nucl_db_fasta))
        save_flag = False
        contigs_saved = set()
        with open(backup_file, 'w') as outfile:
            with open(nucl_db_fasta, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        contig_uid = line[1:].rstrip('\n\r')
                        if contig_uid in contigs:
                            save_flag = True
                            contigs_saved.add(contig_uid)
                        else:
                            save_flag = False
                    if save_flag:
                        outfile.write(line)
        # Sanity check before making chenges in the database or deleting files
        if contigs.difference(contigs_saved):
            raise ValueError('Nucleotide DB file corrupted. Add contigs: ' + str(contigs.difference(contigs_saved)))
        if contigs_saved.difference(contigs):
            raise ValueError('Nucleotide DB file corrupted. Remove contigs: ' + str(contigs_saved.difference(contigs)))
