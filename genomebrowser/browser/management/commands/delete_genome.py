import os
import shutil
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.importer import Importer

class Command(BaseCommand):
    help = 'Deletes a genome with all genes and annotations from the database'
    
    def add_arguments(self, parser):
        parser.add_argument('-g', default='', help='Genome name')

    def handle(self, *args, **options):
        genome_name = options['g']
        if genome_name == '':
            raise CommandError('Genome name required')
        print('Looking for genome', genome_name)
        genome_set = Genome.objects.filter(name=genome_name)
        if genome_set.count() == 0:
            print('Genome ' + genome_name + ' not found')
            raise CommandError()
        elif genome_set.count() > 1:
            print('Not unique genome name: ' + genome_name)
            raise CommandError()
        importer = Importer()
        print('Note 1: Genome deletion removes contigs, genes and gene annotations for this genome. It also removes proteins not linked to any gene. But it does not remove a strain associated with this genome if the strain has other genomes.')
        self._clean_nucl_database(genome_name, importer)
        importer.export_proteins()
        importer.delete_search_databases()
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_nucl'])), importer.config['cgcms.search_db_nucl'])
        shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_prot'])), importer.config['cgcms.search_db_prot'])
        print('Deleting genome...')
        genome_set.delete()
        print('Deleting proteins not linked to genes...')
        Protein.objects.filter(gene=None).delete()
        print('Deleting strains not linked to genomes...')
        Strain.objects.filter(genome=None).delete()
        print('Deleting samples not linked to genomes...')
        Sample.objects.filter(genome=None).delete()
        importer.create_search_databases()
        shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'], genome_name))
        # importer.cleanup()
        print('Done!')
        
        
    def _clean_nucl_database(self, genome_name, importer):
        """Copies nucleotide database from static directory for appending. 
        Returns unique sequence IDs.
        """
        result = set()
        contigs_remove = set()
        for item in Contig.objects.filter(genome__name = genome_name).values('contig_id'):
            contigs_remove.add(item['contig_id'] + '|' + genome_name)
        nucl_db_fasta = importer.config['DEFAULT']['search_db_nucl']
        backup_file = os.path.join(importer.config['DEFAULT']['temp_dir'], os.path.basename(nucl_db_fasta))
        save_flag = True
        contigs_saved = set()
        with open(backup_file, 'w') as outfile:
            with open(nucl_db_fasta, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        contig_uid = line[1:].rstrip('\n\r')
                        if contig_uid in contigs_remove:
                            save_flag = False
                        else:
                            save_flag = True
                            contigs_saved.add(contig_uid)
                    if save_flag:
                        outfile.write(line)
        # Sanity check before making chenges in the database or deleting files
        contig_data = Contig.objects.values('contig_id', 'genome__name')
        contigs = set()
        for item in contig_data:
            if item['genome__name'] != genome_name:
                contig_uid = item['contig_id'] + '|' + item['genome__name']
                contigs.add(contig_uid)
        if contigs.difference(contigs_saved):
            raise ValueError('Nucleotide DB file corrupted. Add contigs: ' + str(contigs.difference(contigs_saved)))
        if contigs_saved.difference(contigs):
            raise ValueError('Nucleotide DB file corrupted. Remove contigs: ' + str(contigs_saved.difference(contigs)))
        return result
