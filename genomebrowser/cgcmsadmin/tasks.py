# Code implementing admin tasks
# Functions should be called from async_tasks.py
import os
import shutil
from browser.dataimport.importer import Importer
from browser.dataimport.annotator import Annotator
from browser.models import Strain, Sample, Genome, Protein

def test_task_impl(request, genome_names):
    print(request)
    print(genome_names)
    return 'Genomes:' + genome_names

def import_genomes_impl(lines):
    print ('Asynchronous task received. Starting import.')
    print(lines)
    importer = Importer()
    result = importer.import_genomes(lines)
    return result

def delete_genomes_impl(genome_names):
    importer = Importer()
    for genome_name in genome_names:
        Genome.objects.filter(name=genome_name).delete()
    Protein.objects.filter(gene=None).delete()
    Strain.objects.filter(genome=None).delete()
    Sample.objects.filter(genome=None).delete()
    importer.export_proteins()
    importer.export_contigs()
    importer.delete_search_databases()
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_nucl'])), importer.config['cgcms.search_db_nucl'])
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_prot'])), importer.config['cgcms.search_db_prot'])
    importer.create_search_databases()
    for genome in genome_names:
        shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'], genome))
    os.remove(importer.config['cgcms.search_db_nucl'])
    os.remove(importer.config['cgcms.search_db_prot'])
    
def import_sample_metadata_impl(lines):
    print ('Asynchronous task received. Starting import.')
    annotator = Annotator()
    annotator.add_sample_metadata(lines)

def import_sample_descriptions_impl(lines):
    print ('Asynchronous task received. Starting import.')
    annotator = Annotator()
    annotator.update_sample_descriptions(lines)

def update_strain_metadata_impl(xlsx_file):
    print ('Asynchronous task received. Starting import.')
    annotator = Annotator()
    annotator.update_strain_metadata(xlsx_path=None, xlsx_file=xlsx_file)
    
def import_annotations_impl(lines):
    print ('Asynchronous task received. Starting import.')
    annotator = Annotator()
    annotator.import_annotations(lines)

def import_regulon_impl(lines):
    print ('Asynchronous task received. Starting import.')
    annotator = Annotator()
    annotator.add_regulons(lines)
