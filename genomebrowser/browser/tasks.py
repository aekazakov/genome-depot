# Code implementing admin tasks
# Functions should be called from async_tasks.py
import os
import shutil
import uuid
import zipfile
from browser.dataimport.importer import Importer, download_ncbi_assembly
from browser.dataimport.annotator import Annotator
from browser.models import Strain, Sample, Genome, Protein, Config

def test_task_impl(request, genome_names):
    print(request)
    print(genome_names)
    return 'Genomes:' + genome_names

def import_genomes_impl(args):
    lines, email, zip_file = args
    print ('Asynchronous task received. Starting import.')
    lines_import = []
    temp_dir = Config.objects.get(param='cgcms.temp_dir').value
    upload_dir = os.path.join(temp_dir, str(uuid.uuid4()))
    if zip_file is None:
        zip_content = {}
    else:
        zip_content = handle_zip_upload(zip_file, upload_dir)
    print('Zip archive', zip_content)
    for line in lines:
        row = line.split('\t')
        if row[0] in zip_content:
            row[0] = os.path.join(zip_content[row[0]])
        elif row[0] == '' and row[-1].startswith('NCBI:'):
            row[0] = download_ncbi_assembly(row[-1][5:].rstrip('\n\r'), email, upload_dir)
        elif not os.path.exists(row[0]):
            print(row[0], 'not found')
        lines_import.append('\t'.join(row))
    
    print(lines)
    importer = Importer()
    result = importer.import_genomes(lines_import)
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

def handle_zip_upload(zipf, upload_dir):
    result = {}
    if not os.path.exists(upload_dir):
        os.mkdir(upload_dir)
    with zipfile.ZipFile(zipf, "r", zipfile.ZIP_STORED) as openzip:
        filelist = openzip.infolist()
        for member in filelist:
            if member.is_dir():
                continue
            openzip.extract(member.filename, path=upload_dir)
            result[str(member.filename)] = os.path.join(upload_dir, member.filename)
    return result

