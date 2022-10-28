# Code implementing admin tasks
# Functions should be called from async_tasks.py
import os
import sys
import shutil
import uuid
import zipfile
from django.conf import settings
from django.core.mail import mail_admins, send_mail
from browser.dataimport.importer import Importer, download_ncbi_assembly
from browser.dataimport.annotator import Annotator
from browser.models import Strain, Sample, Genome, Protein, Config

def test_task_impl(request, genome_names):
    print(request)
    print(genome_names)
    try:
        print(genome_names[999])
        subject = 'CGCMS test task finished'
        message = f'Test task finished successfuly at {settings.BASE_URL}'  #f'Hi {user.username}, thank you for registering in geeksforgeeks.'
        mail_admins(subject, message)
        #send_mail(subject, message, settings.EMAIL_HOST_USER, recipient_list)
    except Exception:
        subject = 'CGCMS test task finished with error'
        message = f'CGCMS test task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise
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
    try:
        result = importer.import_genomes(lines_import)
        subject = 'CGCMS task finished successfuly'
        message = f'"Import Genomes" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Import Genomes" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise
    return result

def update_static_files_impl(genome_ids):
    '''
        Deletes and re-creates Jbrowse static files for an input list of genome IDs, then deletes and re-creates search databases.
    '''
    importer = Importer()
    for genome_id in genome_ids:
        genome = Genome.objects.filter(id=genome_id)
        if genome is not None:
            self.inputgenomes[genome.name]['strain'] = genome.strain.strain_id
            self.inputgenomes[genome.name]['sample'] = genome.sample.sample_id
            self.inputgenomes[genome.name]['gbk'] = genome.gbk_filepath
            self.inputgenomes[genome.name]['url'] = genome.external_url
            self.inputgenomes[genome.name]['external_id'] = genome.external_id
    try:
        importer.export_jbrowse_data()
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
        subject = 'CGCMS task finished successfuly'
        message = f'"Update static files" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Update static files" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

def delete_genomes_impl(genome_names):
    '''
        Deletes genomes from an input list of genomes, then deletes and re-creates search databases
    '''
    try:
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
        subject = 'CGCMS task finished successfuly'
        message = f'"Delete genomes" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Delete genomes" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

    
def import_sample_metadata_impl(lines):
    print ('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_sample_metadata(lines)
        subject = 'CGCMS task finished successfuly'
        message = f'"Import sample metadata" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Import sample metadata" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

def import_sample_descriptions_impl(lines):
    print ('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_sample_descriptions(lines)
        subject = 'CGCMS task finished successfuly'
        message = f'"Import sample descriptions" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Import sample descriptions" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

def update_strain_metadata_impl(xlsx_file):
    print ('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_strain_metadata(xlsx_path=None, xlsx_file=xlsx_file)
        subject = 'CGCMS task finished successfuly'
        message = f'"Update strain metadata" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Update strain metadata" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise
    
def import_annotations_impl(lines):
    print ('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.import_annotations(lines)
        subject = 'CGCMS task finished successfuly'
        message = f'"Import annotations" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Import annotations" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

def import_regulon_impl(lines):
    print ('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_regulons(lines)
        subject = 'CGCMS task finished successfuly'
        message = f'"Import regulon" task finished successfuly at {settings.BASE_URL}'
        mail_admins(subject, message)
    except Exception:
        subject = 'CGCMS task finished with error'
        message = f'CGCMS "Import regulon" task at {settings.BASE_URL} finished with error.\nError:{sys.exc_info()[0]}. {sys.exc_info()[1]}, {sys.exc_info()[2].tb_frame.f_code.co_filename}:{sys.exc_info()[2].tb_lineno}'
        mail_admins(subject, message)
        raise

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

