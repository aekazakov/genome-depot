# Code implementing admin tasks
# Functions should be called from async_tasks.py
import os
import gzip
import shutil
import uuid
import logging
import traceback
from Bio import GenBank
from django.conf import settings
from django.core.mail import mail_admins
from browser.util import download_ncbi_assembly
from browser.pipeline.genome_import import Importer
from browser.pipeline.annotate import Annotator
from browser.models import Strain, Sample, Genome, Protein, Config
from browser.admin_utils import is_email_configured

logger = logging.getLogger("GenomeDepot")

def run_annotation_pipeline_impl(args):
    '''
        Implements the Run Annotation Pipeline task
    '''
    logger.debug(
        'Asynchronous task run_annotation_tools received. Starting the pipeline.'
    )
    (genomes, plugins, task_name) = args
    annotator = Annotator()
    messages = [task_name + ' task report',]
    subject = ''
    for plugin_ind, plugin in enumerate(plugins):
        logger.debug(
            f'Starting tool {str(plugin_ind + 1)} of {str(len(plugins))}: {plugin}'
        )
        try:
            plugin_conf_enabled = annotator.config['plugins.' + plugin + '.enabled']
            if plugin_conf_enabled in ('1', 'yes', 'Yes', 'y', 'Y'):
                status = annotator.run_external_tools(genomes, plugin_name=plugin)
                messages.append(
                    'Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) +
                    ': ' + plugin + ' plugin ' + status
                )
            else:
                messages.append(
                    'Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) +
                    ': ' + plugin + ' plugin disabled'
                )
        except Exception:
            messages.append(
                'Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) +
                ': ' + plugin + ' plugin ERROR'
            )
            messages.append(traceback.format_exc())
            subject = f'{task_name} task finished with error'
            logger.exception(
                f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
            )
    if subject == '':
        subject = f'{task_name} task finished successfuly'
        message = '"Run annotation tools" task finished successfuly at ' + \
        f'{settings.BASE_URL}' + '\n\n' + '\n'.join(messages)
        ret = 'Annotation pipeline finished.'
    else:
        subject = f'{task_name} task finished with error'
        message = '"Run annotation tools" task finished with error at ' + \
        f'{settings.BASE_URL}' + '\n\n' + '\n'.join(messages)
        ret = 'Annotation pipeline error.'
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)
    return ret
    

def import_genomes_impl(args):
    '''
        Implements the Import Genomes task
    '''
    lines, email, task_name = args
    logger.debug('Asynchronous task import_genomes received. Starting import.')
    temp_dir = Config.objects.get(param='core.temp_dir').value
    result = ''
    try:
        upload_dir = os.path.join(temp_dir, str(uuid.uuid4()))
        for i, line in enumerate(lines):
            row = line.split('\t')
            if row[0] == '' and row[-1].startswith('NCBI:'):
                assembly_path = download_ncbi_assembly(row[-1][5:].rstrip('\n\r'),
                                                email,
                                                upload_dir
                                                )
                if assembly_path == '':
                    logger.warning(
                        'Genome not found for NCBI identifier ' + row[-1][5:]
                    )
                else:
                    logger.debug('File uploaded to ' + assembly_path)
                    lines[i] = assembly_path + line
        genome_import_batch_size = 0
        genome_import_batch_limit = 1000
        genome_import_batch = []
        genome_batch_count = 0
        for line in lines: 
            genome_import_batch_size += 1
            genome_import_batch.append(line)
            if genome_import_batch_size == genome_import_batch_limit:
                genome_batch_count += 1
                logger.debug('Importing genome batch %d', genome_batch_count)
                importer = Importer()
                result += importer.import_genomes(genome_import_batch) + '\n'
                genome_import_batch_size = 0
                genome_import_batch = []
        if genome_import_batch:
            genome_batch_count += 1
            logger.debug('Importing genome batch %d', genome_batch_count)
            importer = Importer()
            result += importer.import_genomes(genome_import_batch) + '\n'
        if 'error' in result:
            subject = f'{task_name} task finished with error'
            message = '"Import Genomes" task finished with error at ' + \
            f'{settings.BASE_URL}\n\n' + result
        else:
            subject = f'{task_name} task finished successfuly'
            message = (
                f'{task_name} report\n"Import Genomes" task finished successfuly ' +
                f'at {settings.BASE_URL}\n\n{result}'
            )
    except Exception:
        result = 'Error!'
        subject = f'{task_name} task raised an unhandled exception'
        message = (
            f'{task_name} report\nGenomeDepot "Import Genomes" task at ' +
            f'{settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)
    return result


def update_static_files_impl(args):
    '''
        Deletes and re-creates Jbrowse static files for an input list of
        genome IDs, then deletes and re-creates search databases.
    '''
    genome_names, task_name = args
    importer = Importer()
    try:
        for genome_name in genome_names:
            genome = Genome.objects.get(name=genome_name)
            if genome is not None:
                if genome.strain:
                    importer.inputgenomes[genome.name]['strain'] = \
                    genome.strain.strain_id
                else:
                    importer.inputgenomes[genome.name]['strain'] = ''
                if genome.sample:
                    importer.inputgenomes[genome.name]['sample'] = \
                    genome.sample.sample_id
                else:
                    importer.inputgenomes[genome.name]['sample'] = ''
                importer.inputgenomes[genome.name]['gbk'] = genome.gbk_filepath
                importer.inputgenomes[genome.name]['url'] = genome.external_url
                importer.inputgenomes[genome.name]['external_id'] = genome.external_id
            genome_fasta = os.path.join(importer.config['core.temp_dir'],
                                        genome_name + '.contigs.fasta'
                                        )
            gbk_file = genome.gbk_filepath
            if gbk_file.endswith('.gz'):
                gbk_handle = gzip.open(gbk_file, 'rt')
            else:
                gbk_handle = open(gbk_file, 'r')
            parser = GenBank.parse(gbk_handle)
            with open(genome_fasta, 'w') as outfile:
                for gbk_record in parser:
                    outfile.write('>' + gbk_record.locus + '\n' + 
                                  ''.join(gbk_record.sequence) + '\n'
                                  )
            gbk_handle.close()
        importer.export_jbrowse_data()
        importer.export_proteins()
        importer.export_contigs()
        importer.delete_search_databases()
        shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                        os.path.basename(importer.config['core.search_db_nucl'])),
                        importer.config['core.search_db_nucl']
                        )
        shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                        os.path.basename(importer.config['core.search_db_prot'])),
                        importer.config['core.search_db_prot']
                        )
        importer.create_search_databases()
        os.remove(importer.config['core.search_db_nucl'])
        os.remove(importer.config['core.search_db_prot'])
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Update static files" task finished successfuly ' +\
            f'at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\n"Update static files" task at {settings.BASE_URL} ' +
            f'finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)

def delete_genomes_impl(args):
    '''
        Deletes genomes from an input list of genomes, then 
        deletes and re-creates search databases
    '''
    genome_names, task_name = args
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
        shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                        os.path.basename(importer.config['core.search_db_nucl'])),
                        importer.config['core.search_db_nucl']
                        )
        shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                        os.path.basename(importer.config['core.search_db_prot'])),
                        importer.config['core.search_db_prot']
                        )
        importer.create_search_databases()
        for genome in genome_names:
            if os.path.exists(os.path.join(importer.config['core.json_dir'], genome)):
                shutil.rmtree(os.path.join(importer.config['core.json_dir'], genome))
        os.remove(importer.config['core.search_db_nucl'])
        os.remove(importer.config['core.search_db_prot'])
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Delete genomes" task finished ' +
            f'successfuly at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report \nGenomeDepot "Delete genomes" task at ' +
            f'{settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)

    
def import_sample_metadata_impl(args):
    '''
        Implements the Import Sample Metadata task
    '''
    lines, task_name = args
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_sample_metadata(lines)
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Import sample metadata" task finished ' +\
            f'successfuly at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\nGenomeDepot "Import sample metadata" task ' +
            f'at {settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)


def import_sample_descriptions_impl(args):
    '''
        Implements the Import Sample Descriptions task
    '''
    lines, task_name = args
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_sample_descriptions(lines)
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Import sample descriptions" task finished ' +
            f'successfuly at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\nGenomeDepot "Import sample descriptions" task at ' +
            f'{settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)


def update_strain_metadata_impl(args):
    '''
        Implements the Update Strain Metadata task
    '''
    xlsx_file, task_name = args
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_strain_metadata(xlsx_path=None, xlsx_file=xlsx_file)
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Update strain metadata" task finished successfuly ' +
            f'at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\nGenomeDepot "Update strain metadata" task ' +
            f'at {settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)

    
def import_annotations_impl(args):
    '''
        Implements the Import Annotations task
    '''
    lines, task_name = args
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.import_annotations(lines)
        subject = f'{task_name} task finished successfuly'
        message = (
            f'"Import annotations" task finished successfuly at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\nGenomeDepot "Import annotations" task at ' +
            f'{settings.BASE_URL} finished with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            f'{task_name} raised an unhandled exception at {settings.BASE_URL}'
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)


def import_regulon_impl(args):
    '''
        Implements the Import Regulon task
    '''
    lines, task_name = args
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_regulons(lines)
        subject = f'{task_name} task finished successfuly'
        message = (
            f'{task_name} report\n"Import regulon" task finished ' +
            f'successfuly at {settings.BASE_URL}'
        )
    except Exception:
        subject = f'{task_name} task finished with error'
        message = (
            f'{task_name} report\nGenomeDepot "Import regulon" task at ' + \
            f' {settings.BASE_URL} finished  with error.\n{traceback.format_exc()}'
        )
        logger.exception(
            task_name + ' raised an unhandled exception at ' + settings.BASE_URL
        )
    if is_email_configured():
        mail_admins(subject, message)
    else:
        logger.warning('Email backend is not configured. All reports are redirected to the log file.')
        logger.info(subject)
        logger.info(message)
