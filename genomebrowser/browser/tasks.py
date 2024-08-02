# Code implementing admin tasks
# Functions should be called from async_tasks.py
import os
import sys
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

logger = logging.getLogger("GenomeDepot")

def run_annotation_pipeline_impl(args):
    '''
        Implements the Run Annotation Pipeline task
    '''
    logger.debug('Asynchronous task run_annotation_tools received. Starting the pipeline.')
    (genomes, plugins) = args
    annotator = Annotator()
    messages = []
    subject = None
    for plugin_ind, plugin in enumerate(plugins):
        logger.debug('Starting tool ' + str(plugin_ind + 1) + ' of ' + str(len(plugins)) + ':' + plugin)
        try:
            plugin_conf_enabled = annotator.config['plugins.' + plugin + '.enabled']
            if plugin_conf_enabled in ('1', 'yes', 'Yes', 'y', 'Y'):
                status = annotator.run_external_tools(genomes, plugin_name=plugin)
                messages.append('Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) + ': ' + plugin + ' plugin ' + status)
            else:
                messages.append('Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) + ': ' + plugin + ' plugin disabled')
        except Exception as err:
            messages.append('Step ' + str(plugin_ind + 1) + ' out of ' + str(len(plugins)) + ': ' + plugin + ' plugin ERROR')
            messages.append(traceback.format_exc())
            subject = 'GenomeDepot annotation pipeline error'
            logger.exception('GenomeDeport annotation pipeline raised an unhandled exception at ' + settings.BASE_URL)
    if subject is None:
        subject = 'GenomeDepot annotation pipeline finished'
        message = '"Run annotation tools" task finished successfuly at ' + \
        f'{settings.BASE_URL}' + '\n\n' + '\n'.join(messages)
        ret = 'Annotation pipeline finished.'
    else:
        subject = 'GenomeDepot annotation pipeline error'
        message = '"Run annotation tools" task finished with error at ' + \
        f'{settings.BASE_URL}' + '\n\n' + '\n'.join(messages)
        ret = 'Annotation pipeline error.'
    mail_admins(subject, message)
    return ret
    

def import_genomes_impl(args):
    '''
        Implements the Import Genomes task
    '''
    lines, email = args
    logger.debug('Asynchronous task import_genomes received. Starting import.')
    temp_dir = Config.objects.get(param='core.temp_dir').value
    result = ''
    try:
        upload_dir = os.path.join(temp_dir, str(uuid.uuid4()))
        for line in lines:
            row = line.split('\t')
            if row[0] == '' and row[-1].startswith('NCBI:'):
                row[0] = download_ncbi_assembly(row[-1][5:].rstrip('\n\r'),
                                                email,
                                                upload_dir
                                                )
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
        subject = 'GenomeDepot task finished successfuly'
        message = '"Import Genomes" task finished successfuly at ' + \
        f'{settings.BASE_URL}\n\n' + result
    except Exception as e:
        result = 'Error!'
        subject = 'GenomeDeport genome import pipeline raised an unhandled exception'
        message = f'GenomeDepot "Import Genomes" task at {settings.BASE_URL} finished ' +\
        f'with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport import pipeline raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)
    return result


def update_static_files_impl(genome_names):
    '''
        Deletes and re-creates Jbrowse static files for an input list of
        genome IDs, then deletes and re-creates search databases.
    '''
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
        subject = 'GenomeDepot task finished successfuly'
        message = '"Update static files" task finished successfuly ' +\
        f'at {settings.BASE_URL}'
    except Exception:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Update static files" task at {settings.BASE_URL} finished' +\
        f'with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)

def delete_genomes_impl(genome_names):
    '''
        Deletes genomes from an input list of genomes, then 
        deletes and re-creates search databases
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
        subject = 'GenomeDepot task finished successfuly'
        message = f'"Delete genomes" task finished successfuly at {settings.BASE_URL}'
    except Exception:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Delete genomes" task at {settings.BASE_URL} finished ' +\
        f'with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)

    
def import_sample_metadata_impl(lines):
    '''
        Implements the Import Sample Metadata task
    '''
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_sample_metadata(lines)
        subject = 'GenomeDepot task finished successfuly'
        message = '"Import sample metadata" task finished successfuly ' +\
        f'at {settings.BASE_URL}'
    except Exception:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Import sample metadata" task at {settings.BASE_URL} ' +\
        f'finished with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)


def import_sample_descriptions_impl(lines):
    '''
        Implements the Import Sample Descriptions task
    '''
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_sample_descriptions(lines)
        subject = 'GenomeDepot task finished successfuly'
        message = '"Import sample descriptions" task finished successfuly ' +\
        f'at {settings.BASE_URL}'
    except Exception:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Import sample descriptions" task at {settings.BASE_URL} ' +\
        f'finished with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)


def update_strain_metadata_impl(xlsx_file):
    '''
        Implements the Update Strain Metadata task
    '''
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.update_strain_metadata(xlsx_path=None, xlsx_file=xlsx_file)
        subject = 'GenomeDepot task finished successfuly'
        message = '"Update strain metadata" task finished successfuly ' +\
        f'at {settings.BASE_URL}'
    except Exception as e:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Update strain metadata" task at {settings.BASE_URL} ' +\
        f'finished with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)

    
def import_annotations_impl(lines):
    '''
        Implements the Import Annotations task
    '''
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.import_annotations(lines)
        subject = 'GenomeDepot task finished successfuly'
        message = '"Import annotations" task finished successfuly ' +\
        f'at {settings.BASE_URL}'
    except Exception as e:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Import annotations" task at {settings.BASE_URL} finished' +\
        f' with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)


def import_regulon_impl(lines):
    '''
        Implements the Import Regulon task
    '''
    logger.debug('Asynchronous task received. Starting import.')
    try:
        annotator = Annotator()
        annotator.add_regulons(lines)
        subject = 'GenomeDepot task finished successfuly'
        message = f'"Import regulon" task finished successfuly at {settings.BASE_URL}'
    except Exception as e:
        subject = 'GenomeDepot task finished with error'
        message = f'GenomeDepot "Import regulon" task at {settings.BASE_URL} finished ' +\
        f'with error.\n' + traceback.format_exc()
        logger.exception('GenomeDeport raised an unhandled exception at ' + settings.BASE_URL)
    mail_admins(subject, message)
