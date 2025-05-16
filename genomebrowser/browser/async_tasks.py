import logging
from django.utils import timezone
from django.utils import dateformat
from django_q.tasks import async_task
from browser.tasks import import_genomes_impl
from browser.tasks import delete_genomes_impl
from browser.tasks import update_static_files_impl
from browser.tasks import import_sample_metadata_impl
from browser.tasks import import_sample_descriptions_impl
from browser.tasks import update_strain_metadata_impl
from browser.tasks import import_annotations_impl
from browser.tasks import import_regulon_impl
from browser.tasks import run_annotation_pipeline_impl

logger = logging.getLogger("GenomeDepot")

def async_run_annotation_tools(request, queryset, plugins):
    logger.debug('Starting asynchronous task for running annotation tools')
    genomes = {item.name:item.gbk_filepath for item in queryset}
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'run-annotation-pipeline.' + timestamp
    args = (genomes, plugins, task_name)
    task_id = async_task(run_annotation_pipeline_impl,
                         args,
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_import_genomes(lines, email):
    logger.debug('Starting asynchronous task for genome import')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'genome-import.' + timestamp
    for item in lines:
        logger.debug(item)
    task_id = async_task(import_genomes_impl,
                         (lines, email, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_delete_genomes(request, queryset):
    logger.debug('Starting asynchronous task of deleting genomes')
    genomes = [item.name for item in queryset]
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'delete-genomes.' + timestamp
    task_id = async_task(delete_genomes_impl,
                         (genomes, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_update_static_files(request, queryset):
    logger.debug('Starting asynchronous task of static file update')
    genomes = [item.name for item in queryset]
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'update-static.' + timestamp
    task_id = async_task(update_static_files_impl,
                         (genomes, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_import_sample_metadata(lines):
    logger.debug('Starting asynchronous import of sample metadata')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'strain-metadata-import.' + timestamp
    task_id = async_task(import_sample_metadata_impl,
                         (lines, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_import_sample_descriptions(lines):
    logger.debug('Starting asynchronous import of strain descriptions')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'strain-descriptions-import.' + timestamp
    task_id = async_task(import_sample_descriptions_impl,
                         (lines, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name
    
def async_update_strain_metadata(xlsx_file):
    logger.debug('Starting asynchronous import of strain descriptions')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'strain-descriptions-import.' + timestamp
    task_id = async_task(update_strain_metadata_impl,
                         (xlsx_file, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name
    
def async_import_annotations(lines):
    logger.debug('Starting asynchronous task for annotations list:')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'annotation-import.' + timestamp
    for item in lines:
        logger.debug(item)
    task_id = async_task(import_annotations_impl,
                         (lines, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name

def async_import_regulon(lines):
    logger.debug('Starting asynchronous task for regulons list:')
    timestamp = dateformat.format(timezone.localtime(timezone.now()),'Y-m-d_H:i:s',)
    task_name = 'regulon-import.' + timestamp
    for item in lines:
        logger.debug(item)
    task_id = async_task(import_regulon_impl,
                         (lines, task_name),
                         task_name=task_name,
                         sync=False
                         )
    logger.debug('Asynchronous task was started with ID ' + task_id +
        ' and name ' + task_name)
    return task_name
