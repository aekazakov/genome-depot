import time
from django.utils import timezone
from django_q.tasks import async_task, result
from browser.tasks import test_task_impl, import_genomes_impl, delete_genomes_impl, update_static_files_impl, import_sample_metadata_impl, import_sample_descriptions_impl, update_strain_metadata_impl, import_annotations_impl, import_regulon_impl

def test_async_task(request, queryset):
    timeout = 60
    task_id = async_task(test_task_impl, str(request), ';'.join([item.name for item in queryset]), sync=False)
    print(task_id + ' submitted')
    task_result = None
    while timeout > 0 :
        print('Calling ' + task_id + ' time left:' + str(timeout))
        task_result = result(task_id)
        if task_result is not None:
            break
        time.sleep(2)
        timeout -= 2
    if task_result is not None:
        print(task_result)
    else:
        print('Time out! No result returned in 60 seconds.')

def async_import_genomes(lines, email):
    print ('Starting asynchronous task for genome list:')
    timestamp = str(timezone.now())
    task_name = 'genome-import-' + timestamp.replace(' ', '_')
    for item in lines:
        print(item)
    task_id = async_task(import_genomes_impl, (lines, email), task_name=task_name, sync=False)
    return task_id

def async_delete_genomes(request, queryset):
    genomes = [item.name for item in queryset]
    timestamp = str(timezone.now())
    task_name = 'delete-genomes-' + timestamp.replace(' ', '_')
    task_id = async_task(delete_genomes_impl, genomes, task_name=task_name, sync=False)
    return task_name

def async_update_static_files(request, queryset):
    genomes = [item.name for item in queryset]
    timestamp = str(timezone.now())
    task_name = 'update-static-' + timestamp.replace(' ', '_')
    task_id = async_task(update_static_files_impl, genomes, task_name=task_name, sync=False)
    return task_name

def async_import_sample_metadata(lines):
    print ('Starting asynchronous import of sample metadata')
    timestamp = str(timezone.now())
    task_name = 'strain-metadata-import-' + timestamp.replace(' ', '_')
    task_id = async_task(import_sample_metadata_impl, lines, task_name=task_name, sync=False)
    return task_name

def async_import_sample_descriptions(lines):
    print ('Starting asynchronous import of strain descriptions')
    timestamp = str(timezone.now())
    task_name = 'strain-descriptions-import-' + timestamp.replace(' ', '_')
    task_id = async_task(import_sample_descriptions_impl, lines, task_name=task_name, sync=False)
    return task_name
    
def async_update_strain_metadata(xlsx_file):
    print ('Starting asynchronous import of strain descriptions')
    timestamp = str(timezone.now())
    task_name = 'strain-descriptions-import-' + timestamp.replace(' ', '_')
    task_id = async_task(update_strain_metadata_impl, xlsx_file, task_name=task_name, sync=False)
    return task_name
    
def async_import_annotations(lines):
    print ('Starting asynchronous task for annotations list:')
    timestamp = str(timezone.now())
    task_name = 'annotation-import-' + timestamp.replace(' ', '_')
    for item in lines:
        print(item)
    task_id = async_task(import_annotations_impl, lines, task_name=task_name, sync=False)
    return task_id

def async_import_regulon(lines):
    print ('Starting asynchronous task for regulons list:')
    timestamp = str(timezone.now())
    task_name = 'regulon-import-' + timestamp.replace(' ', '_')
    for item in lines:
        print(item)
    task_id = async_task(import_regulon_impl, lines, task_name=task_name, sync=False)
    return task_id
