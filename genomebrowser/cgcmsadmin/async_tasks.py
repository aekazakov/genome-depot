import time
from django.utils import timezone
from django_q.tasks import async_task, result
from cgcmsadmin.tasks import test_task_impl, import_genomes_impl

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


def async_import_genomes(lines):
    print ('Starting asynchronous task for genome list:')
    timestamp = str(timezone.now())
    task_name = 'genome-import-' + timestamp.replace(' ', '_')
    for item in lines:
        print(item)
    task_id = async_task(import_genomes_impl, lines, task_name=task_name, sync=False)
    print(task_id + ' submitted')
    return task_id


