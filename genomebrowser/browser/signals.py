from django.db.models.signals import post_save, pre_delete
from django.dispatch import receiver
from django.core import serializers
from django.utils import timezone
from browser.models import ChangeLog
from browser.models import Genome, Strain, Sample

# this receiver is executed every-time some data is saved in Genome table
@receiver(post_save, sender=Genome)
def genome_saved(sender, instance, created, **kwargs):
    # code to execute before every model save
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Genome created:", instance.name)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Genome saved:", instance.name)


@receiver(pre_delete, sender=Genome)
def genome_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    print("Genome deleted:", instance.name)


# Strain created, updated or deleted
@receiver(post_save, sender=Strain)
def srain_saved(sender, instance, created, **kwargs):
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Strain created:", instance.strain_id)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Strain saved:", instance.strain_id)


@receiver(pre_delete, sender=Strain)
def strain_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    print("Strain deleted:", instance.strain_id)


# Sample created, updated or deleted
@receiver(post_save, sender=Sample)
def sample_saved(sender, instance, created, **kwargs):
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Sample created:", instance.sample_id)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        print("Sample saved:", instance.sample_id)


@receiver(pre_delete, sender=Sample)
def sample_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    print("Sample deleted:", instance.sample_id)
