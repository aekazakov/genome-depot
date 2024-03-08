import logging
from django.db.models.signals import post_save, pre_delete
from django.dispatch import receiver
from django.core import serializers
from django.utils import timezone
from browser.models import ChangeLog
from browser.models import Genome, Strain, Sample

logger = logging.getLogger("GenomeDepot")

# this receiver is executed every-time some data is saved in Genome table
@receiver(post_save, sender=Genome)
def genome_saved(sender, instance, created, **kwargs):
    # code to execute before every model save
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Genome created: %s", instance.name)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Genome saved: %s", instance.name)


@receiver(pre_delete, sender=Genome)
def genome_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    logger.info("Genome deleted: %s", instance.name)


# Strain created, updated or deleted
@receiver(post_save, sender=Strain)
def srain_saved(sender, instance, created, **kwargs):
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Strain created: %s", instance.strain_id)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Strain saved: %s", instance.strain_id)


@receiver(pre_delete, sender=Strain)
def strain_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    logger.info("Strain deleted: %s", instance.strain_id)


# Sample created, updated or deleted
@receiver(post_save, sender=Sample)
def sample_saved(sender, instance, created, **kwargs):
    if created:
        ChangeLog.objects.create(action='created',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Sample created: %s", instance.sample_id)
    else:
        ChangeLog.objects.create(action='saved',
                                 data=serializers.serialize('json', [instance]),
                                 timestamp=timezone.now()
                                 )
        logger.info("Sample saved: %s", instance.sample_id)


@receiver(pre_delete, sender=Sample)
def sample_deleted(sender, instance, **kwargs):
    ChangeLog.objects.create(action='deleted',
                             data=serializers.serialize('json', [instance]),
                             timestamp=timezone.now()
                             )
    logger.info("Sample deleted: %s", instance.sample_id)
