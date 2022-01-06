import os
from django.core.management.base import BaseCommand, CommandError
from browser.models import *
from browser.dataimport.annotator import Annotator

# Register plugins here
from browser.dataimport.plugins.fama import application as fama
from browser.dataimport.plugins.amrfinder import application as amrfinder
from browser.dataimport.plugins.antismash import application as antismash
from browser.dataimport.plugins.phispy import application as phispy
from browser.dataimport.plugins.ecis_screen import application as ecis_screen


class Command(BaseCommand):
    help = 'For protein sequences from selected genomes uploaded to Django database, this program runs hmmsearch with PFAM and TIGRFAM HMMs, and writes domain annotations to database. Existing domain annotations are deleted.'

    def add_arguments(self, parser):
        parser.add_argument('-i', default='genomes.txt', help='Path to input file')

    def handle(self, *args, **options):
        if os.path.exists(options['i']):
            annotator = Annotator()
            genomes = {}
            with open(options['i'], 'r') as infile:
                for line in infile:
                    if line.startswith('#'):
                        continue
                    filepath, genome_name, _, _, _, _ = line.rstrip('\n\r').split('\t')
                    genomes[genome_name] = filepath
            
            plugins_available = set()
            for param in annotator.config:
                if param.startswith('plugins.') and annotator.config[param] != '':
                    plugin = param.split('.')[1]
                    plugins_available.add(plugin)
                    
            print('Available plugins:', ','.join(list(plugins_available)))

            # Run Fama
            if 'fama' in plugins_available:
                fama_annotations = fama(annotator, genomes)
                print(fama_annotations, 'ready for upload')
                annotator.add_custom_annotations(fama_annotations)

            # Run amrfinder
            if 'amrfinder' in plugins_available:
                amrfinder_annotations = amrfinder(annotator, genomes)
                print(amrfinder_annotations, 'ready for upload')
                annotator.add_custom_annotations(amrfinder_annotations)

            # Run antiSMASH
            if 'antismash' in plugins_available:
                antismash_annotations = antismash(annotator, genomes)
                print(antismash_annotations, 'ready for upload')
                annotator.add_custom_annotations(antismash_annotations)
            
            # Run PhiSpy
            if 'phispy' in plugins_available:
                phispy_annotations = phispy(annotator, genomes)
                print(phispy_annotations, 'ready for upload')
                annotator.add_custom_annotations(phispy_annotations)

            # Run eCIS-screen
            if 'ecis_screen' in plugins_available:
                ecis_screen_annotations = ecis_screen(annotator, genomes)
                print(ecis_screen_annotations, 'ready for upload')
                annotator.add_custom_annotations(ecis_screen_annotations)
            
        else:
            raise CommandError('Genomes file ' + options['i'] + ' not found.')
