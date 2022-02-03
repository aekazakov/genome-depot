import os
import configparser
import gzip
import shutil
#import mmh3
import hashlib
from collections import defaultdict, OrderedDict
from pathlib import Path
from Bio import GenBank
from subprocess import Popen, PIPE, CalledProcessError

from django.db import connection
from django.utils import timezone

from browser.models import *
from genomebrowser.settings import STATIC_URL
from browser.dataimport.annotator import Annotator
from browser.dataimport.taxonomy import load_taxonomy

""" Imports genomes into database from GenBank files.
Input file must have five columns:
1. Path to GenBank file (gzipped or not)
2. Genome ID.
3. Strain ID. Can be blank if sample ID is not blank
4. Sample ID. Can be blank if strain ID is not blank
5. URL for the sequence (can be blank)
6. External ID to be used as a label (can be blank)

"""

class Importer(object):

    def __init__(self):
        #print('Config path:', config_path)
        #self.config_path = config_path
        self.config = {}
        self.read_config() # configparser.ConfigParser()
        #self.config.read(config_path)
        #print(self.config['Taxonomy']['names_file'])
        self.inputgenomes = defaultdict(dict)
        self.staticfiles = defaultdict(list)
        self.taxonomy, self.taxonomy_id_lookup, self.eggnog_taxonomy_lookup = load_taxonomy(self.config['ref.taxonomy'], self.config['cgcms.eggnog_taxonomy'])
        #self.taxonomy_id_lookup = {}
        self.protein_hash2id = {}
        self.protein_hash2gene = defaultdict(list)
#        self.load_taxonomy()
#        if not os.path.exists(self.config['genomes.eggnog_output_stored']):
#            with open(self.config['genomes.eggnog_output_stored'], 'w') as outfile:
#                outfile.write('')
        self.eggnog_input_file = os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_input.faa')
        if os.path.exists(self.eggnog_input_file):
            os.remove(self.eggnog_input_file)
        # Data tables
        self.strain_instances = {}
        self.strain_metadata = {}
        self.sample_instances = {}
        self.sample_metadata = {}
        self.taxon_instances = {}
        self.genome_data = {}
        self.contig_data = {}
        self.protein_data = {}
        self.gene_data = {}
        self.mappings = {}
        self.relations = {}
        
    def read_config(self):
        for item in Config.objects.values('param', 'value'):
            self.config[item['param']] = item['value']
    
    def make_dirs(self):
        """
            Try to create directories for static files. 
        """
        Path(self.config['cgcms.temp_dir']).mkdir(parents=True, exist_ok=True)
        Path(self.config['cgcms.static_dir']).mkdir(parents=True, exist_ok=True)
        Path(self.config['cgcms.static_dir'] + '/gbff').mkdir(parents=True, exist_ok=True)
        Path(self.config['cgcms.json_dir']).mkdir(parents=True, exist_ok=True)
        Path(self.config['cgcms.search_db_dir']).mkdir(parents=True, exist_ok=True)
        
    def load_proteins(self):
        """
            Populates protein_hash2id dictionary
        """
        for item in Protein.objects.values('id', 'protein_hash'):
            self.protein_hash2id[item['protein_hash']] = item['id']
        
    def load_genome_list(self, in_file):
        """
            Reads file with genome paths, names and links.
        """
        with open(in_file, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                filepath, genome_id, strain_id, sample_id, url, external_id = line.rstrip('\n\r').split('\t')
                # Sanity check
                if sample_id == '' and strain_id == '':
                    raise ValueError('Either strain ID or sample ID required for genome ' + genome_id)

                self.inputgenomes[genome_id]['strain'] = strain_id
                self.inputgenomes[genome_id]['sample'] = sample_id
                self.inputgenomes[genome_id]['gbk'] = filepath
                self.inputgenomes[genome_id]['url'] = url
                self.inputgenomes[genome_id]['external_id'] = external_id
    
    def check_genomes(self):
        result = {}
        genomes_exist = []
        saved_genomes = set([item['name'] for item in Genome.objects.values('name')])
        for genome_id in self.inputgenomes:
            if genome_id in saved_genomes:
                genomes_exist.append(genome_id)
            else:
                result[genome_id] = self.inputgenomes[genome_id]
        if genomes_exist:
            print('CHECK LIST OF GENOMES. SOME GENOME NAMES ALREADY EXIST IN THE DATABASE:')
            for genome in genomes_exist:
                print(genome)
            raise ValueError('Genome identifiers must be unique. Change genome ID or delete conflicting genome from DB.')
        self.genomefiles = result
        
    def download_genomes(self):
        """
            This function is a placeholder for real downloader.
            Currently, it only checks if a genome file already exists.
        """
        for genome_id in self.inputgenomes:
            genome_file = self.inputgenomes[genome_id]['gbk']
            if not os.path.exists(genome_file):
                print('File not found: ' + genome_file)
                print('Please download genome from ' + self.inputgenomes[genome_id]['url'])
                raise FileNotFoundError('File not found')

    def get_taxonomic_order(self, taxonomy_id):
        result = 'Unknown'
        if taxonomy_id in self.taxonomy:
            if self.taxonomy[taxonomy_id]['rank'] == 'order':
                result = self.taxonomy[taxonomy_id]['name']
            else:
                parent_id = self.taxonomy[taxonomy_id]['parent']
                while parent_id != '1':
                    rank = self.taxonomy[parent_id]['rank']
                    if rank == 'order':
                        result = self.taxonomy[parent_id]['name']
                        break
                    parent_id = self.taxonomy[parent_id]['parent']
        return result

    def get_gbk_organism(self, gbk_file):
        organism_data = {}
        if gbk_file.endswith('.gz'):
            gbk_handle = gzip.open(gbk_file, 'rt')
        else:
            gbk_handle = open(gbk_file, 'r')
        parser = GenBank.parse(gbk_handle)
        organism = ''
        for gbk_record in parser:
            organism = str(gbk_record.organism)
            for feature in gbk_record.features:
                if feature.key == 'source':
                    for qualifier in feature.qualifiers:
                        if qualifier.key == '/organism=':
                            if 'organism' not in organism_data:
                                organism_data['organism'] = qualifier.value[1:-1]
#                            elif organism_data['organism'] != qualifier.value[1:-1]:
#                                raise ValueError('Inconsistent organism names ' + ' ' + organism_data['organism'] + ' ' + qualifier.value[1:-1])
                            else:
                                if 'tax_id' in organism_data:
                                    break  # Only the first occurrence of the 'source' feature is considered
                        elif qualifier.key == '/db_xref=':
                            qualifier_value = qualifier.value[1:-1]
                            if qualifier_value.startswith('taxon:'):
                                taxonomy_id = qualifier_value.split(':')[1]
                                if taxonomy_id == '6666666':  # Fix for fake taxonomy ID of SEED genomes
                                    taxonomy_id = '0'  # 
                                if taxonomy_id not in self.taxonomy:
                                    raise KeyError(gbk_file + ' parsing error. ' + taxonomy_id + ' not found in taxonomy reference file. Stop now and update taxonomy reference data.')
                                organism_data['tax_id'] = taxonomy_id
                    break  # Only the first occurrence of the 'source' feature is considered
            if 'organism' in organism_data and 'tax_id' in organism_data:
                break  # Only the first sequence record with the 'source' feature is considered
        gbk_handle.close()
        # Fill in missing values
        if 'organism' not in organism_data:
            if organism != '':
                organism_data['organism'] = organism  # if 'source' feature is missing in the entire file, let's try organism attribute of sequence records
            else:
                print('Organism name not found in', gbk_file)
                organism_data['organism'] = 'Unknown'
        if 'tax_id' not in organism_data:
            print('Taxonomy ID not found in', gbk_file)
            organism_data['tax_id'] = '0'
        return organism_data
    
    def generate_strain_data(self, gbk_file, strain_id, order):
        """
            Generates strain data, if the strain was not found in ENIGMA strain collection.
            Finds organism name and taxonomy identifier in the gbk_file, and creates new Taxon object.
            Creates and returns Strain object with strain_id identifier
        """
        # returns Strain object
        organism_data = self.get_gbk_organism(gbk_file)

        if 'organism' not in organism_data:
            if strain_id == '':
                organism_data['organism'] = 'Unknown organism'
            else:
                organism_data['organism'] = 'Environmental isolate ' + strain_id

        if 'tax_id' not in organism_data:
            organism_data['tax_id'] = '0'

        if order == '':
            if organism_data['tax_id'] == '0':
                organism_data['order'] = 'Unknown'
            else:
                organism_data['order'] = self.get_taxonomic_order(organism_data['tax_id'])
        else:
            organism_data['order'] = order
            
        if 'strain_id' not in organism_data:
            organism_data['strain_id'] = strain_id
        try:
            taxon = Taxon.objects.get(taxonomy_id = organism_data['tax_id'])
        except Taxon.DoesNotExist:
            taxon = Taxon(taxonomy_id=organism_data['tax_id'],
                eggnog_taxid=self.taxonomy[organism_data['tax_id']]['eggnog_taxid'],
                name=self.taxonomy[organism_data['tax_id']]['name'],
                rank=self.taxonomy[organism_data['tax_id']]['rank'],
                parent_id=self.taxonomy[organism_data['tax_id']]['parent']
                )
            taxon.save()
        strain = Strain(strain_id=strain_id,
                        full_name=organism_data['organism'],
                        order=organism_data['order'],
                        taxon=taxon)
        return strain

    def prepare_strain_data(self):
        """
            This function creates strain entries and strain metadata for new genomes.
        """
        saved_metadata = defaultdict(dict)
        with open(self.config['strains.metadata_file'], 'r') as infile:
            # Strain metadata fields: strain, source, url, key, value
            for line in infile:
                if line.startswith('#'):
                    continue
                row = line.rstrip('\n\r').split('\t')
                saved_metadata[row[0]][row[3]] = [row[1], row[2], row[4]]
            
        saved_strains = {strain_data['strain_id']:strain_data['taxon__taxonomy_id'] for strain_data in Strain.objects.values('strain_id', 'taxon__taxonomy_id')}
        for genome_id in self.inputgenomes:
            strain_id = self.inputgenomes[genome_id]['strain']
            if strain_id == '':
                continue
            if strain_id in saved_strains:
                if saved_strains[strain_id] == '0':
                    # try to update strain taxonomy from GBK file
                    organism_data = self.get_gbk_organism(self.inputgenomes[genome_id]['gbk'])
                    if 'tax_id' in organism_data and organism_data['tax_id'] != '0':
                        try:
                            taxon = Taxon.objects.get(taxonomy_id = organism_data['tax_id'])
                        except Taxon.DoesNotExist:
                            taxon = Taxon(taxonomy_id=organism_data['tax_id'],
                                eggnog_taxid=self.taxonomy[organism_data['tax_id']]['eggnog_taxid'],
                                name=self.taxonomy[organism_data['tax_id']]['name'],
                                rank=self.taxonomy[organism_data['tax_id']]['rank'],
                                parent_id=self.taxonomy[organism_data['tax_id']]['parent']
                                )
                            taxon.save()
                        strain = Strain.objects.get(strain_id=strain_id)
                        strain.taxon = taxon
                        strain.full_name = organism_data['organism']
                        strain.order = self.get_taxonomic_order(organism_data['tax_id'])
                        strain.save()
                        self.inputgenomes[genome_id]['taxon'] = taxon
                    else:
                        self.inputgenomes[genome_id]['taxon'] = saved_strains[strain_id]
                else:
                    self.inputgenomes[genome_id]['taxon'] = saved_strains[strain_id]
            else:
                if strain_id not in self.strain_instances:
                    order = ''
                    if 'Phylogenetic order' in saved_metadata[strain_id]:
                        order = saved_metadata[strain_id]['Phylogenetic order']
                    self.strain_instances[strain_id] = self.generate_strain_data(self.inputgenomes[genome_id]['gbk'], strain_id, order)
                    if strain_id in saved_metadata:
                        # add strain metadata
                        self.strain_metadata[strain_id] = []
                        for key in saved_metadata[strain_id]:
                            self.strain_metadata[strain_id].append({'source':saved_metadata[strain_id][key][0],
                                                                'url':saved_metadata[strain_id][key][1],
                                                                'key':key,
                                                                'value':saved_metadata[strain_id][key][2]
                                                                })
                self.inputgenomes[genome_id]['taxon'] = self.strain_instances[strain_id].taxon.taxonomy_id

    def prepare_sample_data(self):
        """
            This function creates sample entries in the database for new genomes.
        """
        saved_samples = set([item['sample_id'] for item in Sample.objects.values('sample_id')])
        for genome_id in self.inputgenomes:
            sample_id = self.inputgenomes[genome_id]['sample']
            if sample_id == '':
                continue
            if sample_id in saved_samples:
                continue
            else:
                sample = Sample(sample_id=sample_id,
                                full_name = '',
                                description = '')
                self.sample_instances[sample_id] = sample
                
    def prepare_taxonomy_data(self):
        taxonomy_ids = set()
        for strain_id, strain in self.strain_instances.items():
            taxon_id = strain.taxon.taxonomy_id
            taxonomy_ids.add(taxon_id)
            parent_id = self.taxonomy[taxon_id]['parent']
            while parent_id != '1':
                taxonomy_ids.add(parent_id)
                parent_id = self.taxonomy[parent_id]['parent']
        saved_taxa = set([taxon.taxonomy_id for taxon in Taxon.objects.all()])
        for taxon_id in taxonomy_ids:
            if taxon_id not in saved_taxa:
                taxon = Taxon(taxonomy_id=taxon_id,
                    eggnog_taxid=self.taxonomy[taxon_id]['eggnog_taxid'],
                    name=self.taxonomy[taxon_id]['name'],
                    rank=self.taxonomy[taxon_id]['rank'],
                    parent_id=self.taxonomy[taxon_id]['parent']
                    )
                self.taxon_instances[taxon_id] = taxon
    
    def process_protein(self, sequence, protein_name):
        #protein_hash = hex(mmh3.hash128(sequence))[2:]
        if sequence.endswith('*'):
            sequence = sequence[:-1]
        protein_hash = hashlib.md5(sequence.encode('utf-8')).hexdigest()
        if protein_hash not in self.protein_data:
            self.protein_data[protein_hash] = {}
            self.protein_data[protein_hash]['name'] = protein_name
            self.protein_data[protein_hash]['length'] = len(sequence)
            self.protein_data[protein_hash]['sequence'] = sequence
        if protein_hash not in self.protein_hash2id:
            self.protein_hash2id[protein_hash] = protein_name
            with open(self.eggnog_input_file, 'a') as outfile:
                outfile.write('>' + protein_hash + '\n' + sequence + '\n')
        return protein_hash
        
    def parse_location(self, location):
        strand = 1
        if location.startswith('complement('):
            strand = -1
            location = location[11:-1]
        if location.startswith('join('):
            location = location[5:-1]
        location = location.replace('>','')
        location = location.replace('<','')
        location_coords = location.split('..')
        try:
            start = int(location_coords[0])
            end = int(location_coords[-1])
        except ValueError:
            print('Location', location)
            raise
        if start > end:
            # Feature that starts at the end of circular genome and ends at the beginning
            end = start
            start = 1
        return start, end, strand
        
    def process_feature(self, feature, genome_id, contig_id, locus_tags):
        feature_data = {}
        result = ''
        translation = ''
        accepted_features = ['gene', 'CDS', 'rRNA', 'tRNA']
        if feature.key not in accepted_features:
            return '', ''
        start, end, strand = self.parse_location(feature.location)
        locus_tag = ''
        old_locus_tag = ''
        translation = ''
        pseudogene = False
        function = ''
        name = ''
        for qualifier in feature.qualifiers:
            if qualifier.key == '/locus_tag=':
                locus_tag = qualifier.value[1:-1]
            elif qualifier.key == '/old_locus_tag=':
                old_locus_tag = qualifier.value[1:-1]
            elif qualifier.key == '/translation=':
                translation = qualifier.value[1:-1]
            elif qualifier.key == '/pseudo':
                pseudogene = True
            elif qualifier.key == '/pseudogene=':
                pseudogene = True
            elif qualifier.key == '/product=':
                function = qualifier.value[1:-1]
            elif qualifier.key == '/gene=':
                name = qualifier.value[1:-1]
            elif qualifier.key == '/gene_synonym=' and name == '':
                name = qualifier.value[1:-1]

        if locus_tag == '' and old_locus_tag != '':
            locus_tag = old_locus_tag # Just in case
        feature_uid = ' '.join([genome_id, contig_id, str(start), str(end), str(strand)])
        if feature_uid in self.gene_data:
            if feature.key != 'gene' and self.gene_data[feature_uid]['type'] != 'pseudogene':
                if pseudogene:
                    self.gene_data[feature_uid]['type'] = 'pseudogene'
                else:
                    self.gene_data[feature_uid]['type'] = feature.key
            if translation != '' and 'protein_hash' not in self.gene_data[feature_uid]:
                protein_name = locus_tag + '|' + str(len(translation)) + 'aa|' + function
                if len(protein_name) > 100:
                    protein_name = protein_name[:97] + '...'
                protein_hash = self.process_protein(translation, protein_name)
                self.gene_data[feature_uid]['protein_hash'] = protein_hash
            if function != '' and 'function' not in self.gene_data[feature_uid]:
                self.gene_data[feature_uid]['function'] = function
            if self.gene_data[feature_uid]['name'] == '':
                self.gene_data[feature_uid]['name'] = name
        else:
            if locus_tag in locus_tags:
                return '','' # duplicated locus tags for different locations are not allowed
            self.gene_data[feature_uid] = {}
            self.gene_data[feature_uid]['name'] = name
            self.gene_data[feature_uid]['locus_tag'] = locus_tag
            self.gene_data[feature_uid]['contig_id'] = contig_id
            self.gene_data[feature_uid]['start'] = start
            self.gene_data[feature_uid]['end'] = end
            self.gene_data[feature_uid]['strand'] = strand
            # TODO: remove
            # self.gene_data[feature_uid]['taxon'] = self.inputgenomes[genome_id]['taxon']
            self.gene_data[feature_uid]['genome_id'] = genome_id
            if pseudogene:
                self.gene_data[feature_uid]['type'] = 'pseudogene'
            else:
                self.gene_data[feature_uid]['type'] = feature.key
            if translation != '':
                protein_name = locus_tag + '|' + str(len(translation)) + 'aa|' + function
                if len(protein_name) > 100:
                    protein_name = protein_name[:97] + '...'
                protein_hash = self.process_protein(translation, protein_name)
                self.gene_data[feature_uid]['protein_hash'] = protein_hash
                self.protein_hash2gene[protein_hash].append(feature_uid)
            if function != '':
                self.gene_data[feature_uid]['function'] = function
        return feature_uid, locus_tag
            
    def process_gbk(self):
        strored_seq_uids = self.export_contigs()
        nucl_db_file = os.path.join(self.config['cgcms.temp_dir'], os.path.basename(self.config['cgcms.search_db_nucl'])) # 
        self.staticfiles[self.config['cgcms.search_db_dir']].append(os.path.basename(nucl_db_file))
        gbff_dir = os.path.join(self.config['cgcms.static_dir'], 'gbff')
        nucl_db_file_handle = open(nucl_db_file, 'a')
        for genome_id in self.inputgenomes:
            gbk_file = self.inputgenomes[genome_id]['gbk']
            gbk_copy = os.path.join(gbff_dir, genome_id + '.genome.gbff.gz')
            organism_data = self.get_gbk_organism(gbk_file)
            try:
                taxon = Taxon.objects.get(taxonomy_id = organism_data['tax_id'])
            except Taxon.DoesNotExist:
                taxon = Taxon(taxonomy_id=organism_data['tax_id'],
                    eggnog_taxid=self.taxonomy[organism_data['tax_id']]['eggnog_taxid'],
                    name=self.taxonomy[organism_data['tax_id']]['name'],
                    rank=self.taxonomy[organism_data['tax_id']]['rank'],
                    parent_id=self.taxonomy[organism_data['tax_id']]['parent']
                    )
                taxon.save()

            if gbk_file.endswith('.gz'):
                shutil.copy(gbk_file, gbk_copy)
                gbk_handle = gzip.open(gbk_file, 'rt')
            else:
                with open(gbk_file, 'rb') as infile:
                    with gzip.open(gbk_copy, 'wb') as outfile:
                        shutil.copyfileobj(infile, outfile)
                gbk_handle = open(gbk_file, 'r')
            self.inputgenomes[genome_id]['gbk_filepath'] = gbk_copy
            parser = GenBank.parse(gbk_handle)
            genome_size = 0
            contig_count = 0
            contig_sizes = []
            genome_sequence = []
            features = OrderedDict()
            genome_fasta = os.path.join(self.config['cgcms.temp_dir'], genome_id + '.contigs.fasta')
            locus_tags = set()
            with open(genome_fasta, 'w') as outfile:
                for gbk_record in parser:
                    contig_sequence = gbk_record.sequence
                    genome_sequence.append(contig_sequence)
                    contig_id = gbk_record.locus
                    contig_name = gbk_record.accession
                    nucl_seq_uid = '>' + contig_id + '|' + genome_id
                    if nucl_seq_uid not in strored_seq_uids:
                        nucl_db_file_handle.write('>' + contig_id + '|' + genome_id + '\n' + ''.join(contig_sequence) + '\n')
                    outfile.write('>' + contig_id + '\n' + ''.join(contig_sequence) + '\n')
                    contig_size = int(gbk_record.size)
                    contig_sizes.append(contig_size)
                    self.contig_data[(genome_id, contig_id)] = {'contig_id':contig_id,
                        'name':gbk_record.accession,
                        'size':contig_size,
                        'genome':genome_id
                        }
                    for feature in gbk_record.features:
                        feature_location, locus_tag = self.process_feature(feature, genome_id, contig_id, locus_tags)
                        locus_tags.add(locus_tag)
                        if feature_location != '':
                            features[feature_location] = ''
                    genome_size += contig_size
                    contig_count += 1
            self.genome_data[genome_id] = {'name':genome_id,
                'description':genome_id,
                'contigs':contig_count,
                'size':genome_size,
                'genes':len(features),
                'json_url':'../genomes/json/' + genome_id,
                'pub_date': timezone.now(),
                'external_url':self.inputgenomes[genome_id]['url'],
                'external_id':self.inputgenomes[genome_id]['external_id'],
                'gbk_filepath':gbk_copy,
                'taxon':taxon
                }
            if self.inputgenomes[genome_id]['strain'] == '':
                self.genome_data[genome_id]['strain'] = None
            else:
                self.genome_data[genome_id]['strain'] = self.inputgenomes[genome_id]['strain']
            if self.inputgenomes[genome_id]['sample'] == '':
                self.genome_data[genome_id]['sample'] = None
            else:
                self.genome_data[genome_id]['sample'] = self.inputgenomes[genome_id]['sample']
            gbk_handle.close()
        nucl_db_file_handle.close()
    
    def run_eggnog_mapper(self):
        """
        Runs eggnog-mapper for eggnog_mapper_input.faa file in temp directory.
        """
        chunk_size = '200000'
        result = os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_output.emapper.annotations')
        work_dir = self.config['cgcms.eggnog_outdir']
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        orthologs_file = os.path.join(work_dir, 'input_file.emapper.seed_orthologs')
        if os.path.exists(result):
            os.remove(result)
        if os.path.exists(orthologs_file):
            os.remove(orthologs_file)
        eggnog_mapper_script = os.path.join(self.config['cgcms.temp_dir'], 'run_emapper.sh')
        with open(eggnog_mapper_script, 'w') as outfile:
            outfile.write('cd ' + work_dir + '\n')
            outfile.write('split -l ' + chunk_size + ' -a 3 -d ' +
                          os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_input.faa') +
                          ' input_file.chunk_ \n')
            outfile.write('for f in input_file.chunk_*; do\n')
            outfile.write(self.config['cgcms.eggnog_command'] + ' -m diamond --no_annot --no_file_comments --cpu 10 -i $f -o $f;\n')
            outfile.write('done\n')
            outfile.write('cat input_file.chunk_*.emapper.seed_orthologs >> ' + orthologs_file + '\n')
            outfile.write(self.config['cgcms.eggnog_command'] + 
                          ' --annotate_hits_table ' +
                          orthologs_file +
                          ' --no_file_comments -o ' +
                          os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_output') +
                          ' --cpu 10\n')
            outfile.write('\n')
        '''
        # split input file into chunks
        split -l 200000 -a 3 -d eggnog_mapper_input.faa input_file.chunk_
        # generate orhtologs for each chunk
        for f in input_file.chunk_*; do
        ./emapper.py -m diamond --no_annot --no_file_comments --cpu 10 -i $f -o $f; 
        done
        # join orthologs files
        cat input_file.chunk_*.emapper.seed_orthologs >> input_file.emapper.seed_orthologs
        # generate annotations
        ./emapper.py --annotate_hits_table input.emapper.seed_orthologs --no_file_comments -o eggnog_mapper_output --cpu 10        
        '''
        '''
        Old command calling emapper.py
        cmd = [self.config['cgcms.eggnog_command'],
              '-m', 'diamond', '--no_file_comments', '--cpu', '12',
              '-i', os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_input.faa'),
              '-o', os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_output')
              ]
        '''
        cmd = ['/bin/bash', eggnog_mapper_script]
        with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
            for line in proc.stdout:
                print(line, end='')
        if proc.returncode != 0:
            # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
            # pylint: disable=no-member
            raise CalledProcessError(proc.returncode, proc.args)
        return result
        
    def parse_eggnog_output(self, eggnog_outfile):
        """Reads eggnog-mapper output. Returns new eggnog annotations to be added into permanent storage"""
        result = []
        existing_eggnog_annotations = set(Protein.objects.exclude(ortholog_groups=None).values_list('protein_hash', flat=True))
        
#        with open(self.config['genomes.eggnog_output_stored'], 'r') as infile:
#            for line in infile:
#                row = line.rstrip('\n\r').split('\t')
#                existing_eggnog_annotations.add(row[0])
            
        data_fields = ['query_name', 'eggNOG_ortholog', 'evalue', 'score', 'taxonomic_group', 
            'name', 'GO_terms', 'EC_number', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
            'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'tax_scope', 'eggNOG_OG', 'bestOG',
            'COG_Functional_Category', 'description']
        list_data_fields = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        with open(eggnog_outfile, 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                protein_hash = row[0]
                if protein_hash not in self.protein_data:
                    continue # Skip proteins, which are not from the genomes of interest
                if protein_hash in existing_eggnog_annotations:
                    continue # Skip proteins that already have eggnog mappings
                result.append(line)
                for row_index, row_value in enumerate(row):
                    if row_index < 2:
                        continue
                    if row_value.startswith('ko:'):
                        row_value = row_value[3:]
                    self.protein_data[protein_hash]['eggnog.' + data_fields[row_index]] = row_value
        return result
        
    def prepare_mappings_data(self, ref_file, mapping_id):
        """Loads file of KEGG reactions, collects IDs from self.gene_data,
        makes list of cross-mappings
        """
        ref_data = defaultdict(dict)
        selected_ref_data = {}
        relations = set()
        with open(ref_file, 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                if ':' in row[0]:
                    row[0] = row[0].split(':')[-1]
                ref_data[row[0]]['id'] = row[1]
                ref_data[row[0]]['description'] = row[2]
        print(len(ref_data), 'reference entries found')
        for protein_hash in self.protein_data:
            if mapping_id not in self.protein_data[protein_hash]:
                continue
            if mapping_id == 'eggnog.COG_Functional_Category':
                ref_items = self.protein_data[protein_hash][mapping_id]
            else:
                ref_items = self.protein_data[protein_hash][mapping_id].split(',')
            for ref_item in ref_items:
                if ref_item == '':
                    continue
                if mapping_id == 'eggnog.KEGG_ko' and ref_item.startswith('ko:'):
                    ref_item = ref_item[3:]
                try:
                    ref_id = ref_data[ref_item]['id']
                    ref_description = ref_data[ref_item]['description']
                except KeyError:
                    if not ref_item.startswith('ko') and mapping_id == 'eggnog.KEGG_Reaction':
                        print(ref_item, 'not found in', ref_file)
                    continue
                selected_ref_data[ref_id] = ref_description
                relations.add((protein_hash, ref_id))
        return selected_ref_data, relations

    def prepare_go_mappings_data(self, ref_file, mapping_id):
        """Loads file of GO terms, collects IDs from self.gene_data,
        makes list of cross-mappings
        """
        ref_data = defaultdict(dict)
        selected_ref_data = {}
        relations = set()
        with open(ref_file, 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                ref_data[row[0]]['id'] = row[1]
                ref_data[row[0]]['description'] = (row[2], row[3])
        print(len(ref_data), 'reference entries found')
        for protein_hash in self.protein_data:
            if mapping_id not in self.protein_data[protein_hash]:
                continue
            ref_items = self.protein_data[protein_hash][mapping_id].split(',')
            for ref_item in ref_items:
                if ref_item == '':
                    continue
                if ref_item in ref_data:
                    ref_id = ref_data[ref_item]['id']
                    ref_description = ref_data[ref_item]['description']
                    selected_ref_data[ref_id] = ref_description
                    relations.add((protein_hash, ref_id))
        return selected_ref_data, relations

    
    def make_mappings(self):
        mappings, relations = self.prepare_mappings_data(self.config['ref.kegg_reactions_file'],
            'eggnog.KEGG_Reaction'
            )
        self.mappings['kegg_reaction'] = mappings
        self.relations['kegg_reaction'] = relations

        mappings, relations = self.prepare_mappings_data(self.config['ref.kegg_orthologs_file'],
            'eggnog.KEGG_ko'
            )
        self.mappings['kegg_ortholog'] = mappings
        self.relations['kegg_ortholog'] = relations

        mappings, relations = self.prepare_mappings_data(self.config['ref.kegg_pathways_file'],
            'eggnog.KEGG_Pathway'
            )
        self.mappings['kegg_pathway'] = mappings
        self.relations['kegg_pathway'] = relations
    
        mappings, relations = self.prepare_mappings_data(self.config['ref.ec_file'],
            'eggnog.EC_number'
            )
        self.mappings['ec'] = mappings
        self.relations['ec'] = relations
    
        mappings, relations = self.prepare_mappings_data(self.config['ref.tc_file'],
            'eggnog.KEGG_TC'
            )
        self.mappings['tc'] = mappings
        self.relations['tc'] = relations

        mappings, relations = self.prepare_mappings_data(self.config['ref.cazy_file'],
            'eggnog.CAZy'
            )
        self.mappings['cazy'] = mappings
        self.relations['cazy'] = relations

        mappings, relations = self.prepare_mappings_data(self.config['ref.cog_codes_file'],
            'eggnog.COG_Functional_Category'
            )
        self.mappings['cog'] = mappings
        self.relations['cog'] = relations

        mappings, relations = self.prepare_go_mappings_data(self.config['ref.go_file'],
            'eggnog.GO_terms'
            )
        self.mappings['go'] = mappings
        self.relations['go'] = relations
    
    def prepare_og_data(self):
        self.mappings['og'] = {}
        self.relations['og'] = set()
        saved_taxonomy_ids = set([taxon.taxonomy_id for taxon in Taxon.objects.all()])
        for protein_hash in self.protein_data:
            if 'eggnog.eggNOG_OG' not in self.protein_data[protein_hash]:
                continue
            orthogroup_ids = self.protein_data[protein_hash]['eggnog.eggNOG_OG'].split(',')
            for orthogroup_id in orthogroup_ids:
                if orthogroup_id == '':
                    continue
                eggnog_id, eggnog_taxonomy_id = orthogroup_id.split('@')
                taxonomy_id = self.eggnog_taxonomy_lookup[eggnog_taxonomy_id]
                if taxonomy_id not in saved_taxonomy_ids and taxonomy_id not in self.taxon_instances:
                    self.taxon_instances[taxonomy_id] = Taxon(taxonomy_id=taxonomy_id,
                        eggnog_taxid = eggnog_taxonomy_id,
                        name=self.taxonomy[taxonomy_id]['name'],
                        rank=self.taxonomy[taxonomy_id]['rank'],
                        parent_id=self.taxonomy[taxonomy_id]['parent']
                        )
                self.mappings['og'][eggnog_id + '@' + taxonomy_id] = (eggnog_id, taxonomy_id)
                self.relations['og'].add((protein_hash, eggnog_id, taxonomy_id))
        for protein_hash in self.protein_data:
            if 'eggnog.taxonomic_group' not in self.protein_data[protein_hash]:
                continue
            taxon_name = self.protein_data[protein_hash]['eggnog.taxonomic_group']
            taxonomy_id = self.taxonomy_id_lookup[taxon_name]
            self.protein_data[protein_hash]['taxonomic_group'] = taxonomy_id
            if taxonomy_id not in saved_taxonomy_ids and taxonomy_id not in self.taxon_instances:
                self.taxon_instances[taxonomy_id] = Taxon(taxonomy_id=taxonomy_id,
                    eggnog_taxid = self.taxonomy[taxonomy_id]['eggnog_taxid'],
                    name=self.taxonomy[taxonomy_id]['name'],
                    rank=self.taxonomy[taxonomy_id]['rank'],
                    parent_id=self.taxonomy[taxonomy_id]['parent']
                    )
        
    def prepare_eggnog_description_data(self):
        self.mappings['description'] = {}
        self.relations['description'] = set()
        saved_descriptions = set([description.fingerprint for description in Eggnog_description.objects.all()])
        for protein_hash in self.protein_data:
            if 'eggnog.description' not in self.protein_data[protein_hash]:
                continue
            description = self.protein_data[protein_hash]['eggnog.description']
            if description == '':
                continue
            fingerprint = hashlib.md5(description.encode('utf-8')).hexdigest()
            if fingerprint not in saved_descriptions and fingerprint not in self.mappings['description']:
                self.mappings['description'][fingerprint] = description
                self.protein_data[protein_hash]['eggnog_description'] = fingerprint


    def write_data(self):
        # write taxonomy data
        Taxon.objects.bulk_create(self.taxon_instances.values(), batch_size=1000)
        print('Taxonomy imported')
        taxon_instances = {taxon.taxonomy_id:taxon for taxon in Taxon.objects.all()}
        
        # write strain data
        Strain.objects.bulk_create(self.strain_instances.values(), batch_size=1000)
        print('Strains imported')

        # write strain metadata
        strain_instances = {strain.strain_id:strain for strain in Strain.objects.all()}
        strain_metadata_instances = []
        for strain_id, metadata_entries in self.strain_metadata.items():
            for metadata_entry in metadata_entries:
                strain_metadata_instances.append(Strain_metadata(
                    strain=strain_instances[strain_id],
                    source=metadata_entry['source'],
                    url=metadata_entry['url'],
                    key=metadata_entry['key'],
                    value=metadata_entry['value']
                    ))
        Strain_metadata.objects.bulk_create(strain_metadata_instances, batch_size=1000)
        print('Strain metadata imported')

        # write sample data
        Sample.objects.bulk_create(self.sample_instances.values(), batch_size=1000)
        print('Samples imported')
        sample_instances = {sample.sample_id:sample for sample in Sample.objects.all()}

        # TODO: create sample metadata
        
        # write genome data
        genome_instances = []
        for genome_id,genome_data in self.genome_data.items():
            if genome_data['strain'] is None:
                strain = None
            else:
                strain = strain_instances[genome_data['strain']]
            if genome_data['sample'] is None:
                sample = None
            else:
                sample = sample_instances[genome_data['sample']]
            genome_instances.append(Genome(name=genome_data['name'],
                description=genome_data['description'],
                strain=strain,
                sample=sample,
                contigs=genome_data['contigs'],
                size=genome_data['size'],
                genes=genome_data['genes'],
                json_url=genome_data['json_url'],
                pub_date=genome_data['pub_date'],
                external_url=genome_data['external_url'],
                external_id=genome_data['external_id'],
                gbk_filepath = genome_data['gbk_filepath'],
                taxon=genome_data['taxon']
                ))
        Genome.objects.bulk_create(genome_instances, batch_size=1000)
        print('Genomes imported')

        # write eggnog descriptions
        eggnog_description_instances = []
        for fingerprint, description in self.mappings['description'].items():
            eggnog_description_instances.append(Eggnog_description(fingerprint=fingerprint,
                description=description
                ))
        Eggnog_description.objects.bulk_create(eggnog_description_instances, batch_size=1000)
        print('Descriptions imported')

        # write protein data
        protein_instances = []
        saved_proteins = set([item['protein_hash'] for item in Protein.objects.values('protein_hash')])
        eggnog_descriptions = {description.fingerprint:description for description in Eggnog_description.objects.all()}
        for protein_hash, protein_data in self.protein_data.items():
            if protein_hash in saved_proteins:
                continue
            protein_instance = Protein(name=protein_data['name'],
                length=protein_data['length'],
                protein_hash=protein_hash,
                sequence=protein_data['sequence']
            )
            if 'taxonomic_group' in self.protein_data[protein_hash]:
                protein_instance.taxonomy_id = taxon_instances[self.protein_data[protein_hash]['taxonomic_group']]
            if 'eggnog_description' in self.protein_data[protein_hash]:
                protein_instance.eggnog_description = eggnog_descriptions[self.protein_data[protein_hash]['eggnog_description']]
            protein_instances.append(protein_instance)
        Protein.objects.bulk_create(protein_instances, batch_size=1000)
        print('Proteins imported')
        
        # write contig data
        contig_instances = []
        genome_ids = {genome.name:genome for genome in Genome.objects.all()}
        for contig_uid,contig_entry in self.contig_data.items():
            contig_instances.append(Contig(contig_id=contig_entry['contig_id'],
                name=contig_entry['name'],
                size=contig_entry['size'],
                genome=genome_ids[contig_entry['genome']]
                ))
        Contig.objects.bulk_create(contig_instances, batch_size=1000)
        print('Contigs imported')
        
        # write gene data
        gene_instances = []
        contig_ids = {(contig.genome.name, contig.contig_id):contig for contig in Contig.objects.all()}
        protein_ids = {protein.protein_hash:protein for protein in Protein.objects.all()}
        for gene_id,gene_entry in self.gene_data.items():
            gene_instance = Gene(name=gene_entry['name'],
                locus_tag=gene_entry['locus_tag'],
                contig=contig_ids[(gene_entry['genome_id'], gene_entry['contig_id'])],
                type=gene_entry['type'],
                start=gene_entry['start'],
                end=gene_entry['end'],
                strand=gene_entry['strand'],
                genome=genome_ids[gene_entry['genome_id']]
                )
            if 'protein_hash' in gene_entry:
                gene_instance.protein = protein_ids[gene_entry['protein_hash']]
            if 'function' in gene_entry:
                gene_instance.function = gene_entry['function']
            gene_instances.append(gene_instance)
        Gene.objects.bulk_create(gene_instances, batch_size=1000)
        print('Genes imported')

        # write kegg reactions data
        existing_ref_data_objects = set([item['kegg_id'] for item in Kegg_reaction.objects.values('kegg_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['kegg_reaction'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Kegg_reaction(kegg_id=reference_id,
                description=reference_name
                ))
        Kegg_reaction.objects.bulk_create(reference_data_instances, batch_size=1000)
        
        relations = []
        reference_data_objects = {ref_object.kegg_id:ref_object for ref_object in Kegg_reaction.objects.all()}
        for (protein_hash, reference_id) in self.relations['kegg_reaction']:
            relations.append(Protein.kegg_reactions.through(protein=protein_ids[protein_hash], kegg_reaction=reference_data_objects[reference_id]))
        Protein.kegg_reactions.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('KEGG reactions imported')

        # write kegg pathways data
        existing_ref_data_objects = set([item['kegg_id'] for item in Kegg_pathway.objects.values('kegg_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['kegg_pathway'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Kegg_pathway(kegg_id=reference_id,
                description=reference_name
                ))
        Kegg_pathway.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.kegg_id:ref_object for ref_object in Kegg_pathway.objects.all()}
        for (protein_hash, reference_id) in self.relations['kegg_pathway']:
            relations.append(Protein.kegg_pathways.through(protein=protein_ids[protein_hash], kegg_pathway=reference_data_objects[reference_id]))
        Protein.kegg_pathways.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('KEGG pathways imported')

        # write kegg orthologs data
        existing_ref_data_objects = set([item['kegg_id'] for item in Kegg_ortholog.objects.values('kegg_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['kegg_ortholog'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Kegg_ortholog(kegg_id=reference_id,
                description=reference_name
                ))
        Kegg_ortholog.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.kegg_id:ref_object for ref_object in Kegg_ortholog.objects.all()}
        for (protein_hash, reference_id) in self.relations['kegg_ortholog']:
            relations.append(Protein.kegg_orthologs.through(protein=protein_ids[protein_hash], kegg_ortholog=reference_data_objects[reference_id]))
        Protein.kegg_orthologs.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('KEGG orthologs imported')

        # write go data
        existing_ref_data_objects = set([item['go_id'] for item in Go_term.objects.values('go_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['go'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Go_term(go_id=reference_id,
                go_namespace=reference_name[1],
                description=reference_name[0]
                ))
        Go_term.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.go_id:ref_object for ref_object in Go_term.objects.all()}
        for (protein_hash, reference_id) in self.relations['go']:
            relations.append(Protein.go_terms.through(protein=protein_ids[protein_hash], go_term=reference_data_objects[reference_id]))
        Protein.go_terms.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('GO terms imported')

        # write og data
        # orthogroup id is eggnog_id + '@' + taxonomy_id
        existing_ref_data_objects = set([item['eggnog_id'] + '@' + item['taxon__taxonomy_id'] for item in Ortholog_group.objects.values('eggnog_id', 'taxon__taxonomy_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['og'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Ortholog_group(eggnog_id=reference_name[0],
                taxon=taxon_instances[reference_name[1]]
                ))
        Ortholog_group.objects.bulk_create(reference_data_instances, batch_size=1000)
        
        og_ids = {(og.eggnog_id, og.taxon.taxonomy_id):og for og in Ortholog_group.objects.all()}
        relations = []
        for (protein_hash, eggnog_id, taxonomy_id) in self.relations['og']:
            relations.append(Protein.ortholog_groups.through(protein=protein_ids[protein_hash], ortholog_group=og_ids[(eggnog_id, taxonomy_id)]))
        Protein.ortholog_groups.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('EggNOG orthologs imported')

        # write ec data
        existing_ref_data_objects = set([item['ec_number'] for item in Ec_number.objects.values('ec_number')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['ec'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Ec_number(ec_number=reference_id,
                description=reference_name
                ))
        Ec_number.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.ec_number:ref_object for ref_object in Ec_number.objects.all()}
        for (protein_hash, reference_id) in self.relations['ec']:
            relations.append(Protein.ec_numbers.through(protein=protein_ids[protein_hash], ec_number=reference_data_objects[reference_id]))
        Protein.ec_numbers.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('EC numbers imported')

        # write tc data
        existing_ref_data_objects = set([item['tc_id'] for item in Tc_family.objects.values('tc_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['tc'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Tc_family(tc_id=reference_id,
                description=reference_name
                ))
        Tc_family.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.tc_id:ref_object for ref_object in Tc_family.objects.all()}
        for (protein_hash, reference_id) in self.relations['tc']:
            relations.append(Protein.tc_families.through(protein=protein_ids[protein_hash], tc_family=reference_data_objects[reference_id]))
        Protein.tc_families.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('TC families imported')

        # write cazy data
        existing_ref_data_objects = set([item['cazy_id'] for item in Cazy_family.objects.values('cazy_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['cazy'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Cazy_family(cazy_id=reference_id,
                description=reference_name
                ))
        Cazy_family.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.cazy_id:ref_object for ref_object in Cazy_family.objects.all()}
        for (protein_hash, reference_id) in self.relations['cazy']:
            relations.append(Protein.cazy_families.through(protein=protein_ids[protein_hash], cazy_family=reference_data_objects[reference_id]))
        Protein.cazy_families.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('CAZy families imported')

        # write cog data
        existing_ref_data_objects = set([item['cog_id'] for item in Cog_class.objects.values('cog_id')])
        reference_data_instances = []
        for reference_id,reference_name in self.mappings['cog'].items():
            if reference_id in existing_ref_data_objects:
                continue
            reference_data_instances.append(Cog_class(cog_id=reference_id,
                description=reference_name
                ))
        Cog_class.objects.bulk_create(reference_data_instances, batch_size=1000)

        relations = []
        reference_data_objects = {ref_object.cog_id:ref_object for ref_object in Cog_class.objects.all()}
        for (protein_hash, reference_id) in self.relations['cog']:
            relations.append(Protein.cog_classes.through(protein=protein_ids[protein_hash], cog_class=reference_data_objects[reference_id]))
        Protein.cog_classes.through.objects.bulk_create(relations, batch_size = 1000, ignore_conflicts=True)
        print('COG classes imported')
    
    def export_gff(self, genome_id):
        result = os.path.join(self.config['cgcms.temp_dir'], genome_id + '.gff3')
        if 'kbase.us/' in self.inputgenomes[genome_id]['url']:
            source = 'KBase'
        elif 'ncbi.nlm.nih.gov/' in self.inputgenomes[genome_id]['url']:
            source = 'NCBI'
        else:
            source = 'Unknown_source'
        with open(result, 'w') as outfile:
            outfile.write('##gff-version 3\n')
            genes = defaultdict(dict)
            for contig in Contig.objects.filter(genome__name=genome_id).order_by('contig_id'):
                outfile.write('##sequence-region ' + contig.contig_id + '\n')
                for gene in Gene.objects.filter(genome__name=genome_id).filter(contig__contig_id=contig.contig_id).order_by('start'):
                    if gene.strand == 1:
                        strand = '+'
                    else:
                        strand = '-'
                    column9 = ['ID=' + gene.locus_tag, 'locus_tag='+ gene.locus_tag, 'internal_id=' + gene.genome.name + '/' + gene.locus_tag]
                    if gene.name != '':
                        column9.append('name=' + gene.name)
                    if gene.function != '':
                        column9.append('product=' + gene.function)
                    row = [contig.contig_id, source, gene.type, str(gene.start), str(gene.end), '.', strand,
                              '0', '; '.join(column9)]
                    outfile.write(('\t').join(row) + '\n')
                for operon in Operon.objects.filter(genome__name=genome_id).filter(contig__contig_id=contig.contig_id).order_by('start'):
                    if operon.strand == 1:
                        strand = '+'
                    else:
                        strand = '-'
                    column9 = ['ID=' + operon.name, 'internal_id=' + operon.genome.name + '/' + operon.name]
                    row = [contig.contig_id, source, 'operon', str(operon.start), str(operon.end), '.', strand,
                              '0', '; '.join(column9)]
                    outfile.write(('\t').join(row) + '\n')
                    
        return result
        
    def export_jbrowse_data(self):
        for genome in self.inputgenomes:
            gff3_file = self.export_gff(genome)
            genome_fasta = os.path.join(self.config['cgcms.temp_dir'], genome + '.contigs.fasta')
            temp_dir = os.path.join(self.config['cgcms.temp_dir'], genome)
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            os.mkdir(temp_dir)
            cmd = [self.config['cgcms.prepare_refseqs_command'], '--fasta', genome_fasta, '--out', temp_dir]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)
            
            cmd = [self.config['cgcms.flatfile_to_json_command'], '--gff', gff3_file, '--out', temp_dir, '--trackLabel', 'CDSs', '--type', 'CDS', '--classname', 'exon' ]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)

            cmd = [self.config['cgcms.flatfile_to_json_command'], '--gff', gff3_file, '--out', temp_dir, '--trackLabel',
                   'Pseudogenes', '--type', 'pseudogene', '--classname', 'feature4' ]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)

            cmd = [self.config['cgcms.flatfile_to_json_command'], '--gff', gff3_file, '--out', temp_dir, '--trackLabel',
                   'tRNAs', '--type', 'tRNA', '--classname', 'transcript-UTR' ]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)

            cmd = [self.config['cgcms.flatfile_to_json_command'], '--gff', gff3_file, '--out', temp_dir, '--trackLabel',
                   'rRNAs', '--type', 'rRNA', '--classname', 'feature3' ]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)

            cmd = [self.config['cgcms.flatfile_to_json_command'], '--gff', gff3_file, '--out', temp_dir, '--trackLabel',
                   'Operons', '--type', 'operon', '--classname', 'feature4' ]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)

            cmd = [self.config['cgcms.generate_names_command'], '--verbose', '--out', temp_dir]
            print(' '.join(cmd))
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            if proc.returncode != 0:
                # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
                # pylint: disable=no-member
                raise CalledProcessError(proc.returncode, proc.args)
            dest_dir = os.path.join(self.config['cgcms.json_dir'], genome)
            if os.path.exists(dest_dir):
                shutil.rmtree(dest_dir)
            shutil.copytree(temp_dir, dest_dir)
            if os.path.exists(os.path.join(self.config['cgcms.temp_dir'], 'trackList.json')):
                shutil.copyfile(os.path.join(self.config['cgcms.temp_dir'], 'trackList.json'), os.path.join(self.config['cgcms.json_dir'], genome, 'trackList.json'))
            shutil.rmtree(temp_dir)

    def copy_static_files(self):
        for directory in self.staticfiles:
            for filename in self.staticfiles[directory]:
                shutil.copyfile(os.path.join(self.config['cgcms.temp_dir'], filename), os.path.join(directory, os.path.basename(filename)))

    def create_search_databases(self):
        # Make blast search databases
        cmd = ['makeblastdb', '-dbtype', 'nucl', '-in', self.config['cgcms.search_db_nucl'], '-out', os.path.join(self.config['cgcms.search_db_dir'], 'blast_nucl')]
        print(' '.join(cmd))
        with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
            for line in proc.stdout:
                print(line, end='')
        if proc.returncode != 0:
            # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
            # pylint: disable=no-member
            raise CalledProcessError(proc.returncode, proc.args)

        cmd = ['makeblastdb', '-dbtype', 'prot', '-in', self.config['cgcms.search_db_prot'], '-out', os.path.join(self.config['cgcms.search_db_dir'], 'blast_prot')]
        print(' '.join(cmd))
        with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
            for line in proc.stdout:
                print(line, end='')
        if proc.returncode != 0:
            # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
            # pylint: disable=no-member
            raise CalledProcessError(proc.returncode, proc.args)

#    def append_eggnog_stored_data(self, new_eggnog_data):
#        with open(self.config['genomes.eggnog_output_stored'], 'a') as outfile:
#            for line in new_eggnog_data:
#                outfile.write(line)
        
    def delete_search_databases(self):
        for filename in os.listdir(self.config['cgcms.search_db_dir']):
            if filename.startswith('mmseqs_') or filename.startswith('blast_'):
                os.remove(os.path.join(self.config['cgcms.search_db_dir'], filename))
        
    def cleanup(self):
        os.remove(os.path.join(self.config['cgcms.temp_dir'], os.path.basename(self.config['cgcms.search_db_nucl'])))
        os.remove(os.path.join(self.config['cgcms.temp_dir'], os.path.basename(self.config['cgcms.search_db_prot'])))
        os.remove(self.config['cgcms.search_db_nucl'])
        os.remove(self.config['cgcms.search_db_prot'])
        for filename in os.listdir(self.config['cgcms.eggnog_outdir']):
            os.remove(os.path.join(self.config['cgcms.eggnog_outdir'], filename))
        
    def import_genomes(self, in_file):
        '''
            This function contains a workflow for genome import.
            
            Input parameters:
              in_file(str): full path to the file with genome list
              
            Input file must contain five fields in each row:
                1. path to gbk file
                2. genome_id: unique genome identifier
                3. strain_id: unique strain identifier or empty string. Either strain_id or sample_id can be empty, but not both of them.
                4. sample_id: unique sample identifier or empty string. Either strain_id or sample_id can be empty, but not both of them.
                5. url: external URL referring to genome source (KBase, NCBI, etc.). Can be empty.
                6. external_id: genome identifier on the external resource, which is used as a text for external link. Should not be empty, even if there is no URL.
            Lines starting with # symbol are ignored

            This function does not return anything.
        '''
        print('Trying to create directories for static files')
        self.make_dirs()
        print('Reading file of genomes', in_file)
        # read genome files list. Populate self.inputgenomes
        self.load_genome_list(in_file)
        
        print('Checking genome names')
        # check if any genomes already exist in the database
        self.check_genomes()
        
        print('Loading genomes')
        # Check if genome files exist. Update self.genomefiles. TODO: download genomes, if needed. 
        self.download_genomes()
        
        print('Preparing strain data')
        # make strain data for upload
        self.prepare_strain_data()
        
        print('Preparing sample data')
        # make sample data for upload
        self.prepare_sample_data()
        
        print('Preparing taxonomy data')
        # make taxonomy data file for upload
        self.prepare_taxonomy_data()
        
        print('Loading proteins from the database')
        # read protein IDs and hashes
        self.load_proteins()

        print('Processing GBK files')
        # make genome, contigs, genes data files
        self.process_gbk()

        print('Running eggnog-mapper')
        # Close MySQL connection before starting eggnog-mapper because it may run for days resulting in "MySQL server has gone away" error
        connection.close()
        # run eggnog-mapper for all proteins
        eggnog_outfile = self.run_eggnog_mapper()
        # TODO: remove mockup and uncomment run_eggnog_mapper call if commented out
        # eggnog_outfile = os.path.join(self.config['cgcms.temp_dir'], 'eggnog_mapper_output.emapper.annotations')
        
        print('Reading eggnog-mapper output')
        # separate eggnog-mapper output by genome?
        new_eggnog_data = self.parse_eggnog_output(eggnog_outfile)
        
        print('Prepare eggnog mappings')
        # make mappings and relations data
        self.make_mappings()

        print('Prepare eggnog ortholog mappings')
        # make og data file for upload
        self.prepare_og_data()
        
        print('Prepare eggnog descriptions')
        # make eggnog description data file for upload
        self.prepare_eggnog_description_data()
        
        print('Populating mysql database')
        # At this point, genes and annotations are actually written to the database
        self.write_data()
        
        print('Predict operons')
        # Export contigs and proteins and run POEM_py3
        operons_data = self.predict_operons()
        # Write operons to the database
        self.create_operons(operons_data)
        
        #export JBrowse data
        print('Creating Jbrowser files')
        self.export_jbrowse_data()
        print('Exporting protein sequences')
        self.export_proteins()
        
        print('At this moment, we are ready to delete search databases')
        self.delete_search_databases()
        
        # copy static files
        print('Copying static files')
        self.copy_static_files()
        print('Generating search databases')
        self.create_search_databases()

#        self.append_eggnog_stored_data(new_eggnog_data)

        # TODO:
        #new_genome_ids = [x['id'] for x in Genome.objects.filter(name__in = self.genome_data.keys()).values('id')]
        # annotator = Annotator(self.config_path)
        # annotator.update_pfam_domains(new_genome_ids)
        # annotator.update_tigrfam_domains(new_genome_ids)

        # delete temp files
        print('Removing temporary files')
        self.cleanup()

    def export_proteins(self):
        prot_db_file = self.config['cgcms.search_db_prot']
        with open(os.path.join(self.config['cgcms.temp_dir'], os.path.basename(prot_db_file)), 'w') as outfile:
            for protein in Protein.objects.values('protein_hash', 'sequence'):
                outfile.write('>' + protein['protein_hash'] + '\n' + protein['sequence'] + '\n')
        self.staticfiles[self.config['cgcms.search_db_dir']].append(os.path.basename(prot_db_file))


    def export_contigs(self):
        result = set()
        nucl_db_file = self.config['cgcms.search_db_nucl']
        with open(os.path.join(self.config['cgcms.temp_dir'], os.path.basename(nucl_db_file)), 'w') as outfile:
            for item in Genome.objects.values('name', 'gbk_filepath'):
                if item['name'] in self.inputgenomes:
                    continue
                gbk_handle = gzip.open(item['gbk_filepath'], 'rt')
                parser = GenBank.parse(gbk_handle)
                for gbk_record in parser:
                    contig_uid = gbk_record.locus + '|' + item['name']
                    outfile.write('>' + contig_uid + '\n' + ''.join(gbk_record.sequence) + '\n')
                    result.add(contig_uid)
                gbk_handle.close()
        self.staticfiles[self.config['cgcms.search_db_dir']].append(os.path.basename(nucl_db_file))
        return result


    def create_operon(self, operon, operon_index, genes):
        result = []
        genome_id = operon[0][0].split('|')[-1]
        operon_id = genome_id + '_operon_' + str(operon_index)
        operon_start = int(operon[0][0].split('|')[4])
        operon_end = operon[-1][1].split('|')[5]
        operon_end, contig_id = operon_end.split('$$')
        operon_end = int(operon_end)
        if operon[0][2] == '+':
            operon_strand = 1
        elif operon[0][2] == '-':
            operon_strand = -1
        else:
            raise ValueError('Unknown strand identifier: ' + operon[0][2])
        genome = Genome.objects.get(name = genome_id)
        operon_instance = Operon.objects.create(name = operon_id,
                                                start = operon_start,
                                                end = operon_end,
                                                strand = operon_strand,
                                                genome = genome,
                                                contig = Contig.objects.get(contig_id=contig_id, genome=genome))
        operon_instance.save()
        for gene_data in [elem[0] for elem in operon] + [operon[-1][1]]:
            locus_tag = genes[gene_data[5:].split('|')[0]]
            gene = Gene.objects.get(genome=genome, locus_tag=locus_tag)
            gene.operon = operon_instance
            gene.save()
            result.append('\t'.join([operon_id, genome_id, contig_id, str(operon_start), str(operon_end), str(operon_strand), locus_tag]))
        return '\n'.join(result) + '\n'
    
    def predict_operons(self):
        working_dir = os.path.join(self.config['cgcms.temp_dir'], 'poem-temp')
        if os.path.exists(working_dir) and os.path.isdir(working_dir):
            shutil.rmtree(working_dir)
        os.mkdir(working_dir)
        
        # Export nucleotide sequences of all new genomes
        with open(os.path.join(working_dir, 'input.fsa'), 'w') as outfile:
            for genome_id in self.inputgenomes:
                gbk_file = self.inputgenomes[genome_id]['gbk']
                if gbk_file.endswith('.gz'):
                    gbk_handle = gzip.open(gbk_file, 'rt')
                else:
                    gbk_handle = open(gbk_file, 'r')
                parser = GenBank.parse(gbk_handle)
                for gbk_record in parser:
                    contig_sequence = gbk_record.sequence
                    contig_id = gbk_record.locus
                    outfile.write('>' + contig_id + '|' + genome_id + '\n' + ''.join(contig_sequence) + '\n')
                gbk_handle.close()
        
        genes = {}
        with open(os.path.join(working_dir, 'input.fsa_prod_aa.fsa'), 'w') as outfile:
            gene_index = 0
            for genome_id in self.inputgenomes:
                for gene in Gene.objects.filter(genome__name = genome_id).select_related('protein', 'genome', 'contig'):
                    if gene.protein is not None:
                        gene_index += 1
                        genes['%08d'%int(gene_index)] = gene.locus_tag
                        outfile.write('>' + ' # '.join([gene.contig.contig_id + '|' + gene.genome.name + '_' + str(gene_index),
                                                        str(gene.start),
                                                        str(gene.end),
                                                        str(gene.strand),
                                                        'ID=' + gene.locus_tag]) + '\n')
                        outfile.write(gene.protein.sequence + '\n')
        poem_script = os.path.join(working_dir, 'run_poem.sh')
        with open(poem_script, 'w') as outfile:
            outfile.write('#!/bin/bash\n')
            outfile.write('source ~/miniconda3/etc/profile.d/conda.sh\n')
            outfile.write('conda activate poem_py3\n')
            outfile.write('bash /home/aekazakov/Soft/POEM_py3k/bin/run_poem_cgcms.sh -f ' + working_dir + ' -a n -p pro\n')
            outfile.write('conda deactivate\n')
            
        cmd = ['/bin/bash', poem_script]
        print(' '.join(cmd))
        # Close MySQL connection before starting external process because it may run for too long resulting in "MySQL server has gone away" error
        connection.close()
        with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
            for line in proc.stdout:
                print(line, end='')
        if proc.returncode != 0:
            # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
            # pylint: disable=no-member
            raise CalledProcessError(proc.returncode, proc.args)
        
        poem_outfile = os.path.join(working_dir, 'input.fsa.adjacency')
        operons_data = defaultdict(list)
        with open (poem_outfile, 'r') as infile:
            operon_index = Operon.objects.count() + 1
            operon = []
            for line in infile:
                j = line[:-1].split('\t')
                if len(j) != 11:
                    continue
                if j[-1] != 'True' and j[-1] != 'False':
                    continue
                qid, strand, sid, label = j[0], j[2], j[5], j[-1]
                if label == 'False':
                    continue
                contig_id, genome_id = j[1].split('|')
                if operon and operon[-1][1] == qid and operon[-1][2] == strand:
                    operon.append([qid, sid, strand])
                else:
                    if operon:
                        operon_index += 1
                        operon_id = genome_id + '_operon_' + str(operon_index)
                        operon_start = int(operon[0][0].split('|')[4])
                        operon_end = operon[-1][1].split('|')[5]
                        operon_end = int(operon_end.split('$$')[0])
                        if operon[0][2] == '+':
                            operon_strand = 1
                        elif operon[0][2] == '-':
                            operon_strand = -1
                        operon_members = []
                        for gene_data in [elem[0] for elem in operon] + [operon[-1][1]]:
                            operon_members.append(genes[gene_data[5:].split('|')[0]])
                        operons_data[operon[0][4]].append([operon_id, operon_start, operon_end, operon_strand, operon[0][3], operon_members])
                        #outfile.write(self.create_operon(operon, operon_index, genes))
                    operon = [[qid, sid, strand, contig_id, genome_id]]
            if operon:
                operon_index += 1
                operon_id = genome_id + '_operon_' + str(operon_index)
                operon_start = int(operon[0][0].split('|')[4])
                operon_end = operon[-1][1].split('|')[5]
                operon_end = int(operon_end.split('$$')[0])
                if operon[0][2] == '+':
                    operon_strand = 1
                elif operon[0][2] == '-':
                    operon_strand = -1
                operon_members = []
                for gene_data in [elem[0] for elem in operon] + [operon[-1][1]]:
                    operon_members.append(genes[gene_data[5:].split('|')[0]])
                operons_data[operon[0][4]].append([operon_id, operon_start, operon_end, operon_strand, operon[0][3], operon_members])
                #outfile.write(self.create_operon(operon, operon_index, genes))
        with open(os.path.join(working_dir, 'poem_output.txt'), 'w') as outfile:
            for genome_id in operons_data:
                for item in operons_data[genome_id]:
                    outfile.write(genome_id + '\t' + '\t'.join([str(x) for x in item]) + '\n')
        return operons_data

    def create_operons(self, operons_data):
        operon_instances = []
        for genome_id in operons_data:
            genome = Genome.objects.get(name = genome_id)
            contigs = {item.contig_id:item for item in Contig.objects.filter(genome__name=genome_id)}
            print('Contigs:')
            print(contigs.keys())
            for operon in operons_data[genome_id]:
                operon_instance = Operon(name = operon[0],
                                         start = operon[1],
                                         end = operon[2],
                                         strand = operon[3],
                                         genome = genome,
                                         contig = contigs[operon[4]])
                operon_instances.append(operon_instance)
        Operon.objects.bulk_create(operon_instances, batch_size=1000)

        changed_genes = []
        for genome_id in operons_data:
            genes = {item.locus_tag:item for item in Gene.objects.filter(genome__name = genome_id)}
            operons = {item.name:item for item in Operon.objects.filter(genome__name = genome_id)}
            for operon in operons_data[genome_id]:
                operon_name = operon[0]
                for locus_tag in operon[5]:
                    gene = genes[locus_tag]
                    gene.operon = operons[operon_name]
                    changed_genes.append(gene)
        Gene.objects.bulk_update(changed_genes, ['operon'], batch_size=1000)
    
    def load_genome_data(self):
        """Populate self.inputgenomes with taxon data"""
        for genome in Genome.objects.all():
            if genome.strain is not None and genome.name in self.inputgenomes:
                self.inputgenomes[genome.name]['taxon'] = genome.strain.taxon.taxonomy_id

    def generate_static_files(self, in_file):
        if os.path.exists(in_file):
            print('Genomes file found', in_file)
        self.load_genome_list(in_file)
        self.load_genome_data()
        self.process_gbk()
            
        self.export_proteins()
        self.export_jbrowse_data()
        self.copy_static_files()
        self.create_search_databases()
        self.cleanup()
        
