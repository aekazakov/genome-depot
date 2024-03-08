""" This file contains functions for taxonomy import and update"""
import os
import tarfile
import logging
import shutil
import urllib.request
from pathlib import Path
from collections import defaultdict
from browser.models import Config, Taxon, Strain, Genome, Protein, Ortholog_group

logger = logging.getLogger("GenomeDepot")

def load_taxonomy(taxonomy_file, eggnog_taxonomy_file):
    logger.info('Loading taxonomy...')
    taxonomy = {}
    taxonomy_lookup = {}
    with open(taxonomy_file, 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            taxonomy_id = row[0]
            name = row[1]
            rank = row[2]
            parent = row[3]
            current = row[4]
            if current == '1':
                taxonomy_lookup[name] = taxonomy_id
            taxonomy[taxonomy_id] = {}
            taxonomy[taxonomy_id]['name'] = name
            taxonomy[taxonomy_id]['eggnog_taxid'] = taxonomy_id
            taxonomy[taxonomy_id]['rank'] = rank
            taxonomy[taxonomy_id]['parent'] = parent
    eggnog_lookup = {}
    with open(eggnog_taxonomy_file, 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            eggnog_lookup[row[0]] = row[1]
            taxonomy_lookup[row[2]] = row[1]
            taxonomy[row[1]]['eggnog_taxid'] = row[0]
    logger.info('Taxonomy loaded')
    return taxonomy, taxonomy_lookup, eggnog_lookup


def update_taxonomy():
    """ 
        Downloads taxonomy dump from NCBI FTP and updates taxonomy 
        data in the GenomeDepot database.
        
        1. Create temporary folder
        2. Download fresh NCBI taxonomy archive.
        3. Extract files.
        4. Read eggNOG taxonomy file from core.eggnog_taxonomy.
        5. Read taxonomy data from database.
        6. Update taxon entries in the database, if needed.
        7. Generate new taxonomy file.
        8. Delete old taxonomy file.
        9. Delete temporary folder.
    """
    ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    config = {}
    for item in Config.objects.values('param', 'value'):
        config[item['param']] = item['value']
    tmp_dir = os.path.join(config['core.temp_dir'], 'taxonomy')

    out_file = os.path.join(tmp_dir, 'ref_taxonomy.txt')
    eggnog_taxa = {}
    taxonomy = defaultdict(dict)
    old_taxonomy_names = {}
    taxonomy_id_lookup = {}
    filenames = ['names.dmp', 'nodes.dmp', 'merged.dmp']
    # create directory
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    ncbi_archive = os.path.join(tmp_dir, 'taxdump.tar.gz')
    with open(config['core.eggnog_taxonomy'], 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            eggnog_taxa[row[0]] = (row[1], row[2])

    with open(config['ref.taxonomy'], 'r') as infile:
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            taxonomy_id = row[0]
            name = row[1]
            parent = row[3]
            old_taxonomy_names[taxonomy_id] = (name, parent)

    logger.info('Downloading ' + ncbi_url)
    _ = urllib.request.urlretrieve(url=ncbi_url, filename=ncbi_archive)
    tar = tarfile.open(ncbi_archive)
    for filename in filenames:
        logger.info('Extracting ' + filename)
        with open(os.path.join(tmp_dir, filename), 'wb') as outfile:
            outfile.write(tar.extractfile(filename).read())

    logger.info('Reading names.dmp')
    with open(os.path.join(tmp_dir, 'names.dmp'), 'r') as f:
        for line in f:
            row = line.rstrip('\n\r').split('\t|\t')
            row[3] = row[3].split('\t')[0]
            if row[3] == 'scientific name':
                taxonomy[row[0]] = {}
                taxonomy[row[0]]['name'] = row[1]
                if row[0] in taxonomy_id_lookup:
                    logger.warn('Warining: non-unique taxon name ' + row[1] +
                          ' with IDs ' + taxonomy_id_lookup[row[1]] + ',' + row[0]
                          )
            if row[3] in ('scientific name', 'synonym', 'includes', 'equivalent name'):
                taxonomy_id_lookup[row[1]] = row[0]
    #initialize nodes
    logger.info('Reading nodes.dmp')
    with open(os.path.join(tmp_dir, 'nodes.dmp'), 'r') as f:
        for line in f:
            row = line.rstrip('\n\r').split('\t|\t')
            taxid = row[0]
            parent = row[1]
            rank = row[2]
            taxonomy[taxid]['parent'] = parent
            taxonomy[taxid]['rank'] = rank
    # inject 'Unknown' entry
    taxonomy['1']['rank'] = 'norank'
    taxonomy['0'] = {}
    taxonomy['0']['name'] = 'Unknown'
    taxonomy['0']['parent'] = '1'
    taxonomy['0']['rank'] = 'norank'
    taxonomy_id_lookup['Unknown'] = '0'
    logger.info('Reading merged.dmp')
    with open(os.path.join(tmp_dir, 'merged.dmp'), 'r') as f:
        for line in f:
            row = line.rstrip('\n\r').split('\t')
            old_id = row[0]
            new_id = row[2]
            if new_id in taxonomy:
                if old_id not in taxonomy:
                    taxonomy[old_id] = {}
                taxonomy[old_id]['name'] = taxonomy[new_id]['name']
                taxonomy[old_id]['parent'] = taxonomy[new_id]['parent']
                taxonomy[old_id]['rank'] = taxonomy[new_id]['rank']
                taxonomy[old_id]['merged'] = 1

    # update eggnog mappings
    for eggnog_id in eggnog_taxa:
        taxonomy_id = eggnog_taxa[eggnog_id][0]
        if taxonomy_id not in taxonomy:
            name = old_taxonomy_names[taxonomy_id][0]
            parent = old_taxonomy_names[taxonomy_id][1]
            if name in taxonomy_id_lookup:
                eggnog_taxa[eggnog_id] = (taxonomy_id_lookup[name], name)
            else:
                logger.warn(name + '(' + taxonomy_id + 
                      ') cannot be located. Transferring mapping to parent ' +
                      parent
                      )
                eggnog_taxa[eggnog_id] = (parent, taxonomy[parent]['name'])
        else:
            # Check if current entry is merged
            if 'merged' in taxonomy[taxonomy_id]:
                name = old_taxonomy_names[taxonomy_id][0]
                eggnog_taxa[eggnog_id] = (taxonomy_id_lookup[name], name)

    # write results
    with open(out_file, 'w') as outfile:
        for taxonomy_id in sorted(taxonomy.keys()):
            taxonomy_data = taxonomy[taxonomy_id]
            if 'merged' in taxonomy_data:
                outfile.write('\t'.join([taxonomy_id,
                                         taxonomy_data['name'],
                                         taxonomy_data['rank'],
                                         taxonomy_data['parent'],
                                         '0']) + '\n')
            else:
                outfile.write('\t'.join([taxonomy_id,
                                         taxonomy_data['name'],
                                         taxonomy_data['rank'],
                                         taxonomy_data['parent'],
                                         '1']) + '\n')
    
    eggnog_taxonomy_temp_file = os.path.join(tmp_dir, 'eggnog_taxonomy_rules.txt')
    with open(eggnog_taxonomy_temp_file, 'w') as outfile:
        for eggnog_id in eggnog_taxa:
            outfile.write('\t'.join([eggnog_id, eggnog_taxa[eggnog_id][0], eggnog_taxa[eggnog_id][1]]) + '\n')
        
    # Update Taxon entries
    logger.info('Updating taxonomy entries in the database')
    delete_taxonomy_ids = set()
    for taxon in Taxon.objects.all():
        # First, check by taxonomy ID in latest taxonomy merged records
        if taxon.taxonomy_id in taxonomy:
            if 'merged' in taxonomy[taxon.taxonomy_id]:
                new_name = taxonomy[taxon.taxonomy_id]['name']
                new_id = taxonomy_id_lookup[new_name]
                logger.info('Changing ' + taxon.name + ':' + taxon.taxonomy_id +
                'to' + new_name + ':' + new_id)
                # Check if new_id already exists in the database
                try:
                    new_taxon = Taxon.objects.get(taxonomy_id=new_id)
                    update_taxonomy_fields(taxon, new_taxon)
                    delete_taxonomy_ids.add(taxon.taxonomy_id)
                except Taxon.DoesNotExist:
                    # Change to new taxonomy ID
                    if taxon.eggnog_taxid not in eggnog_taxa:
                        taxon.eggnog_taxid = new_id
                    taxon.name = new_name
                    taxon.taxonomy_id = new_id
                    taxon.rank = taxonomy[new_id]['rank']
                    taxon.parent_id = taxonomy[new_id]['parent']
                    taxon.save()
    # Delete records
    for taxonomy_id in delete_taxonomy_ids:
        taxon = Taxon.objects.get(taxonomy_id=taxonomy_id)
        taxon.delete()
    
    for taxon in Taxon.objects.all():
        # First, check by taxonomy ID in latest taxonomy
        if taxon.taxonomy_id in taxonomy:
            if taxon.name != taxonomy[taxon.taxonomy_id]['name']:
                logger.info('Renaming ' +
                      taxon.name +
                      ' to ' +
                      taxonomy[taxon.taxonomy_id]['name']
                      )
                taxon.name = taxonomy[taxon.taxonomy_id]['name']
                taxon.save()
            if taxon.rank != taxonomy[taxon.taxonomy_id]['rank']:
                taxon.rank = taxonomy[taxon.taxonomy_id]['rank']
                taxon.save()
            if taxon.parent_id != taxonomy[taxon.taxonomy_id]['parent']:
                taxon.parent_id = taxonomy[taxon.taxonomy_id]['parent']
                taxon.save()
        # Second, check by taxonomy ID in updated entries
        elif taxon.eggnog_taxid in eggnog_taxa:
            new_id = eggnog_taxa[taxon.eggnog_taxid]
            logger.info('Changing ' + old_id + ' to ' + new_id)
            taxon.taxonomy_id = new_id
            taxon.name = taxonomy[new_id]['name']
            taxon.rank = taxonomy[new_id]['rank']
            taxon.parent_id = taxonomy[new_id]['parent']
            taxon.save()
        else:
            logger.error('Taxon ' +
                  taxon.name +
                  '(' +
                  taxon.taxonomy_id +
                  ') cannot be found in NCBI Taxonomy. Check the name and taxonomy ID'+
                  ' and fix the record in the GenomeDepot database.'
                  )
                  
    # Taxonomy consistency check:
    # All parents must exist in the database and be unique
    create_ids = set()
    for taxon in Taxon.objects.all():
        if not Taxon.objects.filter(taxonomy_id=taxon.parent_id).exists():
            create_ids.add(taxon.parent_id)
    
    for taxonomy_id in create_ids:
        if taxonomy_id in taxonomy:
            create_taxonomy_records(taxonomy_id, taxonomy, eggnog_taxa)
        elif taxonomy_id in old_taxonomy_names:
            old_name = old_taxonomy_names[taxonomy_id][0]
            if old_name in taxonomy_id_lookup:
                new_taxonomy_id = taxonomy_id_lookup[old_name]
                create_taxonomy_records(new_taxonomy_id, taxonomy, eggnog_taxa)
        else:
            logger.error('Taxonomy ID' +
                  taxonomy_id +
                  'cannot be found in the downloaded NCBI Taxonomy. Check the NCBI web-site'+
                  ' and fix the record in the GenomeDepot database.'
                  )
            
    os.rename(eggnog_taxonomy_temp_file, config['core.eggnog_taxonomy'])
    os.rename(out_file, config['ref.taxonomy'])
    shutil.rmtree(tmp_dir)

def create_taxonomy_records(taxonomy_id, taxonomy, eggnog_taxa):
    if taxonomy_id not in taxonomy:
        return
    if Taxon.objects.filter(taxonomy_id=taxonomy_id).exists():
        return
    eggnog_taxid = taxonomy_id
    if taxonomy_id in eggnog_taxa:
        eggnog_taxid = eggnog_taxa[taxonomy_id][0]
    logger.info('Creating ' + taxonomy[taxonomy_id]['name'] + '(' + taxonomy_id + ')')
    new_parent_taxon=Taxon.objects.create(
        taxonomy_id=taxonomy_id,
        eggnog_taxid=eggnog_taxid,
        name=taxonomy[taxonomy_id]['name'],
        rank=taxonomy[taxonomy_id]['rank'],
        parent_id=taxonomy[taxonomy_id]['parent']
        )
    new_parent_id = new_parent_taxon.parent_id
    while new_parent_id != '1':
        if new_parent_id not in taxonomy:
            logger.error('Cannot create record for ' + new_parent_id + ': not in NCBI taxonomy')
            break
        if Taxon.objects.filter(taxonomy_id=new_parent_id).exists():
            break
        eggnog_taxid = new_parent_id
        if new_parent_id in eggnog_taxa:
            eggnog_taxid = eggnog_taxa[new_parent_id][0]
        logger.info('Creating ' + taxonomy[new_parent_id]['name'] + '(' + new_parent_id + ')')
        new_parent_taxon=Taxon.objects.create(
            taxonomy_id=new_parent_id,
            eggnog_taxid=eggnog_taxid,
            name=taxonomy[new_parent_id]['name'],
            rank=taxonomy[new_parent_id]['rank'],
            parent_id=taxonomy[new_parent_id]['parent']
            )
        new_parent_id = new_parent_taxon.parent_id
    

def update_taxonomy_fields(old_taxon, new_taxon):
    for strain in Strain.objects.filter(taxon=old_taxon):
        strain.taxon = new_taxon
        strain.save()

    for genome in Genome.objects.filter(taxon=old_taxon):
        genome.taxon = new_taxon
        genome.save()

    for protein in Protein.objects.filter(taxonomy_id=old_taxon):
        protein.taxonomy_id = new_taxon
        protein.save()

    for ogroup in Ortholog_group.objects.filter(taxon=old_taxon):
        ogroup.taxon = new_taxon
        ogroup.save()
