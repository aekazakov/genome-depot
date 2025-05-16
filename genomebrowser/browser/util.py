"""
    Various utility functions
"""
import os
import re
import gzip
import logging
import shutil
import urllib
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqFeature
from Bio import GenBank
from django.core.management.base import CommandError

from browser.models import Config
from browser.models import Contig
from browser.models import Site
from browser.models import Annotation
from browser.models import Eggnog_description
from browser.models import Cog_class
from browser.models import Kegg_reaction
from browser.models import Kegg_pathway
from browser.models import Kegg_ortholog
from browser.models import Go_term
from browser.models import Cazy_family
from browser.models import Ec_number
from browser.models import Tc_family
from browser.models import Ortholog_group
from browser.models import Gene
from browser.models import Protein
from browser.models import Genome
from browser.models import Strain_metadata
from browser.models import Strain
from browser.models import Sample_metadata
from browser.models import Sample
from browser.models import Tag
from browser.models import Taxon
from browser.pipeline.genome_import import Importer
from browser.pipeline.util import export_proteins

logger = logging.getLogger("GenomeDepot")

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def export_genome(genome, output_buffer):
    '''
        Reads one genome from GenBank file,
        adds gene mappings and annotations from the GenomeDepot database,
        writes the genome in GenBank format into the output buffer
    '''
                                        
    # Put genes here genedata[contig_id][gene_id] = list of notes
    genedata = autovivify(2, list)
    # Put operons here operondata[contig_id][gene_id] = 
    # (name, start, end, strand)
    operondata = defaultdict(dict)
    # Put sites here operondata[contig_id][gene_id] = 
    # list of (regulon, regulator, start, end, strand)
    sitedata = autovivify(2, list)
    for contig_id in list(Contig.objects
        .values_list('contig_id',flat=True)
        .filter(genome=genome)
    ):
        operons_seen = {}
        sites_seen = set()
        for gene in Gene.objects.select_related(
            'protein',
            'operon'
        ).prefetch_related(
            'protein__ortholog_groups',
            'protein__ortholog_groups__taxon',
            'protein__kegg_orthologs',
            'protein__kegg_pathways',
            'protein__kegg_reactions',
            'protein__ec_numbers',
            'protein__go_terms',
            'protein__tc_families',
            'protein__cog_classes',
            'protein__cazy_families'
        ).filter(
            genome=genome,
            contig__contig_id=contig_id
        ):
            if gene.protein:
                for item in gene.protein.ortholog_groups.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'EGGNOG5:'
                        + item.eggnog_id + '['
                        + item.taxon.name + '] [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.kegg_orthologs.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.kegg_pathways.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.kegg_reactions.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.ec_numbers.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'EC:'
                        + item.ec_number + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.go_terms.all():
                    genedata[contig_id][gene.locus_tag].append(
                        ''
                        + item.go_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.tc_families.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'TC:'
                        + item.tc_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.cog_classes.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'COG:'
                        + item.cog_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
                for item in gene.protein.cazy_families.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'CAZy:'
                        + item.cazy_id + ' ('
                        + item.description
                        + ') [GD/eggNOG-mapper]'
                    )
            for annotation in Annotation.objects.filter(gene_id=gene):
                genedata[contig_id][gene.locus_tag].append(
                    annotation.source + ':' +
                    ''.join([i for i in annotation.value if ord(i)<128])
                    + ' (' + ''.join([i for i in annotation.note if ord(i)<128])
                    + ') [GD/'
                    + annotation.source + ']'
                )
                
            if gene.operon:
                if gene.operon.id not in operons_seen:
                    operondata[contig_id][gene.locus_tag] = (
                        gene.operon.name,
                        gene.operon.start,
                        gene.operon.end, 
                        gene.operon.strand
                    )
                    operons_seen[gene.operon.id] = gene.locus_tag
        sites = Site.objects.filter(
            genome=genome,
            contig__contig_id=contig_id
        ).prefetch_related(
            'operons',
            'genes',
            'regulon__regulators'
        ).select_related('regulon')
        for site in sites:
            if site.id in sites_seen:
                continue
            if site.genes:
                for gene in site.genes.all()[:1]:
                    for regulator in site.regulon.regulators.all():
                        sitedata[contig_id][gene.locus_tag].append(
                            (site.regulon.name,
                             regulator.locus_tag,
                             site.start,
                             site.end,
                             site.strand
                             )
                        )
                        sites_seen.add(site.id)
            elif site.operons:
                for operon in site.operons.all()[:1]:
                    if operon.id in operons_seen:
                        if site.regulon.regulators:
                            for regulator in site.regulon.regulators.all():
                                sitedata[contig_id][operons_seen[operon.id]]\
                                .append((
                                    site.regulon.name,
                                    regulator.locus_tag,
                                    site.start,
                                    site.end,
                                    site.strand
                                ))
                                sites_seen.add(site.id)
                        else:
                            sitedata[contig_id][operons_seen[operon.id]]\
                            .append((
                                site.regulon.name,
                                'unknown regulator',
                                site.start,
                                site.end,
                                site.strand
                            ))
                            sites_seen.add(site.id)
                                
                             
    genome_file = genome.gbk_filepath
    if genome_file.endswith('.gz'):
        fh = gzip.open(genome_file, 'rt')
    else:
        fh = open(genome_file, 'r')
    for seq_record in SeqIO.parse(fh, "genbank"):
        accession = seq_record.id
        if '.' in accession:
            accession = accession.split('.')[0]
        if accession in genedata:
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    try:
                        gene_id = feature.qualifiers['locus_tag'][0]
                    except KeyError:
                        feature_start = str(feature.location.start + 1)
                        if (
                            feature_start.startswith('>') or 
                            feature_start.startswith('<')
                        ):
                            feature_start = feature_start[1:]
                        gene_id = '_'.join([accession,
                                            feature_start,
                                            str(feature.location.end),
                                            str(feature.location.strand)
                                            ])
                        logger.warning('Using %s instead of locus tag',
                                       gene_id)
                    if gene_id in genedata[accession]:
                        if 'note' not in feature.qualifiers:
                            feature.qualifiers['note'] = []
                        for ref in genedata[accession][gene_id]:
                            feature.qualifiers['note'].append(ref)
        if accession in operondata:
            for gene_id, operon in operondata[accession].items():
                operon_feature = SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(
                        operon[1] - 1, operon[2], strand=operon[3]),
                    type='operon',
                    id=operon[0]
                    )
                operon_feature.strand = operon_feature.location.strand
                operon_feature.qualifiers['operon'] = []
                operon_feature.qualifiers['operon'].append(operon[0])
                seq_record.features.append(operon_feature)
        if accession in sitedata:
            for gene_id, site_items in sitedata[accession].items():
                for site in site_items:
                    site_feature = SeqFeature.SeqFeature(
                        SeqFeature.FeatureLocation(
                            site[2] - 1, site[3],
                            strand=site[4]
                            ),
                        type='site',
                        id=site[0] + '_site_at_' + gene_id
                        )
                site_feature.strand = site_feature.location.strand
                if 'note' not in site_feature.qualifiers:
                    site_feature.qualifiers['note'] = []
                site_feature.qualifiers['note'].append('Regulon ' + 
                                                   site[0] +
                                                   '. Binding site of ' +
                                                   site[1])
                seq_record.features.append(site_feature)
        seq_record.features = sorted(seq_record.features,
                                     key=lambda x: x.location.start
                                     )
        SeqIO.write(seq_record, output_buffer, 'genbank')
        

def export_genomes(out_dir, genome_ids = []):
    '''
        Exports genomes as genbank files
    '''
    if not os.path.exists(out_dir):
        logger.error('Output directory does not exist: ' + out_dir)
        return
    if not genome_ids:
        genome_ids = list(Genome.objects.values_list('name',
                                                     flat=True
                                                     ).distinct())
    for genome_id in genome_ids:
        try:
            genome = Genome.objects.get(name = genome_id)
        except Genome.DoesNotExist:
            logger.error('Genome not found: ' + str(genome_id))
            continue
        with gzip.open(os.path.join(out_dir, 
                                    genome_id + '.gbff.gz'
                                    ), 'wt') as outfile:
            export_genome(genome, outfile)                            
        logger.info(genome_id + ' exported')
    

def download_ncbi_assembly(assembly_id, email, upload_dir):
    """Download genbank assemblies for a given search term.
    Args:
        assembly_id: search term, usually assembly RefSeq of GenBank accession
        email: email address, required for NCBI requests
        upload_dir: folder to save to
    """
    if not os.path.exists(upload_dir):
        try:
            os.mkdir(upload_dir)
        except OSError:
            logger.error(
                "Critical error: genome import pipeline can't create genome " +
                "download directory (%s)!" % (upload_dir)
            )
            raise
    Entrez.email = email
    handle = Entrez.esearch(db="assembly", term=assembly_id, retmax='200')
    record = Entrez.read(handle)
    try:
        record_id = record['IdList'][0]
    except IndexError:
        return ''
    esummary_handle = Entrez.esummary(db="assembly", id=record_id, report="full")
    summary = Entrez.read(esummary_handle)
    if 'FtpPath_RefSeq' in summary['DocumentSummarySet']['DocumentSummary'][0]:
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    elif 'FtpPath_GenBank' in summary['DocumentSummarySet']['DocumentSummary'][0]:
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
    else:
        url = ''
    if url == '':
        return ''
    link = os.path.join(url, os.path.basename(url) + '_genomic.gbff.gz')
    outfile = os.path.join(upload_dir, assembly_id + '_genomic.gbff.gz')
    logger.info('Getting assembly from ' + str(link))
    urllib.request.urlretrieve(link, filename=outfile)
    return outfile

    
def delete_all_data(confirm=True):
    '''
    This function deletes all data from the GenomeDepot database except configs.
    It would be called from CLI only.
    User must explicitely confirm data deletion!
    This function doesn't delete static files.
    '''
    print('')
    while confirm:
        answer = input('All data in the GenomeDepot database will be deleted.'
            'This action cannot be reversed. Continue? (y/n)'
        )
        if answer.lower() in ["y","yes"]:
            break
        elif answer.lower() in ["n","no"]:
            print('Exiting.')
            return
        else:
            print('Please answer yes/y or no/n.')
    # Delete all mappings before deleting genes
    logger.info("Deleting annotations...")
    Annotation.objects.all().delete()
    logger.info("Deleting classifications...")
    Eggnog_description.objects.all().delete()
    Cog_class.objects.all().delete()
    Kegg_reaction.objects.all().delete()
    Kegg_pathway.objects.all().delete()
    Kegg_ortholog.objects.all().delete()
    Go_term.objects.all().delete()
    Cazy_family.objects.all().delete()
    Ec_number.objects.all().delete()
    Tc_family.objects.all().delete()
    logger.info("Deleting ortholog families...")
    Ortholog_group.objects.all().delete()
    # Delete genes before deleting proteins
    logger.info("Deleting genes...")
    Gene.objects.all().delete()
    # Delete proteins
    logger.info("Deleting proteins...")
    Protein.objects.all().delete()
    # Delete genomes
    logger.info("Deleting genomes...")
    Genome.objects.all().delete()
    # Delete strains
    logger.info("Deleting strains...")
    Strain_metadata.objects.all().delete()
    Strain.objects.all().delete()
    # Delete samples
    logger.info("Deleting samples...")
    Sample_metadata.objects.all().delete()
    Sample.objects.all().delete()
    # Delete taxonomy entries
    logger.info("Deleting taxonomy...")
    Taxon.objects.all().delete()
    # Delete tags
    Tag.objects.all().delete()
    logger.info("...done!")
    

def delete_all_genomes(confirm=True):
    '''
    This function deletes all genomes and associated data from the GenomeDepot database.
    It would be called from CLI only.
    User must explicitely confirm data deletion!
    This function doesn't delete static files.
    '''
    print('')
    while confirm:
        answer = input('All genomes, strains and samples in the GenomeDepot database '
            'will be deleted. This action cannot be reversed. Continue? (y/n)'
        )
        if answer.lower() in ["y","yes"]:
            break
        elif answer.lower() in ["n","no"]:
            print('Exiting.')
            return
        else:
            print('Please answer yes/y or no/n.')
    # Delete all mappings before deleting genes
    logger.info("Deleting annotations...")
    Annotation.objects.all().delete()
    # Delete genes before deleting proteins
    logger.info("Deleting genes...")
    Gene.objects.all().delete()
    # Delete proteins
    logger.info("Deleting proteins...")
    Protein.objects.all().delete()
    # Delete genomes
    logger.debug("Deleting genomes...")
    Genome.objects.all().delete()
    # Delete strains
    logger.info("Deleting strains...")
    Strain_metadata.objects.all().delete()
    Strain.objects.all().delete()
    # Delete samples
    logger.info("Deleting samples...")
    Sample_metadata.objects.all().delete()
    Sample.objects.all().delete()
    logger.info("...done.")


def delete_genome(genome_name):
    '''
        Deletes one genome with all contigs, genes and annotations 
        from the database.
        
        Also removes proteins not linked to any gene. But it does not remove a strain 
        associated with this genome if the strain has other genomes linked to.
    '''
    if genome_name == '':
        logger.error('Genome name required')
        return
    logger.info('Looking for genome ' + genome_name)
    genome_set = Genome.objects.filter(name=genome_name)
    if genome_set.count() == 0:
        logger.error('Genome ' + genome_name + ' not found')
        return
    elif genome_set.count() > 1:
        logger.error('Genome name is not unique: ' + genome_name)
        return
    importer = Importer()
    logger.info('Deleting genome...')
    genome_set.delete()
    logger.info('Deleting proteins not linked to genes...')
    Protein.objects.filter(gene=None).delete()
    logger.info('Deleting strains not linked to genomes...')
    Strain.objects.filter(genome=None).delete()
    logger.info('Deleting samples not linked to genomes...')
    Sample.objects.filter(genome=None).delete()
    logger.info('Re-creating search databases...')
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
    logger.info('Removing temporary files...')
    json_dir = os.path.join(importer.config['core.json_dir'], genome_name)
    if os.path.exists(json_dir):
        shutil.rmtree(json_dir)
    os.remove(importer.config['core.search_db_nucl'])
    os.remove(importer.config['core.search_db_prot'])
    # importer.cleanup()
    logger.debug('Done!')

def delete_genomes(genomes_file):
    '''
        Deletes one or more genomes listed in a text file
        with all genes and annotations from the database

        Also removes proteins not linked to any gene. But it does not remove strains
        associated with the genomes if a strain has other genomes linked to.
    '''
    if not os.path.exists(genomes_file):
        logger.error(genomes_file + ' not found')
        return
    logger.debug('Deleting genomes...')
    importer = Importer()
    with open(genomes_file, 'r') as infile:
        for line in infile:
            genome_name = line.rstrip('\n\r')
            if genome_name == '':
                continue
            logger.debug('Looking for genome' + genome_name)
            genome_set = Genome.objects.filter(name=genome_name)
            if genome_set.count() == 0:
                logger.warning('Genome ' + genome_name + ' not found. Skipped.')
                continue
            elif genome_set.count() > 1:
                logger.warning('Non-unique genome name: ' + genome_name + '. Skipped.')
                continue
            if os.path.exists(os.path.join(importer.config['core.json_dir'],
                                           genome_name
                                           )):
                shutil.rmtree(os.path.join(importer.config['core.json_dir'],
                                           genome_name
                                           ))
            genome_set.delete()
    logger.info('Deleting proteins not linked to genes...')
    Protein.objects.filter(gene=None).delete()
    logger.info('Deleting strains not linked to genomes...')
    Strain.objects.filter(genome=None).delete()
    logger.info('Deleting samples not linked to genomes...')
    Sample.objects.filter(genome=None).delete()
    logger.info('Re-creating search databases...')
    importer.export_proteins()
    importer.export_contigs()
    importer.delete_search_databases()
    shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                    os.path.basename(importer.config['core.search_db_nucl'])
                    ),
                    importer.config['core.search_db_nucl']
                    )
    shutil.copyfile(os.path.join(importer.config['core.temp_dir'],
                    os.path.basename(importer.config['core.search_db_prot'])),
                    importer.config['core.search_db_prot']
                    )
    importer.create_search_databases()
    os.remove(importer.config['core.search_db_nucl'])
    os.remove(importer.config['core.search_db_prot'])
    # importer.cleanup()
    logger.info('Done!')


def generate_static_files(genomes_file):
    '''
        For genomes listed in a text file, this function 
        generates static files
    '''
    if os.path.exists(genomes_file):
        importer = Importer()
        logger.debug('Genomes file ' + genomes_file)
        importer.generate_static_files(genomes_file)
    else:
        logger.error('Genomes file not found: ' + genomes_file)
    

def import_config(config_file, overwrite=False):
    '''Imports settings from a config file
    (text file with key/value entries separated by "=" symbol)
    '''
    if config_file is not None and os.path.exists(config_file):
        logger.debug('Importing parameters from' + config_file)
        configs = {}
        with open(config_file, 'r') as infile:
            for line in infile:
                key, val = line.rstrip('\n\r').split('=')
                configs[key.strip()] = val.strip()
        configs_saved = set()
        for item in Config.objects.all():
            if item.param in configs:
                configs_saved.add(item.param)
                if overwrite:
                    if item.value != configs[item.param]:
                        item.value = configs[item.param]
                        item.save()
        for param in configs:
            if param not in configs_saved:
                _ = Config.objects.create(param=param, value=configs[param])
    else:
        logger.error('Input file not found:' + config_file)


def export_config(out_file):
    if os.path.exists(out_file):
        logger.error('Output file already exists: ' + out_file)
        return
    with open(out_file, 'w') as outfile:
        for item in Config.objects.values_list('param', 'value'):
            outfile.write('='.join((str(x) for x in item)) + '\n')


def regenerate_jbrowse_files(genome_id):
    '''
    Deletes JBrowse files for a genome and generates new files
    '''
    # Check genome ID
    if genome_id == '':
        logger.error('Genome name required')
        return
    logger.debug('Looking for genome' + genome_id)
    genome_set = Genome.objects.filter(name=genome_id)
    if genome_set.count() == 0:
        logger.error('Genome ' + genome_id + ' not found')
        return
    elif genome_set.count() > 1:
        logger.error('Non-unique genome name: ' + genome_id)
        return
    logger.debug('Genome found:' + genome_id)
    genome = genome_set[0]
    # Configure importer
    importer = Importer()
    importer.inputgenomes[genome_id]['gbk'] = genome.gbk_filepath
    importer.inputgenomes[genome_id]['url'] = genome.external_url
    importer.inputgenomes[genome_id]['external_id'] = genome.external_id
    if genome.strain is None:
        importer.inputgenomes[genome_id]['strain'] = ''
    else:
        importer.inputgenomes[genome_id]['strain'] = genome.strain.strain_id
    if genome.sample is None:
        importer.inputgenomes[genome_id]['sample'] = ''
    else:
        importer.inputgenomes[genome_id]['sample'] = genome.sample.sample_id
    importer.export_jbrowse_genome_data(genome_id)


def recreate_search_databases():
    '''Deletes and re-creates nucleotide and protein search
    databases. Use the function if the files are missing or corrupted,
    or if the genome import pipeline crashed before creating the 
    search database files.
    '''
    importer = Importer()
    nucl_db_fasta = os.path.join(importer.config['core.temp_dir'],
        os.path.basename(importer.config['core.search_db_nucl'])
        )
    if os.path.exists(nucl_db_fasta):
        os.remove(nucl_db_fasta)
    with open(nucl_db_fasta, 'w') as outfile:
        for genome_data in Genome.objects.values_list('name', 'gbk_filepath'):
            if genome_data[1].endswith('.gz'):
                gbk_handle = gzip.open(genome_data[1], 'rt')
            else:
                gbk_handle = open(genome_data[1], 'r')
            parser = GenBank.parse(gbk_handle)
            for gbk_record in parser:
                contig_sequence = gbk_record.sequence
                contig_id = gbk_record.locus
                outfile.write('>' + contig_id + '|' +
                              genome_data[0] + '\n' +
                              ''.join(contig_sequence) + '\n'
                              )

    prot_db_fasta = os.path.join(importer.config['core.temp_dir'],
        os.path.basename(importer.config['core.search_db_prot'])
        )
    if os.path.exists(prot_db_fasta):
        os.remove(prot_db_fasta)
    export_proteins(None, prot_db_fasta)
    if os.path.exists(importer.config['core.search_db_nucl']):
        os.remove(importer.config['core.search_db_nucl'])
    if os.path.exists(importer.config['core.search_db_prot']):
        os.remove(importer.config['core.search_db_prot'])
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
    

def update_tags(genome_file, tag_names):
    '''
    Assigns one or more tags to genomes listed
    in a text or tab-separated file.
    '''
    if not os.path.exists(genome_file):
        logger.error('Genomes file ' + genome_file + ' not found.')
        return
    # create tags if missing
    tags = {}
    for tag_name in tag_names.split(','):
        tag_name = tag_name.strip()
        if tag_name == '':
            continue
        if not bool(re.match('[A-Za-z0-9\-\_]+$', 'tag_name')):
            raise CommandError(
                'Tag must contain only alphabet letters (a-z), ' +
                'numbers (0-9), hyphen or underscore. Correct the tag ' + tag_name
            )
        try:
            tag = Tag.objects.get(name=tag_name)
            tags[tag_name] = tag
        except Tag.DoesNotExist:
            tag = Tag(name = tag_name, description = '')
            tag.save()
            tags[tag_name] = tag
    # add tags to genome
    with open(genome_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n\r').split('\t')
            if len(row) == 1:
                genome_name = row[0]
            else:
                filepath, genome_name, _, _, _, _ = line.rstrip('\n\r').split('\t')
            try:
                genome = Genome.objects.get(name = genome_name)
                for tag_name, tag in tags.items():
                    genome.tags.add(tag)
                genome.save()
            except Genome.DoesNotExist:
                logger.warning('Genome not found:' + genome_name)
                continue
