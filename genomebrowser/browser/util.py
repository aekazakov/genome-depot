"""
    Various utility functions
"""
import os
import gzip
import logging
import urllib
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqFeature
from django.core.management.base import CommandError

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


logger = logging.getLogger("CGCMS")

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def export_genome(genome, output_buffer):
    '''
        Reads genome from GenBank file,
        adds gene mappings and annotations from CGCMS,
        writes genome in GenBank format into output buffer
    '''
                                        
    # Put genes here genedata[contig_id][gene_id] = list of notes
    genedata = autovivify(2, list)
    # Put operons here operondata[contig_id][gene_id] = 
    # (name, start, end, strand)
    operondata = defaultdict(dict)
    # Put sites here operondata[contig_id][gene_id] = 
    # list of (regulon, regulator, start, end, strand)
    sitedata = autovivify(2, list)
    for contig_id in list(Contig.objects.values_list('contig_id',
                                                 flat=True
                                                 ).filter(
                                                 genome=genome
                                                 )):
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
                    genedata[contig_id][gene.locus_tag]\
                        .append(
                            'EGGNOG5:'
                            + item.eggnog_id + '['
                            + item.taxon.name + '] [CGCMS/eggNOG-mapper]'
                        )
                for item in gene.protein.kegg_orthologs.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.kegg_pathways.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.kegg_reactions.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'KEGG:'
                        + item.kegg_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.ec_numbers.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'EC:'
                        + item.ec_number + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.go_terms.all():
                    genedata[contig_id][gene.locus_tag].append(
                        ''
                        + item.go_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.tc_families.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'TC:'
                        + item.tc_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.cog_classes.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'COG:'
                        + item.cog_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
                for item in gene.protein.cazy_families.all():
                    genedata[contig_id][gene.locus_tag].append(
                        'CAZy:'
                        + item.cazy_id + ' ('
                        + item.description
                        + ') [CGCMS/eggNOG-mapper]'
                    )
            for annotation in Annotation.objects.filter(gene_id=gene):
                genedata[contig_id][gene.locus_tag].append(
                    annotation.source + ':' + annotation.value 
                    + ' (' + annotation.note + ') [CGCMS/'
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
        logger.info(genome_id, 'exported')
    

def download_ncbi_assembly(assembly_id, email, upload_dir):
    """Download genbank assemblies for a given search term.
    Args:
        assembly_id: search term, usually assembly RefSeq of GenBank accession
        email: email address, requierd for NCBI requests
        upload_dir: folder to save to
    """
    if not os.path.exists(upload_dir):
        os.mkdir(upload_dir)
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

    
def delete_all_data():
    # Delete all mappings before deleting genes
    Annotation.objects.all().delete()
    Eggnog_description.objects.all().delete()
    Cog_class.objects.all().delete()
    Kegg_reaction.objects.all().delete()
    Kegg_pathway.objects.all().delete()
    Kegg_ortholog.objects.all().delete()
    Go_term.objects.all().delete()
    Cazy_family.objects.all().delete()
    Ec_number.objects.all().delete()
    Tc_family.objects.all().delete()
    Ortholog_group.objects.all().delete()
    # Delete genes before deleting proteins
    Gene.objects.all().delete()
    # Delete proteins
    Protein.objects.all().delete()
    # Delete genomes
    Genome.objects.all().delete()
    # Delete strains
    Strain_metadata.objects.all().delete()
    Strain.objects.all().delete()
    # Delete samples
    Sample_metadata.objects.all().delete()
    Sample.objects.all().delete()
    # Delete taxonomy entries
    Taxon.objects.all().delete()
    # Delete tags
    Tag.objects.all().delete()
    

def delete_all_genomes():
    # Delete all mappings before deleting genes
    print("Deleting annotations...")
    Annotation.objects.all().delete()
    # Delete genes before deleting proteins
    print("Deleting genes...")
    Gene.objects.all().delete()
    # Delete proteins
    print("Deleting proteins...")
    Protein.objects.all().delete()
    # Delete genomes
    print("Deleting genomes...")
    Genome.objects.all().delete()
    # Delete strains
    print("Deleting strains...")
    Strain_metadata.objects.all().delete()
    Strain.objects.all().delete()
    # Delete samples
    print("Deleting samples...")
    Sample_metadata.objects.all().delete()
    Sample.objects.all().delete()
    print("...done.")


def delete_genome(genome_name):
    if genome_name == '':
        raise CommandError('Genome name required')
    print('Looking for genome', genome_name)
    genome_set = Genome.objects.filter(name=genome_name)
    if genome_set.count() == 0:
        print('Genome ' + genome_name + ' not found')
        raise CommandError()
    elif genome_set.count() > 1:
        print('Not unique genome name: ' + genome_name)
        raise CommandError()
    importer = Importer()
    print('''Note 1: Genome deletion removes contigs, genes and 
    gene annotations for this genome. It also removes proteins 
    not linked to any gene. But it does not remove a strain 
    associated with this genome if the strain has other genomes.''')
    print('Deleting genome...')
    genome_set.delete()
    print('Deleting proteins not linked to genes...')
    Protein.objects.filter(gene=None).delete()
    print('Deleting strains not linked to genomes...')
    Strain.objects.filter(genome=None).delete()
    print('Deleting samples not linked to genomes...')
    Sample.objects.filter(genome=None).delete()
    importer.export_proteins()
    importer.export_contigs()
    importer.delete_search_databases()
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                    os.path.basename(importer.config['cgcms.search_db_nucl'])),
                    importer.config['cgcms.search_db_nucl']
                    )
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'],
                    os.path.basename(importer.config['cgcms.search_db_prot'])),
                    importer.config['cgcms.search_db_prot']
                    )
    importer.create_search_databases()
    shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'],
                               genome_name
                               ))
    os.remove(importer.config['cgcms.search_db_nucl'])
    os.remove(importer.config['cgcms.search_db_prot'])
    # importer.cleanup()
    print('Done!')
