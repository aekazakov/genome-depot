import binascii
import csv
import logging
import gzip
from io import BytesIO

from django.http import HttpResponse
from django.core.exceptions import SuspiciousOperation
from django.db.models import Q
from browser.models import Annotation
from browser.models import Gene
from browser.models import Genome
from browser.models import Protein
from browser.models import Kegg_ortholog
from browser.models import Kegg_pathway
from browser.models import Kegg_reaction
from browser.models import Ec_number
from browser.models import Tc_family
from browser.models import Cog_class
from browser.models import Cazy_family
from browser.models import Go_term
from browser.util import export_genome

logger = logging.getLogger("CGCMS")


def export_csv(request):
    '''
        Main gateway function for CSV export.
    '''
    query_type = request.GET.get('type')
    if query_type == 'genome':
        return _export_genomes_csv(request)
    elif query_type == 'annotation':
        return _export_annotations_csv(request)
    else:
        return _export_genes_csv(request)
    
    
def _export_annotations_csv(request):
    '''
        Returns list of annotations in tab-separated text format
        
        Template gene_list_subpage.html contains the Ajax JS calling this function 
        
    '''

    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="export.tab"'
    writer = csv.writer(response, delimiter='\t')
    search_context = ('','')
    annotation_query = request.GET.get('annotation_query')
    genome = request.GET.get('genome')
    search_context = ('Gene annotation query', annotation_query)
    # get objects
    if genome:
        object_list = Annotation.objects.filter(
            (
                Q(source__icontains=annotation_query) |
                Q(value__icontains=annotation_query) |
                Q(note__icontains=annotation_query)
            ) &
            Q(gene_id__genome__name=genome)
        ).order_by(
            'gene_id__locus_tag'
        ).select_related(
            'gene_id', 'gene_id__genome', 'gene_id__genome__taxon'
        )
    else:
        object_list = Annotation.objects.filter(
            Q(source__icontains=annotation_query) |
            Q(value__icontains=annotation_query) |
            Q(note__icontains=annotation_query)
        ).order_by(
            'gene_id__locus_tag'
        ).select_related(
            'gene_id', 'gene_id__genome', 'gene_id__genome__taxon'
        )
    # write output
    writer.writerow(['Locus tag',
                     'Name',
                     'Organism',
                     'Genome',
                     'Contig',
                     'Start',
                     'End',
                     'Strand',
                     'Type',
                     'Function',
                     'Annotation_source',
                     'Annotation_type',
                     'Annotation_note'
                     ])
    for obj in object_list:
        writer.writerow([obj.gene_id.locus_tag,
                         obj.gene_id.name,
                         obj.gene_id.genome.taxon.name,
                         obj.gene_id.genome.name,
                         obj.gene_id.contig.contig_id,
                         str(obj.gene_id.start),
                         str(obj.gene_id.end),
                         str(obj.gene_id.strand),
                         obj.gene_id.type,
                         obj.gene_id.function,
                         obj.source,
                         obj.key,
                         obj.value,
                         obj.note
                         ])
    return response
    

def _export_genes_csv(request):
    '''
        Returns list of genes in tab-separated text format
        
        Template gene_list_subpage.html contains the Ajax JS calling this function 
        
    '''
    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="export_genes.tab"'
    writer = csv.writer(response, delimiter='\t')
    search_context = ('','')
    annotation_query = request.GET.get('annotation_query')
    query_type = request.GET.get('type')
    query = request.GET.get('query')
    genome = request.GET.get('genome')

    if query_type == 'gene':
        if query:
            search_context = ('Query', query)
            if genome:
                object_list = Gene.objects.filter(
                    genome__name=genome
                ).filter(
                    Q(name__icontains=query) |
                    Q(locus_tag__icontains=query) |
                    Q(function__icontains=query)
                ).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) |
                    Q(locus_tag__icontains=query)
                ).order_by('locus_tag')
        elif genome:
            # Generate gene list with all mappings and annotations
            search_context = ('Genome', genome)
            object_list = Gene.objects.filter(
                genome__name=genome
            ).order_by(
                'locus_tag'
            ).select_related(
                'genome',
                'contig',
                'protein',
                'protein__eggnog_description'
            ).prefetch_related(
                'protein__kegg_orthologs',
                'protein__kegg_reactions',
                'protein__kegg_pathways',
                'protein__go_terms',
                'protein__cazy_families',
                'protein__ec_numbers',
                'protein__tc_families',
                'protein__ortholog_groups',
                'protein__cog_classes',
            )
            writer.writerow([
                'Locus tag',
                'Name',
                'Organism',
                'Genome',
                'Contig',
                'Start',
                'End',
                'Strand',
                'Type',
                'Function',
                'EggNOG description'
                'EggNOG families'
                'KEGG Orthologs',
                'KEGG Pathways',
                'KEGG Reactions',
                'GO Terms',
                'CAZy Families',
                'EC',
                'TC',
                'COG classes',
                'Annotations'
            ])
            for gene in object_list:
                eggnog_description = ''
                eggnogs = ''
                ko = ''
                kp = ''
                kr = ''
                go = ''
                cazy = ''
                ec = ''
                tc = ''
                cog = ''
                if gene.protein:
                    if gene.protein.eggnog_description:
                        eggnog_description = gene.protein.eggnog_description.description
                    if gene.protein.ortholog_groups:
                        eggnogs = ';'.join([item.eggnog_id + '[' + item.taxon.name + ']' for item in gene.protein.ortholog_groups.all()])
                    if gene.protein.kegg_orthologs:
                        ko = ';'.join([item.kegg_id for item in gene.protein.kegg_orthologs.all()])
                    if gene.protein.kegg_pathways:
                        kp = ';'.join([item.kegg_id for item in gene.protein.kegg_pathways.all()])
                    if gene.protein.kegg_reactions:
                        kr = ';'.join([item.kegg_id for item in gene.protein.kegg_reactions.all()])
                    if gene.protein.go_terms:
                        go = ';'.join([item.go_id for item in gene.protein.go_terms.all()])
                    if gene.protein.cazy_families:
                        cazy = ';'.join([item.cazy_id for item in gene.protein.cazy_families.all()])
                    if gene.protein.ec_numbers:
                        ec = ';'.join([item.ec_number for item in gene.protein.ec_numbers.all()])
                    if gene.protein.tc_families:
                        tc = ';'.join([item.tc_id for item in gene.protein.tc_families.all()])
                    if gene.protein.cog_classes:
                        cog = ';'.join([item.cog_id for item in gene.protein.cog_classes.all()])
                annotations = ';'.join([item.key + '=' + item.value + '[' + item.source + ']' for item in gene.annotation_set.all()])

                writer.writerow([
                    gene.locus_tag,
                    gene.name,
                    gene.genome.taxon.name,
                    gene.genome.name,
                    gene.contig.contig_id,
                    str(gene.start),
                    str(gene.end),
                    str(gene.strand),
                    gene.type,
                    gene.function,
                    eggnog_description,
                    eggnogs,
                    ko,
                    kp,
                    kr,
                    go,
                    cazy,
                    ec,
                    tc,
                    cog,
                    annotations
                ])
            return response
            
        else:
            object_list = Gene.objects.none()
    else:
        if query_type == 'og':
            search_context = ('Ortholog group', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        ortholog_groups__id=query).values('protein_hash')
                        ]
        elif query_type == 'ko_id':
            search_context = ('KEGG ortholog', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_orthologs__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'kp_id':
            search_context = ('KEGG pathway', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_pathways__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'kr_id':
            search_context = ('KEGG reaction', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_reactions__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'ec_id':
            search_context = ('EC number', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        ec_numbers__ec_number=query).values('protein_hash')
                        ]
        elif query_type == 'tc_id':
            search_context = ('TCDB family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        tc_families__tc_id=query).values('protein_hash')
                        ]
        elif query_type == 'cazy_id':
            search_context = ('CAZy family', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cazy_families__cazy_id=query).values('protein_hash')
                        ]
        elif query_type == 'cog_id':
            search_context = ('COG class', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cog_classes__cog_id=query).values('protein_hash')
                        ]
        elif query_type == 'go_id':
            search_context = ('GO term', query)
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        go_terms__go_id=query).values('protein_hash')
                        ]
        elif query_type == 'ko':
            search_context = ('KEGG ortholog query', query)
            ko_ids = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')
                        ]
        elif query_type == 'kp':
            search_context = ('KEGG pathway query', query)
            kp_ids = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')
                        ]
        elif query_type == 'kr':
            search_context = ('KEGG reaction query', query)
            kr_ids = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')
                        ]
        elif query_type == 'ec':
            search_context = ('EC number query', query)
            ec_ids = Ec_number.objects.filter(
                Q(ec_number__icontains=query) |
                Q(description__icontains=query)
            ).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                ec_numbers__ec_number__in=ec_ids
            ).values('protein_hash')]
        elif query_type == 'tc':
            search_context = ('TCDB family query', query)
            tc_ids = Tc_family.objects.filter(
                Q(tc_id__icontains=query) |
                Q(description__icontains=query)
            ).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        tc_families__tc_id__in=tc_ids).values('protein_hash')
                        ]
        elif query_type == 'cazy':
            search_context = ('CAZy family query', query)
            cazy_ids = Cazy_family.objects.filter(
                Q(cazy_id__icontains=query) |
                Q(description__icontains=query)
            ).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cazy_families__cazy_id__in=cazy_ids).values('protein_hash')
                        ]
        elif query_type == 'cog':
            search_context = ('COG class query', query)
            cog_ids = Cog_class.objects.filter(
                Q(cog_id__icontains=query) |
                Q(description__icontains=query)
            ).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cog_classes__cog_id__in=cog_ids).values('protein_hash')
                        ]
        elif query_type == 'go':
            search_context = ('GO term query', query)
            go_ids = Go_term.objects.filter(
                Q(go_id__icontains=query) |
                Q(description__icontains=query)
            ).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        go_terms__go_id__in=go_ids).values('protein_hash')
                        ]
        else:
            search_context = ('Unknown', query)
            proteins = []

        if genome:
            if not query_type:
                object_list = Gene.objects.filter(
                    genome__name=genome
                ).order_by('locus_tag')
            else:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by('locus_tag')
        else:
            if not query_type:
                object_list = Gene.objects.none()
            else:
                object_list = Gene.objects.filter(
                    protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'contig'
                )
    if search_context[0] == 'Genome':
        writer.writerow(['Locus tag',
                         'Name',
                         'Organism',
                         'Genome',
                         'Contig',
                         'Start',
                         'End',
                         'Strand',
                         'Type',
                         'Function'
                         ])
        for gene in object_list:
            writer.writerow([gene.locus_tag,
                             gene.name,
                             gene.genome.taxon.name,
                             gene.genome.name,
                             gene.contig.contig_id,
                             str(gene.start),
                             str(gene.end),
                             str(gene.strand),
                             gene.type,
                             gene.function])
    else:
        writer.writerow(['Locus tag',
                         'Name',
                         'Organism',
                         'Genome',
                         'Contig',
                         'Start',
                         'End',
                         'Strand',
                         'Type',
                         'Function',
                         search_context[0]
                         ])
        for gene in object_list:
            writer.writerow([gene.locus_tag,
                             gene.name,
                             gene.genome.taxon.name,
                             gene.genome.name,
                             gene.contig.contig_id,
                             str(gene.start),
                             str(gene.end),
                             str(gene.strand),
                             gene.type,
                             gene.function,
                             search_context[1]
                             ])

    return response


def _export_genomes_csv(request):
    '''
        Returns list of genes in tab-separated text format
        
        Template gene_list_subpage.html contains the Ajax JS calling this function 
        
    '''
    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="genomes.tab"'
    writer = csv.writer(response, delimiter='\t')
    search_context = ('','')
    annotation_query = request.GET.get('annotation_query')
    query_type = request.GET.get('type')
    query = request.GET.get('query')
    if query:
        object_list = Genome.objects.filter(
            name__icontains=query
        ).order_by(
            'name'
        ).select_related(
            'strain', 'sample', 'taxon'
        ).prefetch_related(
            'tags'
        )
    else:
        object_list = Genome.objects.order_by(
            'name'
        ).select_related(
            'strain', 'sample', 'taxon'
        ).prefetch_related(
            'tags'
        )
    writer.writerow([
                     'Name',
                     'Tags',
                     'Source',
                     'Taxonomy',
                     'Size, bp',
                     'Contigs',
                     'Genes'
                     ])
    for genome in object_list:
        if genome.strain:
            source = genome.strain.full_name
        else:
            source = genome.sample.full_name
        writer.writerow([genome.name,
                         ';'.join([item.name for item in genome.tags.all()]),
                         source,
                         genome.taxon.name,
                         str(genome.size),
                         str(genome.contigs),
                         str(genome.genes)
                         ])
    return response


def export_fasta(request):
    '''
        Returns protein sequences in FASTA format
    '''
    response = HttpResponse(content_type='text/plain')
    query_type = request.GET.get('type')
    query = request.GET.get('query')
    genome = request.GET.get('genome')
    if genome:
        response['Content-Disposition'] = 'attachment; filename="' + str(genome) + '_proteins.faa"'
    else:
        response['Content-Disposition'] = 'attachment; filename="exported_proteins.faa"'
    annotation_query = request.GET.get('annotation_query')

    if annotation_query:
        if genome:
            object_list = Annotation.objects.filter(
                (
                    Q(source__icontains=annotation_query) |
                    Q(value__icontains=annotation_query) |
                    Q(note__icontains=annotation_query)
                ) &
                Q(gene_id__genome__name=genome)
            ).order_by(
                'gene_id__locus_tag'
            ).select_related(
                'gene_id',
                'gene_id__genome',
                'gene_id__genome__taxon',
                'gene_id__protein'
            )
        else:
            object_list = Annotation.objects.filter(
                Q(source__icontains=annotation_query) |
                Q(value__icontains=annotation_query) |
                Q(note__icontains=annotation_query)
            ).order_by(
                'gene_id__locus_tag'
            ).select_related(
                'gene_id',
                'gene_id__genome',
                'gene_id__genome__taxon',
                'gene_id__protein'
            )
    elif query_type == 'gene':
        if query:
            if genome:
                object_list = Gene.objects.filter(genome__name=genome).filter(
                    Q(name__icontains=query) |
                    Q(locus_tag__icontains=query) |
                    Q(function__icontains=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                )
            else:
                object_list = Gene.objects.filter(
                    Q(name__icontains=query) |
                    Q(locus_tag__icontains=query)
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                )
        elif genome:
            object_list = Gene.objects.filter(
                genome__name=genome
            ).order_by(
                'locus_tag'
            ).select_related(
                'genome', 'genome__taxon', 'protein'
            )
        else:
            object_list = Gene.objects.none()
    else:
        if query_type == 'og':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        ortholog_groups__id=query).values('protein_hash')
                        ]
        elif query_type == 'ko_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_orthologs__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'kp_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_pathways__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'kr_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_reactions__kegg_id=query).values('protein_hash')
                        ]
        elif query_type == 'ec_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        ec_numbers__ec_number=query).values('protein_hash')
                        ]
        elif query_type == 'tc_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        tc_families__tc_id=query).values('protein_hash')
                        ]
        elif query_type == 'cazy_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cazy_families__cazy_id=query).values('protein_hash')
                        ]
        elif query_type == 'cog_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cog_classes__cog_id=query).values('protein_hash')
                        ]
        elif query_type == 'go_id':
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        go_terms__go_id=query).values('protein_hash')
                        ]
        elif query_type == 'ko':
            ko_ids = Kegg_ortholog.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_orthologs__kegg_id__in=ko_ids).values('protein_hash')
                        ]
        elif query_type == 'kp':
            kp_ids = Kegg_pathway.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_pathways__kegg_id__in=kp_ids).values('protein_hash')
                        ]
        elif query_type == 'kr':
            kr_ids = Kegg_reaction.objects.filter(
                Q(kegg_id__icontains=query) |
                Q(description__icontains=query)
            ).values('kegg_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        kegg_reactions__kegg_id__in=kr_ids).values('protein_hash')
                        ]
        elif query_type == 'ec':
            ec_ids = Ec_number.objects.filter(
                Q(ec_number__icontains=query) |
                Q(description__icontains=query)
            ).values('ec_number')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        ec_numbers__ec_number__in=ec_ids).values('protein_hash')]
        elif query_type == 'tc':
            tc_ids = Tc_family.objects.filter(
                Q(tc_id__icontains=query) |
                Q(description__icontains=query)
            ).values('tc_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        tc_families__tc_id__in=tc_ids).values('protein_hash')
                        ]
        elif query_type == 'cazy':
            cazy_ids = Cazy_family.objects.filter(
                Q(cazy_id__icontains=query) |
                Q(description__icontains=query)
            ).values('cazy_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cazy_families__cazy_id__in=cazy_ids).values('protein_hash')
                        ]
        elif query_type == 'cog':
            cog_ids = Cog_class.objects.filter(
                Q(cog_id__icontains=query) |
                Q(description__icontains=query)
            ).values('cog_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        cog_classes__cog_id__in=cog_ids).values('protein_hash')
                        ]
        elif query_type == 'go':
            go_ids = Go_term.objects.filter(
                Q(go_id__icontains=query) |
                Q(description__icontains=query)
            ).values('go_id')
            proteins = [item['protein_hash'] for item in Protein.objects.filter(
                        go_terms__go_id__in=go_ids).values('protein_hash')
                        ]
        else:
            proteins = []

        if genome:
            if not query_type:
                object_list = Gene.objects.filter(
                    genome__name=genome
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                )
            else:
                object_list = Gene.objects.filter(
                    genome__name=genome, protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                )
        else:
            if not query_type:
                object_list = Gene.objects.none()
            else:
                object_list = Gene.objects.filter(
                    protein__protein_hash__in=proteins
                ).order_by(
                    'locus_tag'
                ).select_related(
                    'genome', 'genome__taxon', 'protein'
                )

    if annotation_query:
        for obj in object_list:
            if obj.gene_id.protein:
                response.write('>' + obj.gene_id.locus_tag + '|' + \
                               obj.gene_id.genome.name + ' [' + \
                               obj.gene_id.genome.taxon.name + ']\n' + \
                               obj.gene_id.protein.sequence + '\n'
                               )
    else:
        for gene in object_list:
            if gene.protein:
                response.write('>' + gene.locus_tag + '|' + gene.genome.name + \
                               ' [' + gene.genome.taxon.name + ']\n' + \
                               gene.protein.sequence + '\n'
                               )

    return response
    

def export_gbk(request, name):
    response = HttpResponse(content_type='application/x-gzip')
    response['Content-Disposition'] = 'attachment; filename=exported_' + name + '_genome.gbk'
    response['Content-Encoding'] = 'gzip'
    try:
        genome = Genome.objects.get(name = name)
    except Genome.DoesNotExist:
        logger.error('Genome not found: ' + str(name))
        raise
    gzip_buffer = BytesIO()
    with gzip.open(gzip_buffer, 'wt') as outfile:
        export_genome(genome, outfile)
    response.write(gzip_buffer.getvalue())
    return response
