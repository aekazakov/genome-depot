""" This file contains functions treemap generation"""
import logging
import plotly.express as px
import pandas as pd
from io import StringIO
from django.db.models import Count
from browser.util import autovivify
from browser.models import Gene
from browser.models import Protein
from browser.models import Annotation

logger = logging.getLogger("GenomeDepot")

def generate_og_treemap(ortholog_group):
    logger.debug('Generating treemap for ' + ortholog_group.eggnog_id)
    gene_ids = Gene.objects.filter(
        protein__ortholog_groups = ortholog_group
    #).select_related(
    #    'protein',
    #    'operon'
    ).values_list('id', flat=True)
    '''
    prefetch_related(
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
    )
    '''
    if not gene_ids:
        treemap = (
            '<h5>Genes for ' + ortholog_group.eggnog_id + '[' + 
            ortholog_group.taxon.name + '] not found in the database.</h5>'
        )
        return treemap, '', []
    #gene_ids = [item.id for item in genes.all()]
    treemap, result_table = generate_annotations_treemap(
        gene_ids, ortholog_group=ortholog_group
    )
    return treemap, result_table, list(gene_ids)


def generate_genes_treemap(gene_ids):
    '''
    genes = Gene.objects.filter(
        id__in = gene_ids
    ).select_related(
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
    )
    '''
    return generate_annotations_treemap(gene_ids)

    
def generate_annotations_treemap(gene_ids, ortholog_group=None):
    logger.debug('Filling genedata')
    genedata = autovivify(2, int)
    root_node = 'Annotations'
    annotations = Annotation.objects.filter(
        gene_id__in=gene_ids
    ).values_list(
        'source', 'key', 'value', 'note'
    )
    for item in annotations:
        label = item[3]
        if 'E-value' in label:
            label = label.split('E-value')[0]
        genedata[item[0]][(item[2], item[1] + ':' + label)] += 1

    for item in Gene.objects.filter(id__in=gene_ids).values_list('function', flat=True):
        genedata['Product'][(item, item)] += 1
    
    if ortholog_group is not None:
        for protein in Protein.objects.filter(
            ortholog_groups=ortholog_group, gene__in=gene_ids
        ).annotate(
            num_genes=Count('gene')
        ).prefetch_related(
            'kegg_orthologs',
            'kegg_pathways',
            'kegg_reactions',
            'ec_numbers',
            'tc_families',
            'cog_classes',
            'cazy_families'
        ):
            for item in protein.kegg_orthologs.all():
                genedata['KEGG Orthologs'][(item.kegg_id, item.description)] += \
                    protein.num_genes
            for item in protein.kegg_pathways.all():
                genedata['KEGG Pathways'][(item.kegg_id, item.description)] += \
                    protein.num_genes
            for item in protein.kegg_reactions.all():
                genedata['KEGG Reactions'][(item.kegg_id, item.description)] += \
                    protein.num_genes
            for item in protein.ec_numbers.all():
                genedata['EC'][(item.ec_number, item.description)] += \
                    protein.num_genes
            for item in protein.tc_families.all():
                genedata['TC'][(item.tc_id, item.description)] += \
                    protein.num_genes
            for item in protein.cog_classes.all():
                genedata['COG classes'][(item.cog_id, item.description)] += \
                    protein.num_genes
            for item in protein.cazy_families.all():
                genedata['CAZY'][(item.cazy_id, item.description)] += \
                    protein.num_genes
    else:
        for gene in Gene.objects.filter(
            id__in=gene_ids, protein__isnull=False
        ).prefetch_related(
            'protein__kegg_orthologs',
            'protein__kegg_pathways',
            'protein__kegg_reactions',
            'protein__ec_numbers',
            'protein__tc_families',
            'protein__cog_classes',
            'protein__cazy_families'
        ):
            for item in gene.protein.kegg_orthologs.all():
                genedata['KEGG Orthologs'][(item.kegg_id, item.description)] += 1
            for item in gene.protein.kegg_pathways.all():
                genedata['KEGG Pathways'][(item.kegg_id, item.description)] += 1
            for item in gene.protein.kegg_reactions.all():
                genedata['KEGG Reactions'][(item.kegg_id, item.description)] += 1
            for item in gene.protein.ec_numbers.all():
                genedata['EC'][(item.ec_number, item.description)] += 1
            for item in gene.protein.tc_families.all():
                genedata['TC'][(item.tc_id, item.description)] += 1
            for item in gene.protein.cog_classes.all():
                genedata['COG classes'][(item.cog_id, item.description)] += 1
            for item in gene.protein.cazy_families.all():
                genedata['CAZY'][(item.cazy_id, item.description)] += 1
    data = {'names':[root_node,], 'labels':[root_node,], 'parents':['',], 'values':[0,]}
    root_node = 'Annotations'
    result_table = []
    logger.debug('Building data structure')
    for category, category_data in genedata.items():
        data['names'].append(category)
        data['labels'].append(category)
        data['parents'].append(root_node)
        data['values'].append(0)
        for item, item_value in category_data.items():
            data['names'].append(item[0])
            data['labels'].append(item[1])
            data['parents'].append(category)
            data['values'].append(item_value)
            result_table.append((category, item[0], item[1], item_value))
    logger.debug('Building treemap from dataframe')
    fig = px.treemap(
        data_frame = pd.DataFrame.from_dict(data),
        names = 'names',
        parents = 'parents',
        values = 'values',
        hover_name = 'labels'
    )
    fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=0, l=0, r=0, b=0),
        uniformtext=dict(minsize=12,),
        autosize=True
    )
    logger.debug('Writing treemap HTML')
    html = StringIO()
    fig.write_html(html, include_plotlyjs=False, full_html=False)
    html = html.getvalue()
    result_table = sorted(result_table, key=lambda x: x[3], reverse=True)
    result_table = '\n'.join(
        ['\t'.join([str(item) for item in row]) for row in result_table]
    )
    result_table = 'Category\tFunction\tDescription\tGene count\n' + result_table
    return html, result_table

