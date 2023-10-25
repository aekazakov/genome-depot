""" This file contains functions treemap generation"""
#import plotly
import plotly.express as px
import pandas as pd
#import plotly.graph_objects as go
#import plotly.io as plotlyio
from io import StringIO
from collections import defaultdict
from browser.util import autovivify
from browser.models import Gene
from browser.models import Annotation


def generate_og_treemap(ortholog_group):
    print('Generating treemap for', ortholog_group.eggnog_id)
    genes = Gene.objects.filter(
        protein__ortholog_groups = ortholog_group
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
    if not genes:
        return '<h5>Genes for ' + ortholog_group.eggnog_id + '[' + ortholog_group.taxon.name + '] not found in the database.</h5>'
    gene_ids = [item.id for item in genes.all()]
    return generate_annotations_treemap(genes), gene_ids


def generate_genes_treemap(gene_ids):
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
    return generate_annotations_treemap(genes)

    
def generate_annotations_treemap(genes):
    genedata = autovivify(2, int)
    root_node = 'Annotations'
    root_value = 0

    for gene in genes:
        genedata['Function'][(gene.function, gene.function)] += 1
        for annotation in Annotation.objects.filter(gene_id=gene):
            label = annotation.note
            if 'E-value' in label:
                label = label.split('E-value')[0]
            genedata[annotation.source][(annotation.value, annotation.key + ':' + label)] += 1
        if not gene.protein:
            continue
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
    for category, category_data in genedata.items():
        data['names'].append(category)
        data['labels'].append(category)
        data['parents'].append(root_node)
        data['values'].append(0)
        category_value = 0
        for item, item_value in category_data.items():
            data['names'].append(item[0])
            data['labels'].append(item[1])
            data['parents'].append(category)
            data['values'].append(item_value)
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
    html = StringIO()
    fig.write_html(html, include_plotlyjs=False, full_html=False)
    fig.write_html("fig1.2.html", include_plotlyjs=True, full_html=True)
    html = html.getvalue()
    #print(html)
    return html

