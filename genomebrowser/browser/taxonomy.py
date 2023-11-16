""" This file contains functions for search in taxonomy"""
import plotly.graph_objects as go
from io import StringIO
from collections import defaultdict, Counter
from django.urls import reverse
from browser.models import Gene
from browser.models import Genome
from browser.models import Operon
from browser.models import Taxon

COLORS = ['#FFAA00',
          '#2D5F91',
          '#819FBD',
          '#819FBD',
          '#91D4D2',
          '#96BEE6',
          '#C0D8F0',
          '#E8655F',
          '#F1A39F',
          '#48B7B4'
          ]

def generate_genome_sunburst(taxon_id=None, children = []):
    '''
    If taxon_id is None, it will choose root node as 
    the highest node that have more than one child.
    '''

    if children:
        labels, parents, values, customdata = get_genomes_taxonomy(taxon_id, children)
    else:
        labels, parents, values, customdata = get_genomes_taxonomy(taxon_id)
        
    if not labels:
        return ''
    
    maxdepth = -1
    if taxon_id is None:
        taxon_id = get_root_id(labels, parents)
        maxdepth = 4
    
    colors = []
    for ind, label in enumerate(labels):
        colors.append(COLORS[ind%10])
    
    fig = go.Figure(go.Sunburst(
        ids = labels,
        labels = labels,
        parents = parents,
        values = values,
        customdata = customdata,
        marker=dict(colors=colors),
        texttemplate = '<a href="%{customdata[0]}" style="color:black;">%{label}</a>',
        hovertemplate = '%{customdata[1]}<extra></extra>',
        branchvalues = 'total',
        level = '',
        maxdepth = maxdepth,
    ))
    fig.update_layout(margin = dict(t = 0, l = 0, r = 0, b = 0),
                      height=800,
                      uniformtext=dict(minsize=10, mode='hide'),
                      paper_bgcolor='rgba(120,120,120,0.1)',
                      plot_bgcolor='rgba(120,120,120,0.1)'
                      )
    html = StringIO()
    fig.write_html(html, include_plotlyjs=False, full_html=False)
    html = html.getvalue()
    return html


def generate_genes_sunburst(gene_ids):
    labels, parents, values, customdata = get_genes_taxonomy(gene_ids)
    
    if not labels:
        return ''
    
    taxon_id = get_root_id(labels, parents)
    maxdepth = -1
    
    colors = []
    for ind, label in enumerate(labels):
        colors.append(COLORS[ind%10])
    
    fig = go.Figure(go.Sunburst(
        ids = labels,
        labels = labels,
        parents = parents,
        values = values,
        customdata = customdata,
        marker=dict(colors=colors),
        texttemplate = '<a href="%{customdata[0]}" style="color:black;">%{label}</a>',
        hovertemplate = '%{customdata[1]}<extra></extra>',
        branchvalues = 'total',
        level = '',
        maxdepth = maxdepth,
    ))
    fig.update_layout(margin = dict(t = 0, l = 0, r = 0, b = 0),
                      height=800,
                      uniformtext=dict(minsize=10, mode='hide'),
                      paper_bgcolor='rgba(120,120,120,0.1)',
                      plot_bgcolor='rgba(120,120,120,0.1)'
                      )
    html = StringIO()
    fig.write_html(html, include_plotlyjs=False, full_html=False)
    html = html.getvalue()
    return html

    
def generate_operons_sunburst(operon_ids, gene_ids):
    '''
        Generates sunburst plot for selected operons and 
        MONOCISTRONIC (not assigned to operon) genes 
    '''
    labels, parents, values, customdata = get_operons_taxonomy(operon_ids, gene_ids)
    print(labels)
    print(parents)
    print(customdata)
    
    if not labels:
        return ''
    
    taxon_id = get_root_id(labels, parents)
    maxdepth = -1
    
    colors = []
    for ind, label in enumerate(labels):
        colors.append(COLORS[ind%10])
    
    fig = go.Figure(go.Sunburst(
        ids = labels,
        labels = labels,
        parents = parents,
        values = values,
        customdata = customdata,
        marker=dict(colors=colors),
        texttemplate = '<a href="%{customdata[0]}" style="color:black;">%{label}</a>',
        hovertemplate = '%{customdata[1]}<extra></extra>',
        branchvalues = 'total',
        level = '',
        maxdepth = maxdepth,
    ))
    fig.update_layout(margin = dict(t = 0, l = 0, r = 0, b = 0),
                      height=800,
                      uniformtext=dict(minsize=10, mode='hide'),
                      paper_bgcolor='rgba(120,120,120,0.1)',
                      plot_bgcolor='rgba(120,120,120,0.1)'
                      )
    html = StringIO()
    fig.write_html(html, include_plotlyjs=False, full_html=False)
    html = html.getvalue()
    return html

    
def get_taxon_children(taxonomy_id):
    '''
    Returns list of taxonomy ids in a subtree of a given taxon
    '''
    result = []
    try:
        _ = Taxon.objects.get(taxonomy_id=taxonomy_id)
    except Taxon.DoesNotExist:
        return result
    result.append(taxonomy_id)
    parent_ids = [taxonomy_id]
    while True:
        children_ids = list(Taxon.objects.values_list(
                                                      'taxonomy_id',
                                                      flat=True
                                                      ).filter(
                                                      parent_id__in=parent_ids
                                                      ))
        if len(children_ids) == 0:
            break
        result += children_ids
        parent_ids = children_ids
    return result

    
def get_genomes_taxonomy(target_taxon_id, target_children = []):
    '''
    Returns four lists for sunburst graph generation: list of genomes and taxa names,
    list of parent taxa names, list of genome counts for each taxon, list of 
    HTML-formatted strings containing a link to the current element
    
    Parameters:
    target_taxon_id (str or None): taxonomy ID of root element. If target_taxon_id
    is None, place root at the highest level that have more than one child.

    target_children (list of str): list of taxonomy IDs that must be shown on 
    the sunburst graph
    
    Steps:
    1. Get list of all genomes with all their taxonomy IDs
    2. Get list of all Taxon objects with their parent ID
    3. Calculate number of genomes
    
    '''
    labels = []
    parents = []
    values = []
    customdata = []
    parent_lookup = {}
    # key: taxonomy id, value: parent id
        
    taxon_counts = Counter()
    # key: name, value: count of children genomes
    children = defaultdict(dict)
    taxon_lookup = {}
    move_root_down = False
    
    if target_taxon_id is None:
        genomes = Genome.objects.values_list('name',
                                             'taxon__name',
                                             'taxon__taxonomy_id'
                                             )
        for item in Taxon.objects.values_list(
            'name',
            'taxonomy_id',
            'parent_id',
            'rank'
        ):
            taxon_lookup[item[1]] = item
        children_taxonomy_ids = list(taxon_lookup.keys())
        move_root_down = True
        # Set target_taxon_id to '1' (root node), and move it down later
        target_taxon_id = '1'
    else:
        if target_children:
            children_taxonomy_ids = target_children
        else:
            children_taxonomy_ids = get_taxon_children(target_taxon_id)
        genomes = Genome.objects.values_list(
            'name',
            'taxon__name',
            'taxon__taxonomy_id'
         ).filter(
            taxon__taxonomy_id__in=children_taxonomy_ids
         )
        for item in Taxon.objects.values_list(
            'name',
            'taxonomy_id',
            'parent_id',
            'rank'
        ).filter(
            taxonomy_id__in=children_taxonomy_ids
        ):
            taxon_lookup[item[1]] = item

    for genome in genomes:
        # Add genome to graph
        labels.append(genome[0])
        parents.append(genome[1])
        values.append(1)
        customdata.append([reverse('genomedetails', args=(genome[0],)),
                           genome[0] + ' genome'
                           ])
        # Add taxon
        genome_taxid = genome[2]
        parent_id = taxon_lookup[genome_taxid][2]
        if parent_id in children_taxonomy_ids:
            children[parent_id][genome_taxid] = 1
        taxon_counts[genome_taxid] += 1
        parent_lookup[genome_taxid] = parent_id
        # Add all parents of the genome up to the root or target_taxon_id
        while True:
            if parent_id in children_taxonomy_ids:
                children[parent_id][genome_taxid] = 1
                taxon_counts[parent_id] += 1
            if parent_id == target_taxon_id \
            or parent_id == '1' \
            or parent_id not in taxon_lookup:
                break
            parent_lookup[parent_id] = taxon_lookup[parent_id][2]
            if parent_id != taxon_lookup[parent_id][2]:
                children[taxon_lookup[parent_id][2]][parent_id] = 1
            parent_id = taxon_lookup[parent_id][2]
    
    if len(taxon_counts) == 0:
        # No taxa, nothing to return
        return [], [], [], []
        
    # It's time to move root node
    if move_root_down:
        remove_nodes= set()
        new_root_id = target_taxon_id
        '''
        Going down from current root node, find taxon that 
        has more than one child and call it new root node
        '''
        while True:
            # Count children that have non-zero genome count, excluding root taxon
            children_nodes = [item for item in taxon_lookup.values() 
                              if item[2] == new_root_id 
                              and item[1] in taxon_counts
                              and item[1] != item[2]
                              ]
            if len(children_nodes) == 1:
                remove_nodes.add(new_root_id)
                new_root_id = children_nodes[0][1]
            else:
                target_taxon_id = new_root_id
                break
        for taxon in remove_nodes:
            del(taxon_counts[taxon])
        
    for taxon in taxon_counts:
        # If the root taxon with taxonomy_id '1' is in the list, skip it
        if taxon == '1':
            continue
        labels.append(taxon_lookup[taxon][0])
        if taxon_counts[taxon] == 1:
            customdata.append([reverse('taxondetails', args=(taxon,)),
                           '<br>' + taxon_lookup[taxon][0] + ' [' + \
                           taxon_lookup[taxon][3] + ']<br>' + \
                           str(taxon_counts[taxon]) + ' genome'
                           ])
        else:
            customdata.append([reverse('taxondetails', args=(taxon,)),
                           '<br>' + taxon_lookup[taxon][0] + ' [' + \
                           taxon_lookup[taxon][3] + ']<br>' + \
                           str(taxon_counts[taxon]) + ' genomes'
                           ])
        if taxon == target_taxon_id or taxon_lookup[taxon][2] == '1':
            # Parent of the root node must be empty string
            parents.append('')
        else:
            parents.append(taxon_lookup[parent_lookup[taxon]][0])
        values.append(taxon_counts[taxon])
    
    return labels, parents, values, customdata

    
def get_genes_taxonomy(gene_ids = []):
    '''
    Returns four lists for sunburst graph generation: list of gene/genome/taxon names,
    list of parent genome/taxon names, list of gene counts for each genome and taxon, list of 
    HTML-formatted strings containing a link to the current element
    
    Parameters:
    gene_ids (list of int): list of numerical gene ids
    
    Steps:
    1. Get list of all genes and associated genomes
    2. Get list of all genomes with all their taxonomy IDs
    3. Get list of all Taxon objects with their parent ID
    4. Calculate number of genomes
    
    '''
    labels = []
    parents = []
    values = []
    customdata = []
    parent_lookup = {}
    # key: genome name or taxonomy id, value: parent id
        
    taxon_counts = Counter()
    # key: name, value: count of children genes
    children = defaultdict(dict)
    taxon_lookup = {}
    move_root_down = True
    
    genes = Gene.objects.filter(
        id__in=gene_ids
    ).values_list(
        'locus_tag',
        'genome__name',
        'genome__taxon__name',
        'genome__taxon__taxonomy_id'
    )

    for item in Taxon.objects.values_list(
        'name',
        'taxonomy_id',
        'parent_id',
        'rank'
    ):
        taxon_lookup[item[1]] = item
    
    children_taxonomy_ids = list(taxon_lookup.keys())
    # children_counts holds counts of all children: genes, genomes and taxa
    children_nodes = defaultdict(set)
    # Set target_taxon_id to '1' (root node), and move it down later
    target_taxon_id = '1'
    
    
    genome_data = Counter()
    for gene_data in genes:
        labels.append(gene_data[0])
        parents.append(gene_data[1])
        children_nodes[gene_data[1]].add(gene_data[0])
        genome_data[gene_data[1]] += 1
        values.append(1)
        customdata.append([reverse('genedetails', args=(gene_data[1],gene_data[0])),
                           gene_data[0] + ' gene'
                           ])
    genome_names = list(genome_data.keys())
    genomes = Genome.objects.filter(
        name__in=genome_names
    ).values_list(
        'name',
        'taxon__name',
        'taxon__taxonomy_id'
    )
    
    for genome in genomes:
        # Add genome to graph
        labels.append(genome[0])
        parents.append(genome[1])
        children_nodes[genome[2]].add(genome[0])
        values.append(genome_data[genome[0]])
        customdata.append([reverse('genomedetails', args=(genome[0],)),
                           genome[0] + ' genome'
                           ])
        # Add taxon
        genome_taxid = genome[2]
        parent_id = taxon_lookup[genome_taxid][2]
        if parent_id in children_taxonomy_ids:
            children[parent_id][genome_taxid] = 1
        taxon_counts[genome_taxid] += genome_data[genome[0]]
        parent_lookup[genome_taxid] = parent_id
        children_nodes[parent_id].add(genome_taxid)
        # Add all parents of the genome up to the root or target_taxon_id
        while True:
            if parent_id in children_taxonomy_ids:
                children[parent_id][genome_taxid] = 1
                taxon_counts[parent_id] += genome_data[genome[0]]
            if parent_id == target_taxon_id \
            or parent_id == '1' \
            or parent_id not in taxon_lookup:
                break
            parent_lookup[parent_id] = taxon_lookup[parent_id][2]
            if parent_id != taxon_lookup[parent_id][2]:
                children[taxon_lookup[parent_id][2]][parent_id] = 1
            children_nodes[taxon_lookup[parent_id][2]].add(parent_id)
            parent_id = taxon_lookup[parent_id][2]

                           
    if len(taxon_counts) == 0:
        # No taxa, nothing to return
        return [], [], [], []
        
    # It's time to move root node
    if move_root_down:
        print(children_nodes)
        remove_nodes= set()
        new_root_id = target_taxon_id
        '''
        Going down from current root node, find taxon that 
        has more than one child and call it new root node
        '''
        while True:
            # Count children that have non-zero genome count, excluding root taxon
            '''
            children_nodes = [item for item in taxon_lookup.values() 
                              if item[2] == new_root_id 
                              and item[1] in children_counts
                              and item[1] != item[2]
                              ]
            '''
            if len(children_nodes[new_root_id]) == 1:
                remove_nodes.add(new_root_id)
                new_root_id = min(children_nodes[new_root_id])
            else:
                target_taxon_id = new_root_id
                break
        for taxon in remove_nodes:
            del(taxon_counts[taxon])
        
    for taxon in taxon_counts:
        # If the root taxon with taxonomy_id '1' is in the list, skip it
        if taxon == '1':
            continue
        labels.append(taxon_lookup[taxon][0])
        customdata.append([reverse('taxondetails', args=(taxon,)),
                           '<br>' + taxon_lookup[taxon][0] + ' [' + \
                           taxon_lookup[taxon][3] + ']<br>' + \
                           str(taxon_counts[taxon]) + ' genomes'
                           ])
        if taxon == target_taxon_id or taxon_lookup[taxon][2] == '1':
            # Parent of the root node must be empty string
            parents.append('')
        else:
            parents.append(taxon_lookup[parent_lookup[taxon]][0])
        values.append(taxon_counts[taxon])
    
    return labels, parents, values, customdata


def get_operons_taxonomy(operon_ids = [], gene_ids = []):
    '''
    Returns four lists for sunburst graph generation: list of gene/operon/genome/taxon names,
    list of parent genome/taxon names, list of gene/operon counts for each genome and taxon, list of 
    HTML-formatted strings containing a link to the current element
    
    Parameters:
    operon_ids (list of int): list of numerical operon ids
    gene_ids (list of int): list of numerical gene ids
    
    Steps:
    0.  Get list of all operons and associated genomes 
    1. Get list of all genes and associated genomes
    2. Get list of all genomes with all their taxonomy IDs
    3. Get list of all Taxon objects with their parent ID
    4. Calculate number of genomes
    
    '''
    labels = []
    parents = []
    values = []
    customdata = []
    parent_lookup = {}
    # key: genome name or taxonomy id, value: parent id
        
    taxon_counts = Counter()
    # key: name, value: count of children genes
    children = defaultdict(dict)
    taxon_lookup = {}
    move_root_down = True
    print('operon_ids:', operon_ids)
    operons = Operon.objects.filter(
        id__in=operon_ids
    ).values_list(
        'name',
        'genome__name',
        'genome__taxon__name',
        'genome__taxon__taxonomy_id'
    )
    
    genes = Gene.objects.filter(
        id__in=gene_ids
    ).values_list(
        'locus_tag',
        'genome__name',
        'genome__taxon__name',
        'genome__taxon__taxonomy_id'
    )

    for item in Taxon.objects.values_list(
        'name',
        'taxonomy_id',
        'parent_id',
        'rank'
    ):
        taxon_lookup[item[1]] = item
    
    children_taxonomy_ids = list(taxon_lookup.keys())
    # children_counts holds counts of all children: genes, genomes and taxa
    children_nodes = defaultdict(set)
    # Set target_taxon_id to '1' (root node), and move it down later
    target_taxon_id = '1'
    
    
    genome_data = Counter()
    for operon_data in operons:
        labels.append(operon_data[0])
        parents.append(operon_data[1])
        children_nodes[operon_data[1]].add(operon_data[0])
        genome_data[operon_data[1]] += 1
        values.append(1)
        customdata.append([reverse('operondetails', args=(operon_data[1],operon_data[0])),
                           operon_data[0] + ' operon'
                           ])

    for gene_data in genes:
        labels.append(gene_data[0])
        parents.append(gene_data[1])
        children_nodes[gene_data[1]].add(gene_data[0])
        genome_data[gene_data[1]] += 1
        values.append(1)
        customdata.append([reverse('genedetails', args=(gene_data[1],gene_data[0])),
                           gene_data[0] + ' gene'
                           ])
    genome_names = list(genome_data.keys())
    genomes = Genome.objects.filter(
        name__in=genome_names
    ).values_list(
        'name',
        'taxon__name',
        'taxon__taxonomy_id'
    )
    
    for genome in genomes:
        # Add genome to graph
        labels.append(genome[0])
        parents.append(genome[1])
        children_nodes[genome[2]].add(genome[0])
        values.append(genome_data[genome[0]])
        customdata.append([reverse('genomedetails', args=(genome[0],)),
                           genome[0] + ' genome'
                           ])
        # Add taxon
        genome_taxid = genome[2]
        parent_id = taxon_lookup[genome_taxid][2]
        if parent_id in children_taxonomy_ids:
            children[parent_id][genome_taxid] = 1
        taxon_counts[genome_taxid] += genome_data[genome[0]]
        parent_lookup[genome_taxid] = parent_id
        children_nodes[parent_id].add(genome_taxid)
        # Add all parents of the genome up to the root or target_taxon_id
        while True:
            if parent_id in children_taxonomy_ids:
                children[parent_id][genome_taxid] = 1
                taxon_counts[parent_id] += genome_data[genome[0]]
            if parent_id == target_taxon_id \
            or parent_id == '1' \
            or parent_id not in taxon_lookup:
                break
            parent_lookup[parent_id] = taxon_lookup[parent_id][2]
            if parent_id != taxon_lookup[parent_id][2]:
                children[taxon_lookup[parent_id][2]][parent_id] = 1
            children_nodes[taxon_lookup[parent_id][2]].add(parent_id)
            parent_id = taxon_lookup[parent_id][2]

                           
    if len(taxon_counts) == 0:
        # No taxa, nothing to return
        return [], [], [], []
        
    # It's time to move root node
    if move_root_down:
        print(children_nodes)
        remove_nodes= set()
        new_root_id = target_taxon_id
        '''
        Going down from current root node, find taxon that 
        has more than one child and call it new root node
        '''
        while True:
            # Count children that have non-zero genome count, excluding root taxon
            '''
            children_nodes = [item for item in taxon_lookup.values() 
                              if item[2] == new_root_id 
                              and item[1] in children_counts
                              and item[1] != item[2]
                              ]
            '''
            if len(children_nodes[new_root_id]) == 1:
                remove_nodes.add(new_root_id)
                new_root_id = min(children_nodes[new_root_id])
            else:
                target_taxon_id = new_root_id
                break
        for taxon in remove_nodes:
            del(taxon_counts[taxon])
        
    for taxon in taxon_counts:
        # If the root taxon with taxonomy_id '1' is in the list, skip it
        if taxon == '1':
            continue
        labels.append(taxon_lookup[taxon][0])
        customdata.append([reverse('taxondetails', args=(taxon,)),
                           '<br>' + taxon_lookup[taxon][0] + ' [' + \
                           taxon_lookup[taxon][3] + ']<br>' + \
                           str(taxon_counts[taxon]) + ' genomes'
                           ])
        if taxon == target_taxon_id or taxon_lookup[taxon][2] == '1':
            # Parent of the root node must be empty string
            parents.append('')
        else:
            parents.append(taxon_lookup[parent_lookup[taxon]][0])
        values.append(taxon_counts[taxon])
    
    return labels, parents, values, customdata

    
def get_root_id(labels, parents):
    '''
        Finds root element for sunburst graph.
        Returns the label for the first element, which parent ID is empty string
        If no such element, returns '1' (potentially unsafe operation)
    '''
    result = '1'
    for parent_ind, parent in enumerate(parents):
        if parent == '':
            result = labels[parent_ind]
            break
    return result
