""" This file contains functions for search in taxonomy"""
from collections import defaultdict, Counter
from django.urls import reverse
from browser.models import Config, Taxon, Genome

def generate_sunburst(taxon_id='1'):
    import plotly.graph_objects as go
    from io import StringIO
    
    labels, parents, values, customdata = get_genomes_taxonomy(taxon_id)
    if not labels:
        return ''
    
    maxdepth = -1
    if taxon_id == '1':
        maxdepth = 4
    
    colors = []
    color_sequence = ['', '#FFAA00', '#2D5F91','#819FBD','#819FBD','#91D4D2', '#96BEE6', '#C0D8F0','#E8655F','#F1A39F','#48B7B4']
    for ind, label in enumerate(labels):
        colors.append(color_sequence[ind%10])
    
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
    fig.update_layout(margin = dict(t = 0, l = 0, r = 0, b = 0), height=800, uniformtext=dict(minsize=10, mode='hide'), paper_bgcolor='rgba(120,120,120,0.1)', plot_bgcolor='rgba(120,120,120,0.1)')
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
        parent = Taxon.objects.get(taxonomy_id=taxonomy_id)
    except Taxon.DoesNotExist:
        return result
    result.append(taxonomy_id)
    parent_ids = [taxonomy_id]
    while True:
        children_ids = list(Taxon.objects.values_list('taxonomy_id', flat=True).filter(parent_id__in=parent_ids))
        if len(children_ids) == 0:
            break
        result += children_ids
        parent_ids = children_ids
    return result

    
def get_genomes_taxonomy(target_taxon_id):
    '''
    Returns three lists for sunburst graph generation: list of genome and taxon names, list of parent taxon names, list of genome numbers for each taxon
    
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
    # key: name, value: count
    children = defaultdict(dict)
    taxon_lookup = {}
    taxon_ranks = {}
    
    if target_taxon_id == '1':
        genomes = Genome.objects.values_list('name', 'taxon__name', 'taxon__taxonomy_id')
        for item in Taxon.objects.values_list('name', 'taxonomy_id', 'parent_id', 'rank'):
            taxon_lookup[item[1]] = item
        children_taxonomy_ids = list(taxon_lookup.keys())
    else:
        children_taxonomy_ids = get_taxon_children(target_taxon_id)
        genomes = Genome.objects.values_list('name', 'taxon__name', 'taxon__taxonomy_id').filter(taxon__taxonomy_id__in=children_taxonomy_ids)
        for item in Taxon.objects.values_list('name', 'taxonomy_id', 'parent_id', 'rank').filter(taxonomy_id__in=children_taxonomy_ids):
            taxon_lookup[item[1]] = item

    data = defaultdict(dict)
    for genome in genomes:
        # Add genome to graph
        labels.append(genome[0])
        parents.append(genome[1])
        values.append(1)
        customdata.append([reverse('genomedetails', args=(genome[0],)), genome[0] + ' genome'])
        # Add taxon
        genome_taxid = genome[2]
        parent_id = taxon_lookup[genome_taxid][2]
        if parent_id in children_taxonomy_ids:
            children[parent_id][genome_taxid] = 1
        taxon_counts[genome_taxid] += 1
        parent_lookup[genome_taxid] = parent_id
        # Add parents
        while True:
            if parent_id in children_taxonomy_ids:
                children[parent_id][genome_taxid] = 1
                taxon_counts[parent_id] += 1
            if parent_id == target_taxon_id or parent_id == '1' or parent_id not in taxon_lookup:
                break
            parent_lookup[parent_id] = taxon_lookup[parent_id][2]
            if parent_id != taxon_lookup[parent_id][2]:
                children[taxon_lookup[parent_id][2]][parent_id] = 1
            parent_id = taxon_lookup[parent_id][2]
    
    if len(taxon_counts) == 0:
        # No taxa, nothing to return
        return [], [], [], []
        
    for taxon in taxon_counts:
        # Hide root node
        if taxon == '1':
            continue
        
        labels.append(taxon_lookup[taxon][0])
        customdata.append([reverse('taxondetails', args=(taxon,)), '<br>' + taxon_lookup[taxon][0] + ' [' + taxon_lookup[taxon][3] + ']<br>' + str(taxon_counts[taxon]) + ' genomes'])
        if taxon == target_taxon_id or taxon_lookup[taxon][2] == '1':
            # Parent of the root node must be empty string
            parents.append('')
        else:
            parents.append(taxon_lookup[parent_lookup[taxon]][0])
        values.append(taxon_counts[taxon])
    
    return labels, parents, values, customdata

