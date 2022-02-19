import io
import os
import sys
import parasail
from datetime import datetime
from collections import defaultdict, OrderedDict
from subprocess import Popen, PIPE, CalledProcessError, STDOUT
# Importing necessary libraries from BioPython
from Bio import Phylo, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import toytree
import toyplot  #.html

from django.db.models import Q
from .models import *
from django.urls import reverse


COLORS = [(51, 34, 136),
    (136, 204, 238),
    (68, 170, 153),
    (17, 119, 51),
    (153, 153, 51),
    (221, 204, 119),
    (204, 102, 119),
    (136, 34, 85),
    (170, 68, 153),
    (119, 170, 221),
    (153, 221, 255),
    (68, 187, 153),
    (187, 204, 51),
    (170, 170, 0),
    (238, 221, 136),
    (238, 136, 102),
    (255, 170, 187),

    (209, 187, 215),
    (174, 118, 163),
    (136, 46, 114),
    (25, 101, 176),
    (82, 137, 199),
    (123, 175, 222),
    (78, 178, 101),
    (144, 201, 135),
    (202, 224, 171),
    (247, 240, 86),
    (246, 193, 65),
    (241, 147, 45),
    (232, 96, 28),
    (221, 221, 221),
    (232, 236, 251),
    (217, 204, 227),
    (153, 79, 136),
    (97, 149, 207),
    (247, 203, 69),
    (165, 23, 14),
    (114, 25, 14),
    (66, 21, 10)
]
DARK_GREY = '.setColorGradient(\'rgb(119, 119, 119)\', \'rgb(199, 199, 199)\')'
RED = '.setColorGradient(\'rgb(255, 85, 92)\', \'rgb(220, 5, 12)\')'


def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def get_color(index):
    color_index = index%len(COLORS)
    light_color = []
    for rgb_val in COLORS[color_index]:
        if rgb_val + 80 > 255:
            light_color.append(255)
        else:
            light_color.append(rgb_val + 80)
    # color = '\'rgb(' + ', '.join([str(x) for x in COLORS[color_index]]) + ')\''
    color = '.setColorGradient(\'rgb(' + ', '.join([str(x) for x in light_color]) + ')\', \'rgb(' + ', '.join([str(x) for x in COLORS[color_index]]) + ')\')'
    return color
    
    
def get_gradients():
    gradients = ['setColorGradient(\'rgb(204, 193, 208)\', \'rgb(96, 74, 123)\')',
        'setColorGradient(\'rgb(122, 145, 198)\', \'rgb(144, 171, 234)\',\'rgb(228,230,255)\',\'rgb(144, 171, 234)\', \'rgb(122, 145, 198)\')',
        'setColorGradient(\'rgb(215, 228, 189)\', \'rgb(119, 147, 60)\')',
        'setColorGradient(\'rgb(252, 213, 181)\', \'rgb(228, 108, 10)\')',
        'setColorGradient(\'rgb(255, 189, 192)\', \'rgb(255, 0, 0)\')',
        'setColorGradient(\'rgb(232, 227, 155)\', \'rgb(255, 0, 0)\')',
        'setColorGradient(\'rgb(217, 217, 217)\', \'rgb(38, 38, 38)\')'
        ]
    
    return '\'rgb(204, 193, 208\')'


def add_gene(gene, gene_uid, track_uid, offset, reverse_gene, gene_color, group, request):
    result = []
    locus_size = 10000
    if gene.strand == 1:
        if reverse_gene:
            strand = '-'
        else:
            strand = '+'
    else:
        if reverse_gene:
            strand = '+'
        else:
            strand = '-'
    
    size = gene.end - gene.start + 1
    if reverse_gene:
        start = locus_size + offset - gene.end
    else:
        start = gene.start - offset

    if start < 0:
        size = size + start
        start = 0
    elif size + start > locus_size:
        size = locus_size - start
    
    result.append('\n\t\t// Add gene ' + gene.locus_tag + ' to lane ' + str(track_uid) + '; offset ' + str(offset) + ' coords ' + str(gene.start) + ' ' + str(gene.end))
    result.append('\t\tgene' + str(gene_uid) + ' = lane' + str(track_uid) + '.addGene( ' + str(start) + ', ' + str(size) + ' , \'' + strand + '\');')
    result.append('\t\tgene' + str(gene_uid) + gene_color + ';')
    result.append('\t\tgene' + str(gene_uid) + '.text.style = "Arial";')
    if gene.name == '' or gene.name is None:
        result.append('\t\tgene' + str(gene_uid) + '.name = "' + gene.locus_tag + '";')
    else:
        result.append('\t\tgene' + str(gene_uid) + '.name = "' + gene.name + '";')
    if group == '':
        if gene.type == 'CDS':
            result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' + gene.locus_tag + ': ' + gene.function + '";')
        else:
            result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' + gene.locus_tag + ': ' + gene.function + ' [' + gene.type + ']";')
    else:
        result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' + gene.locus_tag + ': ' + gene.function + ' [' + group + ']";')
    gene_url = request.build_absolute_uri(reverse('genedetails', args=[gene.genome.name, gene.locus_tag]))
    result.append('\t\tgene' + str(gene_uid) + '.onClick = function() {window.open("' + gene_url + '", "_blank");};')
    return result


def make_taxonomy_tree(taxonomy_ids, start_taxon, top_taxon):
    print('Starting tree creation')
    tree = defaultdict(dict)
    root_id = '1'
    all_taxa = Taxon.objects.all().values('taxonomy_id', 'parent_id')
    parents = {item['taxonomy_id']:item['parent_id'] for item in all_taxa}
    children = defaultdict(list)
    for taxon in all_taxa:
        if taxon['parent_id'] == '1' and taxon['taxonomy_id'] == '1':
            continue
        children[taxon['parent_id']].append(taxon['taxonomy_id'])
    tree[root_id]['parent'] = parents[root_id]
    tree[root_id]['children'] = children[root_id]
    tree[top_taxon]['parent'] = parents[top_taxon]
    tree[top_taxon]['children'] = children[top_taxon]
    for taxonomy_id in taxonomy_ids:
        if taxonomy_id in tree:
            continue
        parent_id = parents[taxonomy_id]
        tree[taxonomy_id]['parent'] = parent_id
        if taxonomy_id in children:
            tree[taxonomy_id]['children'] = children[taxonomy_id]
        while parent_id != root_id:
            if parent_id in tree:
                break
            tree[parent_id]['parent'] = parents[parent_id]
            tree[parent_id]['children'] = children[parent_id]
            parent_id = parents[parent_id]
    print('Tree created')
    return tree


def make_clustal_alignment(infasta):
    result = ''
    args = ['clustalo', '-i', '-']
    print('Running clustalo')
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        outfasta, err = p.communicate(infasta)
    print('clustalo finished')
    if p.returncode != 0:
        print('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Clustal omega finished with error:\n'+ ' '.join(err) + '\n')
        return result
    result = outfasta
    return result


def make_muscle_alignment(infasta):
    result = None
    muscle_cline = MuscleCommandline()
    print(muscle_cline)
    
    print('Running MUSCLE')
    with Popen(str(muscle_cline), stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        outfasta, err = p.communicate(infasta)
    print('MUSCLE finished')
    if p.returncode != 0:
        print('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] MUSCLE finished with error:\n'+ ' '.join(err) + '\n')
        return result
    result = outfasta
    return result


def sort_nodes(tree, target):
    result = [target]
    path = [tree.root] + tree.get_path(target)
    for node in reversed(path):
        for leaf in node.get_terminals(order='level'):
            if leaf.name not in result:
                result.append(leaf.name)
    return result


def make_protein_tree(proteins):
    """
        Input:
        proteins: list of tuples. The first element of tuples is gene ID, the second is protein sequence, the third is gene label for the tree
    """
    infasta = []
    gene_labels = {}
    for protein in proteins:
        infasta.append('>' + str(protein[2]))
        gene_labels[protein[2]] = protein[0]
        infasta.append(protein[1])
    #outfasta = make_clustal_alignment('\n'.join(infasta))
    outfasta = make_muscle_alignment('\n'.join(infasta))
    if outfasta is None:
        return [proteins[0][2]], '', ''
    align = AlignIO.read(io.StringIO(outfasta), 'fasta')

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)
    #print(distMatrix)
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using NJ algorithm
    NJTree = constructor.nj(distMatrix)
    NJTree.root_at_midpoint()
    node_ids = sort_nodes(NJTree, proteins[0][2])
    nodes = [gene_labels[item] for item in node_ids]
    # Draw the phlyogenetic tree using terminal
    #Phylo.draw_ascii(NJTree)
    newick = io.StringIO()
    Phylo.write(NJTree, newick, 'newick')
    tree = toytree.tree(newick.getvalue(), tree_format=1)
    node_ids.reverse()
    print('Reversed node IDs', node_ids)
    height = 70 + 56 * len(nodes)
    width = 200
    canvas = toyplot.Canvas(width=width, height=height)
    canvas.style['background-color'] = 'white'
    axes = canvas.cartesian(bounds=(10, width - 20, 82, height - 24), padding=0)  #, ymin=0, ymax=20)
    axes.x.spine.position = 'high'
    tree.draw(axes=axes, width=width + 10, height=height, fixed_order=node_ids, scalebar=True, shrink=10)
    tree_canvas = toyplot.html.tostring(canvas)
    return nodes, tree_canvas, newick.getvalue()


def sort_proteins(proteins, query_hash):
    """Sorts proteins by semi-global alignment score"""
    result = [query_hash]
    profile = parasail.profile_create_stats_sat(proteins[query_hash], parasail.blosum62)
    scores = []
    for protein_hash, protein_seq in proteins.items():
        if protein_hash == query_hash:
            continue
        #print(ref['id'])
        search_result = parasail.sg_stats_scan_profile_sat(profile,protein_seq,11,1)
        #print(ref['id'], result.score, result.matches, result.similar, result.end_query, result.end_ref)
        scores.append((protein_hash,search_result.score))
    result += [x[0] for x in reversed(sorted(scores, key=lambda x: x[1]))]
    return result


def get_sorted_orthologs(eggnog_og, pivot_gene, genelist_size=50):
    ''' Returns list of genes sorted by protein siilarityto the pivot genee'''
    ret_genes = [pivot_gene]
#    genes = {item.id:item for item in Gene.objects.filter(protein__ortholog_groups__id=eggnog_og.id).select_related('protein', 'genome', 'genome__strain', 'genome__sample', 'contig')}
    genes = {item.id:item for item in Gene.objects.filter(protein__ortholog_groups__id=eggnog_og.id).select_related('protein', 'genome', 'genome__taxon', 'contig')}
    if len(genes) == 1:
        return ret_genes, len(genes), '', ''
    proteins = {gene.protein.protein_hash:gene.protein.sequence for gene in genes.values()}
    sorted_proteins = sort_proteins(proteins, pivot_gene.protein.protein_hash)
    gene2protein = defaultdict(list)
    gene_count = 0
    tree_proteins = [(pivot_gene.id, pivot_gene.protein.sequence, pivot_gene.locus_tag)]
    for protein_hash in sorted_proteins:
        if gene_count >= genelist_size:
            break
        prot_genes = [x for x in genes.values() if x.protein.protein_hash == protein_hash]
        for gene in prot_genes:
            if gene.id == pivot_gene.id:
                continue
            #ret_genes.append(gene)
            tree_proteins.append((gene.id, proteins[protein_hash], gene.locus_tag))
            gene_count += 1
            if gene_count >= genelist_size:
                break
    tree_nodes, tree_svg, tree_newick = make_protein_tree(tree_proteins)
    #print(tree_svg)
    #print(tree_newick)
    if tree_nodes[0] == pivot_gene.id:
        print('First gene is ', pivot_gene.locus_tag)
    else:
        print('First gene is ', genes[tree_nodes[0]], 'instead of', pivot_gene.locus_tag)
    for gene_id in tree_nodes:
        if gene_id == pivot_gene.id:
            continue
        ret_genes.append(genes[int(gene_id)])
    return ret_genes, len(genes), tree_svg, tree_newick


def get_orthologs(eggnog_og, pivot_gene):
    ''' Returns prioritized list of genes'''
    result = [pivot_gene]
    genes = Gene.objects.filter(protein__ortholog_groups__id=eggnog_og.id).select_related('protein', 'genome', 'genome__strain__taxon')
    protein_ids = list(set([gene.protein.protein_hash for gene in genes]))
    
    tree = make_taxonomy_tree([gene.genome.strain.taxon.taxonomy_id for gene in genes], pivot_gene.genome.strain.taxon.taxonomy_id, eggnog_og.taxon.taxonomy_id)

    #gene2taxon = {}
    taxon2genes = autovivify(2, dict)
    for gene in genes:
        genome_id = gene.genome.id
        #gene2taxon[gene.id] = gene.genome.strain.taxon
        taxon2genes[gene.genome.strain.taxon.taxonomy_id][genome_id][gene.id] = gene
    eggnog_top_taxon = eggnog_og.taxon
    
    start_genome = pivot_gene.genome.id
    start_taxon = pivot_gene.genome.strain.taxon.taxonomy_id
    visited_taxa = set()
    # Get genes from start genome
    for gene_id in sorted(taxon2genes[start_taxon][start_genome].keys()):
        if gene_id == pivot_gene.id:
            continue
        result.append(taxon2genes[start_taxon][start_genome][gene_id])
    # Get genes from other genomes of the start taxon
    for genome_id in sorted(taxon2genes[start_taxon].keys()):
        if genome_id == start_genome:
            continue
        for gene_id in sorted(taxon2genes[start_taxon][genome_id].keys()):
            if taxon2genes[start_taxon][genome_id][gene_id] not in result:
                result.append(taxon2genes[start_taxon][genome_id][gene_id])
    visited_taxa.add(start_taxon)
    # Traverse the tree starting from start taxon:
    if 'children' in tree[start_taxon]:
        for child_taxon in sorted(tree[start_taxon]['children']):
            if child_taxon in visited_taxa:
                continue
            child_ordered_genes, child_visited_taxa = collect_genes(tree, taxon2genes, visited_taxa, child_taxon)
            visited_taxa = visited_taxa.union(child_visited_taxa)
            result += child_ordered_genes

    parent_id = tree[start_taxon]['parent']
    while True:
        child_ordered_genes, child_visited_taxa = collect_genes(tree, taxon2genes, visited_taxa, parent_id)
        visited_taxa = visited_taxa.union(child_visited_taxa)
        result += child_ordered_genes
        try:
            parent_id = tree[parent_id]['parent']
        except KeyError:
            print('Parent node not found ', parent_id)
            parent_id = '1'
        if parent_id == '1':
            break

    if parent_id == '1':
        child_ordered_genes, child_visited_taxa = collect_genes(tree, taxon2genes, visited_taxa, parent_id)
        result += child_ordered_genes
    return result


def collect_genes(tree, taxon2genes, visited_taxa, taxon):
    print('Visiting ' + taxon)
    if taxon in visited_taxa:
        print(taxon + ' has been visited before')
        return [], visited_taxa
    ordered_genes = []
    if taxon in taxon2genes and taxon not in visited_taxa:
        for genome_id in sorted(taxon2genes[taxon].keys()):
            for gene_id in sorted(taxon2genes[taxon][genome_id].keys()):
                ordered_genes.append(taxon2genes[taxon][genome_id][gene_id])
    visited_taxa.add(taxon)
    if 'children' in tree[taxon]:
        for child_taxon in sorted(tree[taxon]['children']):
            if child_taxon in visited_taxa:
                continue
            child_ordered_genes, child_visited_taxa = collect_genes(tree, taxon2genes, visited_taxa, child_taxon)
            ordered_genes += child_ordered_genes
            visited_taxa = visited_taxa.union(child_visited_taxa)
    return ordered_genes, visited_taxa


def get_scribl(start_gene, eggnog_og, request):
    locus_size = 10000
    genelist_size = 50
    gene_uid = 0
    track_uid = 0
    scribl = []
    eggnog2color = {eggnog_og.eggnog_id:RED, '':DARK_GREY}
    #ordered_orthologs = get_orthologs(eggnog_og, start_gene)
    ordered_orthologs, og_gene_count, tree_canvas, tree_newick = get_sorted_orthologs(eggnog_og, start_gene, genelist_size)
    if og_gene_count == 1:
        return '', '', '', og_gene_count, 1
    plot_gene_count = len(ordered_orthologs)
    #if len(ordered_orthologs) > genelist_size:
    #    ordered_orthologs = ordered_orthologs[:genelist_size]

    scribl.append('\t\tcanvas.width = parent.offsetWidth - 201;')
    scribl.append('\t\tcanvas.height = ' + str(50 + 56 * len(ordered_orthologs)) + ';')
    #scribl.append('parent.offsetHeight = ' + str(50 + 55 * len(ordered_orthologs)) + ';')
    scribl.append('\t\tctx.font = "12px Arial";')
    scribl.append('\t\tctx.fillStyle = "#dbdbdb";')
    scribl.append('\t\tchart = new Scribl(canvas, canvas.width - 25);')
    scribl.append('\t\tchart.scale.pretty = false;')
    scribl.append('\t\tchart.scale.min = -50;')
    scribl.append('\t\tchart.scale.max = ' + str(locus_size) +';')
    scribl.append('\t\tchart.tick.major.size = 500;')
    scribl.append('\t\tchart.tick.minor.size = 50;')
    scribl.append('\t\tchart.tick.major.color = "#c7c7c7";')
    scribl.append('\t\tchart.tick.minor.color = "white";')
    scribl.append('\t\tchart.scale.font.color = "#c7c7c7";')
    scribl.append('\t\tchart.laneSizes = 16;')
    scribl.append('\t\tchart.scale.font.size = 12;')
    scribl.append('\t\tchart.scale.size = 20;')

    scribl.append('\t\tchart.trackBuffer = 40;')
    
    # Add empty track to make space below scale
    scribl.append('\t\ttrack' + str(track_uid) + ' = chart.addTrack();')
    
    for gene in ordered_orthologs:
        # Create first track
        # Crete track for organism name
        reverse = True
        if gene.strand == 1:
            reverse = False
        track_uid += 1
        scribl.append('\t\ttrack' + str(track_uid) + ' = chart.addTrack();')
        scribl.append('\t\tlane' + str(track_uid) + ' = track' + str(track_uid) + '.addLane();')
        middle_point = (gene.end + gene.start) / 2
        offset = int(middle_point - locus_size / 2)
        if reverse:
            display_start = offset - locus_size
            if display_start < 1:
                display_start = 1
            display_end = offset
            if display_end > gene.contig.size:
                display_end = gene.contig.size
            scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' [' + gene.genome.taxon.name + '] ' +
                gene.contig.contig_id + ': complement(' + str(display_start) + '..' + str(display_end) + ')\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
#            if gene.genome.strain is not None:
#                scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' genome [' + gene.genome.strain.full_name + '] ' +
#                    gene.contig.contig_id + ': complement(' + str(display_start) + '..' + str(display_end) + ')\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
#            else:
#                scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' genome [' + gene.genome.sample.full_name + '] ' +
#                    gene.contig.contig_id + ': complement(' + str(display_start) + '..' + str(display_end) + ')\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
        else:
            display_start = offset
            if display_start < 1:
                display_start = 1
            display_end = offset + locus_size
            if display_end > gene.contig.size:
                display_end = gene.contig.size
            scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' [' + gene.genome.taxon.name + '] ' +
                gene.contig.contig_id + ': ' + str(display_start) + '..' + str(display_end) + '\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
#            if gene.genome.strain is not None:
#                scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' genome [' + gene.genome.strain.full_name + '] ' +
#                    gene.contig.contig_id + ': ' + str(display_start) + '..' + str(display_end) + '\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
#            else:
#                scribl.append('\t\tctx.fillText(\'' + gene.genome.name + ' genome [' + gene.genome.sample.full_name + '] ' +
#                    gene.contig.contig_id + ': ' + str(display_start) + '..' + str(display_end) + '\', 10, track' + str(track_uid) + '.getPixelPositionY() - 8);')
        # Crete track for gene glyphs
        #track_uid += 1
        middle_point = (gene.end + gene.start) / 2
        offset = int(middle_point - locus_size / 2)
        range_start = int(middle_point - locus_size / 2)
        range_end = int(middle_point + locus_size / 2)
        neigbor_genes = Gene.objects.filter(genome=gene.genome, contig=gene.contig).filter(Q(end__range=[range_start, range_end]) |Q(start__range=[range_start, range_end])).select_related('protein', 'genome').prefetch_related('protein__ortholog_groups')
        for locus_member in neigbor_genes:
            group = ''
            if locus_member in ordered_orthologs:
                group = eggnog_og.eggnog_id
            elif locus_member.protein is None:
                group = ''
            else:
                eggnogs = sorted([item.eggnog_id for item in locus_member.protein.ortholog_groups.all()])
                for og_id in eggnogs:
                    if og_id.startswith('COG'):
                        group = og_id
                        break
                if group == '' and eggnogs:
                    group = eggnogs[0]
            if group not in eggnog2color:
                eggnog2color[group] = get_color(len(eggnog2color))
            gene_color = eggnog2color[group]
            gene_uid += 1
            scribl += add_gene(locus_member, gene_uid, track_uid, offset, reverse, gene_color, group, request)

    return '\n'.join(scribl), tree_canvas, tree_newick, og_gene_count, plot_gene_count


