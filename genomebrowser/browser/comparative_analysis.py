import io
import parasail
import logging
from datetime import datetime
from collections import defaultdict
from subprocess import Popen, PIPE, STDOUT
# Importing necessary libraries from BioPython
from Bio import Phylo, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import toytree
import toyplot

from django.db.models import Q
from .models import Gene
from django.urls import reverse

logger = logging.getLogger("CGCMS")


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

def _get_color(index):
    color_index = index%len(COLORS)
    light_color = []
    for rgb_val in COLORS[color_index]:
        if rgb_val + 80 > 255:
            light_color.append(255)
        else:
            light_color.append(rgb_val + 80)
    color = '.setColorGradient(\'rgb(' + ', '.join(
                [str(x) for x in light_color]
            ) + ')\', \'rgb(' + ', '.join(
                [str(x) for x in COLORS[color_index]]
            ) + ')\')'
    return color
    
def add_gene(gene, gene_uid, track_uid, offset, reverse_gene,
             gene_color, group, request, locus_size):
    result = []
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
    
    result.append('\n\t\t// Add gene ' + gene.locus_tag + ' to lane ' +
                  str(track_uid) + '; offset ' + str(offset) + ' coords ' +
                  str(gene.start) + ' ' + str(gene.end)
                  )
    result.append('\t\tgene' + str(gene_uid) + ' = lane' + str(track_uid) +
                  '.addGene( ' + str(start) + ', ' + str(size) +
                  ' , \'' + strand + '\');'
                  )
    result.append('\t\tgene' + str(gene_uid) + gene_color + ';')
    result.append('\t\tgene' + str(gene_uid) + '.text.style = "Arial";')
    if gene.name == '' or gene.name is None:
        result.append('\t\tgene' + str(gene_uid) + '.name = "' + gene.locus_tag + '";')
    else:
        result.append('\t\tgene' + str(gene_uid) + '.name = "' + gene.name + '";')
    if group == '':
        if gene.type == 'CDS':
            result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' +
                          gene.locus_tag + ': ' + gene.function + '";'
                          )
        else:
            result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' +
                          gene.locus_tag + ': ' + gene.function + ' [' +
                          gene.type + ']";'
                          )
    else:
        result.append('\t\tgene' + str(gene_uid) + '.onMouseover = "' +
                      gene.locus_tag + ': ' + gene.function + ' [' + group + ']";'
                      )
    gene_url = request.build_absolute_uri(reverse('genedetails',
                                          args=[gene.genome.name, gene.locus_tag])
                                          )
    result.append('\t\tgene' +
                  str(gene_uid) +
                  '.onClick = function() {window.open("' +
                  gene_url +
                  '", "_blank");};'
                  )
    return result

def make_muscle_alignment(infasta):
    result = None
    muscle_cline = MuscleCommandline()
    logger.debug(str(muscle_cline))
    
    logger.info('Running MUSCLE')
    with Popen(str(muscle_cline),
               stdin=PIPE,
               stdout=PIPE,
               stderr=STDOUT,
               bufsize=1,
               universal_newlines=True) as p:
        outfasta, err = p.communicate(infasta)
    logger.info('MUSCLE finished')
    if p.returncode != 0:
        logger.error('[' + 
              datetime.now().strftime("%d/%m/%Y %H:%M:%S") +
              '] MUSCLE finished with error:\n' +
              ' '.join(err) +
              '\n'
              )
        return result
    result = outfasta
    return result

def make_protein_tree(proteins):
    """
        Input:
        proteins: list of tuples. 
        The first element of each tuple is gene ID,
        the second is protein sequence,
        the third is gene label to show on the tree
    """
    infasta = []
    gene_labels = {}
    for protein in proteins:
        infasta.append('>' + str(protein[2]))
        gene_labels[protein[2]] = protein[0]
        infasta.append(protein[1])
    outfasta = make_muscle_alignment('\n'.join(infasta))
    if outfasta is None:
        return [proteins[0][2]], '', ''
    align = AlignIO.read(io.StringIO(outfasta), 'fasta')

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using NJ algorithm
    NJTree = constructor.nj(distMatrix)
    NJTree.root_with_outgroup(proteins[0][2])
    nodes = [gene_labels[term.name] for term
             in NJTree.get_terminals(order='postorder')
             ]
    node_ids = [term.name for term in NJTree.get_terminals(order='postorder')]
    newick = io.StringIO()
    Phylo.write(NJTree, newick, 'newick')
    tree = toytree.mtree(newick.getvalue(), tree_format=1)
    tip_labels = [item.split('|')[0] for item in tree.treelist[0].get_tip_labels()]
    nodes.reverse()
    logger.debug('Node IDs %s', str(node_ids))
    canvas, axes, marks = tree.draw(width=200,
                                    height=70 + 56 * len(nodes),
                                    fixed_order=node_ids,
                                    tip_labels=tip_labels,
                                    scalebar=True
                                    )
    canvas.style['background-color'] = 'white'
    canvas.style['padding-top'] = '60px'
    tree_canvas = toyplot.html.tostring(canvas)
    return nodes, tree_canvas, newick.getvalue()

def sort_proteins(proteins, query_hash):
    """Sorts proteins by semi-global alignment score"""
    result = [query_hash]
    profile = parasail.profile_create_stats_sat(proteins[query_hash],
                                                parasail.blosum62
                                                )
    scores = []
    for protein_hash, protein_seq in proteins.items():
        if protein_hash == query_hash:
            continue
        search_result = parasail.sg_stats_scan_profile_sat(profile,protein_seq,11,1)
        scores.append((protein_hash,search_result.score))
    result += [x[0] for x in reversed(sorted(scores, key=lambda x: x[1]))]
    return result

def get_sorted_orthologs(eggnog_og, pivot_gene, genelist_size=50):
    ''' Returns list of genes sorted by protein siilarityto the pivot genee'''
    ret_genes = [pivot_gene]
    genes = {item.id:item
             for item
             in Gene.objects.filter(
                protein__ortholog_groups__id=eggnog_og.id).select_related(
                'protein', 'genome', 'genome__taxon', 'contig'
                )
            }
    if len(genes) == 1:
        return ret_genes, len(genes), '', ''
    proteins = {gene.protein.protein_hash:gene.protein.sequence
                for gene
                in genes.values()
                }
    sorted_proteins = sort_proteins(proteins, pivot_gene.protein.protein_hash)
    gene_count = 0
    tree_proteins = [(pivot_gene.id,
                      pivot_gene.protein.sequence,
                      pivot_gene.locus_tag + '|' + pivot_gene.genome.name
                      )]
    for protein_hash in sorted_proteins:
        if gene_count >= genelist_size:
            break
        prot_genes = [x for x
                      in genes.values()
                      if x.protein.protein_hash == protein_hash
                      ]
        for gene in prot_genes:
            if gene.id == pivot_gene.id:
                continue
            tree_proteins.append((gene.id,
                                  proteins[protein_hash],
                                  gene.locus_tag + '|' + gene.genome.name
                                  ))
            gene_count += 1
            if gene_count >= genelist_size:
                break
    tree_nodes, tree_svg, tree_newick = make_protein_tree(tree_proteins)
    if tree_nodes[0] == pivot_gene.id:
        logger.debug('First gene is %s', pivot_gene.locus_tag)
    else:
        logger.debug('First gene is %s instead of %s',
                     genes[tree_nodes[0]],
                     pivot_gene.locus_tag
                     )
    for gene_id in tree_nodes:
        if gene_id == pivot_gene.id:
            continue
        ret_genes.append(genes[int(gene_id)])
    return ret_genes, len(genes), tree_svg, tree_newick

def collect_genes(tree, taxon2genes, visited_taxa, taxon):
    logger.debug('Visiting ' + taxon)
    if taxon in visited_taxa:
        logger.debug(taxon + ' has been visited before')
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
            child_ordered_genes, child_visited_taxa = collect_genes(tree,
                                                                    taxon2genes,
                                                                    visited_taxa,
                                                                    child_taxon
                                                                    )
            ordered_genes += child_ordered_genes
            visited_taxa = visited_taxa.union(child_visited_taxa)
    return ordered_genes, visited_taxa

def get_scribl(start_gene, eggnog_og, request):
    locus_size = int(request.GET.get('size')) * 1000 #10000
    genelist_size = int(request.GET.get('lines'))  #10
    gene_uid = 0
    track_uid = 0
    scribl = []
    eggnog2color = {eggnog_og.eggnog_id:RED, '':DARK_GREY}
    ordered_orthologs, og_gene_count, tree_canvas, tree_newick = get_sorted_orthologs(
        eggnog_og, start_gene, genelist_size
        )
    if og_gene_count == 1:
        return '', '', '', og_gene_count, 1
    plot_gene_count = len(ordered_orthologs)
    scribl.append('\t\tcanvas.width = parent.offsetWidth - 201;')
    scribl.append('\t\tcanvas.height = ' +
                  str(50 + 56 * len(ordered_orthologs)) + ';'
                  )
    scribl.append('\t\tctx.font = "12px Arial";')
    scribl.append('\t\tctx.fillStyle = "#dbdbdb";')
    scribl.append('\t\tchart = new Scribl(canvas, canvas.width - 25);')
    scribl.append('\t\tchart.scale.pretty = false;')
    scribl.append('\t\tchart.scale.min = -50;')
    scribl.append('\t\tchart.scale.max = ' + str(locus_size) +';')
    scribl.append('\t\tchart.tick.major.size = ' + str(locus_size / 10) + ';')
    scribl.append('\t\tchart.tick.minor.size = ' + str(locus_size / 100) + ';')
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
        scribl.append('\t\ttrack' + str(track_uid) +
                      ' = chart.addTrack();'
                      )
        scribl.append('\t\tlane' + str(track_uid) +
                      ' = track' + str(track_uid) +
                      '.addLane();'
                      )
        middle_point = (gene.end + gene.start) / 2
        offset = int(middle_point - locus_size / 2)
        taxon_name = gene.genome.taxon.name
        taxon_name = taxon_name.replace("'", "")
        if reverse:
            display_start = offset - locus_size
            if display_start < 1:
                display_start = 1
            display_end = offset
            if display_end > gene.contig.size:
                display_end = gene.contig.size
            scribl.append('\t\tctx.fillText(\'' + 
                          gene.genome.name + ' [' +
                          taxon_name + '] ' +
                          gene.contig.contig_id + 
                          ': complement(' +
                          str(display_start) +
                          '..' + str(display_end) +
                          ')\', 10, track' +
                          str(track_uid) +
                          '.getPixelPositionY() - 8);'
                          )
        else:
            display_start = offset
            if display_start < 1:
                display_start = 1
            display_end = offset + locus_size
            if display_end > gene.contig.size:
                display_end = gene.contig.size
            scribl.append('\t\tctx.fillText(\'' + 
                          gene.genome.name +
                          ' [' + taxon_name + '] ' +
                          gene.contig.contig_id +
                          ': ' + str(display_start) +
                          '..' + str(display_end) +
                          '\', 10, track' +
                          str(track_uid) +
                          '.getPixelPositionY() - 8);'
                          )
        # Crete track for gene glyphs
        middle_point = (gene.end + gene.start) / 2
        offset = int(middle_point - locus_size / 2)
        range_start = int(middle_point - locus_size / 2)
        range_end = int(middle_point + locus_size / 2)
        neigbor_genes = Gene.objects.filter(
                        genome=gene.genome,
                        contig=gene.contig).filter(
                        Q(end__range=[range_start, range_end])|
                        Q(start__range=[range_start, range_end])).select_related(
                        'protein', 'genome'
                        ).prefetch_related(
                        'protein__ortholog_groups'
                        )
        for locus_member in neigbor_genes:
            group = ''
            if locus_member in ordered_orthologs:
                group = eggnog_og.eggnog_id
            elif locus_member.protein is None:
                group = ''
            else:
                eggnogs = sorted([item.eggnog_id
                                  for item 
                                  in locus_member.protein.ortholog_groups.all()])
                for og_id in eggnogs:
                    if og_id.startswith('COG'):
                        group = og_id
                        break
                if group == '' and eggnogs:
                    group = eggnogs[0]
            if group not in eggnog2color:
                eggnog2color[group] = _get_color(len(eggnog2color))
            gene_color = eggnog2color[group]
            gene_uid += 1
            scribl += add_gene(locus_member,
                               gene_uid,
                               track_uid,
                               offset,
                               reverse,
                               gene_color,
                               group,
                               request,
                               locus_size
                               )

    return '\n'.join(scribl), tree_canvas, tree_newick, og_gene_count, plot_gene_count
