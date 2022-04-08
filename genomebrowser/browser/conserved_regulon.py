from .models import *
from collections import defaultdict
from browser.dataimport.annotator import autovivify


def get_lowest_level_og(eggnog_groups):
    result = None
    RANKS = ('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'norank')
    eggnog_ranks = {}
    for eggnog_group in eggnog_groups.all():
        eggnog_ranks[eggnog_group.taxon.rank] = eggnog_group
    for rank in RANKS:
        if rank in eggnog_ranks:
            result = eggnog_group
            break
    return result
    

def build_conserved_regulon(og_id):
    context = {}
    try:
        og = Ortholog_group.objects.get(id=og_id)
    except Ortholog_group.DoesNotExist:
        raise Http404('Ortholog group not found: ' + og_id)
    
    regulators = Gene.objects.filter(protein__ortholog_groups__id=og.id).select_related('protein', 'genome', 'genome__strain__taxon')
    regulons = Regulon.objects.filter(regulators__id__in=[item.id for item in regulators]).select_related('genome')
    sites = Site.objects.filter(regulon__id__in=[item.id for item in regulons]).select_related('genome', 'contig', 'regulon').prefetch_related('genes', 'operons')
    context['regulator_og'] = og
    context['regulators'] = regulators
    context['regulons'] = regulons
    context['sites'] = sites
    
    eggnogs = autovivify(2, list)
    regulon_labels = set()
    eggnog_labels = defaultdict(dict)
    for site in sites:
        regulon_label = site.regulon.name + '[' + site.genome.name + ']'
        regulon_labels.add(regulon_label)
        for gene in site.genes.all():
            if not gene.protein:
                continue
            eggnog_group = get_lowest_level_og(gene.protein.ortholog_groups.select_related('taxon'))
            if eggnog_group is not None:
                eggnogs[eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'][regulon_label].append(gene.locus_tag)
                eggnog_labels[eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'][gene.function] = ''
        for operon in site.operons.all():
            for gene in operon.genes.all():
                if not gene.protein:
                    continue
                eggnog_group = get_lowest_level_og(gene.protein.ortholog_groups.select_related('taxon'))
                if eggnog_group is not None:
                    eggnogs[eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'][regulon_label].append(gene.locus_tag)
                    eggnog_labels[eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'][gene.function] = ''
    
    context['gene_table_header'] = sorted(list(regulon_labels))
    
    gene_table = []
    gene_table.append('			<table>')
    gene_table.append('				<thead>')
    gene_table.append('					<tr>')
    gene_table.append('						<th>EggNOG family</th>')
    for regulon_label in sorted(list(regulon_labels)):
        gene_table.append('						<th>' + regulon_label + '</th>')
    gene_table.append('						<th>Function</th>')
    gene_table.append('					</tr>')
    gene_table.append('				</thead>')
    gene_table.append('				<tbody>')
    for eggnog_label in eggnogs:
        gene_table.append('					<tr>')
        gene_table.append('					    <td>' + eggnog_label + '</td>')
        gene_functions = set()
        for regulon_label in sorted(list(regulon_labels)):
            if regulon_label in eggnogs[eggnog_label]:
                gene_table.append('					    <td>' + ';'.join(eggnogs[eggnog_label][regulon_label]) + '</td>')
            else:
                gene_table.append('					    <td>&nbsp;</td>')
        gene_table.append('					    <td>' + ';'.join(sorted(eggnog_labels[eggnog_label].keys())) + '</td>')
        gene_table.append('					</tr>')
    gene_table.append('				</tbody>')
    gene_table.append('			</table>')
    context['gene_table'] = '\n'.join(gene_table)
    return context
