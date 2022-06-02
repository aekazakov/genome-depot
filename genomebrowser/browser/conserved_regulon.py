from .models import *
from collections import defaultdict
from browser.dataimport.annotator import autovivify
from django.urls import reverse


def get_lowest_level_og(protein):
    result = None
    RANKS = ('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'norank')
    eggnog_ranks = {}
    for eggnog_group in protein.ortholog_groups.all():
        eggnog_ranks[eggnog_group.taxon.rank] = eggnog_group
    for rank in RANKS:
        if rank in eggnog_ranks:
            result = eggnog_ranks[rank]
            break
    return result
    

def build_conserved_regulon(og_id):
    context = {}
    try:
        og = Ortholog_group.objects.get(id=og_id)
    except Ortholog_group.DoesNotExist:
        raise Http404('Ortholog group not found: ' + og_id)
    
    regulators = Gene.objects.filter(protein__ortholog_groups__id=og.id).select_related('protein', 'genome')
    regulons = Regulon.objects.filter(regulators__id__in=[item.id for item in regulators]).select_related('genome', 'genome__taxon').prefetch_related('regulators', 'regulators__genome')
    sites = Site.objects.filter(regulon__id__in=[item.id for item in regulons]).select_related('genome', 'contig', 'regulon').prefetch_related('genes', 'genes__genome', 'operons', 'operons__genome', 'operons__contig', 'operons__genes', 'operons__genes__genome', 'operons__genes__protein', 'operons__genes__protein__ortholog_groups', 'operons__genes__protein__ortholog_groups__taxon', 'genes__protein', 'genes__protein__ortholog_groups', 'genes__protein__ortholog_groups__taxon')
    context['regulator_og'] = og
    context['regulators'] = regulators
    context['regulons'] = regulons
    context['sites'] = sites
    
    eggnogs = autovivify(2, list)
    regulon_labels = set()
    eggnog_labels = {}
    eggnog_functions = defaultdict(dict)
    genes_added = set()
    operons_added = set()
    fake_id = 0
    
    for site in sites:
        regulon_label = site.regulon.name + '[' + site.genome.name + ']'
        regulon_labels.add(regulon_label)
        for gene in site.genes.all():
            if gene.id in genes_added:
                continue
            genes_added.add(gene.id)
            if gene.protein is None:
                current_id = str(fake_id)
                fake_id += 1
                eggnogs[current_id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                eggnog_functions[current_id][gene.function] = ''
                eggnog_labels[current_id] = 'No EggNOG ID'
            else:
                eggnog_group = get_lowest_level_og(gene.protein)
                if eggnog_group is not None:
                    eggnogs[eggnog_group.id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                    eggnog_functions[eggnog_group.id][gene.function] = ''
                    eggnog_labels[eggnog_group.id] = '<a href="' + reverse('searchgene') + "?og=" + str(eggnog_group.id) + '">' + eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']</a>'
                else:
                    current_id = str(fake_id)
                    fake_id += 1
                    eggnogs[current_id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                    eggnog_functions[current_id][gene.function] = ''
                    eggnog_labels[current_id] = 'No EggNOG ID'
        for operon in site.operons.all():
            if operon.id in operons_added:
                continue
            operons_added.add(operon.id)
            for gene in operon.genes.all():
                if gene.id in genes_added:
                    continue
                genes_added.add(gene.id)
                if not gene.protein:
                    current_id = str(fake_id)
                    fake_id += 1
                    eggnogs[current_id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                    eggnog_functions[current_id][gene.function] = ''
                    eggnog_labels[current_id] = 'No EggNOG ID'
                else:
                    eggnog_group = get_lowest_level_og(gene.protein)
                    if eggnog_group is not None:
                        eggnogs[eggnog_group.id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                        eggnog_functions[eggnog_group.id][gene.function] = ''
                        eggnog_labels[eggnog_group.id] = '<a href="' + reverse('searchgene') + "?og=" + str(eggnog_group.id) + '">' + eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']</a>'
                    else:
                        current_id = str(fake_id)
                        fake_id += 1
                        eggnogs[current_id][regulon_label].append('<a href="' + reverse('genedetails', args=(gene.genome.name, gene.locus_tag)) + '" title="' + gene.locus_tag + '">')
                        eggnog_functions[current_id][gene.function] = ''
                        eggnog_labels[current_id] = 'No EggNOG ID'

    context['gene_table_header'] = sorted(list(regulon_labels))
    gene_table_rows = []
    
    for eggnog_id in eggnogs:
        table_row = ['<td>' + eggnog_labels[eggnog_id] + '</td>']
        for regulon_label in sorted(list(regulon_labels)):
            if regulon_label in eggnogs[eggnog_id]:
                table_row.append('<td class="filledcell">' + '<br/>'.join([item + str(index + 1) + '</a>' for index, item in enumerate(eggnogs[eggnog_id][regulon_label])]) + '</td>')
            else:
                table_row.append('<td>&nbsp;</td>')
        table_row.append('<td>' + '; '.join(sorted(eggnog_functions[eggnog_id].keys())) + '</td>')
        gene_table_rows.append(table_row)
    if gene_table_rows:
        context['gene_table_rows'] = gene_table_rows

    return context
