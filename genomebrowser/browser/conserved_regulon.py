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
    regulator_functions = set()
    eggnog_labels = {}
    eggnog_functions = defaultdict(dict)
    genes_added = set()
    operons_added = set()
    fake_id = 0
    
    for item in regulators:
        regulator_functions.add(item.function)
    context['regulator_functions'] = sorted(list(regulator_functions))
    
    for site in sites:
        regulon_label = (site.regulon.name, site.genome.name)
        regulon_labels.add(regulon_label)
        for gene in site.genes.all():
            if gene.id in genes_added:
                continue
            genes_added.add(gene.id)
            if gene.protein is None:
                current_id = str(fake_id)
                fake_id -= 1
                eggnogs[current_id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                if gene.function not in eggnog_functions[current_id]:
                    eggnog_functions[current_id][gene.function] = 0
                eggnog_functions[current_id][gene.function] += 1
            else:
                eggnog_group = get_lowest_level_og(gene.protein)
                if eggnog_group is not None:
                    eggnogs[eggnog_group.id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                    if gene.function not in eggnog_functions[eggnog_group.id]:
                        eggnog_functions[eggnog_group.id][gene.function] = 0
                    eggnog_functions[eggnog_group.id][gene.function] += 1
                    eggnog_labels[eggnog_group.id] = eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'
                else:
                    current_id = str(fake_id)
                    fake_id -= 1
                    eggnogs[current_id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                    if gene.function not in eggnog_functions[current_id]:
                        eggnog_functions[current_id][gene.function] = 0
                    eggnog_functions[current_id][gene.function] += 1
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
                    fake_id -= 1
                    eggnogs[current_id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                    if gene.function not in eggnog_functions[current_id]:
                        eggnog_functions[current_id][gene.function] = 0
                    eggnog_functions[current_id][gene.function] += 1
                else:
                    eggnog_group = get_lowest_level_og(gene.protein)
                    if eggnog_group is not None:
                        eggnogs[eggnog_group.id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                        if gene.function not in eggnog_functions[eggnog_group.id]:
                            eggnog_functions[eggnog_group.id][gene.function] = 0
                        eggnog_functions[eggnog_group.id][gene.function] += 1
                        eggnog_labels[eggnog_group.id] = eggnog_group.eggnog_id + '[' + eggnog_group.taxon.name + ']'
                    else:
                        current_id = str(fake_id)
                        fake_id -= 1
                        eggnogs[current_id][regulon_label].append((gene.locus_tag, gene.genome.name, gene.function))
                        if gene.function not in eggnog_functions[current_id]:
                            eggnog_functions[current_id][gene.function] = 0
                        eggnog_functions[current_id][gene.function] += 1

    context['gene_table_header'] = sorted(list(regulon_labels))
    gene_table_rows = []
    
    for eggnog_id in eggnogs:
        # Add first cell in a row
        max_function_count = 0
        max_function = 'Unknown function'
        for gene_function in eggnog_functions[eggnog_id]:
            if eggnog_functions[eggnog_id][gene_function] > max_function_count:
                max_function_count = eggnog_functions[eggnog_id][gene_function]
                max_function = gene_function

        if eggnog_id in eggnog_labels:
            table_row = ['<td class="sticky-col style3"><a href="' + reverse('searchgene') + "?type=og&query=" + str(eggnog_id) + '" title="' + eggnog_labels[eggnog_id] + '">' + max_function + '</a>' + '</td>']
        else:
            # Do not show eggnog_id in the page, it is fake
            table_row = ['<td class="sticky-col style3"><a href="javascript:alert(\'No EggNOG ID.\');" title="No EggNOG ID">' + max_function + '</a>' + '</td>']
        
        for regulon_label in sorted(list(regulon_labels)):
            if regulon_label in eggnogs[eggnog_id]:
                tooltip = []
                for index,item in enumerate(eggnogs[eggnog_id][regulon_label]):
                    tooltip.append(str(index + 1) + '. Gene: <a href="' + reverse('genedetails', args=(item[1], item[0])) + '" title="' + item[0] + ': ' + item[2] + '">' + item[0] + '</a> [' + item[1] + ']')
                    tooltip.append(item[2])
                table_row.append('<td class="filledcell"><div class="tooltip">' + str(len(eggnogs[eggnog_id][regulon_label])) + '<span>' + '<br/>'.join(tooltip) + '</span></div></td>')
            else:
                table_row.append('<td>&nbsp;</td>')
        #table_row.append('<td>' + '; '.join(sorted(eggnog_functions[eggnog_id].keys())) + '</td>')
        gene_table_rows.append(table_row)
    if gene_table_rows:
        context['gene_table_rows'] = gene_table_rows

    return context
