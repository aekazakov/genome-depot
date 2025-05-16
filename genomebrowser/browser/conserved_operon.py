from browser.models import Gene
from browser.conserved_regulon import RANKS
from browser.treemap import generate_genes_treemap
from browser.taxonomy import generate_operons_sunburst

def build_conserved_operon(operon_id):
    result = ''
    operon_members = Gene.objects.filter(operon__id=operon_id)
    ortholog_groups = set()
    for gene in operon_members.all():
        print(gene.locus_tag)
        if gene.protein:
            gene_ogs = {og.taxon.rank:og.id for og
                in gene.protein.ortholog_groups.all()
                }
            print(gene_ogs)
            for rank in RANKS:
                if rank in gene_ogs:
                    ortholog_groups.add(gene_ogs[rank])
                    print(rank, gene_ogs[rank])
                    break
    if not ortholog_groups:
        return result
    else:
        ortholog_groups = list(ortholog_groups)
    operons = set()
    monocistronic_genes = set()
    for gene in Gene.objects.filter(
        protein__ortholog_groups__in = ortholog_groups
    ).select_related(
        'protein',
        'operon'
    ):
        if gene.operon:
            operons.add(gene.operon.id)
        else:
            monocistronic_genes.add(gene.id)
        
    all_gene_ids = list(monocistronic_genes) + \
        list(Gene.objects.filter(
            operon__id__in = operons
        ).values_list('id', flat=True))

    print(all_gene_ids)
    treemap, functional_profile_tsv = generate_genes_treemap(all_gene_ids)
    
    sunburst = generate_operons_sunburst(list(operons), list(monocistronic_genes))
    print('sunburst:', sunburst)
    return treemap, sunburst, list(operons), functional_profile_tsv
