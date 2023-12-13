import os
import gzip
from collections import defaultdict
from Bio import SeqIO
from browser.models import Genome
from browser.models import Gene


def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def export_proteins_bygenome(genomes, out_dir):
    """
    Makes protein FASTA files, one file for each input genome

    genomes(dict<str:str>): dictionary with genome name
    as key and GBK path as value
    out_dir(str): output directory
    
    """
    ret = {}
    for genome_name in genomes.keys():
        genome_id = Genome.objects.filter(name=genome_name).values('id')[0]['id']
        outfasta = os.path.join(out_dir, str(genome_id) + '.faa') 
        with open(outfasta, 'w') as outfile:
            target_genes = Gene.objects.filter(
                genome__id = genome_id
                ).select_related('protein')
            for gene in target_genes:
                if gene.protein is not None:
                    outfile.write('>' + gene.locus_tag + '\n')
                    outfile.write(gene.protein.sequence + '\n')
        ret[genome_name] = outfasta
    return ret
                    

def export_nucl_bygenome(genomes, out_dir):
    """
    Makes nucleotide FASTA files, one file for each input genome
    
    genomes(dict<str:str>): dictionary with genome name
    as key and GBK path as value
    out_dir(str): output directory

    """
    ret = {}
    for genome_name, gbkfile in genomes.items():
        genome_id = Genome.objects.filter(name=genome_name).values('id')[0]['id']
        outfasta = os.path.join(out_dir, str(genome_id) + '.faa')
        with open(outfasta, 'w') as outfile:
            if gbkfile.endswith('.gz'):
                fh = gzip.open(gbkfile, 'rt')
            else:
                fh = open(gbkfile, 'r')
            for seq_record in SeqIO.parse(fh, "genbank"):
                outfasta.write('>' + seq_record.id + '\n')
                outfasta.write(str(seq_record.seq) + '\n')
            fh.close()
        ret[genome_name] = outfasta
    return ret
