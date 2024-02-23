import os
import gzip
from collections import defaultdict
from Bio import SeqIO
from browser.models import Genome
from browser.models import Gene
from browser.models import Protein


def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def export_proteins_bygenome(genomes, out_dir):
    """
    Makes protein FASTA files, one file for each input genome

    Parameters:
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
    
    Parameters:
    genomes(dict<str:str>): dictionary with genome name
    as key and GBK path as value
    out_dir(str): output directory

    """
    ret = {}
    for genome_name, gbkfile in genomes.items():
        genome_id = Genome.objects.filter(name=genome_name).values('id')[0]['id']
        outfasta = os.path.join(out_dir, str(genome_id) + '.fna')
        with open(outfasta, 'w') as outfile:
            if gbkfile.endswith('.gz'):
                fh = gzip.open(gbkfile, 'rt')
            else:
                fh = open(gbkfile, 'r')
            for seq_record in SeqIO.parse(fh, "genbank"):
                outfile.write('>' + seq_record.id + '\n')
                outfile.write(str(seq_record.seq) + '\n')
            fh.close()
        ret[genome_name] = outfasta
    return ret


def export_proteins(genome_ids, out_filename):
    """
    Populates proteins dictionary and creates protein FASTA file
    The resulting FASTA file is non-redundant: each protein sequence
    is written only once. Sequence headers are protein hashes
    
    Parameters:
    genome_names(list of int): list of genome ids
    out_filename(str): output file path
    """
    with open(out_filename, 'w') as outfile:
        if genome_ids is None:
            for protein in Protein.objects.all():
                outfile.write('>' + protein.protein_hash + '\n')
                outfile.write(protein.sequence + '\n')
        else:
            proteins = {}
            target_genes = Gene.objects.filter(
                                               genome__id__in = genome_ids
                                               ).select_related('protein')
            for gene in target_genes:
                if gene.protein is not None:
                    proteins[gene.protein.protein_hash] = gene.protein
        
            for protein_hash, protein in proteins.items():
                outfile.write('>' + protein_hash + '\n')
                outfile.write(protein.sequence + '\n')
