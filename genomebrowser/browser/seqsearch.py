import os
import uuid
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from Bio import Alphabet
#from Bio import SeqIO
from subprocess import Popen, PIPE, CalledProcessError, STDOUT
from genomebrowser.settings import STATICFILES_DIRS
from browser.models import Config
#log_file = os.path.join(STATICFILES_DIRS[0], 'genomes', 'tmp', 'search.log')


def _verify_alphabet(sequence, alphabet):
    alphabet = set(alphabet) 
    return all(letter in alphabet for letter in sequence)
    
#def validate_nucl(query):
#    query_lines = query.split('\n')
#    seq_record = SeqRecord(Seq(''.join([x.rstrip('\n\r') for x in query_lines[1:]])), id=query_lines[0][1:].rstrip('\r\n'))
#    _verify_alphabet(seq_record.seq)
#    return str(seq_record.seq), str(seq_record.id)


def run_protein_search(query):
    search_dir = Config.objects.get(param='cgcms.search_db_dir').value
    PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
    log_file = os.path.join(search_dir, 'search.log')
    blast_db = os.path.join(search_dir, 'blast_prot')# os.path.join(STATICFILES_DIRS[0], 'genomes', 'search', 'blast_prot')
    result = []
    searchcontext = ''
    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] New BLASTP search started.\n')
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Query sequence:\"'+ query + '\"\n')
    query_lines = query.split('\n')
    seq_record = SeqRecord(Seq(''.join([x.rstrip('\n\r') for x in query_lines[1:]])), id=query_lines[0][1:].rstrip('\r\n'))
    if not _verify_alphabet(seq_record.seq.upper(), PROTEIN_ALPHABET):
        searchcontext = 'Wrong protein sequence format. FASTA header and valid sequence required.'
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)
    if not sequence:
        searchcontext = 'Wrong sequence format. FASTA header and valid sequence required.'
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext
    args = [
        'blastp',
        '-db',
        blast_db,
        '-max_target_seqs',
        '50',
        '-evalue',
        '0.00001',
        '-matrix=PAM30',
        '-outfmt',
        '6'
        ]
    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Run BLASTP with args:\"'+ ' '.join(args) + '\"\n')
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        blastoutput, err = p.communicate(query.strip())
    if p.returncode != 0:
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] BLASTP finished with error:\n'+ ' '.join(err) + '\n')
        searchcontext = 'BLASTP finished with error:\n' + '\n'.join(err)
        return result, searchcontext
    for line in blastoutput.split('\n'):
        if line.startswith('#'):
            continue
        row = line.rstrip('\n\r').split('\t')
        if len(row) < 12:
            continue
        if not sequence_id.startswith(row[0]):
            continue
        result.append('\t'.join(row))
    if not result:
        searchcontext = 'No hits found'
    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] BLASTP finished. ' + str(len(result)) + ' hits found.\n')
    if len(result) > 100:
        result = result[:100]
    return result, searchcontext


def run_nucleotide_search(query):
    search_dir = Config.objects.get(param='cgcms.search_db_dir').value
    DNA_ALPHABET = 'GATCRYWSMKHBVDN'
    log_file = os.path.join(search_dir, 'search.log')
    blast_db = os.path.join(search_dir, 'blast_nucl') #os.path.join(STATICFILES_DIRS[0], 'genomes', 'search', 'blast_nucl')
    result = []
    searchcontext = ''

    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] New megablast search started.\n')
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Query sequence:\"'+ query + '\"\n')

    query_lines = query.split('\n')
    seq_record = SeqRecord(Seq(''.join([x.rstrip('\n\r') for x in query_lines[1:]])), id=query_lines[0][1:].rstrip('\r\n'))
    if not _verify_alphabet(seq_record.seq.upper(), DNA_ALPHABET):
        searchcontext = 'Wrong nucleotide sequence format. FASTA header and valid sequence required. Multiple entries not supported.'
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)

    if not sequence:
        searchcontext = 'Wrong sequence format. FASTA header and sequence required. Multiple entries not supported.'
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext
    args = [
        'megablast',
        '-a', '6',
        '-b', '100',
        '-D', '3',
        '-e', '0.001',
        '-f', 'T',
        '-d', blast_db
        ]
    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Run megablast with args:\"'+ ' '.join(args) + '\"\n')
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        blastoutput, err = p.communicate(query.strip())
    if p.returncode != 0:
        with open(log_file, 'a') as log:
            log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] megablast finished with error:\n'+ ' '.join(err) + '\n')
        searchcontext = 'Megablast finished with error:\n' + '\n'.join(err)
        return result, searchcontext
    for line in blastoutput.split('\n'):
        if line.startswith('#'):
            continue
        row = line.rstrip('\n\r').split('\t')
        if len(row) < 12:
            continue
        if not sequence_id.startswith(row[0]):
            continue
        result.append('\t'.join(row))
    if not result:
        searchcontext = 'No hits found'
    with open(log_file, 'a') as log:
        log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] megablast finished. ' + str(len(result)) + ' hits found.\n')
    if len(result) > 100:
        result = result[:100]
    return result, searchcontext
