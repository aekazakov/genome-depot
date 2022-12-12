import os
import uuid
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import Popen, PIPE, CalledProcessError, STDOUT
from django.core.exceptions import SuspiciousOperation
from browser.models import Config


def _verify_alphabet(sequence, alphabet):
    alphabet = set(alphabet) 
    return all(letter in alphabet for letter in sequence)

def validate_params(params):
    result = {}
    if isinstance(params,str):
        result['sequence'] = params
        result['evalue'] = '0.0001'
        result['hitstoshow'] = '100'
    else:
        result['sequence'] = params['sequence']
        try:
            result['evalue'] = str(float(params['evalue']))
        except TypeError:
            raise SuspiciousOperation("Unacceptable value '%s' for e-value parameter." % result['evalue'])
        if result['evalue'] not in ['1e-20', '1e-10', '1e-08', '1e-06', '0.0001', '0.01', '1.0', '10.0']:
            raise SuspiciousOperation("Unacceptable value '%s' for e-value parameter." % result['evalue'])
        result['hitstoshow'] = str(int(params['hitstoshow']))
        if result['hitstoshow'] not in ['10', '20', '50', '100', '500', '1000']:
            raise SuspiciousOperation("Unacceptable value '%s' for hitstoshow parameter." % result['hitstoshow'])
    return result

        
def run_protein_search(params):
    result = []
    try:
        params = validate_params(params)
    except SuspiciousOperation as e:
        searchcontext = str(e)
        return result, searchcontext, 0
    query = params['sequence']

    search_dir = Config.objects.get(param='cgcms.search_db_dir').value
    PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
    #log_file = os.path.join(search_dir, 'search.log')
    blast_db = os.path.join(search_dir, 'blast_prot')
    searchcontext = ''
    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] New BLASTP search started.\n')
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Query sequence:\"'+ query + '\"\n')
    query_lines = query.split('\n')
    seq_record = SeqRecord(Seq(''.join([x.rstrip('\n\r') for x in query_lines[1:]])), id=query_lines[0][1:].rstrip('\r\n'))
    if not _verify_alphabet(seq_record.seq.upper(), PROTEIN_ALPHABET):
        searchcontext = 'Wrong protein sequence format. FASTA header and valid sequence required.'
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext, 0
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)
    if not sequence:
        searchcontext = 'Wrong sequence format. FASTA header and valid sequence required.'
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext, 0
    query_len = len(seq_record)
    args = [
        'blastp',
        '-db',
        blast_db,
        '-max_target_seqs', params['hitstoshow'],
        '-evalue', params['evalue'],
        '-matrix=PAM30',
        '-outfmt',
        '6'
        ]
    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Run BLASTP with args:\"'+ ' '.join(args) + '\"\n')
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        blastoutput, err = p.communicate(query.strip())
    if p.returncode != 0:
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] BLASTP finished with error:\n'+ ' '.join(err) + '\n')
        searchcontext = 'BLASTP finished with error:\n' + '\n'.join(err)
        return result, searchcontext, 0
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
    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] BLASTP finished. ' + str(len(result)) + ' hits found.\n')
    if len(result) > 100:
        result = result[:100]
    return result, searchcontext, query_len

def run_nucleotide_search(params):
    result = []
    try:
        params = validate_params(params)
    except SuspiciousOperation as e:
        searchcontext = str(e)
        return result, searchcontext, 0
    query = params['sequence']
    
    search_dir = Config.objects.get(param='cgcms.search_db_dir').value
    DNA_ALPHABET = 'GATCRYWSMKHBVDN'
    #log_file = os.path.join(search_dir, 'search.log')
    blast_db = os.path.join(search_dir, 'blast_nucl')
    searchcontext = ''

    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] New megablast search started.\n')
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Query sequence:\"'+ query + '\"\n')

    query_lines = query.split('\n')
    seq_record = SeqRecord(Seq(''.join([x.rstrip('\n\r') for x in query_lines[1:]])), id=query_lines[0][1:].rstrip('\r\n'))
    if not _verify_alphabet(seq_record.seq.upper(), DNA_ALPHABET):
        searchcontext = 'Wrong nucleotide sequence format. FASTA header and valid sequence required. Multiple entries not supported.'
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext, 0
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)
    query_len = len(seq_record)

    if not sequence:
        searchcontext = 'Wrong sequence format. FASTA header and sequence required. Multiple entries not supported.'
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Wrong sequence format.\n')
        return result, searchcontext, 0
    args = [
        'megablast',
        '-a', '2',
        '-b', params['hitstoshow'],
        '-D', '3',
        '-e', params['evalue'],
        '-f', 'T',
        '-d', blast_db
        ]
    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] Run megablast with args:\"'+ ' '.join(args) + '\"\n')
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        blastoutput, err = p.communicate(query.strip())
    if p.returncode != 0:
        #with open(log_file, 'a') as log:
        #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] megablast finished with error:\n'+ ' '.join(err) + '\n')
        searchcontext = 'Megablast finished with error:\n' + '\n'.join(err)
        return result, searchcontext, 0
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
    #with open(log_file, 'a') as log:
    #    log.write('[' + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '] megablast finished. ' + str(len(result)) + ' hits found.\n')
    if len(result) > 100:
        result = result[:100]
    return result, searchcontext, query_len
