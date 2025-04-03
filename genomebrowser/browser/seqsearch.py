import os
import re
import logging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import Popen, PIPE, STDOUT
from django.core.exceptions import SuspiciousOperation
from browser.models import Config

logger = logging.getLogger("GenomeDepot")

def _verify_alphabet(sequence, alphabet):
    '''
        Checks if sequence contains illegal symbols
    '''
    alphabet = set(alphabet) 
    return all(letter in alphabet for letter in sequence)

def _sanitize_sequence(query):
    '''
        Remove all symbols except a-zA-Z
    '''
    query_lines = query.split('\n')
    if query.startswith('>'):
        sequence = ''.join([x.rstrip('\n\r') for x in query_lines[1:]])
        query_id = query_lines[0][1:].rstrip('\r\n')
    else:
        sequence = ''.join([x.rstrip('\n\r') for x in query_lines])
        query_id = 'no_sequence_id'
    sequence = ''.join([i if ord(i) < 128 else '' for i in sequence])
    sequence = re.sub(r"[^^a-zA-Z]", '', sequence)
    return query_id, sequence
    
def validate_params(params):
    '''
        Validates search parameters received from user.
        Raises SuspiciousOperation exception if the validation fails.
    '''
    result = {}
    if isinstance(params,str):
        result['sequence'] = params
        result['evalue'] = '0.0001'
        result['hitstoshow'] = '100'
    else:
        result['sequence'] = params['sequence']
        try:
            evalue = params['evalue']
        except KeyError:
            raise SuspiciousOperation("e-value parameter is missing")
        if evalue is None:
            raise SuspiciousOperation("e-value parameter is missing")
        try:
            result['evalue'] = str(float(evalue))
        except ValueError:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for e-value parameter." % result['evalue']
            )
        except TypeError:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for e-value parameter." % result['evalue']
            )
        if result['evalue'] not in ['1e-20',
                                    '1e-10',
                                    '1e-08',
                                    '1e-06',
                                    '0.0001',
                                    '0.01',
                                    '1.0',
                                    '10.0'
                                    ]:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for e-value parameter." % result['evalue']
            )
        try:
            hitstoshow = params['hitstoshow']
        except KeyError:
            raise SuspiciousOperation("hitstoshow parameter is missing")
        if hitstoshow is None:
            raise SuspiciousOperation("e-value parameter is missing")
        try:
            result['hitstoshow'] = str(int(hitstoshow))
        except TypeError:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for hitstoshow parameter." % result['hitstoshow']
            )
        except ValueError:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for hitstoshow parameter." % result['hitstoshow']
            )
        if result['hitstoshow'] not in ['10', '20', '50', '100', '500', '1000']:
            raise SuspiciousOperation(
            "Unacceptable value '%s' for hitstoshow parameter." % result['hitstoshow']
            )
    return result

        
def run_protein_search(params):
    '''
        Runs BLASTP search for protein sequence in params.
    '''
    result = []
    try:
        params = validate_params(params)
    except SuspiciousOperation as e:
        searchcontext = str(e)
        logger.debug('Parameters: %s', str(params))
        return result, searchcontext, 0, ''
    query = params['sequence']
    search_dir = Config.objects.get(param='core.search_db_dir').value
    PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
    blast_db = os.path.join(search_dir, 'blast_prot')
    searchcontext = ''
    query_id, query_sequence = _sanitize_sequence(query)
    seq_record = SeqRecord(Seq(query_sequence), id=query_id)
    if not _verify_alphabet(seq_record.seq.upper(), PROTEIN_ALPHABET):
        searchcontext = 'Wrong protein sequence format. ' +\
                        'FASTA header and valid sequence required.'
        return result, searchcontext, 0, query_id
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)

    if not sequence:
        searchcontext = 'Wrong sequence format. ' +\
                        'FASTA header and valid sequence required.'
        return result, searchcontext, 0, sequence_id
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
    with Popen(args,
               stdin=PIPE,
               stdout=PIPE,
               stderr=PIPE,
               bufsize=1,
               universal_newlines=True
               ) as p:
        blastoutput, err = p.communicate('>' + sequence_id + '\n' + sequence)
    if p.returncode != 0:
        searchcontext = 'BLASTP finished with error:\n' + ''.join(err)
        logger.error('BLASTP finished with error. Parameters: %s \n %s',
                     str(params), ''.join(err)
                     )
        return result, searchcontext, 0, sequence_id
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
    if len(result) > int(params['hitstoshow']):
        result = result[:int(params['hitstoshow'])]
    return result, searchcontext, query_len, sequence_id

def run_nucleotide_search(params):
    '''
        Runs Megablast search for nucleotide sequence in params.
    '''
    result = []
    try:
        params = validate_params(params)
    except SuspiciousOperation as e:
        searchcontext = str(e)
        logger.debug('Parameters: %s', str(params))
        return result, searchcontext, 0, ''
    query = params['sequence']
    
    search_dir = Config.objects.get(param='core.search_db_dir').value
    DNA_ALPHABET = 'GATCRYWSMKHBVDN'
    blast_db = os.path.join(search_dir, 'blast_nucl')
    searchcontext = ''
    query_id, query_sequence = _sanitize_sequence(query)
    seq_record = SeqRecord(Seq(query_sequence), id=query_id)
    if not _verify_alphabet(seq_record.seq.upper(), DNA_ALPHABET):
        searchcontext = 'Wrong nucleotide sequence format. ' +\
                        'FASTA header and valid sequence required. ' +\
                        'Multiple entries not supported.'
        return result, searchcontext, 0, query_id
    sequence = str(seq_record.seq)
    sequence_id = str(seq_record.id)
    query_len = len(seq_record)

    if not sequence:
        searchcontext = 'Wrong sequence format. FASTA header and sequence required.' +\
                        'Multiple entries not supported.'
        return result, searchcontext, 0, sequence_id
    args = [
        'blastn',
        '-task',  'megablast',
        '-num_threads', '2',
        '-max_target_seqs', params['hitstoshow'],
        '-evalue', params['evalue'],
        '-db', blast_db,
        '-outfmt', '6'
        ]
    print(' '.join(args))

    with Popen(args,
               stdin=PIPE,
               stdout=PIPE,
               stderr=PIPE,
               bufsize=1,
               universal_newlines=True
               ) as p:
        blastoutput, err = p.communicate('>' + sequence_id + '\n' + sequence)
    if p.returncode != 0:
        searchcontext = 'Megablast finished with error:\n' + ''.join(err)
        logger.error('Megablast finished with error. Parameters: %s \n %s',
                     str(params), ''.join(err)
                     )
        return result, searchcontext, 0, sequence_id
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
    if len(result) > int(params['hitstoshow']):
        result = result[:int(params['hitstoshow'])]
    return result, searchcontext, query_len, sequence_id
