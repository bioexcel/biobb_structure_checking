""" Command sequences """
import biobb_structure_checking.constants as cts

def check(strcheck):
    if strcheck.strucm.sequence_data.has_canonical:
        print('Canonical sequence')
        can_seq = strcheck.strucm.sequence_data.get_canonical()
        print(can_seq)
    else:
        print(cts.MSGS['NO_CANONICAL'])
        can_seq = ''

    pdb_seq = strcheck.strucm.sequence_data.get_pdbseq()
    print('Structure sequence')
    print(pdb_seq)
    strcheck.summary['FASTA'] = {
        'canonical': can_seq,
        'structure': pdb_seq
    }

    return {}

def fix(strcheck, opts, fix_data=None):
    pass