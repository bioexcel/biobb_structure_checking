""" Command sequences """
import biobb_structure_checking.constants as cts


def check(strcheck):
    fasta = ''
    if strcheck.strucm.sequence_data.has_canonical[0]:
        print('Canonical sequence')
        can_seq = strcheck.strucm.sequence_data.get_canonical()
        print(can_seq)
        fasta += f"{can_seq}\n".replace('\n\n', '\n')
    else:
        print(cts.MSGS['NO_CANONICAL'])
        can_seq = ''

    pdb_seq = strcheck.strucm.sequence_data.get_pdbseq()
    print('Structure sequence')
    print(pdb_seq)
    fasta += f"{pdb_seq}".replace('\n\n', '\n')
    strcheck.summary['FASTA'] = {
        'canonical': can_seq,
        'structure': pdb_seq
    }
    return fasta


def fix(strcheck, opts, fix_data=None):
    if opts['output_fasta']:
        try:
            with open(opts['output_fasta'], 'w') as fasta:
                fasta.write(fix_data)
            print(f"Sequences written on {opts['output_fasta']}")
        except Exception:
            print(
                cts.MSGS['NO_WRITE_FASTA'].format(strcheck.args['output_fasta'])
            )
