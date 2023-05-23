""" Module supporting backbone command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput
from biobb_structure_checking.structure_manager import NotEnoughAtomsError


def check(strcheck):
    if not strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['NO_CA_ONLY_STRUCTURE'])
        return None
    return True

def fix(strcheck, opts, fix_data=None):

    no_int_recheck = opts['fix'] is not None or\
        strcheck.args['non_interactive']

    if isinstance(opts, str):
        fix = opts
    else:
        fix = opts['fix']

    input_line = ParamInput(
        'Build protein',
        strcheck.args['non_interactive'],
        set_none='no'
    )
    input_line.add_option_yes_no()
    input_line.set_default('yes')
    input_option, fix = input_line.run(fix)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], fix

    if input_option == 'yes':
        print("Build full structure")

        # Checking for canonical sequence
        if strcheck.strucm.models_data.has_models():
            sys.err("Fix CA-only not implemented for Models, exiting")
        if not strcheck.strucm.sequence_data.has_canonical[0]:
            read_ok = False
            if not strcheck.args['non_interactive']:
                while not read_ok:
                    input_line = ParamInput(
                        "Enter canonical sequence path (FASTA)",
                        strcheck.args['non_interactive']
                    )
                    strcheck.args['fasta_seq_path'] = input_line.run(
                        strcheck.args['fasta_seq_path']
                    )
                    if not strcheck.args['fasta_seq_path']:
                        print(cts.MSGS['FASTA_MISSING'])

                    read_ok = strcheck.strucm.sequence_data.load_sequence_from_fasta(
                        strcheck.args['fasta_seq_path']
                    )
            if not read_ok:
                strcheck.strucm.sequence_data.fake_canonical_sequence(strcheck.strucm)
            else:
                strcheck.strucm.sequence_data.read_canonical_seqs(
                    strcheck.strucm,
                    False
                )
            strcheck.strucm.sequence_data.match_sequence_numbering(strcheck.strucm)

        return strcheck.strucm.fix_ca_only(strcheck.args['modeller_key'])

