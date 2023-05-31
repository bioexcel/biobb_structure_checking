""" Module supporting amide command"""

import numpy as np
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput
from biobb_structure_checking.structure_manager import NotAValidResidueError


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    amide_check = strcheck.strucm.check_amide_contacts()
    if 'list' not in amide_check:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_AMIDES'])
        return {}
    strcheck.summary['amide']['n_amides'] = len(amide_check['list'])

    if not amide_check['cont_list']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_UNUSUAL_AMIDES'])
        return {}

    print(cts.MSGS['UNUSUAL_AMIDES'].format(len(amide_check['cont_list'])))

    strcheck.summary['amide']['detected'] = []
    for at_pair in sorted(
        amide_check['cont_list'],
        key=mu.key_sort_atom_pairs
    ):
        print(
            f" {mu.atom_id(at_pair[0]):12}"
            f" {mu.atom_id(at_pair[1]):12}"
            f" {np.sqrt(at_pair[2]):8.3f} A"
        )
        strcheck.summary['amide']['detected'].append({
            'at1': mu.atom_id(at_pair[0]),
            'at2': mu.atom_id(at_pair[1]),
            'dist': round(float(np.sqrt(at_pair[2])), 4)
        })
    return amide_check


def fix(strcheck, opts, fix_data=None):
    if not fix_data:
        return False
    if isinstance(opts, str):
        amide_fix = opts
    else:
        amide_fix = opts['fix']
    no_int_recheck = amide_fix is not None or strcheck.args['non_interactive']
    while fix_data:
        input_line = ParamInput(
            'Fix amide atoms',
            strcheck.args['non_interactive']
        )
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_auto()
        input_line.add_option_list(
            'resnum',
            sorted(mu.prep_rnums_list(fix_data['res_to_fix'])),
            case='sensitive',
            multiple=True
        )
        input_line.set_default('All')
        input_option, amide_fix = input_line.run(amide_fix)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], amide_fix

        if input_option == 'none':
            if strcheck.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'auto':
            to_fix = strcheck.strucm.amide_auto_fix(fix_data)
        else:
            to_fix = [
                res
                for res in fix_data['res_to_fix']
                if mu.residue_num(res) in amide_fix.split(',') or input_option == 'all'
            ]
        fix_num = 0
        done = set()
        for res in to_fix:
            if res not in done:
                try:
                    strcheck.strucm.invert_amide_atoms(res)
                except NotAValidResidueError as err:
                    return [err.message]
                done.add(res)  # To avoid double change
                fix_num += 1

        print(cts.MSGS['AMIDES_FIXED'].format(amide_fix, fix_num))
        strcheck.summary['amide']['fixed'] = amide_fix

        strcheck.strucm.modified = True
        fix_data = {}
        if not opts['no_recheck']:
            if not strcheck.args['quiet']:
                print(cts.MSGS['AMIDES_RECHECK'])
            fix_data = check(strcheck)
            amide_fix = ''
            if no_int_recheck:
                fix_data = {}
    return False
