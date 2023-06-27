""" Module supporting chiral command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    chiral_check = strcheck.strucm.check_chiral_sides()

    if 'list' not in chiral_check:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_CHIRALS'])
        return {}
    strcheck.summary['chiral']['n_chirals'] = len(chiral_check['list'])

    if not chiral_check['res_to_fix']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_WRONG_CHIRAL_SIDE'])
        return {}

    print(cts.MSGS['WRONG_CHIRAL_SIDE'].format(
        len(chiral_check['res_to_fix'])
    ))
    strcheck.summary['chiral']['detected'] = []
    for res in chiral_check['res_to_fix']:
        print(f" {mu.residue_id(res):10}")
        strcheck.summary['chiral']['detected'].append(mu.residue_id(res))

    return chiral_check

def fix(strcheck, opts, fix_data=None):
    if not fix_data:
        return False
    if isinstance(opts, str):
        chiral_fix = opts
    else:
        chiral_fix = opts['fix']
    input_line = ParamInput(
        'Fix chiralities',
        strcheck.args['non_interactive']
    )
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'resnum',
        mu.prep_rnums_list(fix_data['res_to_fix']),
        case='sensitive',
        multiple=True
    )
    input_line.set_default('All')
    input_option, chiral_fix = input_line.run(chiral_fix)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], chiral_fix

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    if input_option == 'all':
        to_fix = fix_data['res_to_fix']
    else:
        to_fix = [
            res
            for res in fix_data['res_to_fix']
            if mu.residue_num(res) in chiral_fix.split(',')
        ]
    fix_num = 0
    for res in to_fix:
        strcheck.strucm.fix_chiral_chains(res)
        fix_num += 1
    print(cts.MSGS['CHIRAL_SIDE_FIXED'].format(chiral_fix, fix_num))

    if not opts['no_check_clashes']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['CHECKING_CLASHES'])

        strcheck.summary['chiral_clashes'] = strcheck.check_report_clashes(to_fix)

    strcheck.strucm.modified = True
    return False
