""" Module supporting getss command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):

    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    SS_bonds = strcheck.strucm.get_SS_bonds()
    if not SS_bonds:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_SS'])
        return {}
    print(cts.MSGS['POSSIBLE_SS'].format(len(SS_bonds)))
    strcheck.summary['getss'] = {'found':[]}
    for ssb in SS_bonds:
        print(
            f" {mu.atom_id(ssb[0]):12}"
            f" {mu.atom_id(ssb[1]):12}"
            f" {ssb[2]:8.3f}"
        )
        strcheck.summary['getss']['found'].append({
            'at1': mu.atom_id(ssb[0]),
            'at2': mu.atom_id(ssb[1]),
            'dist': round(float(ssb[2]), 4)
        })
    return SS_bonds


def fix(strcheck, opts, fix_data=None):
    if not fix_data:
        return False
    if isinstance(opts, str):
        getss_mark = opts
    else:
        getss_mark = opts['mark']

    pairs_list = [
        f"{mu.residue_num(a[0].get_parent())}-{mu.residue_num(a[1].get_parent())}"
        for a in fix_data
    ]
    input_line = ParamInput('Mark SS', strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'bypair', pairs_list, multiple=True
    )
    input_line.set_default('All')
    input_option, getss_mark = input_line.run(getss_mark)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], getss_mark

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    cys_to_mark = []

    for pair in fix_data:
        pair_id = f"{mu.residue_num(pair[0].get_parent())}-{mu.residue_num(pair[1].get_parent())}"
        if input_option == 'all' or pair_id in getss_mark.split(','):
            cys_to_mark.append(pair[0].get_parent())
            cys_to_mark.append(pair[1].get_parent())
    strcheck.summary['getss']['marked'] = [
        mu.residue_id(a)
        for a in cys_to_mark
    ]
    strcheck.strucm.mark_ssbonds(cys_to_mark)
    strcheck.strucm.update_internals()
    return False
