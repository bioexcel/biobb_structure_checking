""" Module supporting fixside command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    miss_at_list = strcheck.strucm.get_missing_atoms('side')
    extra_at_list = strcheck.strucm.check_extra_atoms()

    if not miss_at_list and not extra_at_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_MISSING_SIDE_ATOMS'])
        return {}

    fix_data = {
        'miss_at_list': miss_at_list,
        'extra_at_list': extra_at_list
    }
    if miss_at_list:
        strcheck.summary['fixside']['detected_missing'] = {}
        print(cts.MSGS['MISSING_SIDE_ATOMS'].format(len(miss_at_list)))
        for r_at in miss_at_list:
            res, at_list = r_at
            print(f" {mu.residue_id(res):10} {','.join(at_list)}")
            strcheck.summary['fixside']['detected_missing'][mu.residue_id(res)] = at_list
    if extra_at_list:
        strcheck.summary['fixside']['detected_unknown'] = {}
        print(cts.MSGS['UNKNOWN_SIDE_ATOMS'].format(len(extra_at_list)))
        for r_at in extra_at_list:
            res, at_list = r_at
            print(f" {mu.residue_id(res):10} {','.join(at_list)}")
            strcheck.summary['fixside']['detected_unknown'][mu.residue_id(res)] = at_list

    return fix_data


def fix(strcheck, opts, fix_data=None):

    if not fix_data:
        return False
    if isinstance(opts, str):
        fix_side = opts
    else:
        fix_side = opts['fix']

    fixside_rnums = [
        mu.residue_num(r_at[0])
        for r_at in fix_data['miss_at_list']
    ]

    input_line = ParamInput('fixside', strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'resnum', fixside_rnums, case='sensitive', multiple=True
    )
    input_line.set_default('All')
    input_option_fix, fix_side = input_line.run(fix_side)

    if input_option_fix == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], fix_side

    strcheck.summary['fixside']['fix'] = fix_side

    if input_option_fix == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    if input_option_fix == 'all':
        to_add = fix_data['miss_at_list']
        to_remove = fix_data['extra_at_list']
    else:
        to_add = [
            r_at
            for r_at in fix_data['miss_at_list']
            if mu.residue_num(r_at[0]) in fix_side.split(',')
        ]

        to_remove = [
            r_at
            for r_at in fix_data['extra_at_list']
            if mu.residue_num(r_at[0]) in fix_side.split(',')
        ]

    if not strcheck.args['quiet']:
        print(cts.MSGS['FIXING_SIDE_CHAINS'])

    fix_num = 0
    rem_num = 0
    strcheck.summary['fixside']['fixed'] = []
    strcheck.summary['fixside']['removed'] = []
    fixed_res = []
    if not opts['no_rem_extra']:
        for r_at in to_remove:
            print(mu.residue_id(r_at[0]))
            for at_id in r_at[1]:
                print("  Removing", at_id)
                mu.remove_atom_from_res(r_at[0], at_id)
                rem_num += 1
            strcheck.summary['fixside']['removed'].append(
                mu.residue_id(r_at[0])
            )
            fixed_res.append(r_at[0])
    if not opts['rebuild']:
        for r_at in to_add:
            mu.remove_H_from_r(r_at[0], verbose=True)
            strcheck.strucm.fix_side_chain(r_at)
            fix_num += 1
            strcheck.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))
            fixed_res.append(r_at[0])
    else:
        strcheck.strucm.rebuild_side_chains(to_add)
        fixed_res = [r_at[0] for r_at in to_add]

    print(cts.MSGS['SIDE_CHAIN_FIXED'].format(fix_num))
    strcheck.strucm.st_data.fixed_side = True
    strcheck.strucm.modified = True
    # Checking new clashes
    if not opts['no_check_clashes']:
        print(cts.MSGS['CHECKING_CLASHES'])
        strcheck.summary['fixside_clashes'] = strcheck.check_report_clashes(fixed_res)
    return False
