""" Module supporting altloc command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    alt_loc_res = strcheck.strucm.get_altloc_residues()
    if not alt_loc_res:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_ALTLOC_FOUND'])
        return {}

    print(cts.MSGS['ALTLOC_FOUND'].format(len(alt_loc_res)))

    strcheck.summary['altloc'] = {}

    fix_data = {
        'alt_loc_res': alt_loc_res,
        'altlocs': {}
    }

    for res in sorted(alt_loc_res, key=lambda x: x.index):
        rid = mu.residue_id(res)
        print(rid)
        strcheck.summary['altloc'][rid] = {}
        fix_data['altlocs'][res] = sorted(alt_loc_res[res][0].child_dict)
        for atm in alt_loc_res[res]:
            strcheck.summary['altloc'][rid][atm.id] = []
            alt_str = f"  {atm.id:4}"
            for alt in sorted(atm.child_dict):
                alt_str += f" {alt} ({atm.child_dict[alt].occupancy:4.2f})"
                strcheck.summary['altloc'][rid][atm.id].append({
                    'loc_label': alt,
                    'occupancy': atm.child_dict[alt].occupancy
                })
            print(alt_str)

    return fix_data


def fix(strcheck, opts, fix_data=None):

    if isinstance(opts, str):
        select_altloc = opts
    else:
        select_altloc = opts['select']

    # Prepare the longest possible list of alternatives
    altlocs = []
    max_al_len = 0
    for res in fix_data['altlocs']:
        if len(fix_data['altlocs'][res]) > max_al_len:
            altlocs = fix_data['altlocs'][res]
            max_al_len = len(fix_data['altlocs'][res])

    input_line = ParamInput(
        'Select alternative',
        strcheck.args['non_interactive'],
        set_none='All'
    )
    input_line.add_option_all()
    input_line.add_option_list('occup', ['occupancy'])
    input_line.add_option_list('altids', altlocs, case='upper')
    input_line.add_option_list(
        'resnum',
        mu.prep_rnums_list(fix_data['altlocs']),
        opt_type='pair_list',
        list2=altlocs,
        case='sensitive',
        multiple=True
    )
    input_line.set_default('occupancy')

    input_option, select_altloc = input_line.run(select_altloc)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], select_altloc

    if input_option != 'all':
        print(f"Selecting location {select_altloc}")
        if input_option in ('occup', 'altids'):
            select_altloc = select_altloc.upper()
            to_fix = {
                res: {
                    'ats': value,
                    'select' : select_altloc
                } for res, value in fix_data['alt_loc_res'].items()
            }

        elif input_option == 'resnum':
            to_fix = {}
            selected_rnums = {}
            for rsel in select_altloc.split(','):
                rnum, alt = rsel.split(':')
                selected_rnums[rnum] = alt
            to_fix = {
                res: {
                    'ats': value,
                    'select': selected_rnums[mu.residue_num(res)]
                }
                for res, value in fix_data['alt_loc_res'].items()
                if mu.residue_num(res) in selected_rnums
            }
        for res in to_fix:
            strcheck.strucm.select_altloc_residue(res, to_fix[res])

    strcheck.summary['altloc']['selected'] = select_altloc

    return False
