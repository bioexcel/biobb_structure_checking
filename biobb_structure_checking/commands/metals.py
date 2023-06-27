""" Module supporting metals command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    met_list = strcheck.strucm.get_metal_atoms()

    if not met_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_METALS_FOUND'])
        return {}

    print(cts.MSGS['METALS_FOUND'].format(len(met_list)))
    strcheck.summary['metals'] = {'detected':[]}
    fix_data = {
        'met_list': met_list,
        'met_rids': [],
        'at_groups': {}
    }
    for atm in sorted(met_list, key=lambda x: x.serial_number):
        print(f" {mu.atom_id(atm):12}")
        res = atm.get_parent()
        fix_data['met_rids'].append(mu.residue_num(res))
        if atm.id not in fix_data['at_groups']:
            fix_data['at_groups'][atm.id] = []
        fix_data['at_groups'][atm.id].append(atm)
        strcheck.summary['metals']['detected'].append(mu.residue_id(res))

    return fix_data


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        remove_metals = opts
    else:
        remove_metals = opts['remove']

    input_line = ParamInput("Remove", strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'atids',
        sorted(fix_data['at_groups']),
        case='sensitive',
        multiple=True
    )
    input_line.add_option_list(
        'resids',
        fix_data['met_rids'],
        case='sensitive',
        multiple=True
    )
    input_line.set_default('All')
    input_option, remove_metals = input_line.run(remove_metals)

    if input_option == "error":
        return cts.MSGS['UNKNOWN_SELECTION'], remove_metals

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    if input_option == 'all':
        to_remove = fix_data['met_list']
    elif input_option == 'resids':
        rid_list = remove_metals.split(',')
        to_remove = [
            atm for atm in fix_data['met_list']
            if mu.residue_num(atm.get_parent()) in rid_list
        ]
    elif input_option == 'atids':
        to_remove = []
        for atid in remove_metals.split(','):
            to_remove.extend(fix_data['at_groups'][atid])

    strcheck.summary['metals']['removed'] = []

    rmm_num = 0
    for atm in to_remove:
        strcheck.summary['metals']['removed'].append(mu.residue_id(atm.get_parent()))
        strcheck.strucm.remove_residue(atm.get_parent(), False)
        rmm_num += 1
    strcheck.strucm.update_internals()
    print(cts.MSGS['METALS_REMOVED'].format(remove_metals, rmm_num))
    strcheck.summary['metals']['n_removed'] = rmm_num
    return False
