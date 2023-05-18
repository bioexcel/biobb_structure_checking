""" Module supporting models command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    wat_list = [
        res
        for res in mu.get_ligands(strcheck.strucm.st, incl_water=True)
        if mu.is_wat(res)
    ]

    if not wat_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_WATERS'])
        return {}

    print(cts.MSGS['WATERS_FOUND'].format(len(wat_list)))
    strcheck.summary['water']['n_detected'] = len(wat_list)

    return {'wat_list': wat_list}


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        remove_wat = opts
    else:
        remove_wat = opts['remove']

    input_line = ParamInput(
        'Remove',
        strcheck.args['non_interactive'],
        set_none='no'
    )
    input_line.add_option_yes_no()
    input_line.set_default('yes')
    input_option, remove_wat = input_line.run(remove_wat)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], remove_wat

    if input_option == 'yes':
        rmw_num = 0
        for res in fix_data['wat_list']:
            strcheck.strucm.remove_residue(res, False)
            rmw_num += 1
        strcheck.strucm.update_internals()
        print(cts.MSGS['WATER_REMOVED'].format(rmw_num))
        strcheck.summary['water']['n_removed'] = rmw_num
    return False
