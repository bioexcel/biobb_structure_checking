""" Module supporting models command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    """ Check for existing hydrogens atoms"""
    remh_list = mu.get_residues_with_H(strcheck.strucm.st)
    if not remh_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_RESIDUES_H_FOUND'])
        return {}

    print(cts.MSGS['RESIDUES_H_FOUND'].format(len(remh_list)))
    strcheck.summary['rem_hydrogen']['n_detected'] = len(remh_list)
    return {'remh_list': remh_list}


def fix(strcheck, opts, fix_data=None):
    """ Removes selected hydrogen atoms"""
    if isinstance(opts, str):
        remove_h = opts
    else:
        remove_h = opts['remove']

    input_line = ParamInput(
        'Remove hydrogen atoms',
        strcheck.args['non_interactive'],
        set_none='no'
    )
    input_line.add_option_yes_no()
    input_line.set_default('yes')

    input_option, remove_h = input_line.run(remove_h)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], remove_h

    if input_option == 'yes':
        rmh_num = 0
        for resh in fix_data['remh_list']:
            mu.remove_H_from_r(resh['r'])
            rmh_num += 1
        strcheck.strucm.revert_can_resnames(canonical=True)
        print(cts.MSGS['REMOVED_H'].format(rmh_num))
        strcheck.strucm.modified = True
        strcheck.summary['rem_hydrogen']['n_removed'] = rmh_num
    return False
