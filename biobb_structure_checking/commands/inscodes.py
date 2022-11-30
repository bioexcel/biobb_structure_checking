""" Inscodes command """
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu

def _check(strcheck):
    ins_codes_list = strcheck.strucm.get_ins_codes()
    if not ins_codes_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_INSCODES_FOUND'])
        return {}

    print(cts.MSGS['INSCODES_FOUND'].format(len(ins_codes_list)))
    strcheck.summary['inscodes'] = []
    for res in ins_codes_list:
        print(mu.residue_id(res))
        strcheck.summary['inscodes'].append(mu.residue_id(res))
    return {'ins_codes_list': ins_codes_list}

def _fix(strcheck, opts, fix_data=None):
    # TODO implement method _inscodes_fix
    if opts['renum']:
        print("--renum option not implemented (yet)")
    return False
