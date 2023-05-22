""" Module supporting models command"""
# from unittest.util import strclass
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    print(cts.MSGS['MODELS_FOUND'].format(strcheck.strucm.models_data.nmodels))
    strcheck.summary['models'] = {
        'nmodels': strcheck.strucm.models_data.nmodels
    }
    if strcheck.strucm.models_data.nmodels == 1:
        if not strcheck.args['quiet']:
            print(cts.MSGS['SINGLE_MODEL'])
        return {}

    strcheck.summary['models']['type'] = strcheck.strucm.models_data.models_type
    supimp = ''
    if not strcheck.strucm.models_data.has_superimp_models():
        supimp = 'do not'
    print(
        cts.MSGS['MODELS_GUESS'].format(
            supimp,
            strcheck.strucm.models_data.models_type['rmsd'],
            mu.MODEL_TYPE_LABELS[
                strcheck.strucm.models_data.models_type['type']
            ]
        )
    )
    return True


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        select_model = opts
    else:
        select_model = opts['select']

    input_line = ParamInput('Select Model Num', strcheck.args['non_interactive'], set_none='All')
    input_line.add_option_all()
    input_line.add_option_numeric(
        'modelno',
        [],
        opt_type='int',
        min_val=1,
        max_val=strcheck.strucm.models_data.nmodels,
        multiple=True
    )
    input_option, select_model = input_line.run(select_model)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], select_model

    print(cts.MSGS['SELECT_MODEL'], select_model)

    if input_option != 'all':
        strcheck.strucm.select_model(select_model)

    if opts['superimpose']:
        strcheck.strucm.superimpose_models()
        print(cts.MSGS['SUPIMP_MODELS'].format(
            strcheck.strucm.models_data.models_type["rmsd"]
        ))
        strcheck.summary['models']['superimp_rmsd'] = strcheck.strucm.models_data.models_type["rmsd"]

    if opts['build_complex']:
        added_chains = strcheck.strucm.build_complex()
        print(cts.MSGS['ADDED_CHAINS_TO_COMPLEX'].format(added_chains))
        strcheck.summary['models']['added_chains'] = added_chains

    if opts['save_split']: # tag as modified to force save
        strcheck.strucm.modified = True

    strcheck.summary['models']['selected'] = select_model
    return False
