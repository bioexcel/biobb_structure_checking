""" Module supporting chains command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.model_utils as mu
from biobb_structure_checking.param_input import ParamInput

def _check(strcheck):
    print(cts.MSGS['CHAINS_DETECTED'].format(len(strcheck.strucm.chain_ids)))
    for ch_id in sorted(strcheck.strucm.chain_ids):
        if isinstance(strcheck.strucm.chain_ids[ch_id], list):
            print(
                cts.MSGS['UNKNOWN_CHAINS'].format(
                    ch_id, s=strcheck.strucm.chain_ids[ch_id]
                )
            )
        else:
            print(f" {ch_id}: {mu.CHAIN_TYPE_LABELS[strcheck.strucm.chain_ids[ch_id]]}")
    strcheck.summary['chains'] = {
        k:mu.CHAIN_TYPE_LABELS[v]
        for k, v in strcheck.strucm.chain_ids.items()
    }
    return len(strcheck.strucm.chain_ids) > 1

def _fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        select_chains = opts
        rename_chains = ''
    else:
        select_chains = opts['select']
        if 'rename' in opts:
            rename_chains = opts['rename']
        else:
            rename_chains = ''

    if strcheck.strucm.has_chains_to_rename:
        strcheck.summary['chains']['selected'] = {}
        input_line = ParamInput(
            'Add missing chain label',
            strcheck.args['non_interactive'],
            set_none='None'
        )
        input_line.add_option_none()
        input_line.add_option_list('auto_chain',['auto'])
        input_line.add_option_free_text('new_chain')
        input_line.set_default("auto")
        input_ok = False
        while not input_ok:
            input_option, rename_chains = input_line.run(rename_chains)
            if input_option == 'error':
                return cts.MSGS['UNKNOWN_SELECTION'], rename_chains
            if input_option in ('none', 'auto_chain') or\
                    (rename_chains not in strcheck.strucm.chain_ids):
                input_ok = True
            else:
                print(f"Chain {rename_chains} is already present")
                rename_chains = ''
        if input_option != 'none':
            strcheck.strucm.rename_empty_chain_label(rename_chains)
            _check(strcheck)

    strcheck.summary['chains']['selected'] = {}
    input_line = ParamInput('Select chain', strcheck.args['non_interactive'], set_none='All')
    input_line.add_option_all()
    input_line.add_option_list(
        'type', ['protein'], multiple=False
    )
    input_line.add_option_list(
        'type', ['na'], multiple=False
    )
    input_line.add_option_list(
        'type', ['dna'], multiple=False
    )
    input_line.add_option_list(
        'type', ['rna'], multiple=False
    )
    input_line.add_option_list(
        'chid', sorted(strcheck.strucm.chain_ids), multiple=True, case="sensitive"
    )
    input_line.set_default('All')

    input_option, select_chains = input_line.run(select_chains)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], select_chains

    if input_option == 'all':
        print(cts.MSGS['SELECT_ALL_CHAINS'])
        return False

    strcheck.strucm.select_chains(select_chains)
    print(cts.MSGS['SELECT_CHAINS'], select_chains)
    strcheck.summary['chains']['selected'] = strcheck.strucm.chain_ids

    return False
