""" Module supporting chains command"""

import logging
import biobb_structure_checking.constants as cts
import biobb_structure_checking.model_utils as mu
from biobb_structure_checking.param_input import ParamInput

def _check(strcheck):
    print(cts.MSGS['CHAINS_DETECTED'].format(len(strcheck.strucm.chains_data.chain_ids)))
    for ch_id in sorted(strcheck.strucm.chains_data.chain_ids):
        if isinstance(strcheck.strucm.chains_data.chain_ids[ch_id], list):
            print(
                cts.MSGS['UNKNOWN_CHAINS'].format(
                    ch_id, s=strcheck.strucm.chains_data.chain_ids[ch_id]
                )
            )
        else:
            print(f" {ch_id}: {mu.CHAIN_TYPE_LABELS[strcheck.strucm.chains_data.chain_ids[ch_id]]}")
    strcheck.summary['chains'] = {
        k:mu.CHAIN_TYPE_LABELS[v]
        for k, v in strcheck.strucm.chains_data.chain_ids.items()
    }
    strcheck.summary['chains']['unlabelled'] = strcheck.strucm.chains_data.has_chains_to_rename
    return len(strcheck.strucm.chains_data.chain_ids) > 1 or strcheck.strucm.chains_data.has_chains_to_rename

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
        if 'renumber' in opts:
            renumber_chains = opts['renumber']

    if strcheck.strucm.chains_data.has_chains_to_rename:
        input_line = ParamInput(
            'Add missing chain label',
            strcheck.args['non_interactive'],
            set_none='None'
        )
        input_line.add_option_none()
        input_line.add_option_list('auto_chain', ['auto'])
        input_line.add_option_free_text('new_chain')
        input_line.set_default("auto")
        input_ok = False
        while not input_ok:
            input_option, rename_chains = input_line.run(rename_chains)
            if input_option == 'error':
                return cts.MSGS['UNKNOWN_SELECTION'], rename_chains
            if input_option in ('none', 'auto_chain') or\
                    (rename_chains not in strcheck.strucm.chains_data.chain_ids):
                input_ok = True
            else:
                logging.warning(f"Chain {rename_chains} is already present")
                rename_chains = ''
        if input_option != 'none':
            new_label = strcheck.strucm.rename_empty_chain_label(rename_chains)
            _check(strcheck)
            strcheck.summary['chains']['renamed'] = new_label

    if renumber_chains:
        if strcheck.strucm.chains_data.has_chains_to_rename:
            print("WARNING: unlabelled chains detected")
        result = strcheck.strucm.renumber_chain_residues(
            renumber_chains,
            opts['allow_merge']
        )
        if result:
            strcheck.summary['chains']['renumbered'] = renumber_chains

    strcheck.summary['chains']['selected'] = {}
    input_line = ParamInput('Select chain', strcheck.args['non_interactive'], set_none='All')
    input_line.add_option_all()
    for typ in 'protein', 'na', 'dna', 'rna':
        input_line.add_option_list('type', [typ], multiple=False)
    input_line.add_option_list(
        'chid', sorted(strcheck.strucm.chains_data.chain_ids), multiple=True, case="sensitive"
    )
    input_line.set_default('All')

    input_option, select_chains = input_line.run(select_chains)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], select_chains

    if input_option == 'all':
        logging.warning(cts.MSGS['SELECT_ALL_CHAINS'])
        return False

    strcheck.strucm.select_chains(select_chains)
    logging.log(15, f"{cts.MSGS['SELECT_CHAINS']} {select_chains}")
    strcheck.summary['chains']['selected'] = strcheck.strucm.chains_data.chain_ids

    return False
