""" Module supporting chains command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    single_chain = len(strcheck.strucm.chains_data.chain_ids[0]) == 1
    for mod_id in strcheck.strucm.chains_data.chain_ids:
        if not strcheck.strucm.models_data.has_models():
            print(cts.MSGS['CHAINS_DETECTED'].format(
                len(strcheck.strucm.chains_data.chain_ids[0])
            ))
        else:
            print(
                f"Model {mod_id}: ",
                cts.MSGS['CHAINS_DETECTED'].format(
                    len(strcheck.strucm.chains_data.chain_ids[mod_id])
                ))

            single_chain = single_chain and\
                (len(strcheck.strucm.chains_data.chain_ids[mod_id]) == 1)

        for ch_id in sorted(strcheck.strucm.chains_data.chain_ids[mod_id]):
            modtxt = ''
            if strcheck.strucm.models_data.has_models():
                modtxt = f"/{mod_id}"
            if isinstance(strcheck.strucm.chains_data.chain_ids[mod_id][ch_id], list):
                print(
                    cts.MSGS['UNKNOWN_CHAINS'].format(
                        ch_id,
                        s=strcheck.strucm.chains_data.chain_ids[mod_id][ch_id]
                    )
                )
            else:
                print(f" {ch_id}{modtxt}: {mu.CHAIN_TYPE_LABELS[strcheck.strucm.chains_data.chain_ids[mod_id][ch_id]]}")

        strcheck.summary['chains'] = {
            k: f"{mu.CHAIN_TYPE_LABELS[v]}/{mod_id}"
            for k, v in strcheck.strucm.chains_data.chain_ids[mod_id].items()
        }
    strcheck.summary['chains']['unlabelled'] = strcheck.strucm.chains_data.has_chains_to_rename

    return not single_chain or\
        strcheck.strucm.chains_data.has_chains_to_rename or\
        '--rebuild' in strcheck.args['options'] or\
        '--renumber' in strcheck.args['options']


def fix(strcheck, opts, fix_data=None):
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
        if 'rebuild' in opts:
            rebuild_chains = opts['rebuild']
        else:
            rebuild_chains = False

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
                    (rename_chains not in strcheck.strucm.chains_data.chain_ids[0]):  #TODO heterogenous models
                input_ok = True
            else:
                print(f"Chain {rename_chains} is already present")
                rename_chains = ''
        if input_option != 'none':
            new_label = strcheck.strucm.rename_empty_chain_label(rename_chains)
            check(strcheck)
            strcheck.summary['chains']['renamed'] = new_label

    if rebuild_chains:
        if strcheck.strucm.models_data.has_models():
            print("WARNING: Rebuild chains not (yet) implemented for Models, skipping")
        else:
            print("Rebuilding chains from backbone connectivity")
            result = strcheck.strucm.rebuild_chains(verbose='verbose' in opts)
            if result:
                strcheck.summary['chains']['rebuild'] = result

    if renumber_chains:
        if strcheck.strucm.models_data.has_models():
            print("WARNING: Renumber chains not (yet) implemented for Models, skipping")
        else:
            if strcheck.strucm.chains_data.has_chains_to_rename:
                print("WARNING: unlabelled chains detected")
            if 'verbose' not in opts and opts['rem_inscodes']:
                opts['verbose'] = False
            else:
                opts['verbose'] = True
            result = strcheck.strucm.renumber_chain_residues(
                renumber_chains,
                rem_inscodes=opts['rem_inscodes'],
                verbose=opts['verbose']
            )
            if result:
                strcheck.summary['chains']['renumbered'] = renumber_chains

    strcheck.summary['chains']['selected'] = {}
    input_line = ParamInput(
        'Select chain',
        strcheck.args['non_interactive'],
        set_none='All'
    )
    input_line.add_option_all()
    for typ in 'protein', 'na', 'dna', 'rna':
        input_line.add_option_list('type', [typ], multiple=False)
    if not strcheck.strucm.models_data.has_models():
        input_line.add_option_list(
            'chid',
            sorted(strcheck.strucm.chains_data.chain_ids[0]),
            multiple=True,
            case="sensitive"
        )
    else:
        input_line.add_option_list(
            'chid',
            [
                sorted(strcheck.strucm.chains_data.chain_ids[mod.id])
                for mod in strcheck.strucm.st
            ],
            multiple=True,
            case="sensitive"
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
    strcheck.summary['chains']['selected'] = strcheck.strucm.chains_data.chain_ids

    return False
