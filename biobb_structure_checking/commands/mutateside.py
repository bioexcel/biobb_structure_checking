""" Module supporting fixside command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput


def check(strcheck):
    # TODO Check _mutateside_check function?
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    return True


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        mut_list = opts
    else:
        mut_list = opts['mut']

    if opts['na_seq']:
        mut_list = strcheck.strucm.prepare_mutations_from_na_seq(opts['na_seq'])

    input_line = ParamInput(
        'Mutation list',
        strcheck.args['non_interactive'],
        set_none=''
    )

    mut_list = input_line.run(mut_list)

    if mut_list == '':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    mutations = strcheck.strucm.prepare_mutations(mut_list)

    print(cts.MSGS['MUTATIONS_TO_DO'])
    for mut in mutations.mutation_list:
        print(mut)
    if opts['rebuild']:
        if strcheck.strucm.chains_data.has_NA():
            print(cts.MSGS['WARN_NOBUILD_NA'])
        mutated_res = strcheck.strucm.rebuild_mutations(mutations)
    else:
        mutated_res = strcheck.strucm.apply_mutations(mutations)
    if not opts['no_check_clashes']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['CHECKING_CLASHES'])

        strcheck.summary['mutateside_clashes'] = strcheck.check_report_clashes(mutated_res)

    strcheck.strucm.modified = True
    return False
