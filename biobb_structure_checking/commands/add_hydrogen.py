""" Module supporting add_hydrogen command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput

def check(strcheck):

    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    ion_res_list = strcheck.strucm.get_ion_res_list()

    if not ion_res_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_SELECT_ADDH'])
        return {'ion_res_list': []}

    fix_data = {
        'ion_res_list': ion_res_list,
    }
    print(cts.MSGS['SELECT_ADDH_RESIDUES'].format(len(ion_res_list)))

    for res_type in strcheck.strucm.data_library.residue_codes['protein']:
        residue_list = [
            mu.residue_num(r_at[0])
            for r_at in ion_res_list
            if r_at[0].get_resname() == res_type
        ]
        if residue_list and strcheck.args['verbose']:
            print(f" {res_type} {','.join(residue_list)}")

    strcheck.summary['add_hydrogen']['ionic_detected'] = [
        mu.residue_id(r_at[0]) for r_at in ion_res_list
    ]
    # print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
    return fix_data


def fix(strcheck, opts, fix_data=None):
    if not fix_data:
        return False

    if not strcheck.strucm.st_data.fixed_side and not opts['no_fix_side']:
        print("WARNING: fixing side chains, override with --no_fix_side")
        strcheck.fixside(['--fix', 'all'])

    # Fixing previously marked N and C terms
    strcheck.strucm.revert_terms()

    if isinstance(opts, str):
        add_h_mode = opts
    else:
        add_h_mode = opts['add_mode']

    input_line = ParamInput('Mode', strcheck.args['non_interactive'])
    input_line.add_option_none()
    sel_options = ['auto']
    if fix_data['ion_res_list']:
        sel_options += ['pH', 'list', 'int', 'int_his']
    input_line.add_option_list('selection', sel_options)
    input_line.set_default('auto')
    input_option, add_h_mode = input_line.run(add_h_mode)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], add_h_mode

    if input_option == "none":
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    std_ion = strcheck.strucm.data_library.std_ion

    add_h_mode = add_h_mode.lower()

    ion_to_fix = {
        r_at[0]: std_ion[r_at[0].get_resname()]['std']
        for r_at in fix_data['ion_res_list']
    }

    if add_h_mode == 'auto':
        if not strcheck.args['quiet']:
            print('Selection: auto')
        strcheck.summary['add_hydrogen']['selection'] = 'auto'
    else:
        if add_h_mode == 'ph':
            ph_value = opts['pH']
            input_line = ParamInput(
                "pH Value",
                strcheck.args['non_interactive'],
                set_none=7.0
            )
            input_line.add_option_numeric(
                "pH",
                [],
                opt_type="float",
                min_val=0.,
                max_val=14.
            )
            input_line.set_default(7.0)
            input_option, ph_value = input_line.run(ph_value)

            if not strcheck.args['quiet']:
                print('Selection: pH', ph_value)

            strcheck.summary['add_hydrogen']['selection'] = f"pH {ph_value}"

            for r_at in fix_data['ion_res_list']:
                res = r_at[0]
                rcode = res.get_resname()
                if ph_value <= std_ion[rcode]['pK']:
                    ion_to_fix[res] = std_ion[rcode]['lowpH']
                else:
                    ion_to_fix[res] = std_ion[rcode]['highpH']
        else:
            if add_h_mode == 'list':
                ions_list = opts['list']
                if not strcheck.args['quiet']:
                    print('Selection: list')

                ions_list = ParamInput(
                    "Enter Forms list as [*:]his22hip",
                    strcheck.args['non_interactive']
                ).run(ions_list)

                strcheck.summary['add_hydrogen']['selection'] = ions_list

                # Changes in tautomeric forms made as mutationts eg.HisXXHip
                mutations = strcheck.strucm.prepare_mutations(ions_list)
                for mut in mutations.mutation_list:
                    for mut_res in mut.mutations:
                        ion_to_fix[mut_res['resobj']] = mut_res['new_id']
            else:
                if add_h_mode == 'int':
                    if not strcheck.args['quiet']:
                        print('Selection: interactive')
                    res_list = fix_data['ion_res_list']
                elif add_h_mode == 'int_his':
                    if not strcheck.args['quiet']:
                        print('Selection: int_his')
                    res_list = [
                        r_at for r_at in fix_data['ion_res_list']
                        if r_at[0].get_resname() == 'HIS'
                    ]
                    strcheck.summary['add_hydrogen']['selection'] = []
                for r_at in res_list:
                    rcode = r_at[0].get_resname()
                    input_line = ParamInput(
                        "Select residue form for " + mu.residue_id(r_at[0]),
                        strcheck.args['non_interactive']
                    )
                    input_line.add_option_list('list', r_at[1].keys())
                    input_line.default = std_ion[rcode]['std']

                    form = None
                    input_option, form = input_line.run(form)
                    ion_to_fix[r_at[0]] = form.upper()
                    strcheck.summary['add_hydrogen']['selection'].append(
                        f"{rcode} {form.upper()}"
                    )

    strcheck.strucm.add_hydrogens(
        ion_to_fix,
        add_charges=opts['add_charges'].upper()
    )
    strcheck.strucm.modified = True
    return False
