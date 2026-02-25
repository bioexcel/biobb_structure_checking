""" Module supporting models command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.pdbio.param_input import ParamInput


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
    water_contacts = {wat: set() for wat in wat_list}

    for contact_type, contacts in strcheck.strucm.check_r_list_clashes(
        wat_list,
        contact_types= mu.HB_CONTACT_TYPES,
        use_wat=True
    ).items():
        if contacts:
            for atom1, atom2, dist in contacts.values():
                if atom1.serial_number > atom2.serial_number:
                    atom2, atom1 = atom1, atom2
                res1 = atom1.get_parent()
                res2 = atom2.get_parent()
                if mu.is_wat(res1) and mu.is_wat(res2):
                    continue
                if mu.is_wat(res1):
                    water_contacts[res1].add(res2)
                if mu.is_wat(res2):
                    water_contacts[res2].add(res1)

    strcheck.summary['water']['contacts'] = {}
    for wat, contacts in water_contacts.items():
        strcheck.summary['water']['contacts'][mu.residue_id(wat)] = sorted(mu.residue_id(res) for res in contacts)

    wats_per_n_contacts = [set() for s in range(0, 5)]

    for wat, contacts in sorted(water_contacts.items()):
        n_contacts = len(contacts)
        if n_contacts < 5:
            wats_per_n_contacts[n_contacts].add(wat)
    for n_contacts, wat_set in enumerate(wats_per_n_contacts):
        if len(wat_set) > 0:
            print(f"Water molecules in contact with {n_contacts} residues: {', '.join(mu.residue_id(wat) for wat in sorted(wat_set))}")
    return {'wat_list': wat_list, 'water_contacts': water_contacts}


def fix(strcheck, opts, fix_data=None):

    if opts['remove'] is not None:
        if isinstance(opts, str):
            remove_wat = opts
        else:
            remove_wat = opts['remove']
    else:
        remove_wat = opts['keep'] if opts['keep'] is not None else None

    input_line = ParamInput(
        'Remove',
        strcheck.args['non_interactive'],
        set_none='no'
    )
    input_line.add_option_yes_no()
    input_line.add_option_numeric(
        'keep',
        [0, 1, 2, 3, 4],
        'int',
        0,
        4,
        multiple=False,
        label_text="Keep contacts with at least N residues"
    )

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

    elif input_option == 'keep':
        strcheck.summary['water']['removed'] = []
        n_removed = 0
        for res in fix_data['water_contacts']:
            if len(fix_data['water_contacts'][res]) < int(remove_wat):
                strcheck.summary['water']['removed'].append(mu.residue_id(res))
                strcheck.strucm.remove_residue(res, False)
                n_removed += 1
        print(cts.MSGS['WATER_KEEP'].format(n_removed, remove_wat))
        strcheck.strucm.update_internals()

    return False
