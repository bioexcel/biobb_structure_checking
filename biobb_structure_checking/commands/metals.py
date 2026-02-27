""" Module supporting metals command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.pdbio.param_input import ParamInput


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    met_list = strcheck.strucm.get_metal_atoms()

    if not met_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_METALS_FOUND'])
        return {}

    print(cts.MSGS['METALS_FOUND'].format(len(met_list)))
    strcheck.summary['metals'] = {'detected':[]}
    fix_data = {
        'met_list': met_list,
        'met_rids': [],
        'at_groups': {}
    }
    met_names = {}
    met_res_list = []
    for atom in met_list:
        res = atom.get_parent()
        met_res_list.append(res)
        if res.get_resname() not in met_names:
            met_names[res.get_resname()] = mu.fetch_residue_name_by_id(res.get_resname())
    
    for atm in sorted(met_list, key=lambda x: x.serial_number):
        print(f" {mu.atom_id(atm):12} ({met_names[atm.get_parent().get_resname()]})")
        res = atm.get_parent()
        fix_data['met_rids'].append(mu.residue_num(res))
        if atm.id not in fix_data['at_groups']:
            fix_data['at_groups'][atm.id] = []
        fix_data['at_groups'][atm.id].append(atm)
        strcheck.summary['metals']['detected'].append(mu.residue_id(res))
    
    met_contacts = {met: set() for met in met_res_list}

    for contact_type, contacts in strcheck.strucm.check_r_list_clashes(
        met_res_list,
        contact_types= ['ligand'],
        use_wat=True
    ).items():
        if contacts:
            for atom1, atom2, dist in contacts.values():
                if atom1.serial_number > atom2.serial_number:
                    atom2, atom1 = atom1, atom2
                res1 = atom1.get_parent()
                res2 = atom2.get_parent()
                if mu.is_hetatm(res1) and res1 in met_contacts:
                    met_contacts[res1].add(atom2)
                if mu.is_hetatm(res2) and res2 in met_contacts:
                    met_contacts[res2].add(atom1)
    
    strcheck.summary['metals']['contacts'] = {}
    for lig in sorted(met_contacts, key=lambda x:(x.get_parent().id, mu.residue_num(x))):
        contacts = sorted(
            met_contacts[lig],
            key=lambda x: (
                x.get_parent().get_parent().id, 
                mu.residue_num(x.get_parent())
            )
        )
        strcheck.summary['metals']['contacts'][mu.residue_id(lig)] = [mu.atom_id(c) for c in contacts]
        print(f"Contacts for {mu.residue_id(lig)}: {', '.join([mu.atom_id(c) for c in contacts])}")

    return fix_data


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        remove_metals = opts
    else:
        remove_metals = opts['remove']

    input_line = ParamInput("Remove", strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'atids',
        sorted(fix_data['at_groups']),
        case='sensitive',
        multiple=True
    )
    input_line.add_option_list(
        'resids',
        fix_data['met_rids'],
        case='sensitive',
        multiple=True
    )
    input_line.set_default('All')
    input_option, remove_metals = input_line.run(remove_metals)

    if input_option == "error":
        return cts.MSGS['UNKNOWN_SELECTION'], remove_metals

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    if input_option == 'all':
        to_remove = fix_data['met_list']
    elif input_option == 'resids':
        rid_list = remove_metals.split(',')
        to_remove = [
            atm for atm in fix_data['met_list']
            if mu.residue_num(atm.get_parent()) in rid_list
        ]
    elif input_option == 'atids':
        to_remove = []
        for atid in remove_metals.split(','):
            to_remove.extend(fix_data['at_groups'][atid])

    strcheck.summary['metals']['removed'] = []

    rmm_num = 0
    for atm in to_remove:
        strcheck.summary['metals']['removed'].append(mu.residue_id(atm.get_parent()))
        strcheck.strucm.remove_residue(atm.get_parent(), False)
        rmm_num += 1
    strcheck.strucm.update_internals()
    print(cts.MSGS['METALS_REMOVED'].format(remove_metals, rmm_num))
    strcheck.summary['metals']['n_removed'] = rmm_num
    return False
