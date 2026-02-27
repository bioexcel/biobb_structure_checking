""" Module supporting ligands command"""
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.pdbio.param_input import ParamInput


def check(strcheck):

    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    lig_list = mu.get_ligands(strcheck.strucm.st, incl_water=False)

    if not lig_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_LIGANDS_FOUND'])
        return {}

    print(cts.MSGS['LIGANDS_DETECTED'].format(len(lig_list)))

    lig_names = {}
    for res in lig_list:
        if res.get_resname() not in lig_names:
            lig_names[res.get_resname()] = mu.fetch_residue_name_by_id(res.get_resname())
    
    fix_data = {
        'lig_list': lig_list,
        'ligand_rids': set(),
        'ligand_rnums': []
    }

    strcheck.summary['ligands'] = {'detected': []}

    for res in sorted(lig_list, key=lambda x: x.index):
        if strcheck.strucm.models_data.has_models():
            if res.get_parent().get_parent().id > 0:
                continue
            print(f" {mu.residue_id(res, False)}/* {lig_names[res.get_resname()]} ")
        else:
            print(f" {mu.residue_id(res, False)} {lig_names[res.get_resname()]}" )

        strcheck.summary['ligands']['detected'].append(mu.residue_id(res))
        fix_data['ligand_rids'].add(res.get_resname())
        fix_data['ligand_rnums'].append(mu.residue_num(res))

    lig_contacts = {lig: set() for lig in lig_list}

    for _, contacts in strcheck.strucm.check_r_list_clashes(
        lig_list,
        contact_types= ['ligand'],
        use_wat=True
    ).items():
        if contacts:         
            for atom1, atom2, dist in contacts.values(): 
                if atom1.serial_number > atom2.serial_number:
                    atom2, atom1 = atom1, atom2
                res1 = atom1.get_parent()
                res2 = atom2.get_parent()
                if mu.is_hetatm(res1) and res1 in lig_contacts:
                    lig_contacts[res1].add(atom2)
                if mu.is_hetatm(res2) and res2 in lig_contacts:
                    lig_contacts[res2].add(atom1)

    strcheck.summary['ligands']['contacts'] = {}
    for lig in sorted(lig_contacts, key=lambda x:(x.get_parent().id, mu.residue_num(x))):
        contacts = sorted(
            lig_contacts[lig],
            key=lambda x: (
                x.get_parent().get_parent().id, 
                mu.residue_num(x.get_parent())
                )
        )
        strcheck.summary['ligands']['contacts'][mu.residue_id(lig)] = [mu.atom_id(c) for c in contacts]
        print(f"Contacts for {mu.residue_id(lig)}: {', '.join([mu.atom_id(c) for c in contacts])}")
    return fix_data


def fix(strcheck, opts, fix_data=None):
    if isinstance(opts, str):
        remove_ligands = opts
    else:
        remove_ligands = opts['remove']
    input_line = ParamInput('Remove', strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'byrids', sorted(fix_data['ligand_rids']), multiple=True
    )
    input_line.add_option_list(
        'byresnum', fix_data['ligand_rnums'], case='sensitive', multiple=True
    )
    input_line.set_default('All')
    input_option, remove_ligands = input_line.run(remove_ligands)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], remove_ligands

    strcheck.summary['ligands']['removed'] = {
        'opt': remove_ligands,
        'lst': []
    }

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return False

    if input_option == 'all':
        to_remove = fix_data['lig_list']
    elif input_option == 'byrids':
        to_remove = [
            res
            for res in fix_data['lig_list']
            if res.get_resname() in remove_ligands.split(',')
        ]
    elif input_option == 'byresnum':
        to_remove = [
            res
            for res in fix_data['lig_list']
            if mu.residue_num(res) in remove_ligands.split(',')
        ]
    rl_num = 0
    for res in to_remove:
        strcheck.summary['ligands']['removed']['lst'].append(mu.residue_id(res))
        strcheck.strucm.remove_residue(res, False)
        rl_num += 1
    strcheck.strucm.update_internals()
    print(cts.MSGS['LIGANDS_REMOVED'].format(remove_ligands, rl_num))
    strcheck.summary['ligands']['n_removed'] = rl_num
    return False
