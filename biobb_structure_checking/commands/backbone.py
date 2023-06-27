""" Module supporting backbone command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.io.param_input import ParamInput
from biobb_structure_checking.structure_manager import NotEnoughAtomsError


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    fix_data = {}
    # Residues with missing backbone
    miss_bck_at_list = strcheck.strucm.get_missing_atoms('backbone')
    fix_data['miss_bck_at_list'] = miss_bck_at_list
    if miss_bck_at_list:
        strcheck.summary['backbone']['missing_atoms'] = {}
        print(cts.MSGS['BCK_MISSING_RESIDUES'].format(len(miss_bck_at_list)))
        for r_at in miss_bck_at_list:
            [res, at_list] = r_at
            print(f" {mu.residue_id(res):10} {','.join(at_list)}")
            strcheck.summary['backbone']['missing_atoms'][mu.residue_id(res)] = at_list
    else:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_BCK_MISSING'])

    # Not bound consecutive residues
    bck_check = strcheck.strucm.get_backbone_breaks()
    if bck_check['bck_breaks_list']:
        print(
            cts.MSGS['BACKBONE_BREAKS'].format(
                len(bck_check['bck_breaks_list'])
            )
        )
        strcheck.summary['backbone']['breaks'] = {'detected': []}
        for brk in bck_check['bck_breaks_list']:
            print(f" {mu.residue_id(brk[0]):10} - {mu.residue_id(brk[1]):10}")
            strcheck.summary['backbone']['breaks']['detected'].append([
                mu.residue_id(brk[0]), mu.residue_id(brk[1])
            ])
        fix_data['bck_breaks_list'] = bck_check['bck_breaks_list']
    else:
        fix_data['bck_breaks_list'] = []
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_BCK_BREAKS'])

    if bck_check['wrong_link_list']:
        print(cts.MSGS['UNEXPECTED_BCK_LINKS'])
        strcheck.summary['backbone']['wrong_links'] = []
        for brk in bck_check['wrong_link_list']:
            print(
                f" {mu.residue_id(brk[0]):10}"
                f" linked to {mu.residue_id(brk[1]):10},"
                f" expected {mu.residue_id(brk[2]):10}"
            )
            strcheck.summary['backbone']['wrong_links'].append([
                mu.residue_id(brk[0]),
                mu.residue_id(brk[1]),
                mu.residue_id(brk[2])
            ])
        fix_data['wrong_link_list'] = True
    else:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_BCK_LINKS'])

    if bck_check['not_link_seq_list']:
        print(cts.MSGS['CONSEC_RES_FAR'])
        strcheck.summary['backbone']['long_links'] = []
        for brk in bck_check['not_link_seq_list']:
            print(
                f" {mu.residue_id(brk[0]):10} -"
                f" {mu.residue_id(brk[1]):10}, "
                f" bond distance {brk[2]:8.3f} "
            )
            strcheck.summary['backbone']['long_links'].append([
                mu.residue_id(brk[0]),
                mu.residue_id(brk[1]),
                round(float(brk[2]), 5)
            ])
        fix_data['not_link_seq_list'] = True

    # TODO move this section to ligands
    if strcheck.strucm.st_data.modified_residue_list:
        print(cts.MSGS['MODIF_RESIDUES'])
        strcheck.summary['backbone']['mod_residues'] = []
        for brk in strcheck.strucm.st_data.modified_residue_list:
            print(f" {mu.residue_id(brk):10}")
            strcheck.summary['backbone']['mod_residues'].append(mu.residue_id(brk))
    # Provisional only missing atoms can be fixed
        fix_data['modified_residue_list'] = True
    return fix_data


def fix(strcheck, opts, fix_data=None):

    no_int_recheck = opts['fix_atoms'] is not None or\
        strcheck.args['non_interactive']

    fix_done = not fix_data['bck_breaks_list']

    fixed_main_res = []
    while not fix_done:
        if opts['extra_gap'] is None:
            opts['extra_gap'] = 0
        fixed_main = _backbone_fix_main_chain(
            strcheck,
            opts['fix_chain'],
            fix_data['bck_breaks_list'],
            strcheck.args['modeller_key'],
            opts['extra_gap']
        )
        if not fixed_main:
            fix_done = True
            continue
        fixed_main_res += fixed_main

        strcheck.summary['backbone']['main_chain_fix'] = [
            mu.residue_id(r)
            for r in fixed_main_res
        ]
        if fixed_main:
            strcheck.strucm.modified = True

        if no_int_recheck or not fixed_main or opts['no_recheck']:
            fix_done = True
            # Force silent re-checking to update modified residue pointers
            fix_data = check(strcheck)
        else:
            if not strcheck.args['quiet']:
                print(cts.MSGS['BACKBONE_RECHECK'])
            fix_data = check(strcheck)
            fix_done = not fix_data['bck_breaks_list']

    # Add CAPS
    fixed_caps = _backbone_add_caps(
        strcheck,
        opts['add_caps'],
        fix_data['bck_breaks_list']
    )

    strcheck.summary['backbone']['added_caps'] = [
        mu.residue_id(r)
        for r in fixed_caps
    ]

    if fixed_caps:
        print(f"Added caps: {', '.join(map(mu.residue_num, fixed_caps))}")
        strcheck.strucm.modified = True
        fix_data = check(strcheck)
    else:
        print('No caps added')

    # Add missing atoms
    fixed_res = []
    if fix_data['miss_bck_at_list']:
        fixed_res = _backbone_fix_missing(
            strcheck,
            opts['fix_atoms'],
            fix_data['miss_bck_at_list']
        )
    if not isinstance(fixed_res, list):
        fixed_res = []
    if not isinstance(fixed_caps, list):
        fixed_caps = []
    if not isinstance(fixed_main_res, list):
        fixed_main_res = []
    res_to_check = fixed_res + fixed_caps + fixed_main_res
    if res_to_check and not opts['no_check_clashes']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['CHECKING_CLASHES'])
        strcheck.summary['backbone']['clashes'] = strcheck.check_report_clashes(res_to_check)

    return False


def _backbone_fix_main_chain(strcheck, fix_main_bck, breaks_list, modeller_key, extra_gap):
    print("Main chain fixes")

    brk_rnums = [
        f"({mu.residue_num(resp[0])}-{mu.residue_num(resp[1])})".replace(' ', '')
        for resp in breaks_list
    ]
    input_line = ParamInput(
        'Fix Main Breaks',
        strcheck.args['non_interactive']
    )
    input_line.add_option_none()
    input_line.add_option_all()
    input_line.add_option_list(
        'brk', brk_rnums, case='sensitive', multiple=True
    )
    input_line.set_default('All')
    input_option, fix_main_bck = input_line.run(fix_main_bck)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], fix_main_bck

    strcheck.summary['backbone']['breaks']['option_selected'] = fix_main_bck

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return []

    # Checking for canonical sequence
    if not strcheck.strucm.sequence_data.has_canonical:
        read_ok = False
        while not read_ok:
            input_line = ParamInput(
                "Enter canonical sequence path (FASTA)",
                strcheck.args['non_interactive']
            )
            strcheck.args['fasta_seq_path'] = input_line.run(
                strcheck.args['fasta_seq_path']
            )
            if not strcheck.args['fasta_seq_path']:
                print(cts.MSGS['FASTA_MISSING'])

            read_ok = strcheck.strucm.sequence_data.load_sequence_from_fasta(
                strcheck.args['fasta_seq_path']
            )
            if not read_ok:
                strcheck.args['fasta_seq_path'] = None
            if strcheck.args['non_interactive'] and not read_ok:
                print(cts.MSGS['FASTA_MISSING'])
                return []

        strcheck.strucm.sequence_data.read_canonical_seqs(
            strcheck.strucm,
            False
        )
        strcheck.strucm.sequence_data.match_sequence_numbering()

    to_fix = [
        rpair
        for rpair in breaks_list
        if f"({mu.residue_num(rpair[0])}-{mu.residue_num(rpair[1])})".replace(' ', '')\
            in fix_main_bck.split(',') or input_option == 'all'
    ]
    return strcheck.strucm.fix_backbone_chain(to_fix, modeller_key, extra_gap)


def _backbone_add_caps(strcheck, add_caps, bck_breaks_list):
    print("Capping terminal ends")
    term_res = strcheck.strucm.get_term_res()
    term_rnums = [mu.residue_num(p[1]) for p in term_res]
    brk_res = set()

    for res_0, res_1 in bck_breaks_list:
        brk_res.add(res_0)
        brk_res.add(res_1)
    true_term = []
    for res in term_res:
        if res[1] not in brk_res:
            true_term.append(res)

    print(
        "True terminal residues: "
        f"{','.join([mu.residue_num(r[1]) for r in true_term])}"
    )
    if bck_breaks_list:
        print(
            "Terminal residues from backbone breaks: ",
            ','.join([
                f"{mu.residue_num(r0)}-{mu.residue_num(r1)}"
                for r0, r1 in bck_breaks_list
            ])
        )

    input_line = ParamInput('Cap residues', strcheck.args['non_interactive'])
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list('terms', ['Terms'])
    if bck_breaks_list:
        input_line.add_option_list('brks', ['Breaks'])
    input_line.add_option_list(
        'resnum', term_rnums, case='sensitive', multiple=True
    )
    input_line.set_default('none')
    input_option, add_caps = input_line.run(add_caps)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], add_caps

    strcheck.summary['backbone']['add_caps'] = add_caps

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return []

    if input_option == 'terms':
        to_fix = true_term
    elif input_option == 'brks':
        to_fix = [
            pair
            for pair in term_res
            if pair[1] in brk_res
        ]
    else:
        to_fix = [
            pair
            for pair in term_res
            if mu.residue_num(pair[1]) in add_caps.split(',') or\
            input_option == 'all'
        ]

    return strcheck.strucm.add_main_chain_caps(to_fix)


def _backbone_fix_missing(strcheck, fix_back, fix_at_list):
    print("Fixing missing backbone atoms")
    fixbck_rnums = [mu.residue_num(r_at[0]) for r_at in fix_at_list]
    input_line = ParamInput(
        'Fix bck missing',
        strcheck.args['non_interactive']
    )
    input_line.add_option_all()
    input_line.add_option_none()
    input_line.add_option_list(
        'resnum', fixbck_rnums, case='sensitive', multiple=True
    )
    input_line.set_default('none')
    input_option, fix_back = input_line.run(fix_back)

    if input_option == 'error':
        return cts.MSGS['UNKNOWN_SELECTION'], fix_back

    strcheck.summary['backbone']['missing_atoms']['selected'] = fix_back

    if input_option == 'none':
        if strcheck.args['verbose']:
            print(cts.MSGS['DO_NOTHING'])
        return []

    to_fix = [
        r_at
        for r_at in fix_at_list
        if mu.residue_num(r_at[0]) in fix_back.split(',') or input_option == 'all'
    ]

    if not strcheck.args['quiet']:
        print(cts.MSGS['ADDING_BCK_ATOMS'])

    fix_num = 0
    strcheck.summary['backbone']['missing_atoms']['fixed'] = []
    fixed_res = []
    for r_at in to_fix:
        try:
            if strcheck.strucm.fix_backbone_O_atoms(r_at):
                fix_num += 1
        except NotEnoughAtomsError as err:
            print(err.message)

        strcheck.summary['backbone']['missing_atoms']['fixed'].append(
            mu.residue_id(r_at[0])
        )
        fixed_res.append(r_at[0])

    print(cts.MSGS['BCK_ATOMS_FIXED'].format(fix_num))

    return fixed_res
