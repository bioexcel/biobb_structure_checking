""" Module supporting chiral_bck command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu


def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None
    check_lst = strcheck.strucm.get_chiral_bck_list()
    if 'list' not in check_lst:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_BCK_CHIRALS'])
        strcheck.summary['chiral_bck']['n_chirals'] = 0
        return {}

    strcheck.summary['chiral_bck']['n_chirals'] = len(check_lst['list'])

    if not check_lst['res_to_fix']:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_CHIRAL_BCK_RESIDUES'])
        return {}

    print(cts.MSGS['CHIRAL_BCK_RESIDUES'].format(len(check_lst['res_to_fix'])))
    strcheck.summary['chiral_bck']['detected'] = []
    for res in check_lst['res_to_fix']:
        print(f" {mu.residue_id(res):10}")
        strcheck.summary['chiral_bck']['detected'].append(mu.residue_id(res))

    return {}


def fix(strcheck, opts, fix_data):
    pass

#    def _chiral_bck_fix(strcheck, chiral_fix, fix_data=None, check_clashes=True):
#        return False
# TODO chiral_bck_fix
#        input_line = ParamInput('Fix CA chiralities', strcheck.args['non_interactive'])
#        input_line.add_option_all()
#        input_line.add_option_none()
#        input_line.add_option_list('resnum', strcheck.chiral_bck_rnums,
#            case='sensitive', multiple=True)
#        [input_option, chiral_fix] = input_line.run(chiral_fix)
#        if input_option == 'error':
#            print('Warning: unknown option {}'.format(amide_fix))
#            strcheck.summary['chiral']['error'] = 'Unknown option'
#            return
#
#        if input_option == 'none':
#            if strcheck.args['verbose']:
#                print("Nothing to do")
#        else:
#            if input_option == 'all':
#                to_fix = strcheck.chiral_bck_res_to_fix
#            else:
#                to_fix = []
#                for res in strcheck.chiral_bck_res_to_fix:
#                    if mu.residue_num(res) in chiral_fix.split(','):
#                        to_fix.append(res)
#            n = 0
#            for res in to_fix:
#                mu.strucm.invert_chiral_CA(res)
#                n += 1
#            print('Quiral residues fixed {} ({})'.format(chiral_fix, n))
#            if check_clashes:
#                if not strcheck.args['quiet']:
#                    print("Checking new steric clashes")
#
#            strcheck.strucm.modified = True
