""" Module supporting cistransbck command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.model_utils as mu

def _check(strcheck):
    (cis_backbone_list, lowtrans_backbone_list) = strcheck.strucm.check_cis_backbone()
    if cis_backbone_list:
        strcheck.summary['cistransbck']['cis'] = []
        print(cts.MSGS['CIS_BONDS'].format(len(cis_backbone_list)))
        for lnk in cis_backbone_list:
            res1, res2, dih = lnk
            print(
                '{:10} {:10} Dihedral: {:8.3f}'.format(
                    mu.residue_id(res1),
                    mu.residue_id(res2),
                    dih
                )
            )
            strcheck.summary['cistransbck']['cis'].append([
                mu.residue_id(res1),
                mu.residue_id(res2),
                round(float(dih), 3)
            ])
    else:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_CIS_BONDS'])

    if lowtrans_backbone_list:
        strcheck.summary['cistransbck']['unusual_trans'] = []
        print(cts.MSGS['LOWTRANS_BONDS'].format(len(lowtrans_backbone_list)))
        for lnk in lowtrans_backbone_list:
            res1, res2, dih = lnk
            print(
                '{:10} {:10} Dihedral: {:8.3f}'.format(
                    mu.residue_id(res1), mu.residue_id(res2), dih
                )
            )
            strcheck.summary['cistransbck']['unusual_trans'].append([
                mu.residue_id(res1), mu.residue_id(res2), round(float(dih), 3)
            ])
    else:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_LOWTRANS_BONDS'])

    return {}

def _fix(strcheck, opts, fix_data):
    pass
