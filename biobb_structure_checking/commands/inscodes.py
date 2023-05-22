""" Inscodes command """
import biobb_structure_checking.constants as cts
import biobb_structure_checking.modelling.utils as mu


def check(strcheck):
    ins_codes_list = strcheck.strucm.get_ins_codes()
    if not ins_codes_list:
        if not strcheck.args['quiet']:
            print(cts.MSGS['NO_INSCODES_FOUND'])
        return {}

    print(cts.MSGS['INSCODES_FOUND'].format(len(ins_codes_list)))
    strcheck.summary['inscodes'] = []
    for res in ins_codes_list:
        print(mu.residue_id(res))
        strcheck.summary['inscodes'].append(mu.residue_id(res))
    return {'ins_codes_list': ins_codes_list}


def fix(strcheck, opts, fix_data=None):
    if opts['renumber']:
        min_res = {}
        for res in fix_data['ins_codes_list']:
            ch_id = res.get_parent().id
            if ch_id not in min_res:
                min_res[ch_id] = res.id[1]
        renum_str = []
        for chn, res_num in min_res.items():
            renum_str.append(
                f"{chn}:{res_num}=z:{res_num},z:{res_num}={chn}:{res_num}"
            )
        strcheck.chains({
            "renumber": ','.join(renum_str),
            "rem_inscodes": True,
            "select": "all",
            "verbose": False
        })
    return False
