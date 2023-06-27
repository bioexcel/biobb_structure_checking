''' Class to manage Chain internal data'''
import sys
from Bio.PDB.Chain import Chain
import biobb_structure_checking.modelling.utils as mu
from biobb_structure_checking.modelling.residue_set import ResidueSetList

class ChainsData():
    ''' Class to manage Chain(s) internal data'''
    def __init__(self, st):
        self.chain_ids = {}
        self.chain_details = {}
        self.has_chains_to_rename = False
        self.st = st

    def stats(self, prefix='', use_models=False) -> None:
        """ Print chains info """
        chids = {}
        num_chains = 0
        chains_str = []
        for mod in self.st:
            chids[mod.id] = []
            for ch_id in sorted(self.chain_ids[mod.id]):
                if not use_models:
                    mod_txt = ""
                    if mod.id > 1:
                        continue
                else:
                    mod_txt = f"/{mod.id}"
                if self.chain_ids[mod.id] == mu.UNKNOWN:
                    chids[mod.id].append(
                        f"{ch_id}{mod_txt}: "
                        f" Unknown (P:{self.chain_details[mod.id][ch_id][0]:.1%}"
                        f" DNA:{self.chain_details[mod.id][ch_id][1]:.1%}"
                        f" RNA:{self.chain_details[mod.id][ch_id][2]:.1%}"
                        f" UNK:{self.chain_details[mod.id][ch_id][3]:.1%}"
                    )
                else:
                    chids[mod.id].append(
                        f"{ch_id}{mod_txt}:"
                        f" {mu.CHAIN_TYPE_LABELS[self.chain_ids[mod.id][ch_id]]}"
                    )

            num_chains += len(self.chain_ids[mod.id])
            chains_str.append(', '.join(chids[mod.id]))

        return f"{prefix} Num. chains: {num_chains} ({'|'.join(chains_str)})"

    def set_chain_ids(self) -> None:
        """
        Identifies and sets the chain ids, guessing its nature
        (protein, dna, rna, ...)
        """
        self.chain_ids = {}
        self.chain_details = {}
        self.has_chains_to_rename = False

        for mod in self.st:
            self.chain_ids[mod.id] = {}
            self.chain_details[mod.id] = {}
            for chn in mod.get_chains():
                guess = mu.guess_chain_type(chn)
                self.chain_ids[mod.id][chn.id] = guess['type']
                self.chain_details[mod.id][chn.id] = guess['details']
                if chn.id == ' ':
                    self.has_chains_to_rename = self.has_chains_to_rename  or True

    def rename_empty_chain_label(self, new_label: str) -> str:
        '''Add chain label to empty ones'''
        if not self.has_chains_to_rename:
            return False
        for mod in self.st:
            if new_label == 'auto':
                new_label_char = 65
                while chr(new_label_char) in self.chain_ids[mod.id]\
                        and new_label_char < ord('z'):
                    new_label_char = +1
                new_label = chr(new_label_char)
            for chn in mod:
                if chn.id == ' ':
                    chn.id = new_label
        self.set_chain_ids()
        return new_label

    def _parse_renumber_str(self, renum_str):
        renum_to_do = []
        renum_str = renum_str.replace(' ', '')

        if ',' not in renum_str:
            tasks = [renum_str]
        else:
            tasks = renum_str.split(',')

        for tsk in tasks:
            if '=' not in tsk:
                print(f"ERROR: wrong syntax {tsk}, use origin=target")
                sys.exit()
            tsk_0, tsk_1 = tsk.split('=')
            chn_0, ini_0, fin_0 = _parse_task_str(tsk_0)
            chn_1, ini_1, fin_1 = _parse_task_str(tsk_1)
            if chn_0 != '*':
                if chn_1 == '*':
                    chn_1 = chn_0
                if chn_0 == chn_1: # chain shift
                    tmp_ch_id = _get_tmp_ch_id(self.st[0])
                    renum_to_do.append([
                        {'mod': 0, 'chn': chn_0, 'ini': ini_0, 'fin': fin_0},
                        {'mod': 0, 'chn': tmp_ch_id, 'ini': ini_1, 'fin': fin_1},
                    ])
                    renum_to_do.append([
                        {'mod': 0, 'chn': tmp_ch_id, 'ini': ini_1, 'fin': fin_1},
                        {'mod': 0, 'chn': chn_1, 'ini': ini_1, 'fin': fin_1},
                    ])
                else:
                    renum_to_do.append([
                        {'mod': 0, 'chn': chn_0, 'ini': ini_0, 'fin': fin_0},
                        {'mod': 0, 'chn': chn_1, 'ini': ini_1, 'fin': fin_1}
                    ])
            else:
                # replicate for all chains
                if chn_1 != '*':
                    print(
                        f"Error, use either wild card or explicit labels"
                        f" in both origin and updated ({chn_0}={chn_1})"
                    )
                    sys.exit()
                for mod in self.st:
                    if mod.id > 0:
                        continue  # Not implemented for multiple models
                    for ch_id in self.chain_ids[mod.id]:
                        renum_to_do.append([
                            {'mod': mod.id, 'chn': ch_id, 'ini': ini_0, 'fin': fin_0},
                            {'mod': mod.id, 'chn': ch_id, 'ini': ini_1, 'fin': fin_1}
                        ])
        return renum_to_do

    def renumber(self, renum_str, rem_inscodes=False, verbose=True):
        """ Renumber residues"""
        modified = False
        renum_to_do = []
        if renum_str.lower() == 'auto':
            tmp_ch_id = _get_tmp_ch_id(self.st[0])
            last_res_num = 0
            for mod in self.st:
                for chn in self.chain_ids:
                    n_term, c_term = mu.get_terms(mod, chn)
                    renum_to_do.append([
                        {'mod': mod.id, 'chn': chn, 'ini': n_term.id[1], 'fin': c_term.id[1]},
                        {'mod': mod.id, 'chn': tmp_ch_id, 'ini': last_res_num + 1, 'fin': 0}
                    ])
                    renum_to_do.append([
                        {
                            'mod': mod.id,
                            'chn':tmp_ch_id,
                            'ini': last_res_num + 1,
                            'fin': last_res_num + c_term.id[1] - n_term.id[1] + 1
                        },
                        {
                            'mod': mod.id,
                            'chn':chn,
                            'ini': last_res_num + 1,
                            'fin': 0
                        }
                    ])
                    last_res_num = last_res_num + c_term.id[1] - n_term.id[1] + 1
        else:
            renum_to_do = self._parse_renumber_str(renum_str)

        for mod in self.st:
            for tsk in renum_to_do:
                org, tgt = tsk
                if org['mod'] != mod.id:
                    continue
                if verbose and not mu.check_residue_id_order(mod[org['chn']]):
                    print(
                        f"WARNING: disordered residue ids found in {org['chn']}, "
                        f"may need explicit order combinations"
                    )

                n_term, c_term = mu.get_terms(mod, org['chn'])
                if not org['ini']:
                    org['ini'] = n_term.id[1]
                if not org['fin']:
                    org['fin'] = c_term.id[1]
                if not tgt['ini']:
                    tgt['ini'] = org['ini']
                tgt['fin'] = org['fin'] - org['ini'] + tgt['ini']
                if verbose:
                    print(
                        f"Renumbering {org['chn']}{org['ini']}-{org['chn']}{org['fin']} as "
                        f"{tgt['chn']}{tgt['ini']}-{tgt['chn']}{tgt['fin']}"
                    )
                # Not needed as checked per residue
                # if tgt['chn'] in mod:
                #     n_term1, c_term1 = mu.get_terms(mod, tgt['chn'])
                #     if n_term1.id[1] <= tgt['ini'] <= c_term1.id[1] or\
                #         n_term1.id[1] <= tgt['fin'] <= c_term1.id[1]:
                #         print(
                #             f"ERROR: new residue numbers ({tgt['ini']}-{tgt['fin']})"
                #             f" overlap with existing ones ({n_term1.id[1]}-{c_term1.id[1]}), "
                #             f" consider using --coords_only")
                #         sys.exit()
                if org['ini'] == tgt['ini'] and\
                        org['fin'] == tgt['fin'] and\
                        org['ini'] == n_term.id[1] and\
                        org['fin'] == c_term.id[1] and\
                        tgt['chn'] not in mod:
                    if verbose:
                        print(
                            f"Whole chains selected, just replacing chain"
                            f" labels from {org['chn']}"
                            f" to {tgt['chn']}"
                        )
                    mod[org['chn']].id = tgt['chn']
                    self.set_chain_ids()
                else:
                    to_move = []
                    for res in mod[org['chn']].get_residues():
                        if res.id[1] >= org['ini'] and res.id[1] <= org['fin']:
                            to_move.append(res)
                    if not to_move:
                        print("WARNING: No residues to move, exiting task")
                        continue

                    if tgt['chn'] not in mod:
                        if verbose:
                            print(f"Creating new chain {tgt['chn']}")
                        new_ch = Chain(tgt['chn'])
                        mod.add(new_ch)
                    else:
                        new_ch = mod[tgt['chn']]
                    inscodes_shift = 0
                    for res in sorted(to_move):
                        new_res = res.copy()
                        new_res_num = res.id[1]
                        new_inscode = res.id[2]
                        if rem_inscodes:
                            if mu.has_ins_code(res):
                                inscodes_shift += 1
                                new_inscode = ' '
                            new_res_num = res.id[1] + inscodes_shift
                        new_res.id = res.id[0], new_res_num - org['ini'] + tgt['ini'], new_inscode
                        new_res.parent = new_ch
                        col_res = _check_collision(new_res, new_ch)
                        if col_res:
                            print(
                                f"ERROR. New residue {mu.residue_id(new_res)}"
                                f" collides with existing {mu.residue_id(col_res)}"
                            )
                            sys.exit()
                        new_ch.add(new_res)
                        mod[org['chn']].detach_child(res.id)

                    if not mod[org['chn']]:
                        if verbose:
                            print(f"Removing empty chain {org['chn']}")
                        mod.detach_child(org['chn'])
                    self.set_chain_ids()
                    modified = True
            self._reorder_chains()
        return modified

    def rebuild(self, pairs_list=None):
        '''Rebuild chains from coordinates'''
        renumber_rules = []
        recover_chain_rules = []
        build_chains = ResidueSetList(pairs_list=pairs_list)
        print(f"{len(build_chains.sets)} chains/fragments found")
        for rset in sorted(build_chains.sets, key=lambda ss: ss.id):
            print(f" {rset}")
            new_chid = rset.id
            renumber_rules.append(
                f"{rset.inir.get_parent().id}:{rset.inir.id[1]}-{rset.finr.id[1]}={new_chid}:"
            )
            fin_chid = chr(ord('A') + int(new_chid) - 1)
            recover_chain_rules.append(f"{new_chid}:={fin_chid}:")
        self.renumber(','.join(renumber_rules))
        self.renumber(','.join(recover_chain_rules))
        self._reorder_chains()
        return build_chains.sets

    def _reorder_chains(self):
        #print("Ordering chains")
        for mod in self.st:
            for chn in sorted(mod):
                new_chn = mod[chn.id]
                mod.detach_child(chn.id)
                mod.add(new_chn)

    def get_chain_type(self, res):
        """ Return type of chain for residue"""
        mod = res.get_parent().get_parent()
        if mu.is_hetatm(res):
            return mu.UNKNOWN
        return self.chain_ids[mod.id][res.get_parent().id]

    def select(self, select_chains: str) -> None:
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids, | protein | dna | rna | na
        """
        if not self.chain_ids:
            self.set_chain_ids()

        for mod in self.st:
            if select_chains.lower() in ('protein', 'dna', 'rna', 'na'):
                if select_chains.lower() == 'na':
                    ch_ok = [mu.DNA, mu.RNA]
                else:
                    ch_ok = [mu.TYPE_LABEL[select_chains.lower()]]
            else:
                ch_ok = select_chains.split(',')
                for chn in ch_ok:
                    if chn not in self.chain_ids[mod.id]:
                        print('Warning: skipping unknown chain', chn)
            for chn in self.chain_ids[mod.id]:
                if chn not in ch_ok and self.chain_ids[mod.id][chn] not in ch_ok:
                    self.st[mod.id].detach_child(chn)
            if not self.st[mod.id]:
                print("ERROR: would remove all chains, exiting")
                sys.exit()

    def has_NA(self):
        """ Checks if any of the chains is NA"""
        has_na = False
        for mod in self.st:
            for ch_type in self.chain_ids[mod.id].values():
                has_na = (has_na or (ch_type > 1))
        return has_na


def _parse_task_str(ts_str):
    # Format [A:]ini[-fin]
    if ':' not in ts_str:
        ts_str = '*:' + ts_str
    chn, rnum = ts_str.split(':')
    if not rnum:
        return ts_str[:-1], 0, 0
    if '-' not in rnum:
        return chn, int(rnum), 0
    ini, fin = rnum.split('-')
    if not fin:
        fin = 0
    return chn, int(ini), int(fin)


def _get_tmp_ch_id(mod):
    tmp_id = ord('a')
    while chr(tmp_id) in mod:
        tmp_id += 1
    tmp_id = chr(tmp_id)
    if tmp_id not in mod:
        return tmp_id
    return ''


def _check_collision(new_res, new_ch):
    for res in new_ch.get_residues():
        if res.id[1] == new_res.id[1]:
            return res
    return None
