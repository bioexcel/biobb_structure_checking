''' Class to manage Chain internal data'''
import sys
import biobb_structure_checking.modelling.utils as mu
from Bio.PDB.Chain import Chain

class ChainsData():
    ''' Class to manage Chain(s) internal data'''
    def __init__(self, st):
        self.chain_ids = {}
        self.chain_details = {}
        self.has_chains_to_rename = False
        self.st = st

    def stats(self, prefix='') -> None:
        """ Print chains info """
        chids = []
        for ch_id in sorted(self.chain_ids):
            if self.chain_ids == mu.UNKNOWN:
                chids.append(
                    f"{ch_id}: "
                    f" Unknown (P:{self.chain_details[ch_id][0]:.1%}"
                    f" DNA:{self.chain_details[ch_id][1]:.1%}"
                    f" RNA:{self.chain_details[ch_id][2]:.1%}"
                    f" UNK:{self.chain_details[ch_id][3]:.1%}"
                )
            else:
                chids.append(f"{ch_id}: {mu.CHAIN_TYPE_LABELS[self.chain_ids[ch_id]]}")

        return f"{prefix} Num. chains: {len(self.chain_ids)} ({', '.join(chids)})"

    def set_chain_ids(self, biounit=False) -> None:
        """
        Identifies and sets the chain ids, guessing its nature (protein, dna, rna, ...)
        """
        self.chain_ids = {}
        self.has_chains_to_rename = False
        for chn in self.st.get_chains():
            if not biounit and chn.get_parent().id > 0:
                continue
            guess = mu.guess_chain_type(chn)
            self.chain_ids[chn.id] = guess['type']
            self.chain_details[chn.id] = guess['details']
            if chn.id == ' ':
                self.has_chains_to_rename = self.has_chains_to_rename  or True

    def rename_empty_chain_label(self, new_label: str) -> str:
        '''Add chain label to empty ones'''
        if not self.has_chains_to_rename:
            return False
        if new_label == 'auto':
            new_label_char = 65
            while chr(new_label_char) in self.chain_ids and new_label_char < ord('z'):
                new_label_char =+1
            new_label = chr(new_label_char)
        for mod in self.st:
            for chn in mod:
                if chn.id == ' ':
                    chn.id = new_label
        self.set_chain_ids()
        return new_label

    def _parse_task_str(self, ts_str):
        #Format [A:]ini[-fin]
        if ':' not in ts_str:
            ts_str = '*:' + ts_str
        chn, rnum = ts_str.split(':')
        if not rnum:
            return ts_str[:-1], 0, 0
        if not '-' in rnum:
            return chn, int(rnum), 0
        ini, fin = rnum.split('-')
        if not fin:
            fin = 0
        return chn, int(ini), int(fin)

    def _parse_renumber_str(self, renum_str):
        renum_to_do = []
        renum_str = renum_str.replace(' ','')

        if ',' not in renum_str:
            tasks = [renum_str]
        else:
            tasks = renum_str.split(',')

        for tsk in tasks:
            if '=' not in tsk:
                print(f"ERROR: wrong syntax {tsk}, use origin=target")
                sys.exit()
            tsk_0, tsk_1 = tsk.split('=')
            chn_0, ini_0, fin_0 = self._parse_task_str(tsk_0)
            chn_1, ini_1, fin_1 = self._parse_task_str(tsk_1)
            if chn_0 != '*':
                if chn_1 == '*':
                    chn_1 = chn_0
                renum_to_do.append([
                    {'chn':chn_0, 'ini': ini_0, 'fin':fin_0},
                    {'chn':chn_1, 'ini': ini_1, 'fin':fin_1}
                ])
            else:
                # replicate for all chains
                if chn_1 != '*':
                    print(
                        f"Error, use either wild card or explicit labels"
                        f" in both origin and updated ({chn_0}={chn_1})"
                    )
                    sys.exit()
                for ch_id in self.chain_ids:
                    renum_to_do.append([
                        {'chn':ch_id, 'ini':ini_0, 'fin':fin_0},
                        {'chn':ch_id, 'ini':ini_1, 'fin':fin_1}
                    ])
        return renum_to_do

    def renumber(self, renum_str, allow_merge=False):
        """ Renumber residues"""
        renum_to_do = []
        if renum_str.lower() == 'auto':
            tmp_ch_id = 'z'
            last_res_num = 0
            for chn in self.chain_ids:
                n_term, c_term = mu.get_terms(self.st[0], chn)
                renum_to_do.append([
                    {'chn': chn, 'ini': n_term.id[1], 'fin': c_term.id[1]},
                    {'chn':tmp_ch_id, 'ini': last_res_num + 1, 'fin': 0}
                ]
                )
                renum_to_do.append([
                    {'chn':tmp_ch_id, 'ini': last_res_num + 1, 'fin': last_res_num + c_term.id[1] - n_term.id[1] + 1},
                    {'chn':chn, 'ini': last_res_num + 1, 'fin': 0}
                ])
                last_res_num = last_res_num + c_term.id[1] - n_term.id[1] + 1
        else:
           renum_to_do = self._parse_renumber_str(renum_str)

        for mod in self.st:
            for tsk in renum_to_do:
                org, tgt = tsk
                if not mu.check_residue_id_order(mod[org['chn']]):
                    print(f"WARNING: disordered residue ids found in {org['chn']}, use explicit order combinations")
                if not allow_merge and org['chn'] != tgt['chn'] and tgt['chn'] in self.chain_ids:
                    print(f"ERROR: chain {tgt['chn']} already exists and --allow_merge not set")
                    sys.exit()
                n_term, c_term = mu.get_terms(mod, org['chn'])
                if not org['ini']:
                    org['ini'] = n_term.id[1]
                if not org['fin']:
                    org['fin'] = c_term.id[1]
                if not tgt['ini']:
                    tgt['ini'] = org['ini']
                tgt['fin'] = org['fin'] - org['ini'] + tgt['ini']
                print(
                    f"Renumbering {org['chn']}{org['ini']}-{org['chn']}{org['fin']} as "
                    f"{tgt['chn']}{tgt['ini']}-{tgt['chn']}{tgt['fin']}"
                )
                if tgt['chn'] in mod:
                    n_term1, c_term1 = mu.get_terms(mod, tgt['chn'])
                    if n_term1.id[1] <= tgt['ini'] <= c_term1.id[1] or\
                        n_term1.id[1] <= tgt['fin'] <= c_term1.id[1]:
                        print("WARNING: new residue numbers overlap with existing ones")
                        sys.exit()
                if org['ini'] == tgt['ini'] and\
                    org['fin'] == tgt['fin'] and\
                    org['ini'] == n_term.id[1] and\
                    org['fin'] == c_term.id[1] and\
                    tgt['chn'] not in mod:
                    print(
                        f"Whole chains selected, just replacing chain labels from {org['chn']}"
                        f" to {tgt['chn']}"
                    )
                    mod[org['chn']].id = tgt['chn']
                    self.set_chain_ids()
                else:
                    if tgt['chn'] not in mod:
                        print(f"Creating new chain {tgt['chn']}")
                        new_ch = Chain(tgt['chn'])
                        mod.add(new_ch)
                    else:
                        new_ch = mod[tgt['chn']]
                    to_move = []
                    for res in mod[org['chn']].get_residues():
                        if res.id[1] >= org['ini'] and res.id[1] <= org['fin']:
                            to_move.append(res)
                    for res in sorted(to_move):
                        new_res = res.copy()
                        new_res.id = res.id[0], res.id[1] - org['ini'] + tgt['ini'], res.id[2]
                        new_res.parent = new_ch
                        new_ch.add(new_res)
                        mod[org['chn']].detach_child(res.id)

                    if len(mod[org['chn']]) == 0:
                        mod.detach_child(org['chn'])
                    self.set_chain_ids()
        return 1

    def get_chain_type(self, res):
        """ Return type of chain for residue"""
        if mu.is_hetatm(res):
            return mu.UNKNOWN
        return self.chain_ids[res.get_parent().id]

    def select(self, select_chains: str) -> None:
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids, | protein | dna | rna | na
        """
        if not self.chain_ids:
            self.set_chain_ids()

        if select_chains.lower() in ('protein', 'dna', 'rna', 'na'):
            if select_chains.lower() == 'na':
                ch_ok = [mu.DNA, mu.RNA]
            else:
                ch_ok = [mu.TYPE_LABEL[select_chains.lower()]]
        else:
            ch_ok = select_chains.split(',')
            for chn in ch_ok:
                if chn not in self.chain_ids:
                    print('Warning: skipping unknown chain', chn)
        for mod in self.st:
            for chn in self.chain_ids:
                if chn not in ch_ok and self.chain_ids[chn] not in ch_ok:
                    self.st[mod.id].detach_child(chn)
            if not self.st[mod.id]:
                print("ERROR: would remove all chains, exiting")
                sys.exit()

    def has_NA(self):
        """ Checks if any of the chains is NA"""
        has_na = False
        for ch_type in self.chain_ids.values():
            has_na = (has_na or (ch_type > 1))
        return has_na