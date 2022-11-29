''' Class to manage Chain internal data'''
import sys
import biobb_structure_checking.model_utils as mu
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
    
    def _parse_renumber_str(self, renum_str):
        renum_to_do = []
        renum_str = renum_str.replace(' ','')
        if ',' not in renum_str:
            tasks = [renum_str]
        else:
            tasks = renum_str.split(',')
        for tsk in tasks:
            tsk_strs = tsk.split('=')
            chn = {}
            ini = {}
            fin = {}
            for i in range(2):
                ts_str = tsk_strs[i]
                if ts_str.endswith(':'):
                    chn[i] = ts_str[:-1]
                    ini[i] = 0
                    fin[i] = 0
                else:
                    if ':' in ts_str:
                        chn[i], seq = ts_str.split(':')
                    else:
                        chn[i], seq = '*', ts_str
                    if '-' in seq:
                        ini[i], fin[i] = seq.split('-')
                        ini[i] = int(ini[i])
                        fin[i] = int(fin[i])
                    else:
                        ini[i] = int(seq)
                        fin[i] = 0
            if chn[0] == '*':
                if chn[1] != '*':
                    print(f"Error, use either wild card or explicit labels in both origin and updated ({tsk['chn']})")
                    sys.exit()
                for ch_id in self.chain_ids:
                    renum_to_do.append(
                        {'chn':[ch_id, ch_id], 'ini':ini.copy(), 'fin':fin.copy()}
                    )
            else:
                if chn[1] == '*':
                    chn[1] = chn[0]
                renum_to_do.append({'chn':chn.copy(), 'ini':ini.copy(), 'fin':fin.copy()})
        return renum_to_do

    def _get_terms(self, mod, ch_id):
        n_term = 0
        c_term = 0
        res_num = 0
        for res in mod[ch_id].get_residues():
            res_num = res.id[1]
            if not n_term:
                n_term = res_num
        c_term = res_num
        return n_term, c_term

    def renumber(self, renum_str, allow_merge=False):
        """ Renumber residues"""
        if renum_str.lower() == 'auto':
            print ("Auto renumbering not implemented")            
            return 0
                        
        for mod in self.st:
            renum_to_do = self._parse_renumber_str(renum_str)
            for tsk in renum_to_do:
                if not allow_merge and tsk['chn'][0] != tsk['chn'][1] and tsk['chn'][1] in self.chain_ids:
                    print(f"ERROR: chain {tsk['chn'][1]} already exists and --allow_merge not set")
                    sys.exit()
                n_term, c_term = self._get_terms(mod, tsk['chn'][0])
                if not tsk['ini'][0]:
                    tsk['ini'][0] = n_term
                if not tsk['fin'][0]:
                    tsk['fin'][0] = c_term
                if not tsk['ini'][1]:
                    tsk['ini'][1] = tsk['ini'][0]
                if not tsk['fin'][1]:
                    tsk['fin'][1] = tsk['fin'][0] - tsk['ini'][0] + tsk['ini'][1]
                print(f"Renumbering {tsk['chn'][0]}{tsk['ini'][0]}-{tsk['chn'][0]}{tsk['fin'][0]} to "
                    f"{tsk['chn'][1]}{tsk['ini'][1]}-{tsk['chn'][1]}{tsk['fin'][1]}"
                    )
                if tsk['chn'][1] in mod:
                    n_term1, c_term1 = self._get_terms(mod, tsk['chn'][1])
                    if n_term1 <= tsk['ini'][1] <= c_term1 or n_term1 <= tsk['fin'][1] <= c_term1:
                        print(f"WARNING: new residue numbers overlap with existing ones")
                        sys.error()
                if tsk['ini'][0] == tsk['ini'][1] and\
                    tsk['fin'][0] == tsk['fin'][1] and\
                    tsk['ini'][0] == n_term and\
                    tsk['fin'][0] == c_term:
                    print(f"Whole chains selected, just replacing chain labels from {tsk['chn'][0]} to {tsk['chn'][1]}")
                    mod[tsk['chn'][0]].id = tsk['chn'][1]
                else:
                    if tsk['chn'][1] not in mod:
                        print(f"Creating new chain {tsk['chn'][1]}")
                        new_ch = Chain(tsk['chn'][1])
                        mod.add(new_ch)
                    else:
                        new_ch = mod[tsk['chn'][1]]
                    to_del = []
                    for res in mod[tsk['chn'][0]].get_residues():
                        if res.id[1] >= tsk['ini'][0] and res.id[1] <= tsk['fin'][0]:
                            new_res = res.copy()
                            new_res.id = res.id[0], res.id[1] - tsk['ini'][0] + tsk['ini'][1], res.id[2]
                            to_del.append(res.id)
                            new_res.parent = new_ch
                            new_ch.add(new_res)
                    for res_id in to_del:
                        mod[tsk['chn'][0]].detach_child(res_id)
                    if not len(mod[tsk['chn'][0]]):
                        mod.detach_child(tsk['chn'][0])
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
