''' Class to manage Chain internal data'''
import sys
import biobb_structure_checking.model_utils as mu

class ChainsData():
    ''' Class to manage Chain(s) internal data'''
    def __init__(self, st):
        self.chain_ids = {}
        self.chain_details = {}
        self.has_chains_to_rename = False
        self.st = st
        self.modified = False

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
        self.modified = True
        return new_label

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
