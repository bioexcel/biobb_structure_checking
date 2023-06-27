''' Class to manage models internal data'''
from Bio.PDB.Superimposer import Superimposer
import biobb_structure_checking.modelling.utils as mu


class ModelsData():
    '''Class to manage models internal data'''
    def __init__(self, st):
        # Checking models type according to RMS among models
        self.st = st
        self.nmodels = len(st)
        self.models_type = mu.guess_models_type(st) if self.nmodels > 1 else 0

    def stats(self, prefix='') -> None:
        """ Print stats """
        if self.nmodels > 1:
            return f"{prefix} Num. models: {self.nmodels} (type: {mu.MODEL_TYPE_LABELS[self.models_type['type']]}, {self.models_type['rmsd']:8.3f} A)"
        return f"{prefix} Num. models: {self.nmodels}"

    def select(self, keep_model: str) -> None:
        """ Selects model(s) and delete the others from the structure.
            Models are renumbered

            Args:
                keep_model: Model number(s) to keep
        """
        models = []
        if '-' in keep_model:
            mod_1, mod_2 = keep_model.split('-')
            for mod_id in range(int(mod_1), int(mod_2) + 1):
                models.append(mod_id)
        elif ',' in keep_model:
            models = [int(mod) for mod in keep_model.split(',')]
        else:
            models = [int(keep_model)]

        ids = [mod.id for mod in self.st.get_models()]
        for md_id in ids:
            if self.st[md_id].serial_num not in models:
                self.st.detach_child(md_id)

        # renumbering models
        for i, mod in enumerate(self.st):
            mod.id = i
            mod.serial_num = i + 1

        self.nmodels = len(self.st)
        self.models_type = mu.guess_models_type(self.st) if self.nmodels > 1 else 0

    def superimpose_models(self):
        '''Superimpose models'''
        spimp = Superimposer()
        if self.nmodels > 1:
            fix_atoms = [
                at
                for at in self.st[0].get_atoms()
                if at.id == 'CA'
            ]
            for mod in self.st.get_models():
                if mod.id == 0:
                    continue
                mov_atoms = [
                    at
                    for at in self.st[mod.id].get_atoms()
                    if at.id == 'CA'
                ]
                spimp.set_atoms(fix_atoms, mov_atoms)
                spimp.apply(self.st[mod.id].get_atoms())
            self.models_type = mu.guess_models_type(self.st) if self.nmodels > 1 else 0
            return True
        return False

    def build_complex(self):
        '''Build a complex with the available models (suitable for biounits)
           Use first available model as base
        '''
        added_chains = 0
        last_chain_id = ''
        for mod in self.st:
            if mod.id == 0:
                for chn in mod:
                    last_chain_id = ord(chn.id)
                continue
            new_chain_ids = [ch.id for ch in mod]
            for ch_id in new_chain_ids:
                new_ch = mod[ch_id].copy()
                last_chain_id += 1
                new_ch.id = chr(last_chain_id)
                new_ch.set_parent(self.st[0])
                self.st[0].add(new_ch)
                added_chains += 1
        self.select('1')
        return added_chains

    def has_models(self) -> bool:
        """ Shortcut method to check whether the structure
            has more than one model

            Returns: Boolean
        """
        return self.nmodels > 1

    def has_superimp_models(self) -> bool:
        """ Shotcut method to check whether the structure has superimposed
            models (i.e. NMR or ensemble)

            Returns: Boolean
        """
        return self.models_type and self.models_type['type'] == mu.ENSM
