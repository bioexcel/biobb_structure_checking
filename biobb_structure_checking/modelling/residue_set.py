''' Class to manage Residue Sets to chain building'''
import re
import sys
import biobb_structure_checking.modelling.utils as mu


class ResidueSetList():
    """ Class to manage Residue sets for chain building"""
    def __init__(self, pairs_list=None, debug=False):
        self.sets = []
        self.n = 0
        if pairs_list:
            self._add_data(pairs_list, debug)
            self._add_meta()

    def _add_meta(self):
        ch_num = 1
        for rset in self._get_sorted_sets():
            rset.set_meta(ch_num)
            ch_num += 1

    def _add_data(self, pairs_list, debug=False):
        """ Build residue set"""
        for atom_pair in pairs_list:
            res1 = atom_pair[0].get_parent()
            res2 = atom_pair[1].get_parent()
            i = self._find(res1)
            j = self._find(res2)
            if i == -1 and j == -1:
                rset = ResidueSet()
                rset.add(res1)
                rset.add(res2)
                self._append_rset(rset)
            elif i != -1 and j != -1 and i != j:
                self.sets[i].union(self.sets[j])
                self._delete(j)
            elif j == -1:
                self.sets[i].add(res2)
            elif i == -1:
                self.sets[j].add(res1)
            if debug:
                for rset in self._get_sorted_sets():
                    print(f"#DEBUG: {rset}")

    def _find(self, item):
        i = 0
        while i < self.n and item not in self.sets[i].items:
            i += 1
        if i == self.n:
            return -1
        else:
            return i

    def _append_rset(self, item):
        self.sets.append(item)
        self.n = len(self.sets)

    def _delete(self, i):
        del self.sets[i]
        self.n = len(self.sets)

    def _get_sorted_sets(self):
        return sorted(self.sets, key=lambda rset: rset.ini)


class ResidueSet():
    """ Class to manage residue sets"""
    def __init__(self):
        self.ini = 999999
        self.fin = 0
        self.inir = None
        self.finr = None
        self.items = set()
        self.type = ''
        self.id = ''

    def set_meta(self, id):
        """ Add molecule type to residue set"""
        if self._is_protein():
            self.type = 'prot'
        else:
            self.type = 'na'
        self.id = id

    def _is_protein(self):
        chain_type = mu.guess_chain_type_list(self._get_residues())['type']
        return chain_type == mu.PROTEIN

    def add(self, res):
        """ Add residue to set"""
        self.items.add(res)
        if res.index < self.ini:
            self.ini = res.index
            self.inir = res
        if res.index > self.fin:
            self.fin = res.index
            self.finr = res

    def union(self, other):
        """ join two residue sets"""
        self.items = self.items.union(other.items)
        if other.ini < self.ini:
            self.ini  = other.ini
            self.inir = other.inir
        if other.fin > self.fin:
            self.fin = other.fin
            self.finr = other.finr

    def get_residue_id_list(self):
        """ Residue list"""
        seq = self._get_residues()
        seql = []
        for i in sorted(seq):
            rid = (
                f"{seq[i].chain}{str(seq[i].resNum)}-{seq[i]._getOneLetterResidueCode()}"
            )
            if seq[i].useModels:
                seql.append(f"{rid}/{str(seq[i].model)}")
            else:
                seql.append(rid)
        return seql

    def _get_residues(self):
        return list(sorted(self.items))

    def get_sequence(self):
        """ get sequence for printint"""
        residue_list = self._get_residues()
        return mu.get_sequence_from_list(
            residue_list,
            mu.guess_chain_type_list(residue_list)['type']
        )

    def __str__(self):
        # residue_list = self._get_residues()
        return f"{self.id} ({mu.residue_id(self.inir)}-{mu.residue_id(self.finr)}) ({self.type}): {self.get_sequence()}"

