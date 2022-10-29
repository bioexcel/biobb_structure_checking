''' Class to manage internal Structure data'''
import biobb_structure_checking.model_utils as mu

class StructureData():
    def __init__(self, st, input_format):
        self.st = st
        self.input_format = input_format
        self.backbone_links = []
        self.modified_residue_list = []
        self.non_canonical_residue_list = []
        
        self.hetatm = {}
        self.stats = {
            'num_res' : 0,
            'num_ats' : 0,
            'res_hetats' : 0,
            'num_wat' : 0,
            'num_h' : 0,
            'res_h' : 0,
            'res_insc' : 0,
            'res_ligands' : 0
        }
        self.ca_only = False
        self.ss_bonds = []
        self.has_charges = False      
        self.total_charge = None
        self.all_residues = []
        self.non_canonical_residue_list = []

    def residue_renumbering(self, data_library):
        """Sets the Bio.PDB.Residue.index attribute to residues for a unique,
        consecutive, residue number, and generates the corresponding list
        in **all_residues**
        """
        i = 1
        self.all_residues = []
        self.non_canonical_residue_list = []
        for res in self.st.get_residues():
            res.index = i
            res.resname = res.resname.strip()
            if type(res).__name__ == 'DisorderedResidue':
                for ch_r in res.child_dict:
                    res.child_dict[ch_r].index = i
            self.all_residues.append(res)
            if res.get_resname() in data_library.canonical_codes:
                self.non_canonical_residue_list.append(
                    {
                        'res':res, 
                        'can_res':data_library.canonical_codes[res.get_resname()], 
                        'new_res':res.resname
                    }
                )
            
            i += 1

    def atom_renumbering(self):
        """ Sets  Bio.PDB.Atom.serial_number for all atoms in the structure,
            overriding original if any.
        """
        if self.input_format == 'pqr' or self.has_charges:
            self.total_charge = 0.
        else:
            self.total_charge = None
        i = 1
        for atm in self.st.get_atoms():
            atm.serial_number = i
            if hasattr(atm, 'selected_child'):
                atm.selected_child.serial_number = i
            if atm.pqr_charge is not None and self.total_charge is not None:
                self.total_charge += atm.pqr_charge
            i += 1

    def check_backbone_connect(self, backbone_atoms, covlnk):
        """
        Determines backbone links usign a distance criterium and produces a dict with
        link relationships in the N-term to C-term direction

        Args:
            backbone_atoms: atoms to be considered as backbone
            covlnk: Threshold distance for a covalent bond
        """
        self.backbone_links = mu.get_backbone_links(
            self.st, backbone_atoms, covlnk
        )
        self.next_residue = {}
        self.prev_residue = {}
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            res1 = at1.get_parent()
            res2 = at2.get_parent()
            self.prev_residue[res2] = res1
            self.next_residue[res1] = res2

    def calc_stats(self):
        """ Calculates general statistics about the structure, and guesses whether
            it is a CA-only structure
        """
        self.stats['num_res'] = 0
        self.stats['num_ats'] = 0
        self.stats['res_hetats'] = 0
        self.stats['num_wat'] = 0
        self.stats['res_insc'] = 0
        self.stats['num_h'] = 0
        self.stats['res_h'] = 0
        for res in self.st.get_residues():
            self.stats['num_res'] += 1
            if mu.is_wat(res):
                self.stats['num_wat'] += 1
            if mu.is_hetatm(res):
                self.stats['res_hetats'] += 1
            if mu.has_ins_code(res):
                self.stats['res_insc'] += 1
            self.stats['num_ats'] += len(res.get_list())
        self.stats['res_ligands'] = self.stats['res_hetats'] - self.stats['num_wat']
        for pair in mu.get_residues_with_H(self.st):
            self.stats['res_h'] += 1
            self.stats['num_h'] += pair['num_h']
        # Detecting whether it is a CA-only structure
        # num_ats should be much larger than num_res
        # waters removed
        # Taking polyGly as a lower limit
        self.ca_only = self.stats['num_ats'] - self.stats['num_wat'] < (self.stats['num_res'] - self.stats['num_wat']) * 4
