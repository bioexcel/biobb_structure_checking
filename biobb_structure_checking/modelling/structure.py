''' Class to manage internal Structure data'''

class StructureData(self):
    def __init__(self, st):
        self.st = st
        self.backbone_links = []
        self.modified_residue_list = []
        self.non_canonical_residue_list = []

        self.hetatm = {}
        self.stats = {
            num_res : 0,
            num_ats : 0,
            res_hetats : 0,
            num_wat : 0,
            num_h : 0,
            res_h : 0,
            res_insc : 0,
            res_ligands : 0
        }
        self.ca_only = False
        self.ss_bonds = []
        self.has_charges = False      
        self.total_charge = None
