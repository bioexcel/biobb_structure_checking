"""
  Module to load data library json file
"""

import json
import sys

class DataLibManager:
    """ 
    | data_lib_manager DataLibManager
    | Manages projects' global data file

    Args:
        file_path (str) : Path to library json file
    """

    def __init__(self, file_path):
        try:
            data_file_h = open(file_path)
            json_map = json.load(data_file_h)
            self.residue_codes = json_map['data_library']['residue_codes']
            self.canonical_codes = json_map['data_library']['canonical_codes']
            self.atom_data = json_map['data_library']['atom_data']
            self.residue_data = json_map['data_library']['residue_data']
            self.distances = json_map['data_library']['distances']
            self.ion_res = json_map['data_library']['addH_check_residues']
            self.std_ion = json_map['data_library']['addH_std_ion']
            self.ff_data = {}
            for ff in json_map['data_library']['atom_type_ff']:
                self.ff_data[ff] = {}

        except IOError:
            print("ERROR: unable to open data library " + file_path, file=sys.stderr)
            sys.exit(2)

    def get_valid_codes(self, mol_type):
        """ DataLibManager.get_valid_codes
        Obtain valid residue codes 
        
        Args:
            mol_type (str - One of valid molecular types) :
                * **na** - dna + rna, 
                * **dna** -
                * **rna** - 
                * **protein** -
        """
        if mol_type == 'na':
            codes = self.residue_codes['dna'] + self.residue_codes['rna']
        else:
            codes = self.residue_codes[mol_type]
        return codes

    def get_all_atom_lists(self):
        """ DataLibManager.get_all_atom_lists
        Obtain lists of atoms per protein residue.
        """
        atom_lists = {
            rcode: {
                'backbone': self.residue_data['*']['bck_atoms'],
                'side': self.residue_data[rcode]['side_atoms']
            }
            for rcode in self.residue_codes['protein']
        }

        for rcode in self.residue_codes['cap_residues']:
            atom_lists[rcode]['backbone'] = self.residue_data[rcode]['bck_atoms']

        return atom_lists

    def get_atom_feature_list(self, feature):
        """ DataLibManager.get_atom_feature_list
        Get a residue list with a specific section of data.
        
        Args:
            feature (str - Options) :
                * **bck_atoms** - Backbone atoms
                * **side_atoms** - Side chain atoms
                * **apolar_atoms** - Apolar atoms
                * **polar_acceptor_atoms** - H-bond acceptor polar atoms
                * **polar_donor_atoms** - H-bond donor polar atoms
                * **positive_atoms** - Atoms bearing formal positive charge
                * **negative_atoms** - Atoms bearing formal negative charge
                * **hydrogen_atoms** - Hydrogen atoms
                * **addH_rules** - Rules for adding Hydrogen atoms
                * **mutation_rules** - Rules to mutate to other residue types

        """
        f_list = {
            rcode: self.residue_data[rcode][feature]
            for rcode in self.residue_data
            if feature in self.residue_data[rcode]
        }

        for rcode in self.canonical_codes:
            if feature in self.residue_data[self.canonical_codes[rcode]]:
                f_list[rcode] = self.residue_data[self.canonical_codes[rcode]][feature]

        if '*' not in f_list:
            f_list['*'] = []
        return f_list

    def get_chiral_data(self):
        """ DataLibManager.get_chiral_data
        Get data related to chiral atoms.
        """
        return self.get_atom_feature_list('chiral_atoms')

    def get_hydrogen_atoms(self):
        """ DataLibManager.get_hydrogen_atoms
        Get list of hydrogen atoms per residue. 
        """
        return self.get_atom_feature_list('hydrogen_atoms')

    def get_add_h_rules(self):
        """ DataLibManager.get_add_h_rules
        Get rules for adding Hydrogen atoms to residues. 
        """
        return self.get_atom_feature_list('addH_rules')

    def get_atom_lists(self, contact_types):
        """ DataLibManager.get_atom_lists
        Get a list of atoms organized per contact types. 
        
        Args:
            contact_types (list(str) - List of types of contacts. Options) :
                * **apolar** - Involving any apolar atom
                * **polar_donor** - Involving two Hbond donors
                * **polar_acceptor** - Involving two HBond acceptors
                * **positive** - Involving two positive atoms
                * **negative** - Involving two negative atoms
        """

        atom_lists = {
            cls_type: self.get_atom_feature_list(cls_type + '_atoms')
            for cls_type in contact_types
            if cls_type != 'severe'
        }

        return atom_lists

    def get_amide_data(self):
        """ DataLibManager.get_amide_data
        Get atoms related to amide issues
        """
        alist = []
        rlist = {}
        for rcode in self.residue_data:
            if 'amide_atoms' in self.residue_data[rcode]:
                alist += self.residue_data[rcode]['amide_atoms']
                rlist[rcode] = self.residue_data[rcode]['amide_atoms']
        return rlist, alist

    def get_mutation_map(self):
        """ DataLibManager.get_mutation_map
        Get the complete map of mutation rules per residue.
        """
        mut_rules = {}
        for rcode in self.residue_data:
            if 'mutation_rules' in self.residue_data[rcode]:
                mut_rules[rcode] = self.residue_data[rcode]['mutation_rules']
                mut_rules[rcode]['side_atoms'] = self.residue_data[rcode]['side_atoms']
        return mut_rules

    def get_canonical_resname(self, rcode):
        """ DataLibManager.get_canonical_resname
        Get parent residue names for modified residues 
        
        Args:
            rcode (str) : Original residue code (3 letter for aminoacids)
        """
        if rcode in self.canonical_codes:
            return self.canonical_codes[rcode]
        return rcode

    def get_ff_data(self, file_path):
        """ DataLibManager.get_ff_data
        Load Forcefield data from external file
       
        Args:
            file_path (str) : Path to library file
        """
        try:
            data_file_h = open(file_path)
            json_map = json.load(data_file_h)
            self.ff_data[json_map['id']] = json_map
        except IOError:
            print("ERROR: unable to open forcefield library " + file_path, file=sys.stderr)
            sys.exit(2)
