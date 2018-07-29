"""
    StructureManager: module to handle structure data.
"""


import sys
import re

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser


import structure_checking.util as util

MODELS_MAXRMS = 5.0    # Threshold value to detect NMR models (angs)

class StructureManager():

    def __init__(self):
        self.model_type = ''
        self.num_ats = 0

    def loadStructure(self, input_pdb_path, use_models, remove_H, debug=False):

        if "pdb:"in input_pdb_path:
            pdbl = PDBList(pdb='tmpPDB')

            try:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path)
                parser = MMCIFParser()
                input_format = 'cif'

            except IOError:
                print ("#ERROR: fetching PDB " + input_pdb_path, file=sys.stderr)
                sys.exit(2)
        else:
            real_pdb_path = input_pdb_path
            if '.pdb' in real_pdb_path:
                parser = PDBParser(PERMISSIVE=1)
                input_format = 'pdb'
            elif '.cif' in real_pdb_path:
                parser = MMCIFParser()
                input_format = 'cif'
            else:
                print ('#ERROR: unknown filetype', file=sys.stderr)
                sys.exit(2)
        try:
            self.st = parser.get_structure('st', real_pdb_path)

        except OSError:
            print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.exit(2)

        #====== Internal residue renumbering =========================================
        i = 1
        for r in self.st.get_residues():
            r.index = i
            i += 1

        #Atom renumbering for mmCIF,
        if input_format == 'cif':
            i = 1
            for at in self.st.get_atoms(): # Check numbering in models
                at.serial_number = i
                if hasattr(at, 'selected_child'):
                    at.selected_child.serial_number = i
                i += 1

        for at in self.st[0].get_atoms():
            self.num_ats += 1

        #checking models type
        if len(self.st) > 1:

            if use_models == 'no':
                print ("WARNING: Input Structure contains models, but using only first one due to use_models settings")
                remove_models = True

            elif use_models == 'auto':
                if debug:
                    print ("DEBUG: Found " + str(len(self.st)) + " models")
                    print ("DEBUG: RMS " + str(util.calcRMSdAll(self.st[0], self.st[1])))

                if util.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    if debug:
                        print ("DEBUG: Models look like alternative conformations, will consider only one")
                    self.model_type = 'alt'
                    print ("WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    remove_models = True

                else:
                    self.model_type = 'traj'
                    remove_models = False

            elif use_models == 'force':
                if util.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    print ('#WARNING: Models found look like NMR models, but using all due to use_models = force')
                remove_models = False

            else:
                print ("#ERROR: Unknown use_models option", file=sys.stderr)
                sys.exit(1)

            if remove_models:
                print ("Removing models")
                ids = []

                for md in self.st.get_models():
                    ids.append(md.id)

                for i in range(1, len(ids)):
                    self.st.detach_child(ids[i])
        # Hydrogens

        if remove_H == 'all':
            print ("Removing H atoms")

            for r in self.st.get_residues():
                util.removeHFromRes(r)
                
    def get_structure(self):
        return self.st

    def saveStructure(self, output_pdb_path):

        pdbio = PDBIO()
        pdbio.set_structure(self.st)

        try:
            pdbio.save(output_pdb_path)

        except OSError:
            print ("#ERROR: unable to save PDB data on " + output_path, file=sys.stderr)

    def get_all_CA_distances(self):
        dist_mat = []
        # get CA Ats
        ca_ats = []
        for at in self.st.get_atoms():
            if at.id == 'CA':
                ca_ats.append(at)
        for i in range(0, len(ca_ats)-1):
            for j in range(i+1, len(ca_ats)):
                dist_mat.append ([ca_ats[i], ca_ats[j], ca_ats[i]-ca_ats[j]])
        return dist_mat
    
    def get_all_distances(self):
        dist_mat = []
        for at1 in self.st.get_atoms():
            for at2 in self.st.get_atoms():
                if at1.serial_number < at2.serial_number:
                    dist_mat.append ([at1, at2, at1-at2])
        return dist_mat
    
    def get_nmodels(self):
        return len(self.st)
    
    def select_model(self,nm):
        self.st = self.st[nm-1]
    
    def get_chain_ids(self):
        chain_ids=[]
        for ch in st.get_chains():
            chain_ids.append(ch.id)
        return chain_ids

    def select_chains(self, select_chain):
        chain_ids = self.get_chain_ids
        ch_ok = select_chains.split(',')
        for ch in ch_ok:
            if not ch in chain_ids:
                print ("Error: request chain not present", ch, file=sys.stderr)
                select_chains=''
        for ch in chain_ids:
            if ch not in ch_ok:
                self.st[0].detach_child(ch)
    def get_altloc_residues(self):
        res_list={}
        for at in self.st.get_atoms():
            r = at.get_parent()
            rid = util.residueid(r)
            if at.get_altloc() != ' ':
                if rid not in res_list:
                    res_list[rid]=[]
                res_list[rid].append(at)
        return res_list

    def select_altloc_residues(self, r, select_altloc):
        alt_loc_res = self.get_altloc_residues()
        for at in alt_loc_res[r]:
            res = at.get_parent()
            if select_altloc == 'occupancy':
                newat = at.selected_child
            else:
                newat = at.child_dict[select_altloc]
                newat.disordered_flag=0
                newat.altloc=' '
                res.detach_child(at.id)
                res.add (newat)
            res.disordered=0
    def get_metals(self, metal_ats):
        met_list=[]
        met_rid=[]
        at_groups={}
        for at in self.st.get_atoms():
            if not re.match('H_', at.get_parent().id[0]):
                continue
            if at.id in metal_ats:
                met_list.append(at)
        
        for at in met_list:
            r = at.get_parent()
            rid = r.get_parent().id + str(r.id[1])
            met_rid.append(rid)
            if at.id not in at_groups.keys():
                at_groups[at.id]=[]
            at_groups[at.id].append(rid)
        
        return [met_list, met_rid, at_groups]
    
    #TODO refine this as remove_residue or remove_atom_set
    def remove_metals (self, met_list, to_remove ):
        for at in met_list:
            r = at.get_parent()
            rid = r.get_parent().id + str(r.id[1])
            if rid in to_remove:
                r.detach_child(at.id)
            