# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from Bio.PDB.PDBIO import PDBIO

class PDBIO_extended (PDBIO):
    
    def _get_atom_line( 
        self, 
        atom, 
        hetfield, 
        segid, 
        atom_number, 
        resname, 
        resseq, 
        icode, 
        chain_id, 
        charge="   "):
        
        if hasattr(atom, "charge"):
            charge = atom.charge
            
        
        line = super()._get_atom_line(
            atom, 
            hetfield, 
            segid, 
            atom_number, 
            resname, 
            resseq, 
            icode, 
            chain_id, 
            charge
        )
        # Support for 4 letter residue names (terminals)
        if len(resname) == 4:
            line = line[0:21] + line[22:]
        
        return line
    


