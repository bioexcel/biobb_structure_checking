# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import sys

from Bio.PDB.PDBIO import PDBIO
from Bio.Data.IUPACData import atom_weights 

_ATOM_FORMAT_STRING = ("%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s\n") 
_ATOM_FORMAT_STRING_PDBQT = ("%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f    %6.3f %s\n") 

    
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

        """
            Return an ATOM PDB string (PRIVATE).
            Extended to allow for variable formatting
        """
        
        if hetfield != " ": 
            record_type = "HETATM" 
        else: 
            record_type = "ATOM  " 
        if atom.element: 
            element = atom.element.strip().upper() 
            if element.capitalize() not in atom_weights: 
                raise ValueError("Unrecognised element %r" % atom.element) 
            element = element.rjust(2) 
        else: 
            element = "  " 
        
        name = atom.get_fullname().strip() 
        # Pad atom name if: 
        #     - smaller than 4 characters 
        # AND - is not C, N, O, S, H, F, P, ..., one letter elements 
        # AND - first character is NOT numeric (funky hydrogen naming rules) 
        
        if len(name) < 4 and name[:1].isalpha() and len(element.strip()) < 2: 
            name = " " + name 
   
        altloc = atom.get_altloc() 
        x, y, z = atom.get_coord() 
        bfactor = atom.get_bfactor() 
        occupancy = atom.get_occupancy() 
        
        try: 
            occupancy_str = "%6.2f" % occupancy 
        except TypeError: 
            if occupancy is None: 
                occupancy_str = " " * 6 
                import warnings 
                from Bio import BiopythonWarning 
   
                warnings.warn( 
                   "Missing occupancy in atom %s written as blank" 
                    % repr(atom.get_full_id()), 
                    BiopythonWarning, 
                ) 
            else: 
                raise TypeError( 
                   "Invalid occupancy %r in atom %r" % (occupancy, atom.get_full_id()) 
                ) 
        # Added charges (from res_library and atom_types from data_library)
        # Format PDBQT for Autodock
        if hasattr(atom, "charge"):
            charge = atom.charge
            element = atom.ADT_type
            format = _ATOM_FORMAT_STRING_PDBQT
        else:
            format = _ATOM_FORMAT_STRING
   
        args = ( 
              record_type, 
              atom_number, 
              name, 
              altloc, 
              resname, 
              chain_id, 
              resseq, 
              icode, 
              x, 
              y, 
              z, 
              occupancy_str, 
              bfactor, 
              charge, 
              element
        ) 
        line = format % args 
        
        # Support for 4 letter residue names (terminals)
        if len(resname) == 4:
            line = line[0:21] + line[22:]
        return line
    
    
