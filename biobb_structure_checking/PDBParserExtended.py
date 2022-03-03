""" PDBParser extended with PQR and PDBQT formats """

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.StructureBuilder import StructureBuilder

class StructureBuilderExtended(StructureBuilder):
    def init_atom(self,
        name,
        coord,
        b_factor,
        occupancy,
        altloc,
        fullname,
        serial_number=None,
        element=None,
        pqr_charge=None,
        radius=None,
        is_pqr=False,
        file_type='PDB'):
        residue = self.residue
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_id(name):
            duplicate_atom = residue[name]
            # atom name with spaces of duplicate atom
            duplicate_fullname = duplicate_atom.get_fullname()
            if duplicate_fullname != fullname:
                # name of current atom now includes spaces
                name = fullname
                warnings.warn(
                    "Atom names %r and %r differ only in spaces at line %i."
                    % (duplicate_fullname, fullname, self.line_counter),
                    PDBConstructionWarning,
                )
        if file_type == 'pdb':
            self.atom = Atom(
                name,
                coord,
                b_factor,
                occupancy,
                altloc,
                fullname,
                serial_number,
                element,
            )
        elif file_type == 'pqr':
            self.atom = Atom(
                name,
                coord,
                None,
                None,
                altloc,
                fullname,
                serial_number,
                element,
                pqr_charge,
                radius,
            )
        elif file_type == 'pdbqt':
            self.atom = Atom(
                name,
                coord,
                b_factor,
                occupancy,
                altloc,
                fullname,
                serial_number,
                element,
                pqr_charge,
                pqr_type,
            )            

        if altloc != " ":
            # The atom is disordered
            if residue.has_id(name):
                # Residue already contains this atom
                duplicate_atom = residue[name]
                if duplicate_atom.is_disordered() == 2:
                    duplicate_atom.disordered_add(self.atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.detach_child(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add(disordered_atom)
                    disordered_atom.disordered_add(self.atom)
                    disordered_atom.disordered_add(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn(
                        "WARNING: disordered atom found with blank altloc before "
                        "line %i.\n" % self.line_counter,
                        PDBConstructionWarning,
                    )
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.disordered_add(self.atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add(self.atom)
        
    
    
class PDBParserExtended(PDBParser):
    """ Extends PDB Parser to consider additional variants as PDBQT and PQR """
    def __init__(
        self,
        PERMISSIVE=True,
        get_header=False,
        structure_builder=None,
        QUIET=False,
        is_pqr=False,
        file_type='PDB'
    ):
        if structure_builder in not None:
            self.structure_builder = structure_builder
        else:
            self.structure_builder = StructureBuildExtended()
            
        super().__init__(
            PERMISSIVE=PERMISSIVE,
            get_header=get_header,
            structure_builder=structure_builder,
            QUIET=QUIET,
            is_pqr= file_type is not 'PDB',
        )
        self.file_type=file_type
    
        