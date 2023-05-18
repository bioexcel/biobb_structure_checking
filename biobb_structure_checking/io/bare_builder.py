""" Consumer class that builds a Structure object,
    ignoring Chain labels and Residue Ids.
    User to recover coordinates from incorrect PDB/mmCif files
    This is used by the PDBParser and MMCIFparser classes.
"""

from Bio.PDB.StructureBuilder import StructureBuilder

# SMCRA hierarchy
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue


class BareStructureBuilder(StructureBuilder):
    """ Deals with constructing the Structure object.
        Discards original chain labels and residue ids.
        Used to rebuild faulty structure files
    """

    def __init__(self):
        """Initialize the class."""
        StructureBuilder.__init__(self)
        self.residue_counter = 0

    def init_chain(self, chain_id):
        """Create a new Chain object with given id.

        Arguments:
         - chain_id - string
        """
        # All residues loaded in chains ATOM:"0", HETATM:"h", or WAT:"w"
        if chain_id not in ['0', 'h', 'w']:
            chain_id = '0'
        if self.model.has_id(chain_id):
            self.chain = self.model[chain_id]
        else:
            self.chain = Chain(chain_id)
            self.model.add(self.chain)

    def init_residue(self, resname, field, resseq, icode):
        """Create a new Residue object.

        Arguments:
         - resname - string, e.g. "ASN"
         - field - hetero flag, "W" for waters, "H" for
           hetero residues, otherwise blank.
         - resseq - int, sequence identifier (ignored)
         - icode - string, insertion code

        """
        if field == "H":
            # The hetero field consists of H_ + the residue name (e.g. H_FUC)
            field = "H_" + resname
            self.init_chain('h')
        elif field == "W":
            self.init_chain('w')
        else:
            self.chain = self.model['0']

        # Set consecutive residue ids
        self.residue_counter += 1
        res_id = (field, self.residue_counter, icode)
        self.residue = Residue(res_id, resname, self.segid)
        self.chain.add(self.residue)
