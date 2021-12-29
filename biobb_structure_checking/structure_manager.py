"""Module to manage structure, based on BioPython Bio.PDB
"""
import warnings
import os
import sys
#import re
from typing import List, Dict, Tuple, Iterable, Mapping, Union, Set


from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio import BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionException

from biobb_structure_checking.mmb_server import MMBPDBList
from biobb_structure_checking.mutation_manager import MutationManager, MutationSet
from biobb_structure_checking.data_lib_manager import DataLibManager
from biobb_structure_checking.residue_lib_manager import ResidueLib
from biobb_structure_checking.sequence_manager import SequenceData
from biobb_structure_checking.PDBIO_extended import PDBIO_extended
import biobb_structure_checking.model_utils as mu

CISTHRES = 20  # TODO check vaules with pdb checking
TRANSTHRES = 160

ALL_CONTACT_TYPES = [
    'severe',
    'apolar',
    'polar_acceptor',
    'polar_donor',
    'positive',
    'negative'
]
AMIDE_CONTACT_TYPES = [
    'polar_acceptor',
    'polar_donor',
]


class StructureManager:
    """Main Class wrapping Bio.PDB structure object
    """

    def __init__(
            self,
            input_pdb_path: str,
            data_library_path: str,
            res_library_path: str,
            pdb_server: str = 'ftp://ftp.wwpdb.org',
            cache_dir: str = 'tmpPDB',
            file_format: str = 'mmCif',
            fasta_sequence_path: str = ''
        ) -> None:
        """
            Class constructor. Sets an empty object and loads a structure
            according to parameters

        Args:
            **input_pdb_path** (str): path to input structure either in pdb or mmCIF format. Format is taken from file extension. Alternatively **pdb:pdbId** fetches the mmCIF file from RCSB
            **data_library_path** (str): Path to json data library
            **res_library_path** (str): Path to residue library
            **pdb_server** (str): **default** for Bio.PDB defaults (RCSB), **mmb** for MMB PDB API
            **cache_dir** (str): path to temporary dir to store downloaded structures
            **file_format** (str): structure file format to use
            **fasta_sequence_path** (str): path to canonical sequence file (needed for PDB input)

        Object structure:
            {
                "input_format" (str): Format of input file pdb|cif
                "models_type" (int): Guessed model type when num. models > 1 NMR|Bunit
                "num_ats" (int): Total Number of atoms
                "nmodels" (int): Number of Models
                "chain_ids" (List): Chain composition as [chain_id:{}]
                "modified" (Boolean): Flag to indicated that structure has been modified
                "all_residues" (list): List of pointer to Bio.PDB.Residue objects, ordered acording to input file
                "num_res" (int): Number of residues
                "res_insc" (int): Number of residues with insertion codes
                "res_hetats" (int): Number of residues flagged as HETATM
                "res_ligands" (int): Number of non water residues flagged as HETATM
                "num_wat" (int): Number of water residues
                "ca_only" (boolead): Flag to indicate a possible CA-only structure
                "modified_residue_list" []: List of residues being connected HETATM as PDB.Bio.Residue
                "backbone_links" []: List of found canonical backbone links as [at1, at2] tuples, according to a distance criterium
            TODO Update and complete
        """

        self.chain_ids = {}
        self.chain_details = {}

        self.backbone_links = []
        self.modified_residue_list = []

        self.hetatm = {}
        self.num_res = 0
        self.num_ats = 0
        self.res_hetats = 0
        self.num_wat = 0
        self.res_insc = 0
        self.res_ligands = 0
        self.ca_only = False

        self.ss_bonds = []

        self.meta = {}

        self.all_residues = []
        self.next_residue = {}
        self.prev_residue = {}

        self.sequence_data = SequenceData()

        self.modified = False
        self.biounit = False
        self.fixed_side = False
        self.file_format = file_format
        self.has_charges = False

        self.data_library = DataLibManager(data_library_path)
        for ff in self.data_library.ff_data:
            self.data_library.get_ff_data(os.path.dirname(data_library_path) + '/' + ff  + '_prm.json')

        self.res_library = ResidueLib(res_library_path)


        self.input_format = self._load_structure_file(
            input_pdb_path, cache_dir, pdb_server, file_format
        )

        if fasta_sequence_path:
            self.sequence_data.load_sequence_from_fasta(fasta_sequence_path)

        # Checking models type according to RMS among models
        self.nmodels = len(self.st)
        self.models_type = mu.guess_models_type(self.st) if self.nmodels > 1 else 0

        # Calc internal data
        self.update_internals(cif_warn=True)

    def _load_structure_file(self, input_pdb_path, cache_dir, pdb_server, file_format):
        """ Load structure file """
        if "pdb:" in input_pdb_path:
            # MMBPDBList child defaults to Bio.PDB.PDBList if MMB server is not selected
            pdbl = MMBPDBList(pdb=cache_dir, server=pdb_server)
            if '.' in input_pdb_path:
                [pdbid, biounit] = input_pdb_path.split('.')
                input_pdb_path = pdbid[4:].upper()
                if pdb_server != 'mmb':
                    raise WrongServerError
                real_pdb_path = pdbl.retrieve_pdb_file(
                    input_pdb_path, file_format='pdb', biounit=biounit
                )
                self.biounit = biounit
            else:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(
                    input_pdb_path, file_format=self.file_format
                )
                if file_format == 'pdb':
                    # change file name to id.pdb
                    os.rename(real_pdb_path, input_pdb_path + ".pdb")
                    real_pdb_path = input_pdb_path + ".pdb"
        else:
            real_pdb_path = input_pdb_path

        if '.pdb' in real_pdb_path: 
            parser = PDBParser(PERMISSIVE=1, is_pqr=False)
            input_format = 'pdb'
        elif '.pqr' in real_pdb_path: 
            parser = PDBParser(PERMISSIVE=1, is_pqr=True)
            input_format = 'pqr'
        elif '.cif' in real_pdb_path:
            parser = MMCIFParser()
            input_format = 'cif'
        else:
            raise UnknownFileTypeError(input_pdb_path)

        warnings.simplefilter('ignore', BiopythonWarning)

        try:
            self.st = parser.get_structure('st', real_pdb_path)
        except ValueError as err:
            raise ParseError('ValueError', err)
        except PDBConstructionException as err:
            raise ParseError('PDBBuildError', err)
        if input_format in ['pdb', 'pqr']:
            self.headers = parse_pdb_header(real_pdb_path)
        else:
            self.headers = MMCIF2Dict(real_pdb_path)
            
        return input_format

    def update_internals(self, cif_warn: bool = False):
        """ Update internal data when structure is modified """
        # Add .index field for correlative, unique numbering of residues
        self.residue_renumbering()
        # Atom renumbering for mmCIF, PDB uses atom number in file
        self.atom_renumbering()
        self.set_chain_ids()
        self.calc_stats()
        self.guess_hetatm()

        self.rr_dist = self.get_all_r2r_distances('all', join_models=False)

        # Precalc backbone . TODO Nucleic Acids
        self.check_backbone_connect(
            ('N', 'C'),
            self.data_library.distances['COVLNK']
        )
        # get canonical and structure sequences
        self.sequence_data.read_sequences(self, clean=True, cif_warn=cif_warn)

    def residue_renumbering(self):
        """Sets the Bio.PDB.Residue.index attribute to residues for a unique,
        consecutive, residue number, and generates the corresponding list
        in **all_residues**
        """
        i = 1
        self.all_residues = []
        for res in self.st.get_residues():
            res.index = i
            res.resname = res.resname.strip()
            if type(res).__name__ == 'DisorderedResidue':
                for ch_r in res.child_dict:
                    res.child_dict[ch_r].index = i
            self.all_residues.append(res)
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
        #self.has_charges = (self.total_charge is not None)

    def update_atom_charges(self, ff):
        """ Update atom charges and types from data library """
        
        print("Updating partial charges and atom types")
        
        self.total_charge = 0.
        
        if ff not in self.data_library.ff_data:
            raise UnknownFFError(ff)
        ff_data = self.data_library.ff_data[ff]

        for res in self.st.get_residues():
            rcode = res.get_resname()
            if len(rcode) == 4:
                rcode3 = rcode[1:]
            else:
                rcode3 = rcode
            can_rcode3 = self.data_library.get_canonical_resname(rcode3)
            if rcode not in self.res_library.residues:
                print("Warning: {} not found in residue library atom charges set to 0.".format(rcode))
                for atm in res.get_atoms():
                    atm.pqr_charge = 0.
                    atm.radius = 0.
                    if atm.id in self.data_library.atom_data['metal_atoms']:
                        atm.xtra['atom_type'] = atm.id.lower().capitalize()
                    else:
                        atm.xtra['atom_type'] = atm.element
            else:
                oxt_ok = rcode[0] != 'C' or len(rcode) != 4
                res_chr = 0.
                
                for atm in res.get_atoms():
                    atm.pqr_charge = self.res_library.get_atom_def(rcode, atm.id).chrg
                    if atm.id in ff_data['residue_data'][can_rcode3]:
                        atm.xtra['atom_type'] = ff_data['residue_data'][can_rcode3][atm.id]
                    elif atm.id in ff_data['residue_data']['*']:
                        atm.xtra['atom_type'] = ff_data['residue_data']['*'][atm.id]
                    else:
                        atm.xtra['atom_type'] = atm.element
                    atm.radius = ff_data['rvdw'][atm.xtra['atom_type']]
                    res_chr += atm.pqr_charge
                    self.total_charge += atm.pqr_charge
                    if atm.id == 'OXT':
                        oxt_ok = True
                if not oxt_ok:
                    print("Warning: OXT atom missing in {}. Run backbone --fix_atoms first".format(mu.residue_id(res)))
        print("Total assigned charge: {:10.2f}".format(self.total_charge))

        self.has_charges = True


    def guess_hetatm(self):
        """ Guesses HETATM type as modified res, metal, wat, organic
        """
        self.hetatm = {}
        for typ in [mu.UNKNOWN, mu.MODRES, mu.METAL, mu.ORGANIC, mu.COVORGANIC, mu.WAT]:
            self.hetatm[typ] = []
        for res in self.st.get_residues():
            if not mu.is_hetatm(res):
                continue
            if mu.is_wat(res):
                self.hetatm[mu.WAT].append(res)
            elif len(res) == 1:
                self.hetatm[mu.METAL].append(res)
            elif 'N' in res or 'C' in res:
                # modified aminoacid candidate, TODO check connectivity with n-1 or n+1
                self.hetatm[mu.MODRES].append(res)
                # TODO check modified nucleotides
            else:
                self.hetatm[mu.ORGANIC].append(res)

    def calc_stats(self):
        """ Calculates general statistics about the structure, and guesses whether
            it is a CA-only structure
        """
        self.num_res = 0
        self.num_ats = 0
        self.res_hetats = 0
        self.num_wat = 0
        self.res_insc = 0
        for res in self.st.get_residues():
            self.num_res += 1
            if mu.is_wat(res):
                self.num_wat += 1
            if mu.is_hetatm(res):
                self.res_hetats += 1
            if mu.has_ins_code(res):
                self.res_insc += 1
            self.num_ats += len(res.get_list())
        self.res_ligands = self.res_hetats - self.num_wat
        # Detecting whether it is a CA-only structure
        # num_ats should be much larger than num_res
        # waters removed
        # Taking polyGly as a lower limit
        self.ca_only = self.num_ats - self.num_wat < (self.num_res - self.num_wat)*4

    def get_ins_codes(self) -> List[Residue]:
        """Makes a list with residues having insertion codes"""
        return [
            res
            for res in self.st.get_residues()
            if mu.has_ins_code(res)
        ]

    def get_metal_atoms(self) -> List[Atom]:
        """ Makes a list of possible metal atoms"""
        return mu.get_metal_atoms(self.st, self.data_library.atom_data['metal_atoms'])

    def get_SS_bonds(self) -> List[Union[Atom, Atom, float]]:
        """ Stores and returns possible SS Bonds by distance"""
        self.ss_bonds = mu.get_all_at2at_distances(
            self.st,
            'SG',
            self.data_library.distances['SS_DIST'],
            not self.has_superimp_models()
        )
        return self.ss_bonds

    def check_chiral_sides(self) -> Dict[List[Residue], List[Residue]]:
        """ Returns a list of wrong chiral side chains"""
        chiral_res = self.data_library.get_chiral_data()
        chiral_list = [
            res
            for res in self.st.get_residues()
            if res.get_resname() in chiral_res
        ]

        if not chiral_list:
            return {}

        return {
            'list': chiral_list,
            'res_to_fix': [
                res
                for res in chiral_list
                if not mu.check_chiral_residue(res, chiral_res)
            ]
        }

    def get_chiral_bck_list(self) -> Dict[List[Residue], List[Residue]]:
        """ Returns a list of residues with chiral CAs"""
        prot_chains = 0
        chiral_bck_list = []
        for chn in self.st.get_chains():
            if self.chain_ids[chn.id] == mu.PROTEIN:
                prot_chains += 1
                for res in chn.get_residues():
                    if res.get_resname() != 'GLY' and not mu.is_hetatm(res):
                        chiral_bck_list.append(res)

        if not prot_chains:
            print("No protein chains detected, skipping")
            return {}

        if not chiral_bck_list:
            print("No residues with chiral CA found, skipping")
            return {}

        return {
            'list': chiral_bck_list,
            'res_to_fix': [res for res in chiral_bck_list if not mu.check_chiral_ca(res)]
        }

    def check_r_list_clashes(self, residue_list: Iterable[Residue], contact_types: Iterable[str]) -> Dict[str, Dict[str, Tuple[Residue, Residue, float]]]:
        """ Checks clashes originated by a list of residues"""
        return mu.check_r_list_clashes(
            residue_list,
            self.rr_dist,
            self.data_library.distances['CLASH_DIST'],
            self.data_library.get_atom_lists(contact_types),
            not self.has_superimp_models(),
            severe='severe' in contact_types
        )

    def check_missing_atoms(self) -> List[Tuple[Residue, Dict[str, List[str]]]]:
        """ Makes a **list of missing atoms** in the structure

            Returns:
                List of residues with missing atoms, as a tuples
                ["r",{"backbone":[atoms_list],"side":[atoms_list]}]
        """
        # TODO Nucleic acids
        valid_codes = self.data_library.get_valid_codes('protein')
        residue_data = self.data_library.get_all_atom_lists()
        miss_at_list = []
        for res in self.st.get_residues():
            if res.get_resname() in valid_codes and not mu.is_hetatm(res):
                miss_at = mu.check_all_at_in_r(
                    res, residue_data[res.get_resname().replace(' ', '')]
                )
                if self.is_C_term(res) and res.get_resname() != 'NME' and 'OXT' not in res:
                    if 'backbone' not in miss_at:
                        miss_at['backbone'] = []
                    miss_at['backbone'].append('OXT')
                if miss_at:
                    miss_at_list.append((res, miss_at))
        return miss_at_list

    def check_extra_atoms(self) -> List[Tuple[Residue, Atom]]:
        """ Makes a **list of extra atoms** in the structure

            Returns:
                List of residues with extra atoms, as a tuples
                ["r",atoms_list]
        """
        # TODO Nucleic acids
        valid_codes = self.data_library.get_valid_codes('protein')
        residue_data = self.data_library.get_all_atom_lists()

        extra_at_list = []
        for res in self.st.get_residues():
            if res.get_resname() in valid_codes and not mu.is_hetatm(res):
                extra_ats = mu.check_unk_at_in_r(
                    res, residue_data[res.get_resname().replace(' ', '')]
                )
                if extra_ats:
                    extra_at_list.append((res, extra_ats))
        return extra_at_list

    def get_missing_atoms(self, fragment: str) -> List[Tuple[Residue, List[str]]]:
        """ Returns list of residues with missing atoms

            Returns:
                List of residues with missing atoms, as a tuples like
                ["r",[atoms_list]]
        """
        miss_ats = []
        for res_at in self.check_missing_atoms():
            res, at_list = res_at
            if fragment == 'side':
                if 'side' in at_list and 'N' in res and 'CA' in res and 'C' in res:
                    miss_ats.append((res, at_list['side']))
            else:
                if at_list['backbone']:
                    miss_ats.append((res, at_list['backbone']))
        return miss_ats

    def get_ion_res_list(self) -> List[Tuple[Residue, List[str]]]:
        """
            returns list of residues with potencial selection on adding H

            Returns:
                List of residues that require selection on adding H
                ["r",[atom_list]]
        """
        ion_res = self.data_library.ion_res
        hydrogen_lists = self.data_library.get_hydrogen_atoms()
       
        ion_res_list = []
        for res in self.all_residues:
            rcode = res.get_resname()
            if len(rcode) == 4:
                rcode = rcode[1:]
            if rcode in ion_res:
                ion_res_list.append((res, hydrogen_lists[rcode]))

        return ion_res_list

    def get_backbone_breaks(self) -> Dict[str, List[List[Residue]]]:
        """
            Determines several backbone anomalies
            1. Modified residues
            2. Backbone breaks
            3. Backbone links not corresponding to sequence
            4. Too long residue links
        """
        bck_breaks_list = []
        wrong_link_list = []
        not_link_seq_list = []
        self.modified_residue_list = []
        for i in range(0, len(self.all_residues)-1):
            res1 = self.all_residues[i]
            res2 = self.all_residues[i+1]
            if not mu.same_chain(res1, res2):
                continue
            if mu.is_hetatm(res1) or mu.is_hetatm(res2):
                if res1 in self.next_residue:
                    if self.next_residue[res1] == res2:
                        if mu.is_hetatm(res1):
                            self.modified_residue_list.append(res1)
                else:
                    continue
            # Skip NA
            # TODO include NA Backbone
            if self.chain_ids[res1.get_parent().id] != mu.PROTEIN:
                continue

            if res1 not in self.next_residue:
                bck_breaks_list.append([res1, res2])
                if mu.seq_consecutive(res1, res2):
                    dist = 0.
                    if 'N' in res1 and 'C' in res2:
                        dist = res1['N'] - res2['C']
                    not_link_seq_list.append([res1, res2, dist])

            else:
                if res2 != self.next_residue[res1]:
                    wrong_link_list.append([res1, self.next_residue[res1], res2])
        return {
            'bck_breaks_list': bck_breaks_list,
            'wrong_link_list': wrong_link_list,
            'not_link_seq_list': not_link_seq_list
        }

    def check_backbone_connect(self, backbone_atoms: Iterable[Atom], covlnk: float):
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

#    def check_cis_backbone(self) -> Tuple[List[Tuple[Residue, Residue, float]], List[Residue, Residue, float]]:
    def check_cis_backbone(self):

        """
        Determines omega dihedrals for two bound residues and classifies them
        as normal trans, low angle trans, and cis
        """
        cis_backbone_list = []
        lowtrans_backbone_list = []
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            res1 = at1.get_parent()
            res2 = at2.get_parent()
            if 'CA' in res1 and 'C' in res1 and 'CA' in res2 and 'N' in res2:
                dih = mu.calc_bond_dihedral(res1['CA'], res1['C'], res2['N'], res2['CA'])
                if abs(dih) < CISTHRES:
                    cis_backbone_list.append((res1, res2, dih))
                elif abs(dih) < TRANSTHRES:
                    lowtrans_backbone_list.append((res1, res2, dih))
        return cis_backbone_list, lowtrans_backbone_list

    def check_amide_contacts(self) -> Dict[str, Union[Set[Residue], List[Residue]]]:
        """ Check close contacts involving amide atoms """
        amide_res, amide_atoms = self.data_library.get_amide_data()

        amide_list = set(
            res for res in self.st.get_residues()
            if res.get_resname() in amide_res
        )

        if not amide_list:
            return {}

        c_list = self.check_r_list_clashes(
            amide_list,
            AMIDE_CONTACT_TYPES
        )
        amide_res_to_fix = []
        amide_cont_list = []
        for cls in c_list:
            for rkey in c_list[cls]:
                [at1, at2] = c_list[cls][rkey][0:2]
                res1 = at1.get_parent()
                res2 = at2.get_parent()
                add_pair = False
                if at1.id in amide_atoms and res1 in amide_list:
                    amide_res_to_fix.append(res1)
                    add_pair = True
                if at2.id in amide_atoms and res2 in amide_list:
                    amide_res_to_fix.append(res2)
                    add_pair = True
                if add_pair:
                    amide_cont_list.append(c_list[cls][rkey])
        return {
            'list':  amide_list,
            'res_to_fix': amide_res_to_fix,
            'cont_list': amide_cont_list,
        }

    def get_stats(self) -> Dict[str, Union[int, Union[int, Dict[str, Union[int, float]]], Dict[str, Residue], bool]]:
        """
         Returns a dict with calculated statistics

         Returns:
            Dict as {}
        """
        return {
            'nmodels': self.nmodels,
            'models_type': self.models_type,
            'nchains': len(self.chain_ids),
            'chain_ids': self.chain_ids,
            'chain_details' : self.chain_details,
            'num_res': self.num_res,
            'num_ats': self.num_ats,
            'res_insc': self.res_insc,
            'res_hetats': self.res_hetats,
            'res_ligands': self.res_ligands,
            'num_wat': self.num_wat,
            'ca_only': self.ca_only,
            'biounit': self.biounit,
            'total_charge': self.total_charge
        }

    def get_term_res(self) -> List[Tuple[str, Residue]]:
        """ Get terminal residues """
        term_res = []
        for res in self.all_residues:
            if mu.is_hetatm(res):
                continue
            #TODO NA
            if res.get_resname() not in self.data_library.get_valid_codes('protein'):
                continue
            if self.is_N_term(res):
                term_res.append(('N', res))
            if self.is_C_term(res):
                term_res.append(('C', res))
        return term_res

    def print_headers(self) -> None:
        """
        Prints selected components from structure headers
        """
        self.get_headers()
        if 'entry_id' in self.meta:
            print(' PDB id: {}'.format(self.meta['entry_id']))
        print(' Title: {}'.format(self.meta['title']))
        print(' Experimental method: {}'.format(self.meta['method']))
        if 'keywords' in self.meta:
            print(' Keywords: {}'.format(self.meta['keywords']))
        if 'resolution' in self.meta:
            print(' Resolution (A): {}'.format(self.meta['resolution']))
        if self.biounit:
            print(' Biounit no. {}'. format(self.meta['biounit']))

    def get_headers(self) -> None:
        """
        Extract selected components from structure headers
        """
        self.meta = {}
        if self.input_format == 'cif':
            self.meta['entry_id'] = ', '.join(self.headers['_entry.id'])
            self.meta['title'] = ', '.join(self.headers['_struct.title'])
            self.meta['method'] = ', '.join(self.headers['_exptl.method'])
            self.meta['keywords'] = ', '.join(self.headers['_struct_keywords.pdbx_keywords'])
            if '_refine_hist.d_res_high' in self.headers:
                self.meta['resolution'] = ', '.join(self.headers['_refine_hist.d_res_high'])
        else:
            self.meta['title'] = self.headers['name']
            self.meta['method'] = self.headers['structure_method']
            if 'keywords' in self.headers:
                self.meta['keywords'] = self.headers['keywords']
            if 'resolution' not in self.headers or not self.headers['resolution']:
                self.meta['resolution'] = 'N.A.'
            else:
                self.meta['resolution'] = self.headers['resolution']
        if self.biounit:
            self.meta['biounit'] = self.biounit

    def print_model_stats(self, prefix='') -> None:
        """ Print stats """
        if self.nmodels > 1:
            print(
                '{} Num. models: {} (type: {}, {:8.3f} A)'.format(
                    prefix,
                    self.nmodels,
                    mu.MODEL_TYPE_LABELS[self.models_type['type']],
                    self.models_type['rmsd']
                )
            )
        else:
            print('{} Num. models: {}'.format(prefix, self.nmodels))

    def print_chain_stats(self, prefix='') -> None:
        """ Print chains info """
        chids = []
        for ch_id in sorted(self.chain_ids):
            if ch_id in self.chain_details:
                chids.append('{}: Unknown (P:{g[0]:.1%} DNA:{g[1]:.1%} RNA:{g[2]:.1%} UNK:{g[3]:.1%})'.format(
                    ch_id, g=self.chain_details[ch_id]))
            else:
                chids.append(
                    '{}: {}'.format(
                        ch_id, mu.CHAIN_TYPE_LABELS[self.chain_ids[ch_id]]
                    )
                )
        print('{} Num. chains: {} ({})'.format(prefix, len(self.chain_ids), ', '.join(chids)))

    def print_stats(self, prefix='') -> None:
        """
        Prints statistics to stdout

        Args:
            prefix: Text prefix to prepend to printed data
        """
        stats = self.get_stats()
        self.print_model_stats(prefix)
        self.print_chain_stats(prefix)

        print('{} Num. residues:  {}'.format(prefix, stats['num_res']))
        print('{} Num. residues with ins. codes:  {}'.format(prefix, stats['res_insc']))
        print('{} Num. HETATM residues:  {}'.format(prefix, stats['res_hetats']))
        print('{} Num. ligands or modified residues:  {}'.format(prefix, stats['res_ligands']))
        print('{} Num. water mol.:  {}'.format(prefix, stats['num_wat']))
        print('{} Num. atoms:  {}'.format(prefix, stats['num_ats']))
        if stats['ca_only']:
            print('Possible CA-Only structure')
        if self.hetatm[mu.MODRES]:
            print('Modified residues found')
            for res in self.hetatm[mu.MODRES]:
                print(mu.residue_id(res))
        if self.hetatm[mu.METAL]:
            print('Metal/Ion residues found')
            for res in self.hetatm[mu.METAL]:
                print(mu.residue_id(res))
        if self.hetatm[mu.ORGANIC]:
            print('Small mol ligands found')
            for res in self.hetatm[mu.ORGANIC]:
                print(mu.residue_id(res))
        if self.has_charges:
            print("Total charge {:6.3f}".format(self.total_charge))

    def save_structure(
            self, 
            output_pdb_path: str, 
            mod_id: str = None, 
            rename_terms: bool = False, 
            output_format='pdb'):
        """
        Saves structure on disk in PDB format

        Args:
            output_pdb_path: OS path to the output file
            mod_id (optional): model to write

        Errors:
            OSError: Error saving the file
        """
        if not output_pdb_path:
            raise OutputPathNotProvidedError
        pdbio = PDBIO_extended(is_pqr=self.has_charges, output_format=output_format)

        if rename_terms:
            self.rename_terms(self.get_term_res())
        else:
            self.revert_terms()

        if mod_id is None:
            pdbio.set_structure(self.st)
            pdbio.save(output_pdb_path)
        else:
            pdbio.set_structure(self.st[mod_id])
            pdbio.save(output_pdb_path)

    def get_all_r2r_distances(self, res_group: Union[str, Iterable[str]], join_models: bool) -> List[Tuple[Residue, Residue, float]]:
        """ Determine residue pairs within a given Cutoff distance
            calculated from the first atom available
            Args:
                res_group: list of residues to check | 'all'
                join_models: consider all models as separated molecules
            Output:
                List of tupes (r1,r2,dist)

        """
        if self.has_NA():
            cutoff = self.data_library.distances['R_R_CUTOFF']['NA']
        else:
            cutoff = self.data_library.distances['R_R_CUTOFF']['PROT']
        return mu.get_all_r2r_distances(
            self.st,
            res_group,
            cutoff,
            join_models=join_models
        )

    def get_altloc_residues(self) -> Dict[Residue, List[Atom]]:
        return mu.get_altloc_residues(self.st)

# Methods to modify structure
    def select_model(self, keep_model: str) -> None:
        """ Selects model(s) and delete the others from the structure. Model are
            renumbered

            Args:
                keep_model: Model number(s) to keep
        """
        models = []
        if '-' in keep_model:
            m1, m2 = keep_model.split('-')
            for v in range(int(m1), int(m2) + 1):
                models.append(v)
        elif ',' in keep_model:
            models = [int(m) for m in keep_model.split(',')]
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
        
        # Update internal data
        self.update_internals()
        self.modified = True

    def superimpose_models(self):
        spimp = Superimposer()
        if self.nmodels > 1:
            fix_atoms = [at for at in self.st[0].get_atoms() if at.id == 'CA']
            for mod in self.st.get_models():
                if mod.id == 0:
                    continue
                mov_atoms = [at for at in self.st[mod.id].get_atoms() if at.id == 'CA']
                spimp.set_atoms(fix_atoms, mov_atoms)
                spimp.apply(self.st[mod.id].get_atoms())
            self.models_type = mu.guess_models_type(self.st) if self.nmodels > 1 else 0
            self.modified = True
    def has_models(self) -> bool:
        """ Shotcut method to check whether the structure has more than one model

            Returns: Boolean
        """
        return self.nmodels > 1


    def has_superimp_models(self) -> bool:
        """ Shotcut method to check whether the structure has superimposed
            models (i.e. NMR or ensemble)

            Returns: Boolean
        """
        return self.models_type and self.models_type['type'] == mu.ENSM

    def set_chain_ids(self) -> None:
        """
        Identifies and sets the chain ids, guessing its nature (protein, dna, rna, ...)
        """
        self.chain_ids = {}
        for chn in self.st.get_chains():
            if not self.biounit and chn.get_parent().id > 0:
                continue
            guess = mu.guess_chain_type(chn)
            if isinstance(guess, list):
                self.chain_ids[chn.id] = mu.UNKNOWN
                if chn.id not in self.chain_details:
                    self.chain_details[chn.id] = guess
            else:
                self.chain_ids[chn.id] = guess

    def has_NA(self):
        """ Checks if any of the chains is NA"""
        has_NA = False
        for k, v in self.chain_ids.items():
            has_NA = (has_NA or (v > 1))
        return has_NA

    def select_chains(self, select_chains: str) -> None:
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids, | protein | dna | rna | na
        """
        if not self.chain_ids:
            self.set_chain_ids()

        if select_chains.lower() in ('protein', 'dna', 'rna', 'na'):
            if select_chains.lower() == 'na':
                ch_ok = [mu.DNA, mu.RNA]
            else:
                ch_ok = [mu.TYPE_LABEL[select_chains.lower()]]
        else:
            ch_ok = select_chains.split(',')
            for chn in ch_ok:
                if chn not in self.chain_ids:
                    print('Warning: skipping unknown chain', chn)
        for mod in self.st:
            for chn in self.chain_ids:
                if chn not in ch_ok and self.chain_ids[chn] not in ch_ok:
                    self.st[mod.id].detach_child(chn)
        if not self.st[mod.id]:
            print("ERROR: would remove all chains, exiting")
            sys.exit()
        # Update internal data
        self.update_internals()
        self.modified = True

    def select_altloc_residue(self, res: Residue, to_fix: Mapping[str, Union[str, Atom]]) -> None:
        """ Selects one alternative conformation when altloc exists. Selection is
            done on the occupancy basis or from the conformation id.
            All relevant atoms in the residue **res** at modified in the same way.
            Triggers **modified** flag.

            Args:
                res: residue (as Bio.PDB.Residue)
                to_fix: atoms in residue to fix as
                    {'ats':[atom_list],'select':'occupancy|conf_id'
        """
        for atm in to_fix['ats']:
            if to_fix['select'].lower() == 'occupancy':
                newat = atm.selected_child
            else:
                if to_fix['select'] in atm.child_dict.keys():
                    newat = atm.child_dict[to_fix['select']]
                else:
                    print(
                        'Warning: unknown alternative {} in {}'.format(
                            to_fix['select'], mu.atom_id(atm)
                        )
                    )
                    continue
            newat.disordered_flag = 0
            newat.altloc = ' '
            res.detach_child(atm.id)
            res.add(newat)
        res.disordered = 0
        self.modified = True

    def remove_residue(self, res: Residue, update_int: bool = True) -> None:
        """
        Removes residue **r** from the structure. Triggers **modified** flag
        and atom and residue renumbering
        """
        mu.remove_residue(res)

        # Update internal data
        if update_int:
            self.update_internals()

        self.modified = True

    def fix_side_chain(self, r_at: Tuple[Residue, Iterable[str]]) -> None:
        """
        Fix missing side chain atoms in given residue. Triggers **modified** flag

        Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
        """
        print(mu.residue_id(r_at[0]))
        for at_id in r_at[1]:
            print("  Adding new atom " + at_id)
            if at_id == 'CB':
                coords = mu.build_coords_CB(r_at[0])
            else:
                coords = mu.build_coords_from_lib(
                    r_at[0],
                    self.res_library,
                    r_at[0].get_resname(),
                    at_id
                )
            mu.add_new_atom_to_residue(r_at[0], at_id, coords)
        self.atom_renumbering()
        self.modified = True

    def rebuild_side_chains(self, r_list: Iterable[str]) -> None:
        """ Rebuild side chain as mutation to same residue using Modeller """
        mut_list = [
            r_at[0].get_parent().id + ':' +\
            r_at[0].get_resname()+\
            str(r_at[0].id[1])+\
            r_at[0].get_resname()
            for r_at in r_list
        ]

        self.rebuild_mutations(self.prepare_mutations(','.join(mut_list)))

        self.atom_renumbering()
        self.modified = True


    def fix_backbone_chain(
            self,
            brk_list: Iterable[Atom],
            modeller_key: str = '',
            extra_gap: int = 0
        ) -> str:
        """ Fixes backbone breaks using Modeller """
        ch_to_fix = set()
        for brk in brk_list:
            ch_to_fix.add(brk[0].get_parent().id)

        modeller_result = self.run_modeller(ch_to_fix, brk_list, modeller_key, extra_gap, extra_NTerm=0)

        self.update_internals()

        return modeller_result

    def run_modeller(
            self,
            ch_to_fix,
            brk_list,
            modeller_key='',
            extra_gap: int = 0,
            extra_NTerm: int = 0,
            sequence_data=None
        ):
        """ Runs modeller """
        # environ var depends on MODELLER version!!! TODO Check usage of this feature by later Modeller versions
        if modeller_key:
            os.environ['KEY_MODELLER9v25'] = modeller_key
        from biobb_structure_checking.modeller_manager import ModellerManager, NoCanSeqError

        mod_mgr = ModellerManager()
        if not sequence_data:
            sequence_data = self.sequence_data

        mod_mgr.sequences = sequence_data

        modif_residues = []

        for mod in self.st:
            if self.has_models():
                print('Processing Model {}'.format(mod.id + 1))
                self.save_structure('{}/templ.pdb'.format(mod_mgr.tmpdir), mod.id)
            else:
                self.save_structure('{}/templ.pdb'.format(mod_mgr.tmpdir))

            for ch_id in self.chain_ids:
                if ch_id not in ch_to_fix:
                    continue
                if sequence_data.data[ch_id]['pdb'][mod.id]['wrong_order']:
                    print("Warning: chain {} has a unusual residue numbering, skipping".format(ch_id))
                print("Fixing chain/model {}/{}".format(ch_id, mod.id))

                try:
                    model_pdb = mod_mgr.build(mod.id, ch_id, extra_NTerm)
                except NoCanSeqError as e:
                    print(e.message)
                    continue

                parser = PDBParser(PERMISSIVE=1)
                model_st = parser.get_structure(
                    'model_st',
                    mod_mgr.tmpdir + "/" + model_pdb['name']
                )

                modif_set_residues = self.merge_structure(
                    sequence_data,
                    model_st,
                    mod.id,
                    ch_id,
                    brk_list,
                    sequence_data.data[ch_id]['pdb'][mod.id]['frgs'][0].features[0].location.start,
                    extra_gap
                ) #TODO consider use canonical numbering instead of defining offset
                modif_residues += modif_set_residues

        return modif_residues

    def merge_structure(
            self,
            sequence_data,
            new_st: Structure,
            mod_id: str,
            ch_id: str,
            brk_list: Iterable[Atom],
            offset: int,
            extra_gap: int = 0
        ) -> str:
        """ Merges the required fragments of Modeller results"""

        spimp = Superimposer()

        list_res = self.st[mod_id][ch_id].get_list()

        modif_residues = []

        for i in range(0, len(sequence_data.data[ch_id]['pdb'][mod_id]['frgs']) - 1):
            loc_i = sequence_data.data[ch_id]['pdb'][mod_id]['frgs'][i].features[0].location
            loc_ii = sequence_data.data[ch_id]['pdb'][mod_id]['frgs'][i + 1].features[0].location
            seq_i = sequence_data.data[ch_id]['pdb'][mod_id]['frgs'][i].features[2].location
            seq_ii = sequence_data.data[ch_id]['pdb'][mod_id]['frgs'][i + 1].features[2].location

            gap_start = loc_i.end
            gap_end = loc_ii.start
            # Gap length taken from sequence gap to avoid PDB numbering issues
            gap_length = seq_ii.start - seq_i.end - 1
            # Offset to account for breaks in PDB residue numbering
            seq_off_i_ii = gap_end - seq_ii.start - gap_start + seq_i.end

            if [self.st[mod_id][ch_id][gap_start], self.st[mod_id][ch_id][gap_end]] not in brk_list:
                #Checking for incomplete gap build needed for fixing side chains with rebuild
                n_br = 0
                while n_br < len(brk_list) - 1 and \
                        (self.st[mod_id][ch_id][gap_start] != brk_list[n_br][0]) and\
                        (self.st[mod_id][ch_id][gap_end] != brk_list[n_br][1]):
                    n_br += 1
                if self.st[mod_id][ch_id][gap_start] == brk_list[n_br][0] or\
                        self.st[mod_id][ch_id][gap_end] == brk_list[n_br][1]:
                    # adjusting gap ends to existint residues
                    gap_start = brk_list[n_br][0].id[1]
                    while gap_start not in self.st[mod_id][ch_id]:
                        gap_start += 1

                    gap_end = brk_list[n_br][1].id[1]
                    while gap_end not in self.st[mod_id][ch_id]:
                        gap_end -= 1

                    extra_gap = 0
                else:
                    continue

            print('Fixing {} - {}'.format(
                mu.residue_id(self.st[mod_id][ch_id][gap_start]),
                mu.residue_id(self.st[mod_id][ch_id][gap_end])))

            # Superimposes structures using fragments at both sides of the gap
            fixed_ats = []
            moving_ats = []

            #checking whether there is a chain id in the model (Support for Modeller >= 10)
            new_ch_id = new_st[0].child_list[0].id

            for nres in range(loc_i.start, loc_i.end):
                mod_nres = nres - offset + 1
                if nres in self.st[mod_id][ch_id] and mod_nres in new_st[0][new_ch_id]:
                    fixed_ats.append(self.st[mod_id][ch_id][nres]['CA'])
                    moving_ats.append(new_st[0][new_ch_id][mod_nres]['CA'])

            for nres in range(loc_ii.start, loc_ii.end):
                mod_nres = nres - offset + 1 - seq_off_i_ii
                if nres in self.st[mod_id][ch_id] and\
                        'CA' in self.st[mod_id][ch_id][nres] and\
                        mod_nres in new_st[0][new_ch_id]:
                    fixed_ats.append(self.st[mod_id][ch_id][nres]['CA'])
                    moving_ats.append(new_st[0][new_ch_id][mod_nres]['CA'])

            moving_ats = moving_ats[:len(fixed_ats)]
            spimp.set_atoms(fixed_ats, moving_ats)
            spimp.apply(new_st.get_atoms())

            # Find position if the 1st residue in the internal residue list
            pos = 0
            while pos < len(list_res) and self.st[mod_id][ch_id].child_list[pos].id[1] != gap_start - extra_gap:
                pos += 1

            res_pairs = []
            for nres in range(gap_start - extra_gap, gap_start + gap_length + 1):
                res_pairs.append([nres, nres - offset + 1])
            for nres in range(gap_end, gap_end + extra_gap + 1):
                res_pairs.append([nres, nres - offset + 1 - seq_off_i_ii])

            for res_pair in res_pairs:
                nres, mod_nres = res_pair
                if nres in self.st[mod_id][ch_id]:
                    self.remove_residue(self.st[mod_id][ch_id][nres], update_int=False)

                res = new_st[0][new_ch_id][mod_nres].copy()
                res.id = (' ', nres, ' ')
                self.st[mod_id][ch_id].insert(pos, res)
                pos += 1
                if nres < gap_start or nres > gap_end:
                    print("  Replacing " + mu.residue_id(res))
                else:
                    print("  Adding " + mu.residue_id(res))

                modif_residues.append(self.st[mod_id][ch_id][nres])

            print()

        return modif_residues

    def add_main_chain_caps(self, caps_list: Iterable[Iterable[str]]) -> List[str]:
        """ Adds ACE and NME caps """
        # print(caps_list)
        fixed = []
        for cap in caps_list:
            if cap[0] == 'N':
                if cap[1] in self.next_residue:
                    mu.add_ACE_cap_at_res(cap[1], self.next_residue[cap[1]])
                else:
                    raise NotEnoughAtomsError
            else:
                if cap[1] in self.prev_residue:
                    mu.add_NME_cap_at_res(cap[1], self.prev_residue[cap[1]])
                else:
                    raise NotEnoughAtomsError
            fixed.append(cap[1])

        self.update_internals()
        return fixed

    def fix_backbone_O_atoms(self, r_at: Tuple[Residue, Iterable[str]]) -> bool:
        """Adding missing backbone atoms not affecting main-chain like O and OXT
                Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
        """
        res, at_list = r_at
        print(mu.residue_id(res))
        if 'C' not in res:
            raise NotEnoughAtomsError
        if len(at_list) == 2 or at_list == ['O']:
            if 'CA' not in res or 'N' not in res or 'C' not in res:
                raise NotEnoughAtomsError
            print("  Adding new atom O")
            mu.add_new_atom_to_residue(res, 'O', mu.build_coords_O(res))
        if 'OXT' in at_list:
            if 'CA' not in res or 'C' not in res or 'O' not in res:
                raise NotEnoughAtomsError
            print("  Adding new atom OXT")
            mu.add_new_atom_to_residue(
                res,
                'OXT',
                mu.build_coords_SP2(mu.OINTERNALS[0], res['C'], res['CA'], res['O'])
            )
            res.resname = 'C' + res.resname

        self.atom_renumbering()
        self.modified = True
        return True

    def add_hydrogens(self, ion_res_list, remove_h: bool = True, add_charges: str = 'ADT'):
        """
        Add hydrogens considering selections in ion_res_list

        Args:
           **r_at_list**: dict as Bio.PDB.Residue: Tauromeric Option
           **remove_h**: Remove Hydrogen atom before adding new ones
           **add_charges**: Add charges and atom types acording to ff
        """
        add_h_rules = self.data_library.get_add_h_rules()

        for res in self.all_residues:
            if mu.is_hetatm(res):
                continue

            if remove_h:
                mu.remove_H_from_r(res, verbose=False)

            if res not in self.prev_residue:
                prev_residue = None
            else:
                prev_residue = self.prev_residue[res]

            if res not in self.next_residue:
                next_residue = None
            else:
                next_residue = self.next_residue[res]

            error_msg = mu.add_hydrogens_backbone(res, prev_residue, next_residue)

            if error_msg:
                print(error_msg, mu.residue_id(res))

            rcode = res.get_resname()

            if len(rcode) == 4:
                rcode = rcode[1:]

            if rcode == 'GLY':
                continue
            # Fixed for modified residues already in the original PDB
        
            if rcode in self.data_library.canonical_codes:
                rcode_can = self.data_library.canonical_codes[rcode]
            else:
                rcode_can = rcode

            
            if rcode_can not in add_h_rules:              
                print(NotAValidResidueError(rcode).message)
                continue
            
            if rcode_can == rcode:
                h_rules = add_h_rules[rcode]
            else:
                h_rules = add_h_rules[rcode_can][rcode]

            if res in ion_res_list:
                if rcode != ion_res_list[res]:
                    print(
                        'Replacing {} by {}'.format(
                            mu.residue_id(res), ion_res_list[res]
                        )
                    )
                error_msg = mu.add_hydrogens_side(
                    res,
                    self.res_library,
                    ion_res_list[res],
                    h_rules[ion_res_list[res]]
                )
                res.resname = ion_res_list[res]
            else:
                error_msg = mu.add_hydrogens_side(res, self.res_library, rcode, h_rules)

            if error_msg:
                print(error_msg, mu.residue_id(res))

        self.residue_renumbering()
        if add_charges:
            self.rename_terms(self.get_term_res())
            self.update_atom_charges(add_charges)
        self.atom_renumbering()
        self.modified = True

    def mark_ssbonds(self, cys_list):
        """ Mark Cys residue in cys_list as CYX for further usage """
        for res in cys_list:
            if 'HG' in res:
                mu.remove_atom_from_res(res, 'HG')
            res.resname = 'CYX'
        self.modified = True

    def rename_terms(self, term_res):
        """ Rename Terminal residues as NXXX or CXXX """
        for t in term_res:
            if t[1].resname not in ('ACE', 'NME'):
                t[1].resname = t[0] + t[1].resname

    def revert_terms(self):
        """ Reverts 4 char len residue names to 3 letter codes"""
        for res in self.st.get_residues():
            if len(res.get_resname()) == 4:
                res.resname = res.resname[1:]

    def is_N_term(self, res: Residue) -> bool:
        """ Detects whether it is N terminal residue."""
        return res not in self.prev_residue

    def is_C_term(self, res: Residue) -> bool:
        """ Detects whether it is C terminal residue."""
        return res not in self.next_residue

    def prepare_mutations(self, mut_list: str) -> List[MutationSet]:
        """ Find residues to mutate from mut_list"""
        mutations = MutationManager(mut_list, self.chain_ids)
        mutations.prepare_mutations(self.st)
        return mutations
    
    def apply_mutations(self, mutations: MutationManager) -> Residue:
        """ Perform mutations """
        mutated_res = mutations.apply_mutations(
            self.data_library.get_mutation_map(),
            self.res_library
        )
        self.residue_renumbering()
        self.atom_renumbering()
        self.modified = True
        return mutated_res

    def rebuild_mutations(self, mutations: MutationManager, modeller_key: str = '') -> Residue:
        """ Perform mutations Rebuilding side chain"""
        ch_to_fix = set()
        brk_list = []
        for mut_set in mutations.mutation_list:
            for mut in mut_set.mutations:
                if self.chain_ids[mut['chain']] > 1:
                    continue
                ch_to_fix.add(mut['chain'])
                start_res = mut['resobj']
                if start_res in self.prev_residue:
                    start_res = self.prev_residue[start_res]
                end_res = mut['resobj']
                if end_res in self.next_residue:
                    end_res = self.next_residue[end_res]
                brk = [start_res, end_res]
                brk_list.append(brk)
        if not ch_to_fix:
            print("No proteins chains left, exiting")
            return []
        mutated_sequence_data = SequenceData()
        mutated_sequence_data.fake_canonical_sequence(self, mutations)
        for mut_set in mutations.mutation_list:
            for mut in mut_set.mutations:
                mu.remove_residue(mut['resobj'])
        mutated_sequence_data.read_structure_seqs(self)
        mutated_sequence_data.match_sequence_numbering()

        #TODO Not tested, to be used on changes in the NTerm residue
        extra_NTerm = 0
        mutated_res = self.run_modeller(ch_to_fix, brk_list, modeller_key, 0, extra_NTerm, mutated_sequence_data)

        self.update_internals()
#        mutated_res = []
#        for frag in (modeller_result):
#            m = re.search('(.*)-(.*)/(.*)', frag)
#            mod, ch, start_res, end_res  = \
#                int(m.group(3)) - 1, \
#                m.group(1)[:1], \
#                int(m.group(1)[1:]), \
#                int(m.group(2)[1:])
#            for i in range(start_res, end_res + 1):
#                mutated_res.append(self.st[mod][ch][i])
        return mutated_res

    def invert_amide_atoms(self, res: Residue):
        """ Fix sidechains with incorrect amide assignments"""
        amide_res = self.data_library.get_amide_data()[0]
        res_type = res.get_resname()
        if res_type not in amide_res:
            raise NotAValidResidueError(res_type)
        mu.swap_atoms(
            res[amide_res[res_type][0]],
            res[amide_res[res_type][1]]
        )

    def fix_chiral_chains(self, res: Residue):
        """ Fix sidechains with chiral errors"""
        chiral_res = self.data_library.get_chiral_data()
        res_type = res.get_resname()
        mu.swap_atoms(
            res[chiral_res[res_type][0]],
            res[chiral_res[res_type][1]]
        )
        if res_type == 'ILE':
            mu.delete_atom(res, 'CD1')
            mu.build_atom(res, 'CD1', self.res_library, 'ILE')

    def prepare_mutations_from_na_seq(self, new_seq):
        """ Prepare mutation list to get new_seq, including complementary chain
            Only for canonical Duplexes

        """
        if new_seq.find(':') != -1:
            mut_seq = new_seq.split(':')
        else:
            mut_seq = [new_seq, mu.rev_complement_na_seq(new_seq)]
        #Prepared for std Duplexes
        mut_list = []
        i = 0
        nch = 0
        for ch_id in self.sequence_data.data:
            chn = self.sequence_data.data[ch_id]
            start = chn['pdb'][0]['frgs'][0].features[0].location.start
            seq = chn['pdb'][0]['frgs'][0].seq
            if len(seq) != len(mut_seq[nch]):
                sys.exit("Sequence lengths do not match")
            prefix = ''
            if chn['pdb'][0]['type'] == mu.DNA:
                prefix = 'D'
            for i, r in enumerate(seq):
                mut_list.append('{}:{}{}{}{}{}'.format(ch_id, prefix, r, start + i, prefix, mut_seq[nch][i]))
            nch += 1
        return ','.join(mut_list)


# ===============================================================================


class Error(Exception):
    """ Base class """
    pass


class WrongServerError(Error):
    def __init__(self):
        self.message = 'ERROR: Biounits supported only on MMB server'


class UnknownFileTypeError(Error):
    def __init__(self, typ):
        self.message = 'ERROR: unknown filetype ({})'.format(typ)


class OutputPathNotProvidedError(Error):
    def __init__(self):
        self.message = 'ERROR: output PDB path not provided'


class NotAValidResidueError(Error):
    def __init__(self, res):
        self.message = 'Warning: {} is not a valid residue in this context'.format(res)


class NotEnoughAtomsError(Error):
    def __init__(self):
        self.message = 'Warning: not enough backbone to build missing atoms'


class ParseError(Error):
    def __init__(self, err_id, err_txt):
        self.message = '{} ({}) found when parsing input structure'.format(err_id, err_txt)

class UnknownFFError(Error):
    def __init__(self, ff):
        self.message = '{} is not a valid ff for assigning atom types'.format(ff)