"""Module to manage structure, based on BioPython Bio.PDB
"""
# from cgitb import reset
# from re import I
import warnings
import os
from os.path import join as opj
import shutil
import sys
import re

from urllib.request import urlretrieve
from urllib.request import urlcleanup

from typing import List, Dict, Tuple, Iterable, Mapping, Union, Set

from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionException

import biobb_structure_checking.constants as cts

from biobb_structure_checking.io.mmb_server import MMBPDBList
from biobb_structure_checking.io.PDBIO_extended import PDBIO_extended
from biobb_structure_checking.io.bare_builder import BareStructureBuilder

from biobb_structure_checking.libs.data_lib_manager import DataLibManager
from biobb_structure_checking.libs.residue_lib_manager import ResidueLib

from biobb_structure_checking.modelling.sequencedata import SequenceData
from biobb_structure_checking.modelling.modelsdata import ModelsData
from biobb_structure_checking.modelling.chainsdata import ChainsData
from biobb_structure_checking.modelling.structuredata import StructureData

import biobb_structure_checking.modelling.utils as mu

from biobb_structure_checking.mutation_manager import MutationManager, MutationSet

MODELLER_ENV_VAR = 'KEY_MODELLER10v4'
ACCEPTED_FORMATS = ['mmCif', 'cif', 'pdb', 'pqr', 'pdbqt']
ACCEPTED_REMOTE_FORMATS = ['mmCif', 'cif', 'pdb', 'xml']


class StructureManager:
    """ Main Class wrapping Bio.PDB structure object
    """
    def __init__(
            self,
            input_pdb_path: str,
            data_library_path: str,
            res_library_path: str,
            pdb_server: str,
            cache_dir: str = 'tmpPDB',
            nocache: bool = False,
            copy_dir: str = './',
            file_format: str = 'mmCif',
            fasta_sequence_path: str = '',
            nowarn: bool = True,
            coords_only: bool = False,
            overwrite: bool = False,
            atom_limit: int = 0
    ) -> None:
        """
            Class constructor. Sets an empty object and loads a structure
            according to parameters

        Args:
            **input_pdb_path** (str): path to input structure either in pdb or mmCIF format.
                Format is taken from file extension. Alternatively **pdb:pdbId** fetches the mmCIF file from RCSB
            **data_library_path** (str): Path to json data library
            **res_library_path** (str): Path to residue library
            **pdb_server** (str): **default** for Bio.PDB defaults (RCSB), **mmb** for MMB PDB API
            **cache_dir** (str): path to temporary dir to store downloaded structures
            **nocache** (bool): Do not cache downloaded structures
            **copy_dir** (str=: Folder to copy input structure
            **file_format** (str): structure file format to use
            **fasta_sequence_path** (str): path to canonical sequence file (needed for PDB input)
            **nowarn** (bool): No warnings on Structure Building

        """
        self.data_library = DataLibManager(data_library_path)
        for ff_name in self.data_library.ff_data:
            self.data_library.get_ff_data(
                opj(os.path.dirname(data_library_path), f"{ff_name.upper()}_prm.json")
            )

        self.res_library = ResidueLib(res_library_path)

        self.sequence_data = SequenceData()

        if fasta_sequence_path:
            self.sequence_data.load_sequence_from_fasta(fasta_sequence_path)

        self.st, headers, input_format, biounit = self._load_structure_file(
            input_pdb_path,
            cache_dir,
            nocache,
            copy_dir,
            pdb_server,
            file_format,
            QUIET=nowarn,
            coords_only=coords_only,
            overwrite=overwrite
        )
        # Check limit of atoms now to avoid delays
        num_ats = len(list(self.st.get_atoms()))
        if atom_limit and num_ats > atom_limit:
            sys.exit(
                cts.MSGS['ATOM_LIMIT'].format(
                    num_ats,
                    atom_limit
                )
            )

        self.models_data = ModelsData(self.st)
        self.chains_data = ChainsData(self.st)
        self.st_data = StructureData(self.st, input_format, headers, biounit)
        self.modified = False

        # Calc internal data
        self.update_internals(cif_warn=True)

    def _load_structure_file(
            self,
            input_pdb_path,
            cache_dir,
            nocache,
            copy_dir,
            pdb_server,
            file_format,
            QUIET=False,
            coords_only=False,
            overwrite=False

    ):
        """ Load structure file """
        biounit = False
        self.pdb_id = 'User'
        if input_pdb_path.startswith('pdb:'):
            input_pdb_path = input_pdb_path[4:]
            # MMBPDBList child defaults to Bio.PDB.PDBList
            # if MMB/BSC server is not selected
            pdbl = MMBPDBList(pdb=cache_dir, server=pdb_server)
            if re.search(r'\.[1-9]+$', input_pdb_path):
                pdbid, biounit = input_pdb_path.split('.')
                input_pdb_path = pdbid.upper()

                # if pdb_server not in ALT_SERVERS:
                #    raise WrongServerError
                if not biounit:
                    real_pdb_path = pdbl.retrieve_pdb_file(
                        input_pdb_path,
                        file_format='pdb',
                        biounit=biounit,
                        nocache=nocache,
                        overwrite=overwrite
                    )
                else:
                    real_pdb_path = pdbl.retrieve_assembly_file(
                        input_pdb_path,
                        biounit,
                        nocache=nocache,
                        overwrite=overwrite
                    )
                self.pdb_id = pdbid
            else:
                if '.' in input_pdb_path:
                    pdbid, file_format = input_pdb_path.split('.')
                    input_pdb_path = pdbid.upper()
                    if file_format not in ACCEPTED_REMOTE_FORMATS:
                        print(
                            f"WARNING: format {file_format} not available "
                            "for downloads, reverting to default"
                        )
                        file_format = 'cif'
                    self.pdb_id = pdbid
                else:
                    input_pdb_path = input_pdb_path.upper()

                # Force mmCif as cif is not accepted by biopython
                if file_format == 'cif':
                    file_format = 'mmCif'
                real_pdb_path = pdbl.retrieve_pdb_file(
                    input_pdb_path,
                    file_format=file_format,
                    nocache=nocache,
                    overwrite=overwrite
                )
                self.pdb_id = input_pdb_path
                if file_format == 'pdb':
                    # change file name to id.pdb
                    new_path = opj(
                        os.path.dirname(real_pdb_path),
                        f"{input_pdb_path.lower()}.pdb"
                    )
                    os.rename(real_pdb_path, new_path)
                    real_pdb_path = new_path
                    # Adding sequence input
                    self.sequence_data.load_sequence_from_fasta(f"pdb:{pdbid}")

        elif input_pdb_path.startswith('http'):
            real_pdb_path = opj(cache_dir, os.path.basename(input_pdb_path))

            if '.' in os.path.basename(input_pdb_path):
                file_format = os.path.splitext(input_pdb_path)[1][1:]
                if file_format not in ACCEPTED_FORMATS:
                    print(f'Error: MMB/BSC Server: File format {file_format} not supported')
                    sys.exit(1)
            print(f"Downloading structure from {input_pdb_path} as {file_format} ...")
            try:
                urlcleanup()
                urlretrieve(input_pdb_path, real_pdb_path)
            except IOError:
                print(f"Download failed")

            nocache = True

        else:
            real_pdb_path = input_pdb_path

        if coords_only:
            builder = BareStructureBuilder()
            print(
                "Loading structure coord_only: Chain and residue ids"
                "will be ignored"
            )
        else:
            builder = None

        if '.pdb' in real_pdb_path:  # accepts .pdbqt
            if '.pdbqt' in real_pdb_path:
                print("Warning: PDBQT file will be loaded as PDB")
            parser = PDBParser(
                PERMISSIVE=1,
                is_pqr=False,
                structure_builder=builder,
                QUIET=QUIET
            )
            input_format = 'pdb'
        elif '.pqr' in real_pdb_path:
            parser = PDBParser(
                PERMISSIVE=1,
                is_pqr=True,
                structure_builder=builder
            )
            input_format = 'pqr'
        elif '.cif' in real_pdb_path:
            parser = MMCIFParser(structure_builder=builder, QUIET=QUIET)
            input_format = 'cif'
        else:
            raise UnknownFileTypeError(input_pdb_path)

        try:
            new_st = parser.get_structure(self.pdb_id, real_pdb_path)
        except ValueError as err:
            raise ParseError('ValueError', err) from err
        except PDBConstructionException as err:
            raise ParseError('PDBBuildError', err) from err

        if input_format in ['pdb', 'pqr']:
            headers = parse_pdb_header(real_pdb_path)
        else:
            headers = MMCIF2Dict(real_pdb_path)

        if copy_dir:
            try:
                shutil.copy(real_pdb_path, copy_dir)
                print(
                    f"Storing a copy of the input structure as "
                    f"{opj(copy_dir, os.path.basename(real_pdb_path))}"
                )
            except Exception:
                print(
                    "WARNING: requested copy will overwrite input file, "
                    "skipping"
                )

        if nocache:
            os.remove(real_pdb_path)

        return new_st, headers, input_format, biounit

    def update_internals(self, cif_warn: bool = False):
        """ Update internal data when structure is modified """
        # Add .index field for correlative, unique numbering of residues
        self.st_data.residue_renumbering(self.data_library)

        # Atom renumbering for mmCIF, PDB uses atom number in file
        self.st_data.atom_renumbering()
        self.chains_data.set_chain_ids()
        self.st_data.calc_stats()
        self.st_data.guess_hetatm()

        self.rr_dist = self.get_all_r2r_distances('all', join_models=False)

        # Precalc backbone
        self.st_data.check_backbone_connect(
            ("N", "C", "P", "O3'"),
            self.data_library.distances['COVLNK']
        )
        # get canonical and structure sequences
        self.sequence_data.read_sequences(self, clean=True, cif_warn=cif_warn)

    def update_atom_charges(self, force_field):
        """ Update atom charges and types from data library """

        print("Updating partial charges and atom types")

        self.st_data.total_charge = 0.

        if force_field not in self.data_library.ff_data:
            raise UnknownFFError(force_field)
        ff_data = self.data_library.ff_data[force_field]

        self.rename_terms(self.get_term_res())
        for res in self.st.get_residues():
            ch_type = self.chains_data.get_chain_type(res)
            ch_type_label = mu.CHAIN_TYPE_LABELS[ch_type].lower()
            rcode = res.get_resname()
            rcode3 = rcode  # Normal residues
            if len(rcode) == 4:  # Protein terms
                rcode3 = rcode[1:]
            elif rcode[-1] in ('3', '5'):  # NA Terms
                rcode3 = rcode[:-1]
            can_rcode3 = self.data_library.get_canonical_resname(rcode3)

            if rcode in self.res_library.residues:
                oxt_ok = rcode[0] != 'C' or len(rcode) != 4
                res_chr = 0.
                for atm in res.get_atoms():
                    atm.pqr_charge = self.res_library.get_atom_def(rcode, atm.id).chrg
                    if atm.id in ff_data['residue_data'][can_rcode3]:
                        atm.xtra['atom_type'] = ff_data['residue_data'][can_rcode3][atm.id]
                    elif atm.id in ff_data['residue_data']['*'][ch_type_label]:
                        atm.xtra['atom_type'] = ff_data['residue_data']['*'][ch_type_label][atm.id]
                    else:
                        atm.xtra['atom_type'] = atm.element
                    atm.radius = ff_data['rvdw'][atm.xtra['atom_type']]
                    res_chr += atm.pqr_charge
                    self.st_data.total_charge += atm.pqr_charge
                    if atm.id == 'OXT':
                        oxt_ok = True
                if not oxt_ok:
                    print(
                        f"Warning: OXT atom missing in {mu.residue_id(res)}. "
                        f"Run backbone --fix_atoms first"
                    )
            else:
                print(
                    f"Warning: {rcode} not found in residue library atom, "
                    "charges set to 0."
                )
                for atm in res.get_atoms():
                    atm.pqr_charge = 0.
                    atm.radius = 0.
                    if atm.id in self.data_library.atom_data['metal_atoms']:
                        atm.xtra['atom_type'] = atm.id.lower().capitalize()
                    else:
                        atm.xtra['atom_type'] = atm.element

        print(f"Total assigned charge: {self.st_data.total_charge:10.2f}")

        self.revert_terms()

        self.st_data.has_charges = True

    def get_ins_codes(self) -> List[Residue]:
        """Makes a list with residues having insertion codes"""
        return [
            res
            for res in self.st.get_residues()
            if mu.has_ins_code(res)
        ]

    def get_metal_atoms(self) -> List[Atom]:
        """ Makes a list of possible metal atoms"""
        return mu.get_metal_atoms(
            self.st,
            self.data_library.atom_data['metal_atoms']
        )

    def get_SS_bonds(self) -> List[Union[Atom, Atom, float]]:
        """ Stores and returns possible SS Bonds by distance"""
        self.st_data.ss_bonds = mu.get_all_at2at_distances(
            self.st,
            'SG',
            self.data_library.distances['SS_DIST'],
            not self.models_data.has_superimp_models()
        )
        return self.st_data.ss_bonds

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
            if self.chains_data.chain_ids[chn.get_parent().id][chn.id] == mu.PROTEIN:
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
            'res_to_fix': [
                res
                for res in chiral_bck_list
                if not mu.check_chiral_ca(res)
            ]
        }

    def check_r_list_clashes(
            self,
            residue_list: Iterable[Residue],
            contact_types: Iterable[str],
            get_all_contacts=False
    ) -> Dict[str, Dict[str, Tuple[Residue, Residue, float]]]:
        """ Checks clashes originated by a list of residues"""
        return mu.check_r_list_clashes(
            residue_list,
            self.rr_dist,
            self.data_library.distances['CLASH_DIST'],
            self.data_library.get_atom_lists(contact_types),
            not self.models_data.has_superimp_models(),
            severe='severe' in contact_types,
            get_all_contacts=get_all_contacts
        )

    def check_missing_atoms(self) -> List[Tuple[Residue, Dict[str, List[str]]]]:
        """ Makes a **list of missing atoms** in the structure

            Returns:
                List of residues with missing atoms, as a tuples
                ["r",{"backbone":[atoms_list],"side":[atoms_list]}]
        """
        valid_codes = {}
        residue_data = {}
        for chain_type in ('protein', 'dna', 'rna'):
            valid_codes[mu.TYPE_LABEL[chain_type]] = self.data_library.get_valid_codes(chain_type)
            if chain_type == 'protein':
                valid_codes[mu.TYPE_LABEL[chain_type]] += list(self.data_library.canonical_codes)
            residue_data[mu.TYPE_LABEL[chain_type]] = self.data_library.get_all_atom_lists(chain_type)
        miss_at_list = []
        for res in self.st.get_residues():
            ch_type = self.chains_data.get_chain_type(res)
            if ch_type not in (mu.PROTEIN, mu.NA, mu.DNA, mu.RNA):
                continue

            if res.get_resname() in valid_codes[ch_type] and \
                    not mu.is_hetatm(res):
                if res.get_resname() in self.data_library.canonical_codes:
                    can_rcode = self.data_library.canonical_codes[res.get_resname()]
                else:
                    can_rcode = res.get_resname()
                miss_at = mu.check_all_at_in_r(
                    res, residue_data[ch_type][can_rcode]
                )
                bck_miss = []
                if ch_type == mu.PROTEIN:
                    if self.is_C_term(res) and \
                            res.get_resname() != 'NME' and \
                            'OXT' not in res:
                        bck_miss.append('OXT')
                else:
                    if not self.is_5_term(res):
                        for at_id in ["P", "OP1", "OP2"]:
                            if at_id not in res:
                                bck_miss.append(at_id)
                if bck_miss:
                    if 'backbone' not in miss_at:
                        miss_at['backbone'] = []
                    miss_at['backbone'] += bck_miss
                if miss_at:
                    miss_at_list.append((res, miss_at))
        return miss_at_list

    def check_extra_atoms(self) -> List[Tuple[Residue, Atom]]:
        """ Makes a **list of extra atoms** in the structure

            Returns:
                List of residues with extra atoms, as a tuples
                ["r",atoms_list]
        """
        valid_codes = {}
        residue_data = {}
        for chain_type in ('protein', 'dna', 'rna'):
            valid_codes[mu.TYPE_LABEL[chain_type]] = self.data_library.get_valid_codes(chain_type)
            if chain_type == 'protein':
                valid_codes[mu.TYPE_LABEL[chain_type]] += list(self.data_library.canonical_codes)
            residue_data[mu.TYPE_LABEL[chain_type]] = self.data_library.get_all_atom_lists(chain_type)
        extra_at_list = []
        for res in self.st.get_residues():
            if mu.is_hetatm(res):
                continue
            if res.get_resname() in self.data_library.canonical_codes:
                can_rcode = self.data_library.canonical_codes[res.get_resname()]
            else:
                can_rcode = res.get_resname()
            ch_type = self.chains_data.get_chain_type(res)
            if ch_type == mu.UNKNOWN:
                continue
            rcode = res.get_resname().replace(' ', '')
            if rcode not in valid_codes[ch_type]:
                print(f"Warning: unknown residue {rcode}")
                continue
            if ch_type == mu.PROTEIN:
                extra_ats = mu.check_unk_at_in_r(
                    res,
                    residue_data[ch_type][can_rcode]
                )
                if extra_ats:
                    extra_at_list.append((res, extra_ats))
            else:
                res_at_list = residue_data[ch_type][can_rcode]
                add_ats = []
                if not self.is_5_term(res):
                    if 'P' not in res_at_list['backbone']:  # fix to avoid add repeated groups
                        add_ats = ["P", "OP1", "OP2"]
                extra_ats = mu.check_unk_at_in_r(
                    res,
                    {
                        'backbone':res_at_list['backbone'] + add_ats,
                        'side': res_at_list['side']
                    }
                )
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
                if 'side' in at_list and not self._missing_bck_atoms(res):
                    miss_ats.append((res, at_list['side']))
            else:
                if at_list['backbone']:
                    miss_ats.append((res, at_list['backbone']))
        return miss_ats

    def _missing_bck_atoms(self, res):
        """ Check whether backbone atoms required
            to build side chains are present
        """
        if self.chains_data.get_chain_type(res) == mu.PROTEIN:
            bck_ats = ("N", "CA", "C")
        else:
            bck_ats = ("C1'", "O4'", "C4'")

        missing = False

        for atm in bck_ats:
            missing = missing or atm not in res

        return missing

    def get_ion_res_list(self) -> List[Tuple[Residue, List[str]]]:
        """returns list of residues with potencial selection on adding H

            Returns:
                List of residues that require selection on adding H
                ["r",[atom_list]]
        """
        ion_res = self.data_library.ion_res
        hydrogen_lists = self.data_library.get_hydrogen_atoms()

        ion_res_list = []
        for res in self.st_data.all_residues:
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
        self.st_data.modified_residue_list = []
        for i in range(0, len(self.st_data.all_residues)-1):
            res1 = self.st_data.all_residues[i]
            res2 = self.st_data.all_residues[i + 1]
            if not mu.same_chain(res1, res2):
                continue
            if mu.is_hetatm(res1) or mu.is_hetatm(res2):
                if res1 in self.st_data.next_residue:
                    if self.st_data.next_residue[res1] == res2:
                        if mu.is_hetatm(res1):
                            self.st_data.modified_residue_list.append(res1)
                else:
                    continue
            # Skip NA
            # TODO include NA Backbone
            chn1 = res1.get_parent()
            if self.chains_data.chain_ids[chn1.get_parent().id][res1.get_parent().id] != mu.PROTEIN:
                continue

            if res1 not in self.st_data.next_residue:
                bck_breaks_list.append([res1, res2])
                if mu.seq_consecutive(res1, res2):
                    dist = 0.
                    if 'N' in res1 and 'C' in res2:
                        dist = res1['N'] - res2['C']
                    else:
                        dist = res1.child_list[0] - res2.child_list[0]
                    not_link_seq_list.append([res1, res2, dist])

            else:
                if res2 != self.st_data.next_residue[res1]:
                    wrong_link_list.append([
                        res1,
                        self.st_data.next_residue[res1],
                        res2
                    ])
        return {
            'bck_breaks_list': bck_breaks_list,
            'wrong_link_list': wrong_link_list,
            'not_link_seq_list': not_link_seq_list
        }

    def check_cis_backbone(self):
        """
        Determines omega dihedrals for two bound residues and classifies them
        as normal trans, low angle trans, and cis
        """
        cis_backbone_list = []
        lowtrans_backbone_list = []
        for lnk in self.st_data.backbone_links:
            [at1, at2] = lnk
            res1 = at1.get_parent()
            res2 = at2.get_parent()
            if 'CA' in res1 and 'C' in res1 and 'CA' in res2 and 'N' in res2:
                dih = mu.calc_bond_dihedral(
                    res1['CA'], res1['C'], res2['N'], res2['CA']
                )
                if abs(dih) < mu.CISTHRES:
                    cis_backbone_list.append((res1, res2, dih))
                elif abs(dih) < mu.TRANSTHRES:
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
            mu.AMIDE_CONTACT_TYPES
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
            'nmodels': self.models_data.nmodels,
            'models_type': self.models_data.models_type,
            'nchains': [
                len(self.chains_data.chain_ids[mod_id])
                for mod_id in self.chains_data.chain_ids
            ],
            'chain_ids': [
                {
                    k: mu.CHAIN_TYPE_LABELS[v]
                    for k, v in self.chains_data.chain_ids[mod_id].items()
                }
                for mod_id in self.chains_data.chain_ids
            ],
            'chain_guess_details': [
                self.chains_data.chain_details[mod_id]
                for mod_id in self.chains_data.chain_details
            ],
            'stats': self.st_data.stats,
            'ca_only': self.st_data.ca_only,
            'biounit': self.st_data.biounit,
            'total_charge': self.st_data.total_charge
        }

    def get_term_res(self) -> List[Tuple[str, Residue]]:
        """ Get terminal residues """
        term_res = []
        for res in self.st_data.all_residues:
            if mu.is_hetatm(res):
                continue
            if self.is_N_term(res):
                term_res.append(('N', res))
            elif self.is_C_term(res):
                term_res.append(('C', res))
            elif self.is_5_term(res):
                term_res.append(('5', res))
            elif self.is_3_term(res):
                term_res.append(('3', res))
        return term_res

    def print_stats(self, prefix='') -> None:
        """
        Prints statistics to stdout

        Args:
            prefix: Text prefix to prepend to printed data
        """

        print(self.models_data.stats(prefix))
        print(
            self.chains_data.stats(
                prefix, use_models=self.models_data.has_models()
            )
        )
        st_stats = self.get_stats()
        print(f"{prefix} Num. residues:  {st_stats['stats']['num_res']}")
        print(
            f"{prefix} Num. residues with ins. codes:  "
            f"{st_stats['stats']['res_insc']}"
        )
        if st_stats['stats']['num_h']:
            print(
                f"{prefix} Num. residues with H atoms: "
                f"{st_stats['stats']['res_h']} "
                f"(total {st_stats['stats']['num_h']} H atoms)"
            )
        else:
            print(
                f"{prefix} Num. residues with H atoms: "
                f"{st_stats['stats']['res_h']}"
            )
        print(
            f"{prefix} Num. HETATM residues:  "
            f"{st_stats['stats']['res_hetats']}"
        )
        print(
            f"{prefix} Num. ligands or modified residues:  "
            f"{st_stats['stats']['res_ligands']}"
        )
        print(
            f"{prefix} Num. water mol.:  "
            f"{st_stats['stats']['num_wat']}"
        )
        print(
            f"{prefix} Num. atoms:  "
            f"{st_stats['stats']['num_ats']}"
        )
        if st_stats['ca_only']:
            print(' CA-ONLY structure')
        self.st_data.print_hetatm_stats()

    def save_structure(
            self,
            output_pdb_path: str,
            mod_id: str = None,
            rename_terms: bool = False,
            output_format: str = 'pdb',
            keep_resnames: bool = False):
        """
        Saves structure on disk in PDB format

        Args:
            output_pdb_path: OS path to the output file
            mod_id (optional): model to write
            rename_terms: rename terminal residues
            output_format: select output format
            keep_resnames: keep canonical residue names
        Errors:
            OSError: Error saving the file
        """
        if not output_pdb_path:
            raise OutputPathNotProvidedError

        chain_labels_chr = True
        for chn in self.st.get_chains():
            chain_labels_chr = chain_labels_chr and len(chn.id) == 1

        if output_format in ('cif', 'mmCif'):
            io = MMCIFIO()
        else:
            if not chain_labels_chr:
                print(
                    "WARNING: current structure cannot be saved as PDB, "
                    "using mmCIF instead"
                )
                io = MMCIFIO()
                output_pdb_path += ".cif"
            else:
                io = PDBIO_extended(
                    is_pqr=self.st_data.has_charges,
                    output_format=output_format
                )

        if rename_terms:
            self.rename_terms(self.get_term_res())
        else:
            self.revert_terms()

        if keep_resnames:
            self.revert_can_resnames(canonical=True)
            print("Warning: reverting residue names to canonical on output")

        if mod_id is None:
            io.set_structure(self.st)
            io.save(output_pdb_path)
        else:
            io.set_structure(self.st[mod_id])
            io.save(output_pdb_path)

        if keep_resnames:
            self.revert_can_resnames(canonical=False)

    def get_all_r2r_distances(
        self,
        res_group: Union[str, Iterable[str]], join_models: bool
    ) -> List[Tuple[Residue, Residue, float]]:
        """ Determine residue pairs within a given Cutoff distance
            calculated from the first atom available
            Args:
                res_group: list of residues to check | 'all'
                join_models: consider all models as separated molecules
            Output:
                List of tupes (r1,r2,dist)

        """
        if self.chains_data.has_NA():
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
        """ Returns list of residues with alternative atom locations"""
        return mu.get_altloc_residues(self.st)

# Methods to modify structure
    def select_model(self, keep_model: str) -> None:
        """ Selects model(s) and delete the others from the structure.
            Model are renumbered

            Args:
                keep_model: Model number(s) to keep
        """
        self.models_data.select(keep_model)
        self.update_internals()
        self.modified = True

    def superimpose_models(self):
        ''' Superimpose models by rmsd'''
        self.modified = self.models_data.superimpose_models()

    def build_complex(self):
        ''' Build a complex from biounit models'''
        if self.models_data.models_type['type'] != mu.BUNIT:
            print(
                f"ERROR: No complex can be built. Models superimose "
                f"RMSd {self.models_data.models_type['rmsd']}"
            )
            return 0
        result = self.models_data.build_complex()
        self.update_internals()
        self.modified = True
        return result

    def select_chains(self, select_chains: str) -> None:
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids, | protein | dna | rna | na
        """
        self.chains_data.select(select_chains)
        self.update_internals()
        self.modified = True

    def rename_empty_chain_label(self, new_label):
        '''Add labels to unlabelled chains'''
        result = self.chains_data.rename_empty_chain_label(new_label)
        self.update_internals()
        self.modified = True
        return result

    def renumber_chain_residues(
        self,
        renum_str,
        rem_inscodes=False,
        verbose=False
    ):
        ''' Allow to relabel chains and residues'''
        result = self.chains_data.renumber(
            renum_str,
            rem_inscodes=rem_inscodes,
            verbose=verbose
        )
        if result:
            self.update_internals()
            self.modified = True
        return result

    def rebuild_chains(self, verbose=False):
        ''' Rebuild chains from coordinates'''
        result = self.chains_data.rebuild(self.st_data.backbone_links)
        if result:
            self.modified = True
        return result

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
                        f"Warning: unknown alternative {to_fix['select']} "
                        f"in {mu.atom_id(atm)}"
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
        if self.chains_data.get_chain_type(r_at[0]) == mu.PROTEIN:
            self._fix_side_chain_protein(r_at)
        else:
            self._fix_side_chain_na(r_at)

    def _fix_side_chain_protein(self, r_at):
        for at_id in r_at[1]:
            print(f"  Adding new atom {at_id}")
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
        self.st_data.atom_renumbering()
        self.modified = True

    def _fix_side_chain_na(self, r_at):
        if mu.is_purine(r_at[0]) and ('N9' in r_at[1] or 'C8' in r_at[1]) or\
            mu.is_pyrimidine(r_at[0]) and ('N1' in r_at[1] or 'C6' in r_at[1]):
            print(
                f"Not enough atoms left on {mu.residue_id(r_at[0])} "
                "to recover base orientation, skipping"
            )
        else:
            for at_id in r_at[1]:
                print(f"  Adding new atom {at_id}")
                coords = mu.build_coords_from_lib(
                    r_at[0],
                    self.res_library,
                    r_at[0].get_resname(),
                    at_id
                )
                mu.add_new_atom_to_residue(r_at[0], at_id, coords)
            self.st_data.atom_renumbering()
            self.modified = True

    def rebuild_side_chains(self, r_list: Iterable[str]) -> None:
        """ Rebuild side chain as mutation to same residue using Modeller """
        mut_list = [
            f"{r_at[0].get_parent().id}:"
            f"{r_at[0].get_resname()}{r_at[0].id[1]}"
            f"{r_at[0].get_resname()}"
            for r_at in r_list
        ]

        self.rebuild_mutations(self.prepare_mutations(','.join(mut_list)))

        self.st_data.atom_renumbering()
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

        modeller_result = self.run_modeller(
            ch_to_fix, brk_list, modeller_key,
            extra_gap, extra_NTerm=0
        )

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
        """ Runs modeller
            Args:
                *ch_to_fix* (list(str)): List of chain ids to be fixed
                *brk_list* (list(residue tuples)): Break to fix
                *modeller_key* (str): Modeller license key (optional). If not used Modeller installation license will be used.
                *extra_gap* (int): Additional residues to be taked either side of the gap. Use when obtained model have too long peptide distances (optional, default:0)
                *extra_NTerm* (int): Additional residues to be modelled on the N Terminus
                *sequende_Data* (SequenceData): SequenceData object containing canonical and structure sequences
                *templates* (list(structures)): Structures to be used as additional templates.
        """
        if modeller_key:
            MODELLER_ENV_VAR, MODELLER_INSTALL_ENV_VAR, modeller_install_dir = _guess_modeller_env()
            if not os.environ.get(MODELLER_ENV_VAR):
                os.environ[MODELLER_ENV_VAR] = modeller_key
            if not os.environ.get(MODELLER_INSTALL_ENV_VAR):
                os.environ[MODELLER_INSTALL_ENV_VAR] = modeller_install_dir

        try:
            from biobb_structure_checking.modeller_manager import ModellerManager, NoCanSeqError
        except ImportError:
            sys.exit("Error importing modeller")

        mod_mgr = ModellerManager()
        if not sequence_data:
            sequence_data = self.sequence_data

        mod_mgr.sequences = sequence_data

        modif_residues = []

        for mod in self.st:
            if self.models_data.has_models():
                print(f"Processing Model {mod.id + 1}")
                self.save_structure(opj(mod_mgr.tmpdir, 'templ.pdb'), mod.id)
            else:
                self.save_structure(opj(mod_mgr.tmpdir, 'templ.pdb'))

            for ch_id in self.chains_data.chain_ids[mod.id]:
                if ch_id not in ch_to_fix:
                    continue
                if sequence_data.data[mod.id][ch_id]['pdb']['wrong_order']:
                    print(f"Warning: chain {ch_id} has a unusual residue numbering, skipping")
                print(f"Fixing chain/model {ch_id}/{mod.id}")

                try:
                    model_pdb = mod_mgr.build(mod.id, ch_id, extra_NTerm)
                except NoCanSeqError as err:
                    print(err.message)
                    continue

                warnings.filterwarnings('ignore', 'BioPythonWarning')

                parser = PDBParser(PERMISSIVE=1, QUIET=True)
                model_st = parser.get_structure(
                    'model_st',
                    opj(mod_mgr.tmpdir, model_pdb['name'])
                )

                modif_set_residues = self.merge_structure(
                    sequence_data,
                    model_st,
                    mod.id,
                    ch_id,
                    brk_list,
                    sequence_data.data[mod.id][ch_id]['pdb']['frgs'][0].features[0].location.start,
                    extra_gap
                )  # TODO consider use canonical numbering instead of defining offset
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

        for i in range(0, len(sequence_data.data[mod_id][ch_id]['pdb']['frgs']) - 1):
            loc_i = sequence_data.data[mod_id][ch_id]['pdb']['frgs'][i].features[0].location
            loc_ii = sequence_data.data[mod_id][ch_id]['pdb']['frgs'][i + 1].features[0].location
            seq_i = sequence_data.data[mod_id][ch_id]['pdb']['frgs'][i].features[2].location
            seq_ii = sequence_data.data[mod_id][ch_id]['pdb']['frgs'][i + 1].features[2].location

            gap_start = loc_i.end
            gap_end = loc_ii.start
            # Gap length taken from sequence gap to avoid PDB numbering issues
            gap_length = seq_ii.start - seq_i.end - 1
            # Offset to account for breaks in PDB residue numbering
            seq_off_i_ii = gap_end - seq_ii.start - gap_start + seq_i.end

            if [self.st[mod_id][ch_id][gap_start], self.st[mod_id][ch_id][gap_end]] not in brk_list:
                # Checking for incomplete gap build needed for fixing side chains with rebuild
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

            print(f"Fixing {mu.residue_id(self.st[mod_id][ch_id][gap_start])}"
                  f" - {mu.residue_id(self.st[mod_id][ch_id][gap_end])}")

            # Superimposes structures using fragments at both sides of the gap
            fixed_ats = []
            moving_ats = []

            # checking whether there is a chain id in the model (Support for Modeller >= 10)
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
            while pos < len(list_res) and\
                    self.st[mod_id][ch_id].child_list[pos].id[1] != gap_start - extra_gap:
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
                    print(f"  Replacing {mu.residue_id(res)}")
                else:
                    print(f"  Adding {mu.residue_id(res)}")

                modif_residues.append(self.st[mod_id][ch_id][nres])

            print()

        return modif_residues

    def add_main_chain_caps(self, caps_list: Iterable[Iterable[str]]) -> List[str]:
        """ Adds ACE and NME caps """
        # print(caps_list)
        fixed = []
        for cap in caps_list:
            if cap[0] == 'N':
                if cap[1] in self.st_data.next_residue:
                    mu.add_ACE_cap_at_res(cap[1], self.st_data.next_residue[cap[1]])
                else:
                    raise NotEnoughAtomsError
            else:
                if cap[1] in self.st_data.prev_residue:
                    mu.add_NME_cap_at_res(cap[1], self.st_data.prev_residue[cap[1]])
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

        self.st_data.atom_renumbering()
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

        for res in self.st_data.all_residues:
            if mu.is_hetatm(res):
                continue
            protein_res = self.chains_data.get_chain_type(res) == mu.PROTEIN

            if remove_h:
                mu.remove_H_from_r(res, verbose=False)

            if res not in self.st_data.prev_residue:
                prev_residue = None
            else:
                prev_residue = self.st_data.prev_residue[res]

            if res not in self.st_data.next_residue:
                next_residue = None
            else:
                next_residue = self.st_data.next_residue[res]

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
                    print(f"Replacing {mu.residue_id(res)} by {ion_res_list[res]}")

                error_msg = mu.add_hydrogens_side(
                    res,
                    self.res_library,
                    ion_res_list[res],
                    h_rules[ion_res_list[res]],
                    protein_res
                )
                res.resname = ion_res_list[res]
                self.st_data.non_canonical_residue_list.append(
                    {'res':res, 'can_res':rcode_can, 'new_res':res.resname}
                )
            else:
                error_msg = mu.add_hydrogens_side(
                    res, self.res_library, rcode, h_rules, protein_res
                )

            if error_msg:
                print(error_msg, mu.residue_id(res))

        self.st_data.residue_renumbering(self.data_library)
        if add_charges:
            self.update_atom_charges(add_charges)
        self.st_data.atom_renumbering()
        self.modified = True

    def mark_ssbonds(self, cys_list):
        """ Mark Cys residue in cys_list as CYX for further usage """
        for res in cys_list:
            if 'HG' in res:
                mu.remove_atom_from_res(res, 'HG')
            res.resname = 'CYX'
            self.st_data.non_canonical_residue_list.append({'res':res, 'can_res':'CYS', 'new_res':res.resname})
        self.modified = True

    def rename_terms(self, term_res):
        """ Rename Terminal residues as NXXX or CXXX in proteins or XX5 XX3 in NA """
        for term in term_res:
            if term[0] in ('N', 'C'):
                if term[1].resname not in ('ACE', 'NME') and len(term[1].resname) == 3:
                    term[1].resname = term[0] + term[1].resname
            elif term[0] in ('5', '3'):
                term[1].resname = term[1].resname + term[0]

    def revert_terms(self):
        """ Reverts special term residue names to canonical ones"""
        for res in self.st.get_residues():
            if mu.is_hetatm(res):
                continue
            if self.chains_data.get_chain_type(res) == mu.PROTEIN:
                if len(res.get_resname()) == 4:
                    res.resname = res.resname[1:]
            elif self.chains_data.get_chain_type(res) in (mu.DNA, mu.RNA, mu.NA):
                if res.get_resname()[-1] in ('5', '3'):
                    res.resname = res.resname[:-1]

    def revert_can_resnames(self, canonical=True):
        """ Revert residue names to canonical ones """
        if canonical:
            for mod_res in self.st_data.non_canonical_residue_list:
                mod_res['res'].resname = mod_res['can_res']
        else:
            for mod_res in self.st_data.non_canonical_residue_list:
                mod_res['res'].resname = mod_res['new_res']

    def is_N_term(self, res: Residue) -> bool:
        """ Detects whether it is N terminal residue."""
        return self.chains_data.get_chain_type(res) == mu.PROTEIN and\
            res not in self.st_data.prev_residue

    def is_C_term(self, res: Residue) -> bool:
        """ Detects whether it is C terminal residue."""
        return self.chains_data.get_chain_type(res) == mu.PROTEIN and\
            res not in self.st_data.next_residue

    def is_5_term(self, res: Residue) -> bool:
        """ Detects whether it is 5' terminal residue."""
        return self.chains_data.get_chain_type(res) in (mu.DNA, mu.RNA) and\
            res not in self.st_data.prev_residue

    def is_3_term(self, res: Residue) -> bool:
        """ Detects whether it is 3' terminal residue."""
        return self.chains_data.get_chain_type(res) in (mu.DNA, mu.RNA) and\
            res not in self.st_data.next_residue

    def prepare_mutations(self, mut_list: str) -> List[MutationSet]:
        """ Find residues to mutate from mut_list"""
        mutations = MutationManager(mut_list, self.chains_data.chain_ids)
        mutations.prepare_mutations(self.st)
        return mutations

    def apply_mutations(self, mutations: MutationManager) -> Residue:
        """ Perform mutations """
        mutated_res = mutations.apply_mutations(
            self.data_library.get_mutation_map(),
            self.res_library
        )
        self.st_data.residue_renumbering(self.data_library)
        self.st_data.atom_renumbering()
        self.modified = True
        return mutated_res

    def rebuild_mutations(
        self,
        mutations: MutationManager,
        modeller_key: str = ''
    ) -> Residue:
        """ Perform mutations Rebuilding side chain"""
        ch_to_fix = {}
        brk_list = {}
        num_fix = 0
        for mod in self.st.get_models():
            ch_to_fix[mod.id] = set()
            brk_list[mod.id] = []
            for mut_set in mutations.mutation_list:
                for mut in mut_set.mutations:
                    # Checking if protein on model 0
                    if self.chains_data.chain_ids[mod.id][mut['chain']] > 1:
                        continue
                    ch_to_fix[mod.id].add(mut['chain'])
                    start_res = mut['resobj']
                    if start_res in self.st_data.prev_residue:
                        start_res = self.st_data.prev_residue[start_res]
                    end_res = mut['resobj']
                    if end_res in self.st_data.next_residue:
                        end_res = self.st_data.next_residue[end_res]
                    brk = [start_res, end_res]
                    brk_list[mod.id].append(brk)
            num_fix += len(ch_to_fix[mod.id])

        if not num_fix:
            print("No protein chains left, exiting")
            return []

        mutated_sequence_data = SequenceData()
        mutated_sequence_data.fake_canonical_sequence(self, mutations)
        for mut_set in mutations.mutation_list:
            for mut in mut_set.mutations:
                mu.remove_residue(mut['resobj'])
        mutated_sequence_data.read_structure_seqs(self)
        mutated_sequence_data.match_sequence_numbering(self)

        # TODO Not tested, to be used on changes in the NTerm residue
        extra_NTerm = 0
        mutated_res = self.run_modeller(
            ch_to_fix,
            brk_list,
            modeller_key,
            0,
            extra_NTerm,
            mutated_sequence_data
        )

        self.update_internals()
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

    def _amide_score(self, matr):
        score = 0.
        for amide_res in sorted(matr):
            for atm in matr[amide_res]['cnts']:
                for atm_cnt in matr[amide_res]['cnts'][atm]:
                    d = atm.element != atm_cnt.element
                    m1 = matr[amide_res]['mod']
                    m2 = atm_cnt.get_parent() in matr and matr[atm_cnt.get_parent()]['mod']
                    # !d !m1 !m2 + !d m1 m2 + d !m1 m2 + d m1 !m2
                    # !d (!m1 !m2 + m1 m2) + d (!m1 m2 + m1 !m2)
                    # !d !(m1^m2) + d (m1^m2)
                    # ! d^(m1^m2)
                    if not d^(m1^m2):
                        score += 1/matr[amide_res]['cnts'][atm][atm_cnt]
        return score

    def _amide_cluster(self, to_fix):
        cluster = {}

        for res in to_fix:
            cluster[res] = set()
            cluster[res].add(res)

        for r_pair in self.rr_dist:
            res1, res2 = r_pair[0:2]
            if res1 in to_fix and res2 in to_fix:
                for res in cluster[res2]:
                    cluster[res1].add(res)
                cluster[res2] = cluster[res1].copy()

        for res in to_fix:
            if res in cluster:
                for res2 in cluster[res]:
                    if res2 == res:
                        continue
                    if res2 in cluster:
                        del cluster[res2]
        return cluster

    def _is_amide_atom(self, amide_res, atm):
        res = atm.get_parent()
        return res.get_resname() in amide_res and atm.id in amide_res[res.get_resname()]

    def amide_auto_fix(self, to_fix):
        ''' Fix amide residues automatically'''
        print("Fixing automatically")
        amide_res = self.data_library.get_amide_data()[0]
        c_list = self.check_r_list_clashes(
            to_fix['res_to_fix'],
            ['polar'],
            get_all_contacts=True
        )
        for res_pair in c_list['polar']:
            for cnt in c_list['polar'][res_pair]:
                at1, at2, dist2 = cnt
                if at1.serial_number > at2.serial_number:
                    at2, at1, dist2 = cnt
        matr = {}
        for res_pair in c_list['polar']:
            for cnt in c_list['polar'][res_pair]:
                at1, at2, dist2 = cnt
                if at1.serial_number > at2.serial_number:
                    at2, at1, dist2 = cnt
                res1 = at1.get_parent()
                res2 = at2.get_parent()
                if self._is_amide_atom(amide_res, at1):
                    if res1 not in matr:
                        matr[res1] = {'mod':False, 'cnts':{}}
                    if at1 not in matr[res1]['cnts']:
                        matr[res1]['cnts'][at1] = {}
                    matr[res1]['cnts'][at1][at2] = dist2
                if self._is_amide_atom(amide_res, at2):
                    if res2 not in matr:
                        matr[res2] = {'mod':False, 'cnts':{}}
                    if at2 not in matr[res2]['cnts']:
                        matr[res2]['cnts'][at2] = {}
                    matr[res2]['cnts'][at2][at1] = dist2

        print(f"Initial contact score: {self._amide_score(matr):.3f}")
        print("Clustering amide residues")
        clusters = self._amide_cluster(to_fix['res_to_fix'])
        print(f"{len(clusters)} cluster(s) found, exploring...")
        to_fix = []
        nclust = 0
        for clust in clusters.values():
            to_fix_part = []
            nclust += 1
            amide_list = []
            mod_vec = max_vec = ''
            for amide_res in sorted(clust):
                amide_list.append(amide_res)
                max_vec += '1'
            print(f"Cluster {nclust}:{', '.join([mu.residue_id(r) for r in amide_list])}")
            conf = 0
            min_score = self._amide_score(matr)
            opt_vec = '0' * len(max_vec)
            while conf <= int(max_vec, 2):
                mod_vec = f"{'0' * len(max_vec)}{bin(conf)[2:]}"[-len(max_vec):]
                for pos, res in enumerate(amide_list):
                    matr[res]['mod'] = mod_vec[pos] == '1'
                score = self._amide_score(matr)
                if score < min_score:
                    min_score = score
                    opt_vec = mod_vec
                conf += 1
            for pos, res in enumerate(amide_list):
                if opt_vec[pos] == '1':
                    to_fix_part.append(res)
                    to_fix.append(res)
                matr[res]['mod'] = opt_vec[pos] == '1'
            if to_fix_part:
                print(
                    f"New score: {min_score:.3f}, fixed residue(s): "
                    f"{', '.join([mu.residue_id(r) for r in to_fix_part])}"
                )
            else:
                print("Score not improved, skipping")
        return to_fix

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
        # Prepared for std Duplexes
        mut_list = []
        i = 0
        nch = 0
        for mod in self.st:
            for ch_id, chn in self.sequence_data.data[mod.id].items():
                start = chn['pdb']['frgs'][0].features[0].location.start
                seq = chn['pdb']['frgs'][0].seq
                if len(seq) != len(mut_seq[nch]):
                    raise SequencesDoNotMatch()
                prefix = ''
                if chn['pdb']['type'] == mu.DNA:
                    prefix = 'D'
                for i, r in enumerate(seq):
                    mut_list.append(
                        f"{ch_id}/{mod.id}:{prefix}{r}{start + i}{prefix}{mut_seq[nch][i]}"
                    )
                nch += 1
        return ','.join(mut_list)


# ===============================================================================


def _guess_modeller_env():
    """ Guessing Modeller version from conda installation if available """
    import subprocess
    conda_info = subprocess.run(['conda', 'list', 'modeller'], stdout=subprocess.PIPE)
    for line in conda_info.stdout.decode('ASCII').split('\n'):
        if 'modeller' in line:
            info = line.split()
    if info[1]:
        print(f"Modeller v{info[1]} detected")
        ver1, ver2 = info[1].split('.')
        return f"KEY_MODELLER{ver1}v{ver2}", f"MODINSTALL{ver1}v{ver2}", f"{os.environ.get('CONDA_PREFIX','')}/lib/modeller-{ver1}.{ver2}"

    print("Modeller version not detected, using default")
    return 'KEY_MODELLER', 'MODINSTALL', 'modeller'
# ===============================================================================

class WrongServerError(Exception):
    def __init__(self):
        self.message = 'ERROR: Biounits supported only on MMB server'
class UnknownFileTypeError(Exception):
    def __init__(self, typ):
        self.message = f'ERROR: unknown filetype ({typ})'
class OutputPathNotProvidedError(Exception):
    def __init__(self):
        self.message = 'ERROR: output PDB path not provided'
class NotAValidResidueError(Exception):
    def __init__(self, res):
        self.message = f'Warning: {res} is not a valid residue in this context'
class NotEnoughAtomsError(Exception):
    def __init__(self):
        self.message = 'Warning: not enough backbone to build missing atoms'

class ParseError(Exception):
    def __init__(self, err_id, err_txt):
        self.message = f'{err_id} ({err_txt}) found when parsing input structure'
class UnknownFFError(Exception):
    def __init__(self, ff):
        self.message = f'{ff} is not a valid ff for assigning atom types'
class SequencesDoNotMatch(Exception):
    def __init__(self):
        self.message = "Sequence lengths do not match"
