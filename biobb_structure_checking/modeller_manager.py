"""
    Module to handle an interface to modeller,
    used to rebuild main and side chains and
    to optimize side chain orientation
"""

import sys
import os
from os.path import join as opj
import uuid
import shutil

import Bio
OLD_ALIGN = Bio.__version__ >= '1.79'
from Bio import SeqIO
if OLD_ALIGN:
    from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# Check for back-compatiblity with biopython < 1.77
try:
    from Bio.Seq import IUPAC
    has_IUPAC = True
except ImportError:
    has_IUPAC = False

try:
    from modeller import Environ, log
    from modeller.automodel import AutoModel, assess
except ImportError:
    sys.exit("Error importing Modeller package")

TMP_BASE_DIR = '/tmp'
DEBUG = False


class ModellerManager():
    """
    | modeller_manager ModellerManager
    | Class to handle Modeller calculations """
    def __init__(self):
        self.tmpdir = opj(TMP_BASE_DIR, "mod" + str(uuid.uuid4()))
        # self.tmpdir = "/tmp/modtest"
        # print("Using temporary working dir " + self.tmpdir)
        self.ch_id = ''
        self.sequences = None
        self.templ_file = 'templ.pdb'
        try:
            os.mkdir(self.tmpdir)
        except IOError as err:
            sys.exit(err)
        self.env = Environ()

        self.env.io.atom_files_directory = [self.tmpdir]
        log.none()

    def build(self, target_model, target_chain, extra_NTerm_res):
        """ ModellerManager.build
        Prepare Modeller input and builds the model

        Args:
            target_model (int) : Model to repair
            target_chain (str) : Chain to repair
            extra_NTerm_res (int) : Number of additional residues
                at NTerm (to fix NTerm, experimental)
        """
        alin_file = opj(self.tmpdir, "alin.pir")

        if target_chain not in self.sequences.has_canonical[target_model]:
            raise NoCanSeqError(target_model, target_chain)

        tgt_seq = self.sequences.data[target_model][target_chain]['can'].seq

        # triming N-term of canonical seq
        pdb_seq = self.sequences.data[target_model][target_chain]['pdb']['frgs'][0].seq
        nt_pos = max(tgt_seq.find(pdb_seq) - extra_NTerm_res, 0)
        tgt_seq = tgt_seq[nt_pos:]

        # TODO trim trailing residues in tgt_seq

        templs = []
        knowns = []
        for ch_id in self.sequences.data[target_model][target_chain]['chains']:
            frgs = self.sequences.data[target_model][ch_id]['pdb']['frgs']
            pdb_seq = frgs[0].seq
            for i in range(1, len(frgs)):
                frag_seq = frgs[i].seq
                pdb_seq += frag_seq
            # tuned to open gaps on missing loops only

            if not OLD_ALIGN:
                alin = self.sequences.aligner.align(tgt_seq, pdb_seq)
            else:
                alin = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)

            if has_IUPAC:
                pdb_seq = Seq(alin[0][1], IUPAC.protein)
            else:
                pdb_seq = Seq(alin[0][1])

            templs.append(
                SeqRecord(
                    pdb_seq,
                    f"templ{ch_id}",
                    "",
                    f"structureX:{self.templ_file}:"
                    f"{frgs[0].features[0].location.start}:"
                    f"{ch_id}:{frgs[-1].features[0].location.end}:"
                    f"{ch_id}:::-1.00: -1.00",
                    annotations={'molecule_type': 'protein'}  # required for writing PIR aligment
                )
            )
            knowns.append(f'templ{ch_id}')

            if ch_id == target_chain:
                tgt_seq = tgt_seq[0:len(pdb_seq)]

        _write_align(tgt_seq, templs, alin_file)

        return self._automodel_run(alin_file, knowns)

    def _automodel_run(self, alin_file, knowns):
        amdl = AutoModel(
            self.env,
            alnfile=alin_file,
            knowns=knowns,
            sequence='target',
            assess_methods=(assess.DOPE, assess.GA341)
        )
        amdl.starting_model = 1
        amdl.ending_model = 1

        # amdl.loop.starting_model = 1
        # amdl.loop.ending_model = 1

        orig_dir = os.getcwd()
        os.chdir(self.tmpdir)
        amdl.make()
        os.chdir(orig_dir)

        return amdl.outputs[0]

    def __del__(self):
        if not DEBUG:
            shutil.rmtree(self.tmpdir)
        else:
            print(f"Using temporary folder: {self.tmpdir}")


def _write_align(tgt_seq, templs, alin_file):
    SeqIO.write(
        [
            SeqRecord(
                tgt_seq,
                'target',
                '',
                'sequence:target:::::::0.00: 0.00',
                annotations={'molecule_type': 'protein'}
            )
        ] + templs,
        alin_file,
        'pir'
    )


class NoCanSeqError(Exception):
    """
    | modeller_manager NoCanSeqError
    | Error raised when no canonical sequence exists
    """
    def __init__(self, mod_id, ch_id):
        self.message = f"No canonical sequence found for chain {ch_id}/{mod_id}"
