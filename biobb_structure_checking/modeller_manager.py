"""
    Module to handle an interface to modeller,
    used to rebuild main and side chains and
    to optimize side chain orientation
"""

import re
import sys
import os
from os.path import join as opj
import uuid
import shutil

import Bio
OLD_ALIGN = Bio.__version__ < '1.79'
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
    from modeller import Environ, Selection, log
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
        self.ch_id = ''
        self.sequences = None
        self.templ_file = 'templ.pdb'
        if not DEBUG:
            self.tmpdir = opj(TMP_BASE_DIR, "mod" + str(uuid.uuid4()))
            try:
                os.mkdir(self.tmpdir)
            except IOError as err:
                sys.exit(err)
        else:
            self.tmpdir = "/tmp/modtest"
            print("Using temporary working dir " + self.tmpdir)

        self.env = Environ()

        self.env.io.atom_files_directory = [self.tmpdir]
        log.none()

    def build(self, target_model, target_chain, extra_NTerm_res, fix_known=False):
        """ ModellerManager.build
        Prepare Modeller input and builds the model

        Args:
            target_model (int) : Model to repair
            target_chain (str) : Chain to repair
            extra_NTerm_res (int) : Number of additional residues
                at NTerm (to fix NTerm, experimental)
            fix_known (bool) : Optimize only the missing internal segments, the
                residues in the structure keep their coordinates
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
        gaps = []
        for ch_id in self.sequences.data[target_model][target_chain]['chains']:
            frgs = self.sequences.data[target_model][ch_id]['pdb']['frgs']
            pdb_seq = frgs[0].seq
            for i in range(1, len(frgs)):
                frag_seq = frgs[i].seq
                pdb_seq += frag_seq
            # tuned to open gaps on missing loops only

            # The chain to fix is placed according to the canonical position of its
            # fragments, an alignment may shift the gaps in low complexity regions
            aligned_seq = None
            if ch_id == target_chain:
                aligned_seq = _place_fragments(frgs, len(tgt_seq), nt_pos)
            if aligned_seq is None:
                if not OLD_ALIGN:
                    alin = self.sequences.aligner.align(tgt_seq, pdb_seq)
                else:
                    alin = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)
                aligned_seq = alin[0][1]

            if has_IUPAC:
                pdb_seq = Seq(aligned_seq, IUPAC.protein)
            else:
                pdb_seq = Seq(aligned_seq)

            if ch_id == target_chain:
                # (first, last) residues of every internal gap
                gaps = [(m.start() + 1, m.end()) for m in re.finditer('-+', str(pdb_seq))
                        if 0 < m.start() and m.end() < len(pdb_seq)]

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

        return self._automodel_run(alin_file, knowns, gaps if fix_known else None)

    def _automodel_run(self, alin_file, knowns, gaps=None):
        model_class = AutoModel if gaps is None else _gap_model_class(gaps)
        amdl = model_class(
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


def _place_fragments(frgs, length, offset):
    """ Place the fragments in their canonical position, filling the gaps with '-'

        Args:
            frgs: fragments of the chain, their features[2] hold the canonical position
            length (int): length of the target sequence
            offset (int): residues trimmed from the start of the canonical sequence
    """
    placed = ['-'] * length
    for frg in frgs:
        first = int(frg.features[2].location.start) - 1 - offset
        seq = str(frg.seq)
        if first < 0 or first + len(seq) > length or \
                set(placed[first:first + len(seq)]) != {'-'}:
            return None
        placed[first:first + len(seq)] = seq
    return ''.join(placed)


# Refining only part of the model: https://salilab.org/modeller/10.8/manual/node23.html
def _gap_model_class(gaps):
    """ AutoModel optimizing only the residues of the gaps, so the rest of the
        model keeps the coordinates of the template """
    class _GapAutoModel(AutoModel):
        def select_atoms(self):
            return Selection(*[self.residues[first - 1:last] for first, last in gaps])
    return _GapAutoModel


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
