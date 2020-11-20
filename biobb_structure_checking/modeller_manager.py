""" Module to handle an interface to modeller """

import sys
import os
import uuid
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Seq import Seq

try:
    from modeller import *
    from modeller.automodel import *
except:
    sys.exit("Error importing modeller")

# Check for back-compatiblity with biopython < 1.77
try:
    from Bio.Seq import IUPAC
    has_IUPAC = True
except ImportError:
    has_IUPAC = False

TMP_BASE_DIR = '/tmp'
DEBUG = False

class ModellerManager():
    """ Class to handle Modeller calculations """
    def __init__(self):
        self.tmpdir = TMP_BASE_DIR + "/mod" + str(uuid.uuid4())
        #self.tmpdir = "/tmp/modtest"
        print("Using temporary working dir " + self.tmpdir)
        self.ch_id = ''
        self.sequences = None
        self.templ_file = 'templ.pdb'
        try:
            os.mkdir(self.tmpdir)
        except IOError as err:
            sys.exit(err)
        self.env = environ()
        self.env.io.atom_files_directory = [self.tmpdir]
        log.none()

    def build(self, target_model, target_chain, extra_NTerm_res):
        """ Prepares Modeller input and builds the model """
        alin_file = self.tmpdir + "/alin.pir"
        tgt_seq = self.sequences.data[target_chain]['can'].seq
        #triming N-term of canonical seq
        pdb_seq = self.sequences.data[target_chain]['pdb'][target_model]['frgs'][0].seq
        nt_pos = max(tgt_seq.find(pdb_seq) - extra_NTerm_res, 0)
        tgt_seq = tgt_seq[nt_pos:]

        #TODO trim trailing residues in tgt_seq

        templs = []
        knowns = []

        for ch_id in self.sequences.data[target_chain]['chains']:
            frgs = self.sequences.data[ch_id]['pdb'][target_model]['frgs']
            pdb_seq = frgs[0].seq
            for i in range(1, len(frgs)):
                frag_seq = frgs[i].seq
                pdb_seq += frag_seq
            # tuned to open gaps on missing loops
            alin = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)

            if has_IUPAC:
                pdb_seq = Seq(alin[0][1], IUPAC.protein)
            else:
                pdb_seq = Seq(alin[0][1])

            templs.append(
                SeqRecord(
                    pdb_seq,
                    'templ' + ch_id,
                    '',
                    'structureX:{}:{}:{}:{}:{}:::-1.00: -1.00'.format(
                        self.templ_file,
                        frgs[0].features[0].location.start,
                        ch_id,
                        frgs[-1].features[0].location.end,
                        ch_id
                    ),
                    annotations = {'molecule_type':'protein'} # required for writing PIR aligment
                )
            )
            knowns.append('templ' + ch_id)

            if ch_id == target_chain:
                tgt_seq = tgt_seq[0:len(pdb_seq)]

        _write_alin(tgt_seq, templs, alin_file)

        return self._automodel_run(alin_file, knowns)


    def _automodel_run(self, alin_file, knowns):
        amdl = automodel(
            self.env,
            alnfile=alin_file,
            knowns=knowns,
            sequence='target',
            assess_methods=(assess.DOPE, assess.GA341)
        )
        amdl.starting_model = 1
        amdl.ending_model = 1

        #amdl.loop.starting_model = 1
        #amdl.loop.ending_model = 1

        orig_dir = os.getcwd()
        os.chdir(self.tmpdir)
        amdl.make()
        os.chdir(orig_dir)

        return amdl.outputs[0]

    def __del__(self):
        if not DEBUG:
            shutil.rmtree(self.tmpdir)

def _write_alin(tgt_seq, templs, alin_file):
    SeqIO.write(
        [
            SeqRecord(
                tgt_seq,
                'target',
                '',
                'sequence:target:::::::0.00: 0.00',
                annotations = {'molecule_type':'protein'}
            )
        ] + templs,
        alin_file,
        'pir'
    )
