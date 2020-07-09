""" Module to handle an interface to modeller """

import sys
import os
import uuid
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq, IUPAC, MutableSeq

try:
    from modeller import *
    from modeller.automodel import *
except:
    sys.exit("Error importing modeller")

TMP_BASE_DIR = '/tmp'

class ModellerManager():
    """ Class to handle Modeller calculations """
    def __init__(self):
        self.tmpdir = TMP_BASE_DIR + "/mod" + str(uuid.uuid4())
        #self.tmpdir = "/tmp/modtest"
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
            pdb_seq = self.sequences.data[ch_id]['pdb'][target_model]['frgs'][0].seq
            prev_pos = len(pdb_seq)
            for i in range(1, len(self.sequences.data[ch_id]['pdb'][target_model]['frgs'])):
                frag_seq = self.sequences.data[ch_id]['pdb'][target_model]['frgs'][i].seq
                pdb_seq += frag_seq
            # tuned to open gaps on missing loops
            alignment = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)

            pdb_seq = Seq(alignment[0][1],IUPAC.protein)

            templs.append(
                SeqRecord(
                    pdb_seq,
                    'templ' + ch_id,
                    '',
                    'structureX:{}:{}:{}:{}:{}:::-1.00: -1.00'.format(
                        self.templ_file,
                        self.sequences.data[ch_id]['pdb'][target_model]['frgs'][0].features[0].location.start,
                        ch_id,
                        self.sequences.data[ch_id]['pdb'][target_model]['frgs'][-1].features[0].location.end,
                        ch_id
                    )
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

        orig_dir = os.getcwd()
        os.chdir(self.tmpdir)
        amdl.make()
        os.chdir(orig_dir)

        return amdl.outputs[0]

    def __del__(self):
        #shutil.rmtree(self.tmpdir)
        print(self.tmpdir)

def _write_alin(tgt_seq, templs, alin_file):
    SeqIO.write(
        [
            SeqRecord(
                tgt_seq,
                'target',
                '',
                'sequence:target:::::::0.00: 0.00'
            )
        ] + templs,
        alin_file,
        'pir'
    )
