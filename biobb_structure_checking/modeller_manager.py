""" Module to handle an interface to modeller """

import sys
import os
import uuid
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    from modeller import environ, log
    from modeller.automodel import automodel, assess
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

    def build(self, target_model, target_chain):
        """ Prepares Modeller input and builds the model """
        alin_file = self.tmpdir + "/alin.pir"
        tgt_seq = self.sequences.data[target_chain]['can'].seq
        templs = []
        knowns = []
        for ch_id in self.sequences.data[target_chain]['chains']:
            pdb_seq = self.sequences.data[ch_id]['pdb'][target_model][0].seq
            # Check N-term
            if ch_id == target_chain:
                nt_pos = tgt_seq.find(pdb_seq)
                tgt_seq = tgt_seq[nt_pos:]
            # Make alignment gaps from breaks
            for i in range(1, len(self.sequences.data[ch_id]['pdb'][target_model])):
                gap_len = self.sequences.data[ch_id]['pdb'][target_model][i].features[0].location.start\
                    - self.sequences.data[ch_id]['pdb'][target_model][i-1].features[0].location.end - 1
                pdb_seq += '-'*gap_len
                pdb_seq += self.sequences.data[ch_id]['pdb'][target_model][i].seq
            templs.append(
                SeqRecord(
                    pdb_seq,
                    'templ' + ch_id,
                    '',
                    'structureX:{}:{}:{}:{}:{}:::-1.00: -1.00'.format(
                        self.templ_file,
                        self.sequences.data[ch_id]['pdb'][target_model][0].features[0].location.start,
                        ch_id,
                        self.sequences.data[ch_id]['pdb'][target_model][-1].features[0].location.end,
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
        shutil.rmtree(self.tmpdir)

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
