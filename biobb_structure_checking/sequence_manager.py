
""" Module to manage sequence information for structures """


import sys
from Bio.PDB.Polypeptide import PPBuilder

from Bio.Seq import Seq, IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from biobb_structure_checking.model_utils import PROTEIN

class SequenceData():
    def __init__(self):
        self.data = {}
        self.has_canonical = False
        self.fasta = []

    def add_empty_chain(self, ch_id):
        self.data[ch_id] = {'can':None, 'chains': [], 'pdb':{}}

    def load_sequence_from_fasta(self, fasta_sequence_path):
        """ Loads canonical sequence from external FASTA file"""
        self.fasta = []
        if fasta_sequence_path:
            try:
                for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                    self.fasta.append(record)
            except IOError:
                sys.exit("Error loading FASTA")

    def read_sequences(self, strucm, clean=True, cif_warn=False):
        """ Extracts sequences"""
        if clean:
            self.data = {}
            self.has_canonical = False

        if not self.has_canonical:
            self.read_canonical_seqs(strucm, cif_warn)

        self.read_structure_seqs(strucm)

    def read_canonical_seqs(self, strucm, cif_warn):
        """ Prepare canonical sequences """

        if not strucm.chain_ids:
            strucm.set_chain_ids()

        if self.fasta:
            chids = []
            seqs = []
            for rec in self.fasta:
                chids.append(rec.id.split('_')[1])
                seqs.append(str(rec.seq))
        else:
            if strucm.input_format != 'cif':
                if cif_warn:
                    print("Warning: sequence features only available in mmCIF" +\
                    " format or with external fasta input")
                return 1
            #TODO check for NA
            if '_entity_poly.pdbx_strand_id' in strucm.headers:
                if not isinstance(strucm.headers['_entity_poly.pdbx_strand_id'], list):
                    chids = [strucm.headers['_entity_poly.pdbx_strand_id']]
                    seqs = [strucm.headers['_entity_poly.pdbx_seq_one_letter_code_can']]
                else:
                    chids = strucm.headers['_entity_poly.pdbx_strand_id']
                    seqs = strucm.headers['_entity_poly.pdbx_seq_one_letter_code_can']
            else:
                if cif_warn:
                    print("Warning: sequence data unavailable on cif data")
                return 1

        for i in range(0, len(chids)):
            for ch_id in chids[i].split(','):
                if ch_id not in strucm.chain_ids:
                    continue
                if strucm.chain_ids[ch_id] != PROTEIN:
                    continue
                if ch_id not in self.data:
                    self.add_empty_chain(ch_id)
                self.data[ch_id]['can'] = SeqRecord(
                    Seq(seqs[i].replace('\n', ''), IUPAC.protein),
                    'csq_' + ch_id,
                    'csq_' + ch_id,
                    'canonical sequence chain ' + ch_id
                )
                self.data[ch_id]['can'].features.append(
                    SeqFeature(FeatureLocation(1, len(seqs[i])))
                )

                for chn in chids[i].split(','):
                    if chn in strucm.chain_ids:
                        self.data[ch_id]['chains'].append(chn)

        self.has_canonical = True
        return 0

    def read_structure_seqs(self, strucm):
        """ Extracts sequences from structure"""
        # PDB extrated sequences
        for mod in strucm.st:
            ppb = PPBuilder()
            for chn in mod.get_chains():
                seqs = []
                #self.sequences[ch_id]['pdb'][mod.id] = [1]
                ch_id = chn.id
                wrong_order = False
                for frag in ppb.build_peptides(chn):
                    start = frag[0].get_id()[1]
                    end = frag[-1].get_id()[1]
                    frid = '{}:{}-{}'.format(ch_id, start, end)
                    sqr = SeqRecord(
                        frag.get_sequence(),
                        'pdbsq_' + frid,
                        'pdbsq_' + frid,
                        'PDB sequence chain ' + frid
                    )
                    if start < end:
                        sqr.features.append(SeqFeature(FeatureLocation(start, end)))
                    else:
                        print("Warning: unusual residue numbering at chain ", ch_id)
                        print("Warning: chain reconstruction may not be available")
                        sqr.features.append(SeqFeature(FeatureLocation(end, start)))
                        wrong_order = True
                    seqs.append(sqr)
                if ch_id not in self.data:
                    self.add_empty_chain(ch_id)
                self.data[ch_id]['pdb'][mod.id] = seqs
                self.data[ch_id]['pdb']['wrong_order'] = wrong_order
