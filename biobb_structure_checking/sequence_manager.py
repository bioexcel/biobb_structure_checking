
""" Module to manage sequence information for structures """

import sys
from Bio.PDB.Polypeptide import PPBuilder
from typing import List, Dict

from Bio.Seq import Seq, MutableSeq
#Deprecated in Biopython 1.78
#from Bio.Seq import IUPAC
from Bio.SeqUtils import IUPACData
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from biobb_structure_checking.model_utils import PROTEIN

# Check for back-compatiblity with biopython < 1.77
try:
    from Bio.Seq import IUPAC
    has_IUPAC = True
except ImportError:
    has_IUPAC = False

class SequenceData():
    """ Class to manage sequence data """
    def __init__(self):
        self.data = {}
        self.has_canonical = False
        self.fasta = []

    def add_empty_chain(self, ch_id:str):
        """ Add base structure for a new chain """
        self.data[ch_id] = {
            'can':None,
            'chains': [],
            'pdb':{}
        }

    def load_sequence_from_fasta(self, fasta_sequence_path):
        """ Loads canonical sequence from external FASTA file"""
        read_ok = True
        self.fasta = []
        if fasta_sequence_path:
            try:
                for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                    self.fasta.append(record)
            except IOError:
                sys.exit("Error loading FASTA")
        if not self.fasta:
            print("WARNING: No valid FASTA formatted sequences found in {} ".format(fasta_sequence_path))
            read_ok = False
        return read_ok
        

    def read_sequences(self, strucm, clean=True, cif_warn=False):
        """ Extracts sequences """
        if clean:
            self.data = {}
            self.has_canonical = {}

        if not self.has_canonical:
            self.read_canonical_seqs(strucm, cif_warn)

        self.read_structure_seqs(strucm)
        self.match_sequence_numbering()

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
                return True
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
                return True



        for i in range(0, len(chids)):
            for ch_id in chids[i].split(','):
                if ch_id not in strucm.chain_ids:
                    continue
                if strucm.chain_ids[ch_id] != PROTEIN:
                    continue
                if ch_id not in self.data:
                    self.add_empty_chain(ch_id)
                if has_IUPAC:
                    new_seq = Seq(seqs[i].replace('\n', ''), IUPAC.protein)
                else:
                    new_seq = Seq(seqs[i].replace('\n', ''))
                self.data[ch_id]['can'] = SeqRecord(
                        new_seq,
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
        
        self.has_canonical = {}
        for ch_id in strucm.chain_ids:
            if strucm.chain_ids[ch_id] != PROTEIN:
                continue
            self.has_canonical[ch_id] = (ch_id in self.data) and hasattr(self.data[ch_id]['can'], 'seq')
            if not self.has_canonical[ch_id]:
                print("Warning, no canonical sequence available for chain {}".format(ch_id))
        return False

    def read_structure_seqs(self, strucm):
        """ Extracts sequences from structure fragments """
        # PDB extrated sequences
        for mod in strucm.st:
            ppb = PPBuilder()
            for chn in mod.get_chains():
                seqs = []
                ch_id = chn.id
                wrong_order = False
                for frag in ppb.build_peptides(chn):
                    start = frag[0].get_id()[1]
                    start_index = frag[0].index
                    end = frag[-1].get_id()[1]
                    end_index = frag[-1].index
                    frid = '{}:{}-{}:{}-{}'.format(ch_id, start, end, start_index, end_index)
                    sqr = SeqRecord(
                        frag.get_sequence(),
                        'pdbsq_' + frid,
                        'pdbsq_' + frid,
                        'PDB sequence chain ' + frid
                    )
                    if start < end:
                        sqr.features.append(SeqFeature(FeatureLocation(start, end)))
                        sqr.features.append(SeqFeature(FeatureLocation(start_index, end_index)))
                    else:
                        print("Warning: unusual residue numbering at chain ", ch_id)
                        print("Warning: chain reconstruction may not be available")
                        sqr.features.append(SeqFeature(FeatureLocation(end, start)))
                        wrong_order = True
                    seqs.append(sqr)
                if ch_id not in self.data:
                    self.add_empty_chain(ch_id)
                self.data[ch_id]['pdb'][mod.id] = {
                    'frgs': seqs,
                    'wrong_order': wrong_order
                }

    def match_sequence_numbering(self):
        """ Assign canonical sequence numbering to structural fragments """
        if not hasattr(self, 'has_canonical'):
            return False
        for ch_id in self.data:
            if ch_id not in self.has_canonical or not self.has_canonical[ch_id]:
                continue
            for mod_id in self.data[ch_id]['pdb']:
                frgs = self.data[ch_id]['pdb'][mod_id]['frgs']
                self.data[ch_id]['pdb'][mod_id]['match_numbering'] = True
                for nfrag in range(0, len(frgs)):
                    inic = self.data[ch_id]['can'].seq.find(frgs[nfrag].seq) + 1
                    fin = inic + len(frgs[nfrag].seq) - 1
                    self.data[ch_id]['pdb'][mod_id]['frgs'][nfrag].features.append(
                        SeqFeature(FeatureLocation(inic, fin))
                    )
                    if inic != frgs[nfrag].features[0].location.start or\
                        fin != frgs[nfrag].features[0].location.end:
                        self.data[ch_id]['pdb'][mod_id]['match_numbering'] = False
        return True

    def fake_canonical_sequence(self, strucm, mutations):
        """ Fakes a canonical sequence to support modeller use
            in fixside and mutateside --rebuild """
        self.read_structure_seqs(strucm)
        self.has_canonical = {}
        for ch_id in strucm.chain_ids:
            if strucm.chain_ids[ch_id] != PROTEIN:
                continue
            #build sequence from frags filling gaps with G
            #Warning IUPAC deprecated in Biopython 1.78
            if has_IUPAC:
                seq = MutableSeq('', IUPAC.protein)
            else:
                seq = MutableSeq('')
            last_pos = 0
            start_pos = 0
            for frag in self.data[ch_id]['pdb'][0]['frgs']:
                if not start_pos:
                    start_pos = frag.features[0].location.start
                if last_pos:
                    seq += 'G'*(frag.features[0].location.start - last_pos-1)

                seq += frag.seq
                last_pos = frag.features[0].location.end
            for mut_set in mutations.mutation_list:
                for mut in mut_set.mutations:
                    if mut['chain'] != ch_id:
                        continue
                    res_num = mut['residue'][1]
                    seq[res_num - int(start_pos)] = IUPACData.protein_letters_3to1[mut['new_id'].capitalize()]
            if ch_id not in self.data:
                self.add_empty_chain(ch_id)
            #Warning IUPAC deprecated in Biopython 1.78
            if has_IUPAC:
                new_seq = Seq(str(seq).replace('\n', ''), IUPAC.protein)
            else:
                new_seq = Seq(str(seq).replace('\n', ''))
            self.data[ch_id]['can'] = SeqRecord(
                new_seq,
                'csq_' + ch_id,
                'csq_' + ch_id,
                'canonical sequence chain ' + ch_id,
                annotations = {'molecule_type':'protein'}
            )
            self.data[ch_id]['can'].features.append(
                SeqFeature(FeatureLocation(start_pos, start_pos + len(seq) - 1))
            )
            self.data[ch_id]['chains'].append(ch_id)
            self.has_canonical[ch_id] = True
