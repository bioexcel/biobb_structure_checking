""" Module to manage sequence information for structures """

import sys
#from typing import List, Dict

from Bio import pairwise2, SeqIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq, MutableSeq
#Deprecated in Biopython 1.78
#from Bio.Seq import IUPAC
from Bio.SeqUtils import IUPACData
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
try:
    from Bio.Seq import IUPAC
    has_IUPAC = True
except ImportError:
    has_IUPAC = False
import biobb_structure_checking.modelling.utils as mu

IDENT_THRES = 0.7
class SequenceData():
    """
    | sequence_manager SequenceData
    | Class to manage sequence data
    """
    def __init__(self):
        self.data = {}
        self.has_canonical = False
        self.fasta = []
        self.raw_pdb_seq = ''

    def add_empty_chain(self, ch_id: str):
        """ SequenceData.add_empty_chain
        Add base structure for a new chain

        Args:
            ch_id (str) : Id of the new chain
        """
        self.data[ch_id] = {
            'can':None,
            'chains': [],
            'pdb':{}
        }

    def load_sequence_from_fasta(self, fasta_sequence_path):
        """ SequenceData.load_sequence_from_fasta
        Load canonical sequence from external FASTA file

        Args:
            fasta_sequence_path (str) : Path to FASTA file
        """
        read_ok = True
        self.fasta = []
        if fasta_sequence_path:
            try:
                for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                    self.fasta.append(record)
            except IOError:
                sys.exit("Error loading FASTA")
        if not self.fasta:
            print(f"WARNING: No valid FASTA formatted sequences found in {fasta_sequence_path}")
            read_ok = False
        return read_ok

    def read_sequences(self, strucm, clean=True, cif_warn=False):
        """ SequenceData.read_sequences
        Extracts sequences from input data

        Args:
            strucm (StructureManager) : Object containing loaded structure
            clean (bool) : (True) Clean existing sequence data
            cif_warn (bool) : (False) Issue a warning when structure is not CIF
        """
        if clean:
            self.data = {}
            self.has_canonical = {}

        self.raw_pdb_seq =  self._get_pack_str_seqs(strucm)

        if not self.has_canonical:
            self.read_canonical_seqs(strucm, cif_warn)

        self.read_structure_seqs(strucm)
        self.match_sequence_numbering()

    def read_canonical_seqs(self, strucm, cif_warn):
        """ SequenceData.reachains_data.d_canonical_seqs
        Prepare canonical sequences from the available input

        Args:
            strucm (StructureManager) : Object containing the loaded structure
            cif_warn (bool) : Issue a warning when structure is not CIF
        """

        if not strucm.chains_data.chain_ids:
            strucm.chains_data.set_chain_ids()

        if self.fasta:
            hits = self.match_chain_seqs()
            chids = [h['ch_id'] for h in hits]
            seqs = [h['seq'] for h in hits]
            print('Getting canonical sequences from matching FASTA input')
            for hit in sorted(hits, key=lambda hit:hit['ch_id']):
                print(f"{hit['ch_id']}: \"{hit['desc']}\", score: {hit['score']} {hit['low']}")
        else:
            if strucm.st_data.input_format != 'cif':
                if cif_warn:
                    print("Warning: sequence features only available in mmCIF" +\
                    " format or with external fasta input")
                return True
            if '_entity_poly.pdbx_strand_id' in strucm.st_data.headers:
                if not isinstance(strucm.st_data.headers['_entity_poly.pdbx_strand_id'], list):
                    chids = [strucm.st_data.headers['_entity_poly.pdbx_strand_id']]
                    seqs = [strucm.st_data.headers['_entity_poly.pdbx_seq_one_letter_code_can']]
                else:
                    chids = strucm.st_data.headers['_entity_poly.pdbx_strand_id']
                    seqs = strucm.st_data.headers['_entity_poly.pdbx_seq_one_letter_code_can']
            else:
                if cif_warn:
                    print("Warning: sequence data unavailable on cif data")
                return True

        for index, chid in enumerate(chids):
            for ch_id in chid.split(','):
                if ch_id not in strucm.chains_data.chain_ids:
                    continue
                if strucm.chains_data.chain_ids[ch_id] == mu.UNKNOWN:
                    continue
                if ch_id not in self.data:
                    self.add_empty_chain(ch_id)
                if has_IUPAC:
                    new_seq = Seq(seqs[index].replace('\n', ''), IUPAC.protein)
                else:
                    new_seq = Seq(seqs[index].replace('\n', ''))
                self.data[ch_id]['can'] = SeqRecord(
                    new_seq,
                    'can_sq_' + ch_id,
                    'can_sq_' + ch_id,
                    'canonical sequence chain ' + ch_id
                )
                self.data[ch_id]['can'].features.append(
                    SeqFeature(FeatureLocation(1, len(seqs[index])))
                )

                # for chn in chids[i].split(','):
                #     if chn in strucm.chain_ids:
                #         self.data[ch_id]['chains'].append(chn)
                seq_matches = self._assign_seq(self.data[ch_id]['can'])
                max_score = max([match[1] for match in seq_matches])
                for match in seq_matches:
                    if match[1] > IDENT_THRES * max_score:
                        self.data[ch_id]['chains'].append(match[0])

        self.has_canonical = {}
        for ch_id in strucm.chains_data.chain_ids:
#            if strucm.chain_ids[ch_id] != PROTEIN:
#                continue
            self.has_canonical[ch_id] = (ch_id in self.data) and\
                hasattr(self.data[ch_id]['can'], 'seq')
            if not self.has_canonical[ch_id]:
                print(f"Warning, no canonical sequence available for chain {ch_id}")
        return False

    def read_structure_seqs(self, strucm):
        """ SequenceData.read_structure_seqs
        Extract sequence from structure fragments

        Args:
            strucm (StructureManager) : Object containing the loaded structure
        """
        # PDB extrated sequences
        if strucm.st_data.non_canonical_residue_list:
            strucm.revert_can_resnames(True)
            can_reverted = True
        else:
            can_reverted = False
        
        for mod in strucm.st:
            ppb = PPBuilder()
            for chn in mod.get_chains():
                seqs = []
                wrong_order = False
                if strucm.chains_data.chain_ids[chn.id] == mu.PROTEIN:
                    frags = ppb.build_peptides(chn)
                    if not frags:
                        frags = [[res for res in chn.get_residues() if not mu.is_hetatm(res)]]
                    #TODO patched for a Weird case where a second model lacks a chain
                    if not frags[0]:
                        print(
                            f"Warning: no protein residues found for chain {chn.id}"
                            f" at model {mod.id}, adding hetatm to avoid empty chain"
                        )
                        frags = [[res for res in chn.get_residues()]]
                elif strucm.chains_data.chain_ids[chn.id] in (mu.DNA, mu.RNA, mu.NA):
                    frags = [[res for res in chn.get_residues() if not mu.is_hetatm(res)]]
                else:
                    self.add_empty_chain(chn.id)
                    frags = []

                for frag in frags:
                    start = frag[0].get_id()[1]
                    start_index = frag[0].index
                    end = frag[-1].get_id()[1]
                    end_index = frag[-1].index
                    frid = f"{chn.id}:{start}-{end}:{start_index}-{end_index}"
                    if hasattr(frag, 'get_sequence'):
                        seq = frag.get_sequence()
                    else:
                        seq = mu.get_sequence_from_list(frag, strucm.chains_data.chain_ids[chn.id])

                    sqr = SeqRecord(
                        seq,
                        'pdbsq_' + frid,
                        'pdbsq_' + frid,
                        'PDB sequence chain ' + frid
                    )
                    if start < end:
                        sqr.features.append(SeqFeature(FeatureLocation(start, end)))
                        sqr.features.append(SeqFeature(FeatureLocation(start_index, end_index)))
                    else:
                        print("Warning: unusual residue numbering at chain ", chn.id)
                        print("Warning: chain reconstruction may not be available")
                        sqr.features.append(SeqFeature(FeatureLocation(end, start)))
                        wrong_order = True
                    seqs.append(sqr)
                    if chn.id not in self.data:
                        self.add_empty_chain(chn.id)
                    self.data[chn.id]['pdb'][mod.id] = {
                        'frgs': seqs,
                        'wrong_order': wrong_order,
                        'type': strucm.chains_data.chain_ids[chn.id]
                    }
            if can_reverted:
                strucm.revert_can_resnames(False)

    def match_sequence_numbering(self):
        """ SequenceData.match_sequence_numbering
        Assign canonical sequence numbering to structural fragments
        """
        if not hasattr(self, 'has_canonical'):
            return False
        for ch_id, ch_data in self.data.items():
            if ch_id not in self.has_canonical or not self.has_canonical[ch_id]:
                continue
            for mod_id in ch_data['pdb']:
                frgs = ch_data['pdb'][mod_id]['frgs']
                ch_data['pdb'][mod_id]['match_numbering'] = True
                for nfrag, frag in enumerate(frgs):
                    inic = ch_data['can'].seq.find(frag.seq) + 1
                    fin = inic + len(frag.seq) - 1
                    ch_data['pdb'][mod_id]['frgs'][nfrag].features.append(
                        SeqFeature(FeatureLocation(inic, fin))
                    )
                    if inic != frag.features[0].location.start or\
                        fin != frag.features[0].location.end:
                        ch_data['pdb'][mod_id]['match_numbering'] = False
        return True

    def fake_canonical_sequence(self, strucm, mutations, ins_res_type='G'):
        """ SequenceData.fake_canonical_sequence
        Fakes a canonical sequence to support modeller use in fixside and mutateside --rebuild

        Args:
            strucm (StructureManager) : Object containing the loaded structure
            mutations (MutationsManager) : Object containing the list of mutations to perform
        """
        self.read_structure_seqs(strucm)
        self.has_canonical = {}
        for ch_id in strucm.chains_data.chain_ids:
            if strucm.chain_ids[ch_id] != mu.PROTEIN:
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
                    seq += ins_res_type*(frag.features[0].location.start - last_pos-1)

                seq += frag.seq
                last_pos = frag.features[0].location.end
            for mut_set in mutations.mutation_list:
                for mut in mut_set.mutations:
                    if mut['chain'] != ch_id:
                        continue
                    res_num = mut['resobj'].id[1]
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
                annotations={'molecule_type':'protein'}
            )
            self.data[ch_id]['can'].features.append(
                SeqFeature(FeatureLocation(start_pos, start_pos + len(seq) - 1))
            )
            self.data[ch_id]['chains'].append(ch_id)
            self.has_canonical[ch_id] = True

    def get_canonical(self):
        """ SequenceData.get_canonical
        Prepares a FASTA string with the canonical sequence
        """
        outseq = ''
        for ch_id in sorted(self.data):
            if self.has_canonical[ch_id]:
                outseq += SeqIO.FastaIO.as_fasta(self.data[ch_id]['can'])
        return outseq

    def get_pdbseq(self):
        """ SequenceData.get_pdbseq
        Prepares a FASTA string with the structure sequence, and fragments
        """
        #TODO re-use this on modeller_manager
        outseq = ''
        for ch_id in sorted(self.data):
            if ch_id in self.has_canonical and self.has_canonical[ch_id]:
                tgt_seq = self.data[ch_id]['can'].seq
                frgs = self.data[ch_id]['pdb'][0]['frgs']
                frags_num = []
                pdb_seq = frgs[0].seq
                for i in range(1, len(frgs)):
                    frag_seq = frgs[i].seq
                    pdb_seq += frag_seq
                    frags_num.append(
                        f"{frgs[i].features[0].location.start}-{frgs[i].features[0].location.end}"
                    )
                # tuned to open gaps on missing loops
                alin = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)

                if has_IUPAC:
                    pdb_seq = Seq(alin[0][1], IUPAC.protein)
                else:
                    pdb_seq = Seq(alin[0][1])

                seq = SeqRecord(
                    pdb_seq,
                    'pdb_sq_' + ch_id,
                    '',
                    'Frags: ' + ','.join(frags_num)
                )
                outseq += SeqIO.FastaIO.as_fasta(seq)
            elif self.data[ch_id]['pdb']:
                last_pos = 0
                start_pos = 0
                sequence = ''
                frags_num = []
                for frag in self.data[ch_id]['pdb'][0]['frgs']:
                    if not start_pos:
                        start_pos = frag.features[0].location.start
                    if last_pos:
                        sequence += '-'*(frag.features[0].location.start - last_pos-1)
                    sequence += frag.seq
                    last_pos = frag.features[0].location.end
                    frags_num.append(
                        f"{frag.features[0].location.start}-{frag.features[0].location.end}"
                    )

                if has_IUPAC:
                    pdb_seq = Seq(str(sequence), IUPAC.protein)
                else:
                    pdb_seq = Seq(str(sequence))

                seq = SeqRecord(
                    pdb_seq,
                    'pdb_sq_' + ch_id,
                    '',
                    'Frags: ' + ','.join(frags_num)
                )
                outseq += SeqIO.FastaIO.as_fasta(seq)
            else:
                print(f"Warning: chain {ch_id} sequence cannot be recongized")
        return outseq

    def match_chain_seqs(self):
        """ Identifies chains in Fasta input """
        all_hits = []
        chids = set()
        for rec in self.fasta:
            for hit in self._assign_seq(rec):
                ch_id, score = hit
                chids.add(ch_id)
                all_hits.append({
                    'ch_id':ch_id,
                    'seq': str(rec.seq),
                    'desc': rec.description,
                    'score': score
                })

        hits = []

        for ch_id in chids:
            max_score = -100
            best_hit=''
            for hit in all_hits:
                if hit['ch_id'] != ch_id:
                    continue
                if hit['score'] > max_score:
                    best_hit = hit
                    max_score = hit['score']
            best_hit['low'] = ''
            if max_score < 0.5 * len(best_hit['seq']):
                best_hit['low'] = "*low"
            hits.append(best_hit)
        return hits

    def _get_pack_str_seqs(self, strucm):
        strucm.revert_can_resnames(canonical=True)
        seqs = {}
        for mod in strucm.st:
            for chn in mod.get_chains():
                if chn.id not in seqs:
                    seqs[chn.id] = []
                seqs[chn.id].append(mu.get_sequence_from_list(
                    [res for res in chn.get_residues() if not mu.is_hetatm(res)],
                    strucm.chains_data.chain_ids[chn.id]
                    )
                )
        strucm.revert_can_resnames(canonical=False)
        return seqs

    def _assign_seq(self, rec):
        if hasattr(rec,'seq'):
            tgt = rec.seq
        else:
            tgt = rec
        matches = []
        for ch_id in self.raw_pdb_seq:
            for seq in self.raw_pdb_seq[ch_id]:
                if not seq: # for not protein/NA chains
                    continue
                score = pairwise2.align.globalxs(tgt, seq, -5, -1, score_only=True)
                if score > 0:
                    matches.append((ch_id, score))
        return matches