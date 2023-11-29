""" Module to manage sequence information for structures """

from biobb_structure_checking.constants import FASTA_DOWNLOAD_PREFIX
import biobb_structure_checking.modelling.utils as mu
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import IUPACData
from Bio.Seq import Seq, MutableSeq
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import sys
import os
# from typing import List, Dict
from urllib.request import urlretrieve, urlcleanup

# pairwise2 to be deprecated, replaced by PairwiseAligner
# But Alignment structure has changed from v1.79
import Bio
OLD_ALIGN = Bio.__version__ < '1.79'
if OLD_ALIGN:
    from Bio import pairwise2
else:
    from Bio.Align import PairwiseAligner
# Deprecated in Biopython 1.78
# from Bio.Seq import IUPAC
try:
    from Bio.Seq import IUPAC
    has_IUPAC = True
except ImportError:
    has_IUPAC = False

IDENT_THRES = 0.7


class SequenceData():
    """
    | sequence_manager SequenceData
    | Class to manage sequence data
    """

    def __init__(self):
        self.data = {}
        self.has_canonical = {}
        self.fasta = []
        self.raw_pdb_seq = ''
        self.aligner = None  # for compatibility biopython <= 1.80
        if not OLD_ALIGN:
            self.aligner = PairwiseAligner()
            self.aligner.mode = 'global'
            self.aligner.open_gap_score = -5
            self.aligner.extend_gap_score = -1

    def add_empty_chain(self, mod, ch_id: str):
        """ SequenceData.add_empty_chain
        Add base structure for a new chain

        Args:
            ch_id (str) : Id of the new chain
        """
        if mod.id not in self.data:
            self.data[mod.id] = {}
        self.data[mod.id][ch_id] = {
            'can': None,
            'chains': [],
            'pdb': {}
        }

    def load_sequence_from_fasta(self, fasta_sequence_path):
        """ SequenceData.load_sequence_from_fasta
        Load canonical sequence from external FASTA file

        Args:
            fasta_sequence_path (str) : Path to FASTA file (pdb: forces download)
        """
        read_ok = True
        self.fasta = []
        if fasta_sequence_path:
            if fasta_sequence_path.startswith('pdb:'):
                fasta_sequence_path = fasta_sequence_path[4:]
                try:
                    tmp_file = f"/tmp/{fasta_sequence_path}.fasta"
                    url = f"{FASTA_DOWNLOAD_PREFIX}/{fasta_sequence_path}"
                    print(f"Retrieving sequence from {url}")
                    urlcleanup()
                    urlretrieve(url, tmp_file)
                    for record in SeqIO.parse(tmp_file, 'fasta'):
                        self.fasta.append(record)
                    os.remove(tmp_file)
                except IOError:
                    sys.exit("Error retrieving FASTA")
            elif fasta_sequence_path.startswith('http'):
                tmp_file = f'/tmp/{os.path.basename(fasta_sequence_path)}'
                print(f"Downloading sequence from {fasta_sequence_path} ...")
                try:
                    urlcleanup()
                    urlretrieve(fasta_sequence_path, tmp_file)
                    for record in SeqIO.parse(tmp_file, 'fasta'):
                        self.fasta.append(record)
                    os.remove(tmp_file)
                except IOError as e:
                    print(e)
                    print(f"Error retrieving FASTA")
            else:
                try:
                    for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                        self.fasta.append(record)
                except IOError:
                    sys.exit("Error loading FASTA")
        if not self.fasta:
            print(
                f"WARNING: No valid FASTA formatted sequences found in"
                f" {fasta_sequence_path}"
            )
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
        self.raw_pdb_seq = {}
        for mod in strucm.st:
            if clean:
                self.data[mod.id] = {}
                self.has_canonical[mod.id] = {}
            self.raw_pdb_seq[mod.id] = _get_pack_str_seqs(strucm)

            if not self.has_canonical[mod.id]:
                self.read_canonical_seqs(strucm, cif_warn)

            self.read_structure_seqs(strucm)
            self.match_sequence_numbering(strucm)

    def read_canonical_seqs(self, strucm, cif_warn):
        """ SequenceData.reachains_data.d_canonical_seqs
        Prepare canonical sequences from the available input

        Args:
            strucm (StructureManager) : Object containing the loaded structure
            cif_warn (bool) : Issue a warning when structure is not CIF
        """
        if not strucm.chains_data.chain_ids:
            strucm.chains_data.set_chain_ids()

        has_models = strucm.models_data.has_models()

        if self.fasta:
            hits = self.match_chain_seqs()
            chids = [h['ch_id'] for h in hits]
            seqs = [h['seq'] for h in hits]
            print('Getting canonical sequences from matching FASTA input')
            for hit in sorted(hits, key=lambda hit: hit['ch_id']):
                modtxt = ''
                if has_models:
                    modtxt = f"/{hit['mod_id']}"
                print(f"{hit['ch_id']}{modtxt}: \"{hit['desc']}\", score: {hit['score']} {hit['low']}")
        else:
            if strucm.st_data.input_format != 'cif':
                if cif_warn:
                    print("Warning: sequence features may not be available, use --sequence for ",
                          "external fasta input")
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
        for mod in strucm.st:
            if mod not in self.data:
                self.data[mod.id] = {}
            for index, chid in enumerate(chids):
                for ch_id in chid.split(','):
                    if ch_id not in strucm.chains_data.chain_ids[mod.id]:
                        continue
                    if strucm.chains_data.chain_ids[mod.id][ch_id] == mu.UNKNOWN:
                        continue
                    if ch_id not in self.data[mod.id]:
                        self.add_empty_chain(mod, ch_id)
                    if has_IUPAC:
                        new_seq = Seq(seqs[index].replace('\n', ''), IUPAC.protein)
                    else:
                        new_seq = Seq(seqs[index].replace('\n', ''))
                    modtxt = ''
                    modtxt2 = ''
                    if has_models:
                        modtxt = f"/{mod.id}"
                        modtxt2 = f" model {mod.id}"
                    self.data[mod.id][ch_id]['can'] = SeqRecord(
                        new_seq,
                        f"can_sq_{ch_id}{modtxt}",
                        f"can_sq_{ch_id}{modtxt}",
                        f"canonical sequence chain {ch_id}{modtxt2}"
                    )
                    self.data[mod.id][ch_id]['can'].features.append(
                        SeqFeature(FeatureLocation(1, len(seqs[index])))
                    )

                    # for chn in chids[i].split(','):
                    #     if chn in strucm.chain_ids:
                    #         self.data[ch_id]['chains'].append(chn)
                    seq_matches = self._assign_seq(self.data[mod.id][ch_id]['can'])
                    if not seq_matches:
                        print("Warning: unable to match sequences with structure")
                        self.data[mod.id][ch_id]['chains'].append(ch_id)
                    else:
                        max_score = max([match[2] for match in seq_matches])
                        for match in seq_matches:
                            if match[2] > IDENT_THRES * max_score:
                                self.data[mod.id][ch_id]['chains'].append(match[1])

            print(f"Canonical sequence for model {mod.id}:")
            self.has_canonical[mod.id] = {}
            for ch_id in strucm.chains_data.chain_ids[mod.id]:
                self.has_canonical[mod.id][ch_id] = (ch_id in self.data[mod.id]) and\
                    hasattr(self.data[mod.id][ch_id]['can'], 'seq')
                if not self.has_canonical[mod.id][ch_id]:
                    print(f"Warning, no canonical sequence available for chain {ch_id}/{mod.id}")
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

        for mod in strucm.st.get_models():
            if mod.id not in self.data:
                self.data[mod.id] = {}
            ppb = PPBuilder()
            for chn in mod.get_chains():
                seqs = []
                wrong_order = False
                if strucm.chains_data.chain_ids[mod.id][chn.id] == mu.PROTEIN:
                    frags = ppb.build_peptides(chn)
                    if not frags:
                        frags = [[
                            res
                            for res in chn.get_residues()
                            if not mu.is_hetatm(res)
                        ]]
                    # TODO patched for a Weird case where a second model lacks a chain
                    if not frags[0]:
                        print(
                            f"Warning: no protein residues found for chain {chn.id}"
                            f" at model {mod.id}, adding hetatm to avoid empty chain"
                        )
                        frags = list(chn.get_residues())
                elif strucm.chains_data.chain_ids[mod.id][chn.id] in\
                        (mu.DNA, mu.RNA, mu.NA):
                    frags = [[
                        res
                        for res in chn.get_residues()
                        if not mu.is_hetatm(res)
                    ]]
                else:
                    self.add_empty_chain(mod, chn.id)
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
                        seq = mu.get_sequence_from_list(
                            frag,
                            strucm.chains_data.chain_ids[mod.id][chn.id]
                        )

                    sqr = SeqRecord(
                        seq,
                        'pdbsq_' + frid,
                        'pdbsq_' + frid,
                        'PDB sequence chain ' + frid
                    )
                    if start < end:
                        sqr.features.append(
                            SeqFeature(FeatureLocation(start, end))
                        )
                        sqr.features.append(
                            SeqFeature(FeatureLocation(start_index, end_index))
                        )
                    else:
                        print("Warning: unusual residue numbering at chain ", chn.id)
                        print("Warning: chain reconstruction may not be available")
                        sqr.features.append(SeqFeature(FeatureLocation(end, start)))
                        wrong_order = True
                    seqs.append(sqr)
                    if chn.id not in self.data[mod.id]:
                        self.add_empty_chain(mod, chn.id)
                    self.data[mod.id][chn.id]['pdb'] = {
                        'frgs': seqs,
                        'wrong_order': wrong_order,
                        'type': strucm.chains_data.chain_ids[mod.id][chn.id]
                    }

            if can_reverted:
                strucm.revert_can_resnames(False)

    def match_sequence_numbering(self, strucm):
        """ SequenceData.match_sequence_numbering
        Assign canonical sequence numbering to structural fragments
        """
        if not hasattr(self, 'has_canonical'):
            return False
        for mod in strucm.st:
            for ch_id, ch_data in self.data[mod.id].items():
                # print(f"{ch_id=} {ch_data=}")
                # print(f"{self.has_canonical=} {mod.id=}")
                if ch_id not in self.has_canonical.get(mod.id, []) or not self.has_canonical.get(mod.id, dict()).get(ch_id, []):
                    continue
                frgs = ch_data['pdb']['frgs']
                ch_data['pdb']['match_numbering'] = True
                for nfrag, frag in enumerate(frgs):
                    inic = ch_data['can'].seq.find(frag.seq) + 1
                    fin = inic + len(frag.seq) - 1
                    ch_data['pdb']['frgs'][nfrag].features.append(
                        SeqFeature(FeatureLocation(inic, fin))
                    )
                    if inic != frag.features[0].location.start or\
                            fin != frag.features[0].location.end:
                        ch_data['pdb']['match_numbering'] = False
        return True

    def fake_canonical_sequence(self, strucm, mutations, ins_res_type='G'):
        """ SequenceData.fake_canonical_sequence
        Fakes a canonical sequence to support modeller use in fixside
            and mutateside --rebuild

        Args:
            strucm (StructureManager) : Object containing the loaded structure
            mutations (MutationsManager) : Object containing the list of
                mutations to perform
        """
        has_models = strucm.models_data.has_models()
        self.read_structure_seqs(strucm)
        self.has_canonical = {}
        for mod in strucm.st:
            self.has_canonical[mod.id] = {}
            for ch_id in strucm.chains_data.chain_ids[mod.id]:
                if strucm.chains_data.chain_ids[mod.id][ch_id] != mu.PROTEIN:
                    continue
                # build sequence from frags filling gaps with G
                # Warning IUPAC deprecated in Biopython 1.78
                if has_IUPAC:
                    seq = MutableSeq('', IUPAC.protein)
                else:
                    seq = MutableSeq('')
                last_pos = 0
                start_pos = 0
                for frag in self.data[mod.id][ch_id]['pdb']['frgs']:
                    if not start_pos:
                        start_pos = frag.features[0].location.start
                    if last_pos:
                        seq += ins_res_type * (
                            frag.features[0].location.start - last_pos - 1
                        )

                    seq += frag.seq
                    last_pos = frag.features[0].location.end
                for mut_set in mutations.mutation_list:
                    for mut in mut_set.mutations:
                        if mut['chain'] != ch_id:
                            continue
                        res_num = mut['resobj'].id[1]
                        seq[res_num - int(start_pos)] = IUPACData.protein_letters_3to1[mut['new_id'].capitalize()]
                if ch_id not in self.data[mod.id]:
                    self.add_empty_chain(mod, ch_id)
                # Warning IUPAC deprecated in Biopython 1.78
                if has_IUPAC:
                    new_seq = Seq(str(seq).replace('\n', ''), IUPAC.protein)
                else:
                    new_seq = Seq(str(seq).replace('\n', ''))
                modtxt = ''
                modtxt2 = ''
                if has_models:
                    modtxt = f"/{mod.id}"
                    modtxt2 = f" Model {mod.id}"
                self.data[mod.id][ch_id]['can'] = SeqRecord(
                    new_seq,
                    f"csq_{ch_id}{modtxt}",
                    f"csq_{ch_id}{modtxt}",
                    f"canonical sequence chain {ch_id}{modtxt2}",
                    annotations={'molecule_type': 'protein'}
                )
                self.data[mod.id][ch_id]['can'].features.append(
                    SeqFeature(
                        FeatureLocation(start_pos, start_pos + len(seq) - 1)
                    )
                )
                self.data[mod.id][ch_id]['chains'].append(ch_id)
                self.has_canonical[mod.id][ch_id] = True

    def get_canonical(self):
        """ SequenceData.get_canonical
        Prepares a FASTA string with the canonical sequence
        """
        outseq = ''
        for mod_id, data in self.data.items():
            for ch_id in sorted(data):
                if self.has_canonical[mod_id][ch_id]:
                    outseq += SeqIO.FastaIO.as_fasta(data[ch_id]['can'])
        return outseq

    def get_pdbseq(self):
        """ SequenceData.get_pdbseq
        Prepares a FASTA string with the structure sequence, and fragments
        """
        has_models = len(self.data) > 1
        # TODO re-use this on modeller_manager
        outseq = ''
        for mod_id, data in self.data.items():
            for ch_id in sorted(data):
                if ch_id in self.has_canonical[mod_id] and\
                        self.has_canonical[mod_id][ch_id]:
                    tgt_seq = data[ch_id]['can'].seq
                    frgs = data[ch_id]['pdb']['frgs']
                    pdb_seq = frgs[0].seq
                    if len(frgs) == 1:
                        frags_num = [
                            f"{frgs[0].features[0].location.start}-{frgs[0].features[0].location.end}"
                        ]
                    else:
                        frags_num = []
                        for i in range(1, len(frgs)):
                            frag_seq = frgs[i].seq
                            pdb_seq += frag_seq
                            frags_num.append(
                                f"{frgs[i].features[0].location.start}-{frgs[i].features[0].location.end}"
                            )
                    # tuned to open gaps on missing loops
                    if OLD_ALIGN:
                        alin = pairwise2.align.globalxs(tgt_seq, pdb_seq, -5, -1)
                    else:
                        alin = self.aligner.align(tgt_seq, pdb_seq)
                    if has_IUPAC:
                        pdb_seq = Seq(alin[0][1], IUPAC.protein)
                    else:
                        pdb_seq = Seq(alin[0][1])
                    modtxt = ''
                    if has_models:
                        modtxt = f"/{mod_id}"
                    seq = SeqRecord(
                        pdb_seq,
                        f"pdb_sq_{ch_id}{modtxt}",
                        '',
                        f"Frags: {','.join(frags_num)}"
                    )
                    outseq += SeqIO.FastaIO.as_fasta(seq)
                elif data[ch_id]['pdb']:
                    last_pos = 0
                    start_pos = 0
                    sequence = ''
                    frags_num = []
                    for frag in data[ch_id]['pdb']['frgs']:
                        if not start_pos:
                            start_pos = frag.features[0].location.start
                        if last_pos:
                            sequence += '-' * (
                                frag.features[0].location.start - last_pos-1
                            )
                        sequence += frag.seq
                        last_pos = frag.features[0].location.end
                        frags_num.append(
                            f"{frag.features[0].location.start}-{frag.features[0].location.end}"
                        )

                    if has_IUPAC:
                        pdb_seq = Seq(str(sequence), IUPAC.protein)
                    else:
                        pdb_seq = Seq(str(sequence))
                    modtxt = ''
                    if has_models:
                        modtxt = f"/{mod_id}"
                    seq = SeqRecord(
                        pdb_seq,
                        f"pdb_sq_{ch_id}{modtxt}",
                        '',
                        f"Frags: {','.join(frags_num)}"
                    )
                    outseq += SeqIO.FastaIO.as_fasta(seq)
                else:
                    print(f"Warning: chain {ch_id}/{mod_id} sequence cannot be recognized")
        return outseq

    def match_chain_seqs(self):
        """ Identifies chains in Fasta input """
        all_hits = []
        chids = set()
        for rec in self.fasta:
            for hit in self._assign_seq(rec):
                mod_id, ch_id, score = hit
                chids.add((mod_id, ch_id))
                all_hits.append({
                    'mod_id': mod_id,
                    'ch_id': ch_id,
                    'seq': str(rec.seq),
                    'desc': rec.description,
                    'score': score
                })

        hits = []

        for chain in chids:
            mod_id, ch_id = chain
            max_score = -100
            best_hit = ''
            for hit in all_hits:
                if hit['ch_id'] != ch_id or hit['mod_id'] != mod_id:
                    continue
                if hit['score'] > max_score:
                    best_hit = hit
                    max_score = hit['score']
            best_hit['low'] = ''
            if max_score < 0.5 * len(best_hit['seq']):
                best_hit['low'] = "*low"
            hits.append(best_hit)
        return hits

    def _assign_seq(self, rec):
        if hasattr(rec, 'seq'):
            tgt = rec.seq
        else:
            tgt = rec
        matches = []
        for mod_id in self.raw_pdb_seq:
            for ch_id in self.raw_pdb_seq[mod_id]:
                for seq in self.raw_pdb_seq[mod_id][ch_id]:
                    if not seq:  # for not protein/NA chains
                        continue
                    if OLD_ALIGN:
                        score = pairwise2.align.globalxs(
                            tgt, seq, -5, -1, score_only=True
                        )
                    else:
                        score = self.aligner.align(tgt, seq).score

                    if score > 0:
                        matches.append((mod_id, ch_id, score))
        return matches


def _get_pack_str_seqs(strucm):
    strucm.revert_can_resnames(canonical=True)
    seqs = {}
    for mod in strucm.st:
        for chn in mod.get_chains():
            if chn.id not in seqs:
                seqs[chn.id] = []
            seqs[chn.id].append(
                mu.get_sequence_from_list(
                    [res for res in chn.get_residues() if not mu.is_hetatm(res)],
                    strucm.chains_data.chain_ids[mod.id][chn.id]
                )
            )
    strucm.revert_can_resnames(canonical=False)
    return seqs
