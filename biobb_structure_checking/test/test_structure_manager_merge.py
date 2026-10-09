""" Unit tests for the residue matching of StructureManager.merge_structure, run with pytest """
from types import SimpleNamespace

from biobb_structure_checking.structure_manager import _missing_res_ids, _model_indices


def fragment(pdb_first, pdb_last, seq_first):
    """ Fragment as in sequencedata, features[0] holds its PDB numbers, features[2] its
        canonical sequence positions """
    def feature(start, end):
        return SimpleNamespace(location=SimpleNamespace(start=start, end=end))
    return SimpleNamespace(features=[feature(pdb_first, pdb_last), None,
                                     feature(seq_first, seq_first + pdb_last - pdb_first)])


def residues(*ids):
    return [SimpleNamespace(id=(' ', num, icode)) for num, icode in ids]


def test_model_indices_insertion_codes():
    """ Kabat numbering (3V6F chain A): 52A and 82A-C are residues of their own """
    ids = [(51, ' '), (52, ' '), (52, 'A'), (53, ' '), (82, ' '), (82, 'A'), (82, 'B'), (82, 'C'), (83, ' ')]
    idx = _model_indices([fragment(51, 83, 1)], residues(*ids))
    assert idx[52] == 2   # 52A does not overwrite 52
    assert idx[53] == 4   # not 3: 52A comes before
    assert idx[83] == 9   # not 6


def test_model_indices_numbering_shift():
    """ 1F45 chain B: the numbers jump 74 -> 79 but the sequence misses 5 residues """
    ids = [(73, ' '), (74, ' '), (79, ' '), (80, ' ')]
    idx = _model_indices([fragment(73, 74, 73), fragment(79, 80, 80)], residues(*ids))
    assert (idx[74], idx[79], idx[80]) == (2, 8, 9)   # gap of 5 in the model: 3-7


def test_missing_res_ids_free_numbers():
    assert _missing_res_ids(74, 79, 4) == [(' ', 75, ' '), (' ', 76, ' '), (' ', 77, ' '), (' ', 78, ' ')]


def test_missing_res_ids_overflow_gets_insertion_code():
    ids = _missing_res_ids(74, 79, 6)
    assert ids[:4] == [(' ', n, ' ') for n in (75, 76, 77, 78)]
    assert ids[4:] == [(' ', 78, 'A'), (' ', 78, 'B')]


def test_missing_res_ids_no_real_gap():
    """ Rebuilding a residue (fixside --rebuild) calls the merge with gap_end == gap_start """
    assert _missing_res_ids(155, 155, 3) == [(' ', 156, ' '), (' ', 157, ' '), (' ', 158, ' ')]
