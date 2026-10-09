""" Unit tests for modeller_manager, run with pytest. Skipped if Modeller is missing """
from types import SimpleNamespace

import pytest
from Bio.Align import PairwiseAligner

pytest.importorskip('modeller')
from biobb_structure_checking.modeller_manager import _place_fragments  # noqa: E402


def fragment(seq, first):
    """ Fragment of the structure with its canonical position, as in sequencedata """
    location = SimpleNamespace(start=first, end=first + len(seq) - 1)
    return SimpleNamespace(seq=seq, features=[None, None, SimpleNamespace(location=location)])


# 3HMW chain H: S133 is followed by three missing residues (S134 K135 S136) and T137
TARGET = 'GPSVFPLAPSSKSTSGGTAALG'
FRAGMENTS = [fragment('GPSVFPLAPS', 1), fragment('TSGGTAALG', 14)]
EXPECTED = 'GPSVFPLAPS---TSGGTAALG'
WRONG    = 'GPSVFPLAP---STSGGTAALG'


def test_alignment_shifts_gap():
    """ The reason for placing the fragments: S133 is aligned to S136 """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    observed = ''.join(str(frg.seq) for frg in FRAGMENTS)
    assert aligner.align(TARGET, observed)[0][1] != EXPECTED
    assert aligner.align(TARGET, observed)[0][1] == WRONG


def test_place_fragments():
    assert _place_fragments(FRAGMENTS, len(TARGET), 0) == EXPECTED


def test_place_fragments_trimmed_nterm():
    """ The first 3 residues of the canonical sequence are trimmed from the target """
    assert _place_fragments(FRAGMENTS, len(TARGET) - 3, 3) is None  # first fragment starts at 1
    assert _place_fragments(FRAGMENTS[1:], len(TARGET) - 3, 3) == '-' * 10 + 'TSGGTAALG'


def test_place_fragments_do_not_fit():
    assert _place_fragments([fragment('GPS', 1), fragment('GPS', 2)], 10, 0) is None  # overlap
    assert _place_fragments([fragment('GPS', 9)], 10, 0) is None  # beyond the target
