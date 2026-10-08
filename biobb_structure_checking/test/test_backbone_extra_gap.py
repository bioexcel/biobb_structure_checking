""" Backbone fix with extra_gap on a structure with close gaps, run with pytest.
    Needs Modeller and its key in the MODELLER_KEY environment variable """
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip('modeller')
if not os.environ.get('MODELLER_KEY'):
    pytest.skip('MODELLER_KEY not set', allow_module_level=True)

TEST_DIR = Path(__file__).parent
BIN = TEST_DIR.parents[1] / 'bin' / 'check_structure'


def run_backbone(tmp_path, extra_gap):
    out = tmp_path / 'fixed.pdb'
    result = subprocess.run(
        [sys.executable, str(BIN), '--modeller_key', os.environ['MODELLER_KEY'],
         '-i', str(TEST_DIR / 'backbone_test' / '1ubq_caps.pdb'), '-o', str(out),
         '--sequence', str(TEST_DIR / 'backbone_test' / '1ubq.fasta'), '--non_interactive',
         'backbone', '--fix_chain', 'all', '--add_caps', 'none', '--fix_atoms', 'all',
         '--extra_gap', str(extra_gap)],
        capture_output=True, text=True, check=False,
        env=dict(os.environ, PYTHONPATH=str(TEST_DIR.parents[1])))
    return result, out


@pytest.mark.parametrize('extra_gap', [0, 4])
def test_close_gaps_with_extra_gap(tmp_path, extra_gap):
    """ The gaps of 1ubq_caps are 2-4 residues apart, so with extra_gap their padding
        overlaps: residues inserted by a gap are replaced by the next one """
    result, out = run_backbone(tmp_path, extra_gap)
    assert 'ERROR' not in result.stderr, result.stderr
    residues = {line[22:27] for line in out.read_text().splitlines() if line.startswith('ATOM')}
    assert len(residues) == 76   # every missing residue is rebuilt
