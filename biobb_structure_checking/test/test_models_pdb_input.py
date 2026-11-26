""" Model selection on a multi-model PDB file (no canonical sequence), run with pytest """
import os
import subprocess
import sys
from pathlib import Path

TEST_DIR = Path(__file__).parent
BIN = TEST_DIR.parents[1] / 'bin' / 'check_structure'


def test_select_model_from_pdb_file(tmp_path):
    """ The sequences of model 1 were matched before its canonical sequence data
        existed: KeyError: 1 in SequenceData.match_sequence_numbering """
    atoms = [line for line in (TEST_DIR / 'chiral_pdb_test' / '1ubq.pdb').read_text().splitlines()
             if line.startswith('ATOM')]
    models = tmp_path / '1ubq_2models.pdb'
    models.write_text(
        ''.join(f"MODEL     {num:4d}\n" + '\n'.join(atoms) + "\nENDMDL\n" for num in (1, 2)) + "END\n"
    )
    out = tmp_path / 'model2.pdb'
    result = subprocess.run(
        [sys.executable, str(BIN), '-i', str(models), '-o', str(out), '--force_save',
         '--non_interactive', 'models', '--select', '2'],
        capture_output=True, text=True, check=False,
        env=dict(os.environ, PYTHONPATH=str(TEST_DIR.parents[1])))
    assert result.returncode == 0, result.stderr
    assert sum(line.startswith('ATOM') for line in out.read_text().splitlines()) == len(atoms)
