# Structure Checking from MDWeb

check_structure performs MDWeb structure checking set as a command line
utility.

It includes some structure manipulation options like selecting models or chains,
removing components of the system, completing missing atoms, and some quality
checking as residue quirality, amide orientation, or vdw clashes.
```
usage: check_structure [-h] [-i INPUT_STRUCTURE_PATH]
                       [--file_format FILE_FORMAT] [--sequence FASTA_SEQ_PATH]
                       [--pdb_server PDB_SERVER] [--cache_dir CACHE_DIR_PATH]
                       [--modeller_key MODELLER_KEY]
                       [--res_lib RES_LIBRARY_PATH]
                       [--data_lib DATA_LIBRARY_PATH]
                       [-o OUTPUT_STRUCTURE_PATH] [--json JSON_OUTPUT_PATH]
                       [--quiet] [--limit ATOM_LIMIT] [--debug] [--force_save]
                       [--rename_terms] [--check_only] [--non_interactive] [--version]
                       command ...


positional arguments:
  command               Command to execute (help: check_structure commands)
  options               Specific command options

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_STRUCTURE_PATH, --input INPUT_STRUCTURE_PATH Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}.
                        Biounits pdb:{pdbid}.{bn}. Biounits require MMB server
  -o OUTPUT_STRUCTURE_PATH, --output OUTPUT_STRUCTURE_PATH
                        Output structure. PDB Format
  --json JSON_OUTPUT_PATH
                        Summary activities on a json file
  --file_format FILE_FORMAT
                        Format for retrieving structures (default=mmCif|pdb|xml)
  --sequence FASTA_SEQ_PATH
                        Canonical sequence in FASTA format, pdb_chain[,chain] in header,
                        may be required for backbone rebuilding
  --pdb_server PDB_SERVER
                        Remote server for retrieving structures (default|MMB)
  --cache_dir CACHE_DIR_PATH
                        Path for structure's cache directory (default: ./tmpPDB)
  --modeller_key MODELLER_KEY
                        User key for modeller, required for backbone rebuilding,
                        register at https://salilab.org/modeller/registration.html

  --res_lib RES_LIBRARY_PATH
                        Override settings default residue library (AMBER prep format)
  --data_lib DATA_LIBRARY_PATH
                        Override settings default data library
  --quiet               Reduces output, removing labels and progress info
  --limit ATOM_LIMIT    Limit on number of atoms, 0: nolimit
  --debug               Add debug information
  --force_save          Force saving an output file even if no modification
  --rename_terms        Label N-term and C-term as NXXX and CXXX residues (for Amber compatibility)
  --check_only          Perform checks only, structure is not modified
  --non_interactive     Do not prompt for missing parameters
  --version             show program's version number and exit
```

## Available commands:

```
commands:     Help on available commands
command_list: Run all tests from conf file
checkall:     Perform all checks without fixes
load:         Stores structure on local cache and provide basic statistics

1. System Configuration
=======================
models [--select model_num]
    Detect/Select Models
chains [--select chain_ids]
    Detect/Select Chains
inscodes [--renum]
    Detects residues with insertion codes. No fix provided.
altloc [--select occupancy| alt_id | list of res_id:alt_id]
    Detect/Select Alternative Locations
metals [--remove All | None | Met_ids_list | Residue_list]
    Detect/Remove Metals
ligands [--remove All | None | Res_type_list | Residue_list]
    Detect/Remove Ligands
hetatm [--remove All | None | Res_type_list | Residue_list] (v3.1)
    Detect/Remove Ligands, revert modified residues
water [--remove Yes|No]
    Remove Water molecules
rem_hydrogen [--remove Yes|No]
    Remove Hydrogen atoms from structure
mutateside [--mut mutation_list] [--no_check_clashes] [-rebuild]
    Mutate side chain with minimal atom replacement. Allows multiple mutations.
    Check generated clashes except --no_check_clashes set
    --rebuild Optimize side chains using Modeller.
add_hydrogen [--add_mode auto | pH | list | interactive | interactive_his] [--no_fix_side] [--keep_h] [--add_charges]
    Add Hydrogen Atoms. Auto: std changes at pH 7.0. His->Hie. pH: set pH value
    list: Explicit list as [*:]HisXXHid, Interactive[_his]: Prompts for all selectable residues
    Fixes missing side chain atoms unless --no_fix_side is set
    Existing hydrogen atoms are remove before adding new ones unless --keep_h set.
    --add_charges adds partial charges (from RES_LIBRARY) and autodock atom types, forces PDBQT output

2. Fix Structure Errors

amide  [--fix All|None|Residue List] [--no_recheck]
    Detect/Fix Amide atoms Assignment
    Amide contacts are recheck unless --no_recheck
chiral [--fix All|None|Residue List] [--no_check_clashes]
    Detect/Fix Improper quirality
    Checks generated clashes unless --no_check_clashes set
fixside [--fix All |None|Residue List] [--no_check_clashes] [--no_rem] [--rebuild]
    Complete side chains
    Checks generated clashes unless --no_check_clashes set
    --no_rem   Do not remove unknown atoms
    --rebuild  Rebuild complete side chain using Modeller
backbone [--fix_atoms All|None|Residue List]
         [--fix_main All|None|Break list]
         [--add_caps All|None|Terms|Breaks|list]
         [--extra_gap]
         [--no_recheck]
         [--no_check_clashes]
    Analyze main chain missing atoms and fragments.
    --fix_atoms O, OXT atoms can be fixed
    --fix_main Missing fragments filled using comparative modelling (Modeller License needed)
    --add_caps Adds ACE and NME residues as necessary, preserving existing atoms
    --extra_gap Recovers additional residues from model either side of the break, helps to fix loop connections
    Rechecks beckbone on each op unless --no_recheck is set.
    Generated clashes are checked unless --no_check_clashes

3. Structure Warnings

cistransbck Analyzes cis-trans dihedrals on backbone atoms
getss      Detect SS Bonds
clashes    Steric clashes (Severe, Apolar, Polar Donors, Polar Acceptors,
           Ionic Positive, Ionic Negative)
```
### Dependencies
* python >=3.6
* biopython
* numpy
* modeller 9.25 (optional)
* psutil (for performance debug, optional)
