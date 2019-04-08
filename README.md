[![Codacy Badge](https://api.codacy.com/project/badge/Grade/376891e43cab4cc591fb78ea43dfd380)](https://www.codacy.com/app/jlgelpi/structureChecking?utm_source=mmb.irbbarcelona.org&amp;utm_medium=referral&amp;utm_content=gitlab/BioExcel/structureChecking&amp;utm_campaign=Badge_Grade)
# Structure Checking from MDWeb

check_structure performs MDWeb structure checking set as a command line
utility.

It includes some structure manipulation options like selecting models or chains,
removing components of the system, completing missing atoms, and some quality
checking as residue quirality, amide orientation, or vdw clashes.

```
usage: checkStruc.py [-h] [-i INPUT_STRUCTURE_PATH] [-o OUTPUT_STRUCTURE_PATH]
                     [--version] [--data_dir DATA_DIR]
                     [--res_lib RES_LIB_PATH] [--data_lib DATA_LIBRARY_PATH]
                     [--json JSON_OUTPUT_PATH] [--quiet] [--check_only]
                     [--non_interactive] [--force_save]
                     [--pdb_server PDB_SERVER]
                     command ...
positional arguments:
  command               Command to execute (help: checkStruc commands)
  options               Specific command options

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_STRUCTURE_PATH, --input INPUT_STRUCTURE_PATH
                        Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}
  -o OUTPUT_STRUCTURE_PATH, --output OUTPUT_STRUCTURE_PATH
                        Output structure. Format PDB
  --version             show program's version number and exit
  --data_dir DATA_DIR   Override settings default data dir
  --res_lib RES_LIB_PATH
                        Override settings default residue library (AMBER prep
                        format)
  --data_lib DATA_LIBRARY_PATH
                        Override settings default data library
  --json JSON_OUTPUT_PATH
                        Cummulated checking results on a json file
  --quiet               Reduces output, removing labels and progress info
  --check_only          Perform checks, structure is not modified
  --non_interactive     Do not prompt for missing parameters
  --force_save          Force saving an output file even if no modification
  --pdb_server PDB_SERVER
                        Server for retrieving structures (default: RCSB|mmb)
  --load                Loads structure from PDB server into the local cache
  --stats               Loads structure and get basic statistics and headers

Available commands:

```

## Available commands:

```
commands:  This help
command_list:      Run all tests from conf file
checkall:   Perform all checks without fixes
load: Stores structure on local cache and provide basic statistics

1. System Configuration
=======================
models [--select_model model_num]     
    Detect/Select Models
chains [--select_chains chain_ids]    
    Detect/Select Chains
inscodes
    Detects residues with insertion codes (no fix)
altloc [--select_altloc occupancy| alt_id | list of res_id:alt_id]
    Detect/Select Alternative Locations
metals [--remove All | None | Met_ids_list | Residue_list]   
    Detect/Remove Metals
ligands [--remove All | None | Res_type_list | Residue_list]
    Detect/Remove Ligands
remwat [--remove Yes|No]
    Remove Water molecules
remh [remh --remove Yes|No]
    Remove Hydrogen atoms from structure
mutateside [--mut mutation_list]
    Mutate side chain with minimal atom replacement. Allows multiple mutations
addH [--mode auto | pH | interactive | interactive_his]
    Add Hydrogen Atoms

2. Fix Structure Errors

amide  [--fix All|None|Residue List]    
    Detect/Fix Amide atoms Assignment
chiral [--fix All|None|Residue List]
    Detect/Fix Improper quirality
fixside [--fix All |None|Residue List]    
    Complete side chains
backbone [--fix All|None|Residue List]   
    Analyze main chain missing atoms and fragments. O, OXT atoms can be fixed

3. Structure Warnings

cistransbck Analyzes cis-trans dihedrals on backbone atoms
getss      Detect SS Bonds
clashes    Steric clashes (Severe, Apolar, Polar Donors, Polar Acceptors,
           Ionic Positive, Ionic Negative)

```
### Dependencies
* python 3.x
* biopython
* numpy
* biobb_model (structure_manager)
