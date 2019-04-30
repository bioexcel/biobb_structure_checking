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
commands:     Help on available commands
command_list: Run all tests from conf file
checkall:     Perform all checks without fixes
fixall:       Performs a default fix (v1.1)
load:         Stores structure on local cache and provide basic statistics

1. System Configuration 
=======================
models [--select model_num]     
    Detect/Select Models
chains [--select chain_ids]    
    Detect/Select Chains 
inscodes [--renum]
    Detects residues with insertion codes. --renum (v1.1) renums residues
altloc [--select occupancy| alt_id | list of res_id:alt_id]
    Detect/Select Alternative Locations 
metals [--remove All | None | Met_ids_list | Residue_list]   
    Detect/Remove Metals 
ligands [--remove All | None | Res_type_list | Residue_list]
    Detect/Remove Ligands 
hetatm [--remove All | None | Res_type_list | Residue_list] (v1.1)
    Detect/Remove Ligands, revert modified residues
water [--remove Yes|No]
    Remove Water molecules
rem_hydrogen [--remove Yes|No]
    Remove Hydrogen atoms from structure 
mutateside [--mut mutation_list] [--no_check_clashes]
    Mutate side chain with minimal atom replacement. Allows multiple mutations.
    Check generated clashes except --no_check_clashes set (v1.1)
add_hydrogen [--add_mode auto | pH | list | interactive | interactive_his] [--no_fix_side] [--keep_h]
    Add Hydrogen Atoms. Auto: std changes at pH 7.0. His->Hie. pH: set pH value
    list: Explicit list as [*:]HisXXHid, Interactive[_his]: Prompts for all selectable residues
    Fixes missing side chain atoms unless --no_fix_side is set
    Existing hydrogen atoms are remove before adding new ones unless --keep_h set. 

2. Fix Structure Errors

amide  [--fix All|None|Residue List] [--no_recheck]    
    Detect/Fix Amide atoms Assignment
    Amide contacts are recheck unless --no_recheck
chiral [--fix All|None|Residue List] [--no_check_clashes]
    Detect/Fix Improper quirality
    Checks generated clashes unless --no_check_clashes set
fixside [--fix All |None|Residue List] [--no_check_clashes]
    Complete side chains 
    Checks generated clashes unless --no_check_clashes set
backbone [--fix_atoms All|None|Residue List]
         [--fill_main All|None|Break List] (v1.1)
         [--add_caps All|None|Break list] (v1.1)
         [--no_recheck]
         [--no_check_clashes]
    Analyze main chain missing atoms and fragments.
    --fix_atoms O, OXT atoms can be fixed
    --fill_main Missing fragments filled using comparative modelling (Modeller License needed)
    --add_caps Adds ACE and NME residues as necessary
    Rechecks beckbone on each op unless --no_recheck is set.
    Generated clashes are checked unless --no_check_clashes
    

3. Structure Warnings

cistransbck Analyzes cis-trans dihedrals on backbone atoms
getss      Detect SS Bonds 
clashes    Steric clashes (Severe, Apolar, Polar Donors, Polar Acceptors, 
           Ionic Positive, Ionic Negative)

```
### Dependencies
* python 3.6
* biopython 1.72
* numpy
* biobb_structure_manager
