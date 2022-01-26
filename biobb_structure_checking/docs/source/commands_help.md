
## commands Help

### Preset commands

* **commands** - _Help on available commands (this file)_
* **command -h** - _Specific help on command_
* **command_list** - _Run tests from conf file or inline string_
* **checkall** - _Perform all checks without fixes_
* **load** - _Stores structure on local cache and provides basic statistics_
***
### System Configuration Commands
Commands to manipulate structure composition. 

**sequences** - _Print canonical and structure sequences in FASTA format_

**models** [--select model_num(s)] [--superimpose] [--save_split]- _Detect/Select Models_
 * Accept List of Models (comma separated) or Model number range 
 * **--superimpose** - Superimposes currently selected models
 * **--save_split** Split models as separated output files. 

**chains** [--select chain_ids | molecule_type] - _Detect/Select Chains_

**inscodes** - _Detects residues with insertion codes. No fix provided (yet)_

**altloc** [--select occupancy| alt_id | list of res_id:alt_id] - _Detect/Select Alternative Locations_

**metals** [--remove All | None | Met_ids_list | Residue_list] - _Detect/Remove Metals_

**ligands** [--remove All | None | Res_type_list | Residue_list] - _Detect/Remove Ligands_

**getss** [--mark All | None | Residue_list] - _Detect SS Bonds_  
 * **--mark** - Replace relevant CYS by CYX to mark SS Bond (HG atom removed if present)

**water** [--remove Yes|No] - _Remove Water molecules_

**rem_hydrogen** [--remove Yes|No] - _Remove Hydrogen atoms from structure_

**mutateside** [--mut mutation_list|file:mutations_file] [--no_check_clashes] [--rebuild] -
_Mutate side chain with minimal atom replacement_  
* Allows multiple mutations (comma separated). 
* mutation_list as file: accepts list of mutations in a external file
* Check generated clashes except **--no_check_clashes** set.
* **--rebuild** optimize side chains using Modeller. 

**add_hydrogen** [--add_mode auto | pH | list | interactive | interactive_his] [--no_fix_side] [--keep_h] [--add_charges FF] - _Add Hydrogen Atoms to the strucure_  
* **--add_mode**
  * **Auto** - std changes at pH 7.0. His->Hie. pH: set pH value  
  * **list** - Explicit list as [*:]HisXXHid
  * **Interactive[_his]**: Prompts for all selectable residues  
* Fixes missing side chain atoms unless **--no_fix_side** is set.  
* Existing hydrogen atoms are removed before adding new ones unless **--keep_h** is set.  
* **--add_charges FF** adds partial charges (from RES_LIBRARY) and atom types from FF forcefield (Accepted: ATD, CMIP).  
* Output format taken from file extension (Accepted: pdb, pdbqt, pqr) or --output_format.

### Fix Structure Errors
Commands to detect and fix possible structure errors. 

**Amide** [--fix All|None|Residue List] [--no_recheck] - _Detect/Fix Amide atoms Assignment_
* Amide contacts are rechecked unless **--no_recheck** That can lead to infinite loops if done non-interactively.
  
**chiral** [--fix All|None|Residue List] [--no_check_clashes] - _Detect/Fix Improper side chain chirality_
* Checks for generated clashes unless **--no_check_clashes** set

**fixside** [--fix All |None|Residue List] [--no_check_clashes] - _Complete side chains (heavy atoms, protein only)_
* Checks generated clashes unless **--no_check_clashes** set
* **--rebuild**  Rebuild complete side chain using Modeller
  
**backbone** [--fix_atoms All|None|Residue List] [--fix_chain All|None|Break list] [--add_caps All|None|Break list] [--extra_gap]        [--no_recheck] [--no_check_clashes] - _Analyze main chain missing atoms and fragments (protein only)_
* **--fix_atoms** Add missing O, OXT backbone atoms.
* **--fix_chain** Missing fragments filled using comparative modelling (Modeller License needed)
* **--add_caps** Add ACE and NME residues as necessary, preserving existing atoms
* **--extra_gap** (Experimental) Recovers additional residues from model structure at either side of the break, helps to fix loop connections.
* Backbone is rechecked on each op unless **--no_recheck** is set. Use on non-interactive.
* Generated clashes are checked unless **--no_check_clashes**

### Structure Warnings
Additional checks on structure quality. No fix available.

**cistransbck** - _Analyzes cis-trans dihedrals on backbone atoms_

**clashes** - _Detect steric clashes in groups: Severe, Apolar, Polar Donors, Polar Acceptors, Ionic Positive, Ionic Negative_

***