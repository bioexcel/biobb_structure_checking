# check_structure

## Structure Checking from MDWeb

check_structure performs [MDWeb](https://mmb.irbbarcelona.org/MDWeb) structure checking set as a command line utility. It is intended to prepare a structure for molecular dynamics simulation.

It includes structure manipulation options like selecting models or chains, removing components of the system, completing side chains and backbone, and quality checking as residue quirality, amide orientation, or vdw clashes.

check_structure can be run interactively. It will prompt for any missing parameter or information, but also can be run standalone providing a list of operations. 



~~~
usage: check_structure [-h] [-i INPUT_STRUCTURE_PATH]
                       [--file_format FILE_FORMAT] [--sequence FASTA_SEQ_PATH]
                       [--pdb_server PDB_SERVER] [--cache_dir CACHE_DIR_PATH]
                       [--modeller_key MODELLER_KEY]
                       [--res_lib RES_LIBRARY_PATH]
                       [--data_lib DATA_LIBRARY_PATH]
                       [-o OUTPUT_STRUCTURE_PATH] [--rename_terms]
                       [--json JSON_OUTPUT_PATH] [-nv] [-v]
                       [--limit ATOM_LIMIT] [--debug] [--force_save]
                       [--check_only] [--non_interactive] [--version]
                       command [command_options]

~~~
### positional arguments:
  **command** - _Command to execute (required)_

  **command_options** - _Specific command options (optional)_
* On a interactive session, parameters required for **command** execution will be prompt as necessary, but can be also introduced as command options in the command line. see [Commands Help](https://biobb-structure-checking.readthedocs.io/en/latest/commands_help.html). 
***
### Arguments for input:

**-i --input** INPUT_STRUCTURE_PATH - _Input structure._
* Formats PDB|mmCIF|PQR. Taken from file extension.
* Remote **pdb:{pdbid}**. See **--file_format** for selecting download format (default: cif)
* Biounits **pdb:{pdbid}.{bn}**. Biounits require MMB server (**--pdb_server MMB**). Format PDB.

**--sequence** FASTA_SEQ_PATH - _Canonical sequence in FASTA format_
* Header should start >pdb_chain[,chain] for backbone rebuild. Required only for PDB/PQR structures.

**--file_format** FILE_FORMAT - _Format for retrieving structures (default=cif|pdb|xml)_

**--pdb_server** PDB_SERVER - _Remote server for retrieving structures (default=rcsb|MMB)_

**--cache_dir** CACHE_DIR_PATH - _Path for structure's cache directory (default: ./tmpPDB)_

**--limit** ATOM_LIMIT - _Limit on number of atoms (0: nolimit)_

***
### Additional data input

**--res_lib** RES_LIBRARY_PATH - _Override default residue library (AMBER prep format)_

**--data_lib** DATA_LIBRARY_PATH - _Override default data library_

***
### Arguments for output

**-o --output** OUTPUT_STRUCTURE_PATH - _Output structure._
* PDB|PDBQT|PQR|CMIP Formats Available (selected from file extension)
  
**--json** JSON_OUTPUT_PATH - _Store a summary of all activities on a json file_

**--force_save** - _Force saving an output file even if no modification_

**--rename_terms** - _Label N-term and C-term as NXXX and CXXX residues in output files_

**-nv --quiet** - _Reduces output, removing labels and progress info_

**-v --verbose** - _Additional progress info_

**--debug** - _Add debug information (Memory usage, timings)_

***
### Configuration arguments
  
**--check_only** - _Perform checks only, structure is not modified_

**--non_interactive** - _Do not prompt for missing parameters_

***
### Miscelanea
**--modeller_key** MODELLER_KEY - _User key for Modeller, required for backbone rebuild unless included in Modeller installation. Register at https://salilab.org/modeller/registration.html_
  

**-h, --help** - _show this help message and exit_

**--version** - _show program's version number and exit_

***