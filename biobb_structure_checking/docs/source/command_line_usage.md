# check_structure

## Structure Checking from MDWeb

check_structure performs [MDWeb](https://mmb.irbbarcelona.org/MDWeb) structure checking set as a command line utility. It is intended to prepare a structure for molecular dynamics simulation.

It includes structure manipulation options like selecting models or chains, removing components of the system, completing side chains and backbone, and quality checking as residue quirality, amide orientation, or vdw clashes.

check_structure can be run interactively. It will prompt for any missing parameter or information, but also can be run standalone providing a list of operations.



~~~
usage: check_structure [-h] [-i INPUT_STRUCTURE_PATH] [--file_format {mmCif,pdb,cif,xml}] [--coords_only]
                       [--sequence FASTA_SEQ_PATH] [--pdb_server PDB_SERVER] [--cache_dir CACHE_DIR_PATH]
                       [--nocache] [--overwrite_cache] [--copy_input COPY_INPUT] [--modeller_key MODELLER_KEY]
                       [--res_lib RES_LIBRARY_PATH] [--data_lib DATA_LIBRARY_PATH] [-o OUTPUT_STRUCTURE_PATH]
                       [--output_format {pdb,pdbqt,pqr,cmip,cif,mmCif}] [--keep_canonical_resnames] [--rename_terms]
                       [--json JSON_OUTPUT_PATH] [-nv] [-v] [--limit ATOM_LIMIT] [--time_limit TIME_LIMIT] [--debug]
                       [--build_warnings] [--force_save] [--check_only] [--non_interactive] [--version]
                       command ...
~~~
### positional arguments:
  **command** - _Command to execute (required)_

  **command_options** - _Specific command options (optional)_
* On a interactive session, parameters required for **command** execution will be prompted as necessary, but can be also introduced as command options in the command line. see [Commands Help](https://biobb-structure-checking.readthedocs.io/en/latest/commands_help.html).
***
### Arguments for input:

**-i --input** INPUT_STRUCTURE_PATH - _Input structure._
* Formats pdb(qt)|cif|pqr. Taken from file extension. pdbqt accepted, but read as pdb.
* Remote **pdb:{pdbid}[.format] | {url}**. See also **--file_format** for selecting download format (default: cif). For remote downloads as **pdb:** in pdb format, sequences are also automatically downloaded.
* Biounits/Assemblies **pdb:{pdbid}.{bn}**. Biounits/assemblies. Default format mmCIF.

**--file_format** {mmCif,cif,pdb,xml} - _Format for retrieving structures (mmCif(default)|cif|pdb|xml)_

**--build_warnings** - _Show structure building warnings (may indicate PDB errors)_

**--coords_only** - _Loads structure coordinates, discards chain labels and residue ids from input_

**--sequence** FASTA_SEQ_PATH - _Input canonical sequence in FASTA format_
* accepted **local files**, **URLs**, or **pdb:{pdbid}**. Required for PDB/PQR structures.

**--pdb_server** PDB_SERVER - _DEPRECATED. pdb: defaultr to wwPDB, use direct urls ofor other servers_

**--cache_dir** CACHE_DIR_PATH - _Path for structure's cache directory (default: ./tmpPDB)_

**--nocache** - _Do not cache downloaded structures_

**--overwrite_cache** - _Overwrite cached file if any_

**--copy_input** DIR - _Copy the downloaded structure in the indicated folder_

**--limit** ATOM_LIMIT - _Limit on number of atoms (0: no limit)_

**--time_limit** TIME_LIMIT - _Set limit on the execution time (sec), 0: nolimit. Default: 3600_

***
### Additional data input

**--res_lib** RES_LIBRARY_PATH - _Override default residue library (AMBER prep format)_

**--data_lib** DATA_LIBRARY_PATH - _Override default data library_

***
### Arguments for output

**-o --output** OUTPUT_STRUCTURE_PATH - _Output structure._
* cif|pdb|pdbqt|pqr|cmip formats available (use file extension or --output_format to set format)

**--output_format** OUTPUT_FORMAT - _Format for the Output._
* cif|pdb|pdbqt|pqr|cmip formats available (if empty file extension is used)

**--keep_canonical_resnames** - _Revert output to canonical residue names when modified by any operation_

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
