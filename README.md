# Structure Checking from MDWeb

Biobb_structure_checking performs 3D structure quality checking intended to facilitate the setup of molecular dynamics simulation of protein or nucleic acids systems.

Biobb_structure_checking package allows to configure the system (selection of model/chains,alternative location, addition of disulfide bonds and hydrogen atoms, side chain mutations), detects and fixes structure errors (missing side chain atoms, backbone breaks, amide assignments, incorrect chirality).

The latest documentation of this package can be found in our readthedocs site:
[latest package documentation](http://biobb_structure_checking.readthedocs.io/en/latest/).


```
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
                       command ...


positional arguments:
  command               Command to execute (help: check_structure commands)
  options               Specific command options

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_STRUCTURE_PATH, --input INPUT_STRUCTURE_PATH Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}.
                        Biounits pdb:{pdbid}.{bn}. Biounits require MMB server
  -o OUTPUT_STRUCTURE_PATH, --output OUTPUT_STRUCTURE_PATH
                        Output structure. PDB|PDBQT|PQR|CMIP Formats
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
  command               Command to execute (help: check_structure commands)
  options               Specific command options

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_STRUCTURE_PATH, --input INPUT_STRUCTURE_PATH
                        Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}
  --file_format FILE_FORMAT
                        Format for retrieving structures
                        (default=mmCif|pdb|xml)
  --sequence FASTA_SEQ_PATH
                        Canonical sequence in FASTA format, pdb_chain[,chain]
                        in header
  --pdb_server PDB_SERVER
                        Server for retrieving structures (default|MMB)
  --cache_dir CACHE_DIR_PATH
                        Path for structure's cache directory (default:
                        ./tmpPDB)
  --modeller_key MODELLER_KEY
                        User key for modeller, required for backbone fix,
                        register at
                        https://salilab.org/modeller/registration.html
  --res_lib RES_LIBRARY_PATH
                        Override settings default residue library (AMBER prep
                        format)
  --data_lib DATA_LIBRARY_PATH
                        Override settings default data library
  -o OUTPUT_STRUCTURE_PATH, --output OUTPUT_STRUCTURE_PATH
                        Output structure. Format PDB|PDBQT|PQR|CMIP
  --rename_terms        Show terminal residues as NXXX, CXXX in output files
  --json JSON_OUTPUT_PATH
                        Summary checking results on a json file
  -nv, --quiet          Minimal output, removing labels and progress info
  -v, --verbose         Add extra progress info
  --limit ATOM_LIMIT    Limit on number of atoms,0:nolimit
  --debug               Add debug information (timings, resources)
  --force_save          Force saving an output file even if no modification
  --check_only          Perform checks only, structure is not modified
  --non_interactive     Do not prompt for missing parameters
  --version             show program's version number and exit

```
### Dependencies
* python >=3.7
* biopython >=1.78
* numpy
* modeller 10.1 (optional)
* psutil (for performance debug, optional)
