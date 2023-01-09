## v3.12.0 (2022.4)
### Extended functions
- amide
  - Added --fix auto option to automatically find best fix combination
- chains
  - Added --rename. fixes empty chain labels
  - Added --renumber. allow to renumber/reorganize chains and residues
  - Added --rem_inscodes. removes insertion codes on renumbering
  - Added --rebuild. rebuild chains labels and residue ids from backbone connectivity
- inscodes
  - Added --renumber. Rebuild residue numbering to remove insertion codes
- models
  - Added --build_complex. Converts biounit's selected models into actual complexes
- sequences
  - Added --output_fasta. Writes sequences found in an external FASTA file
- json
  - Extended information on json summary output
- input format
  - File format for remote download can be defined using extension as in pdb:2ki5.pdb
- output format
  - Added mmCif output (only atom records)
- input management
  - Added --nocache to avoid caching downloaded structures
  - Added --copy_input to recover a copy of the input structure
  - Added --coords_only to discard chain labels and residue ids from input. Used to revover faulty structure files

### Bug Fixes
- Added missing defaults for Notebook execution
- Ionized/tautomeric residue names recognized
- Fixed behaviour of structure headers output with missing entries
## v3.10.1 (2022.3)
### Extended functions
- Structure details
  - Added hydrogen atoms count
- Input files
  - Support for PQR format
- backbone, sequences
  - Input Sequences are automatically assigned to the appropriate chain. No special requirements for FASTA headers
- add_hydrogen
  - Extended to support Nucleic Acids
- sequences
  - Structure sequence is always reported even in the absence of canonical one
- chains
  - Unlabelled chains can be fixed.

### Bug Fixes
- Fixed banner format
- Fixed residue id on metals output
- Removed required usage N and C terms 4-letter residue names except of requested output

## v3.9.11 (2021.4)

### New Functions
- Support for PDBQT files
- add_hydrogen
  - Multiple charge sets accepted (currently ADP, CMIP)

## v3.9.10 (2021.4)

### New Functions

### Extended functions
- Adding new option --keep_canonical_resnames to revert output to canonical residue names when modified by any operation.

- models
  - Select multiple models
  - Superimpose current models
  - Output models as multiple PDB files

### Bug Fixes

- add_hydrogen
  - Extended to support Nucleic Acids
  - added option to keep canonical residue names
### Bug Fixes
- add_hydrogen
  - Forced upper case for FF names
## v3.9.9 (2021.4)

### Bug Fixes
- Remove pinned numpy dependency

***
## v3.9.7 (2021.4)

### Bug Fixes
- Fixed tests
- Fixed behaviour with --non_interactive and missing command options


***
## v3.9.6 (2021.4)

### Extended functions

- Jupyter Notebook support extended
  - Help(command) function
  - Settings accepted also as Python dictionaries-

- add_hydrogen
    - Added support for multiple atom type sets (currently ADP, CMIP)
    - Added support for pqr, pdbqt, cmip output formats following file extension
    - Added support for pqr input format
    - Terminal residue names available on 3 or 4 letter (N-, C-) codes

- command_list
    - Added support for inline command lists

- mutationside
    - Added support for RNA mutations
### Bug Fixes
- Fix charge assignment for modified residues and terminals
- Fixed wrong assigment of RNA residues as protein
- Fixed tests
- Fixed behaviour with --non_interactive and missing command options

***
## V3.8.5 (2021.2)

### New functions
- sequences
  - prints both canonical and structural sequences (protein or NA).

### Extended functions
- mutateside
    - added support for mutation of DNA residues
    - --na_seq allows to set a desired final sequence in a single operation (for DNA duplexes)
- chains
    - allows to select chains according to molecular type (protein | dna | rna | na).
    - Improved guess of chain type.
- getss
    - can now mark CYS residues as part of SSBonds. Further commands like add_hydrogen reacts accordingly
- load
    - can now be used to convert downloaded cif to pdb format using --force_save
- Modeller based commands
    - adapted to support Modeller >= 10.1
- several commands
    - adapted to work with NA chains

- Extended support for verbosity
  - Default verbosity has been reduced to errors and warnigs
  - --nv --quiet removes all progress reports
  - -v  adds extra progress report (original default)

### Bugs fixed
- clashes, add_hydrogen did not work for modified or cap residues
- mutateside crashed when trying to delete already modified atoms
- added error message when chains will remove all chains in the structure
- added error message when mutateside has no available mutations

***
## v3.7.3 (2021.1)

### Extended functionality
- Added support for biopython 1.78

### Bugs fixed
- backbone
    - added error messages for incorrect or missing FASTA sequences
    - backbone reconstruction crashed when not all chains sequences were available even if missing chains where not needed

***
## v3.0.2

### New functions
- --rename_terms
    - Rename N and C terms as NXXX, CXXX

### Extended functions:
- mutateside & fixside --rebuild
  - Use Modeller for building the side chains, include conformational search
- backbone --extra_gap
    - Allows to use more residues from the built model to fill backbone breaks (experimental)
- add_hydrogen --add_charges
   - Adds point charges to atoms and produces a PDBQT output suitable for autodock

### Bug fixes
- --debug
    - Psutils import failed even when --debug was not requested
- backbone
    - Residue internal pointers corrupted when backbone was modified
