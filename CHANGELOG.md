## v3.9.0 (2021.4)

### New functions

### Extended functions

- Jupyter Notebook support extended
  - Help() function
  - Settings accepted as Python dictionaries. 
        
- add_hydrogen
    - Added support for multiple atom type sets (currently ADP, CMIP)
    - Added support for pqr, pdbqt, cmip output formats
    - Added support for pqr input format
    - Terminal residue names available on 3 or 4 letter (N-, C-) codes

- command_list
    - Added support for inline command lists

### Bug Fixes
- Fix charge assignment for modified residues and terminals
- Fixed tests

***
## V3.8.5 (2021.2)

### New functions
- sequences
  - prints both canonical and structural sequences (protein or NA).

### Extended functions
- mutateside
    - added support for mutation of DNA/RNA residues
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
  - Uses Modeller for building the side chains, includes conformational search
- backbone --extra_gap
    - Allow to use more residues from the built model to fill backbone breaks
- add_hydrogen --add_charges
   - Adds point charges to atoms and produces a PDBQT output suitable for autodock

### Bug fixes
- --debug
    - Psutils import failed even when --debug was not requested
- backbone
    - Residue internal pointers corrupted when backbone was modified