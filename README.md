[![Codacy Badge](https://api.codacy.com/project/badge/Grade/376891e43cab4cc591fb78ea43dfd380)](https://www.codacy.com/app/jlgelpi/structureChecking?utm_source=mmb.irbbarcelona.org&amp;utm_medium=referral&amp;utm_content=gitlab/BioExcel/structureChecking&amp;utm_campaign=Badge_Grade)
# Structure Checking from MDWeb

checkStruc.py performs MDWeb structure checking set as a command line
utility.

It includes some structure manipulation options like selecting models or chains,
removing components of the system, completing missing atoms, and some quality
checking as residue quirality, amide orientation, or vdw clashes.

```
Usage:  checkStruc [-h|--help] command help|options 
                   -i input_pdb_path -o input_pdb_path
```

## Available commands:

```
commands:  This help
command_list:      Run all tests from conf file

1. System Configuration 

models     Detect/Select Models
chains     Detect/Select Chains 
altloc     Detect/Select Alt Locations 
metals     Detect/Remove Heavy Metals
ligands    Detect/Remove Ligands 
remwat     Remove Water molecules
remh       Remove Hydrogen atoms 

2. Fix Structure Errors

amide      Detect/Fix Amide Assignment
chiral     Detect/Fix Improper quirality
chiral_bck Detect/Fix Improper quirality on backbone CA's
fixside    Complete side chains 

3. Structure Warnings

getss      Detect SS Bonds 
cisbck     Unusual cis/trans backbone 
nonconres  Non Consecutive residues 
bckgaps    Detect missing residues
clashes    Steric clashes (Severe, Polar Donors, Polar Acceptors, Apolar
           Ionic Positive, Ionic Negative)
```
### Dependencies
* python 3.x
* biopython 
* numpy
* biobb_model (structure_manager)

