# 
# Global data for structure manager
#

#=========================================================================
# Defaults. Don't touch unless knowing how
#=========================================================================

#Atom groups
metal_ats = [
    "MG", "MN", "MO", "Mg2", "Mn2", "Ca2", "ZN", "NI", "FE", 
    "Zn2", "Ni2", "Fe2", "CO", "CU", "HG", "Co2", "Cu2", 
    "CD", "AG", "AU", "Cd2", "CA", "Ca"
]
apolar_elements = ["C","S"]
polar_acceptor = ["O", "S"]
polar_donor = ["N"]
pos_ats = ["NZ", "NE", "NH1", "NH2"]
neg_ats = ["OD1", "OD2", "OE1", "OE2"]

# Relevant Distances
SS_DIST = 2.5
CA_DIST = 3.8
CA_DIST_ERR = 1.
STERIC_CLASH_DIST =1.
MAX_DIST = 4.5;
STERIC_APOLAR_CLASH_DIST = 2.9
STERIC_POLAR_CLASH_DIST = 3.1
IONIC_CLASH_DIST = 3.5
CA_CA_THRESHOLD = 15
