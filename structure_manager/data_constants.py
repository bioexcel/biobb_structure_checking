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
#polar_acceptor = ["O", "SD","SG","OD1","OD2","OE1","OE2","OG","OG1"]
polar_acceptor = ["O", "SD","OD1","OD2","OE1","OE2"]
#polar_donor = ["N","NE","NZ","NH1","NH2","ND1","NE2", "SG", "OG", "OG1"]
polar_donor = ["N","NE","NZ","NH1","NH2","ND1","ND2","NE1","NE2"]
polar_both = ["OG", "OG1", "SG"]
pos_ats = ["NZ", "NE", "NH1", "NH2"]
neg_ats = ["OD1", "OD2", "OE1", "OE2"]

# Relevant Distances
SS_DIST = 2.5

CA_DIST = 3.8
CA_DIST_ERR = 1.

CLASH_DIST = {
    'severe': 1.,
    'apolar': 2.9,
    'donor': 3.1,
    'acceptor': 3.1,
    'positive': 3.5,
    'negative': 3.5
}

MAX_DIST = 4.5;

R_R_CUTOFF = 15
