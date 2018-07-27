#
# Template for structureCheck settings
#
#
# Set this r
base_dir_path = "STRUCTURECHECK_HOME_DIR"

#=========================================================================
# Default paths.override at command line
#=========================================================================

help_dir_path = base_dir_path + "/helpdocs"

data_dir_path = base_dir_path + "/dat"

residue_library = data_dir_path + "all_amino03.lib"

#=========================================================================
# Defaults. Don't touch unless knowing how
#=========================================================================

metal_ats = [
    "MG", "MN", "MO", "Mg2", "Mn2", "Ca2", "ZN", "NI", "FE", 
    "Zn2", "Ni2", "Fe2", "CO", "CU", "HG", "Co2", "Cu2", 
    "CD", "AG", "AU", "Cd2", "CA", "Ca"
]
