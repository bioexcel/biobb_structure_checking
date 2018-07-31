#!/bin/csh
setenv APPDIR `pwd`/..
#models
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_models_5_test.pdb --non_interactive models --select_model 5 > 1ark_models_test_5.log
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_models_1_test.pdb --non_interactive models --select_model 1 > 1ark_models_test_1.log
#chains
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_chains_test_all.pdb --non_interactive chains --select_chain All > 2ki5_chains_test_all.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_chains_test_B.pdb --non_interactive chains --select_chain B > 2ki5_chains_test_B.log
#altloc
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_altloc_test_occ.pdb --non_interactive altloc --select_altloc occupancy > 2ki5_altloc_test_occ.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_altloc_test_B.pdb --non_interactive altloc --select_altloc B > 2ki5_altloc_test_B.log
#ligands
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_SO4.pdb --non_interactive ligands --remove SO4 > 2ki5_ligands_test_SO4.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_All.pdb --non_interactive ligands --remove All > 2ki5_ligands_test_All.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_None.pdb --non_interactive ligands --remove None > 2ki5_ligands_test_None.log
#metals
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_SO4.pdb --non_interactive metals --remove ZN > 1bqo_metals_test_ZN.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_A305.pdb --non_interactive metals --remove A305 > 1bqo_metals_test_A305.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_All.pdb --non_interactive metals --remove All > 1bqo_metals_test_All.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_None.pdb --non_interactive metals --remove None > 1bqo_metals_test_None.log
#remwat
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_remwat_test_None.pdb --non_interactive remwat --remove No > 2ki5_remwat_test_None.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_remwat_test_Yes.pdb --non_interactive remwat --remove Yes > 2ki5_remwat_test_Yes.log

#ss bonds
python3 $APPDIR/checkStruc.py -i pdb:4ku1 --check_only --non_interactive getss > 4ku1_getss_test_log

#all check
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_all_test.pdb --json 1ark_all_test.json --non_interactive command_list --list scripts/all_checks > 1ark_all_test.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_all_test.pdb --json 2ki5_all_test.json --non_interactive command_list --list scripts/all_checks > 2ki5_all_test.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_all_test.pdb --json 1bqo_all_test.json --non_interactive command_list --list scripts/all_checks > 1bqo_all_test.log
python3 $APPDIR/checkStruc.py -i pdb:4ku1 -o 4ku1_all_test.pdb --json 4ku1_all_test.json --non_interactive command_list --list scripts/all_checks > 4ku1_all_test.log
