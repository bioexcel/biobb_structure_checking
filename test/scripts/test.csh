#!/bin/csh
setenv APPDIR `pwd`/..
rm *pdb *log *json *diff 
#models
echo "Running models on 1ark"
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_models_5_test.pdb --non_interactive models --select_model 5 > 1ark_models_test_5.log
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_models_1_test.pdb --non_interactive models --select_model 1 > 1ark_models_test_1.log
#chains
echo "Running chains on 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_chains_test_all.pdb --non_interactive chains --select_chain All > 2ki5_chains_test_all.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_chains_test_B.pdb --non_interactive chains --select_chain B > 2ki5_chains_test_B.log
#altloc
echo "Running altloc on 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_altloc_test_occ.pdb --non_interactive altloc --select_altloc occupancy > 2ki5_altloc_test_occ.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_altloc_test_B.pdb --non_interactive altloc --select_altloc B > 2ki5_altloc_test_B.log
#ligands
echo "Running ligands on 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_SO4.pdb --non_interactive ligands --remove SO4 > 2ki5_ligands_test_SO4.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_All.pdb --non_interactive ligands --remove All > 2ki5_ligands_test_All.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_ligands_test_None.pdb --non_interactive ligands --remove None > 2ki5_ligands_test_None.log
#metals
echo "Running metals on 1bqo"
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_ZN.pdb --non_interactive metals --remove ZN > 1bqo_metals_test_ZN.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_A305.pdb --non_interactive metals --remove A305 > 1bqo_metals_test_A305.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_All.pdb --non_interactive metals --remove All > 1bqo_metals_test_All.log
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_metals_test_None.pdb --non_interactive metals --remove None > 1bqo_metals_test_None.log
#remwat
echo "Running remwat on 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 --non_interactive remwat --remove No > 2ki5_remwat_test_None.log
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_remwat_test_Yes.pdb --non_interactive remwat --remove Yes > 2ki5_remwat_test_Yes.log
#remh
echo "Running remh on 1ark"
python3 $APPDIR/checkStruc.py -i pdb:1ark --non_interactive remh --remove No > 1ark_remh_test_None.log
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_remh_test_Yes.pdb --non_interactive remh --remove Yes > 1ark_remh_test_Yes.log
#ss bonds
echo "Running getss on 4ku1"
python3 $APPDIR/checkStruc.py -i pdb:4ku1 --check_only --non_interactive getss > 4ku1_getss_test.log
#clashes
echo "Running clashes on 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 --check_only --non_interactive clashes > 2ki5_clashes_test.log
#amide
echo "Running amide on 1ubq"
python3 $APPDIR/checkStruc.py -i pdb:1ubq --non_interactive amide --fix None > 1ubq_amide_test_None.log
python3 $APPDIR/checkStruc.py -i pdb:1ubq -o 1ubq_amide_test_All.pdb --non_interactive amide --fix All > 1ubq_amide_test_All.log
#chiral
echo "Running chiral on modified 1ubq"
python3 $APPDIR/checkStruc.py -i chiral_pdb_test/1ubq_chi.pdb -o 1ubq_chi_chiral_test.pdb --non_interactive chiral --fix All > 1ubq_chi_chiral_test_All.log
python3 $APPDIR/checkStruc.py -i chiral_pdb_test/1ubq_chi.pdb --non_interactive chiral --fix None > 1ubq_chi_none_test_None.log
#chiral_bck
echo "Running chiral_bck on modified 1ark"
python3 $APPDIR/checkStruc.py -i chiral_pdb_test/1ark_m1_chica_nh.pdb -o 1ark_chiral_bck_test.pdb --non_interactive chiral_bck > 1ark_chiral_bck_test.log
#All_test
echo -n "Running All checks on 1ark"
python3 $APPDIR/checkStruc.py -i pdb:1ark -o 1ark_all_test.pdb --json 1ark_all_test.json --non_interactive command_list --list scripts/all_checks > 1ark_all_test.log
echo -n " 2ki5"
python3 $APPDIR/checkStruc.py -i pdb:2ki5 -o 2ki5_all_test.pdb --json 2ki5_all_test.json --non_interactive command_list --list scripts/all_checks > 2ki5_all_test.log
echo -n " 1bqo"
python3 $APPDIR/checkStruc.py -i pdb:1bqo -o 1bqo_all_test.pdb --json 1bqo_all_test.json --non_interactive command_list --list scripts/all_checks > 1bqo_all_test.log
echo -n " 4kui"
python3 $APPDIR/checkStruc.py -i pdb:4ku1 -o 4ku1_all_test.pdb --json 4ku1_all_test.json --non_interactive command_list --list scripts/all_checks > 4ku1_all_test.log
echo -n " 1ubq"
python3 $APPDIR/checkStruc.py -i pdb:1ubq -o 1ubq_all_test.pdb --json 1ubq_all_test.json --non_interactive command_list --list scripts/all_checks > 1ubq_all_test.log
echo " 1svc"
python3 $APPDIR/checkStruc.py -i pdb:1svc -o :wq1svc_all_test.pdb --json 1svc_all_test.json --non_interactive command_list --list scripts/all_checks > 1svc_all_test.log
echo "Calculating diffs"
foreach f (*pdb *json *log)
diff $f ref/$f > $f.diff
end
