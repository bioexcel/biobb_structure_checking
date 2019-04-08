#!/usr/bin/env bash

#set -x
shopt -s nullglob


CLI="check_structure"

#Get bash script directory
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
REF=$(dirname "$DIR")/ref
CHI=$(dirname "$DIR")/chiral_pdb_test

#Remove all tmp files inside script directory
rm $DIR/*.pdb $DIR/*.log $DIR/*.json $DIR/*.diff &> /dev/null

#models
#echo "Running models on 1ark"
$CLI -i pdb:1ark -o $DIR/1ark_models_5_test.pdb --non_interactive models --select 5 &> $DIR/1ark_models_test_5.log
$CLI -i pdb:1ark -o $DIR/1ark_models_1_test.pdb --non_interactive models --select 1 &> $DIR/1ark_models_test_1.log
#chains
#echo "Running chains on 2ki5"
$CLI -i pdb:2ki5 -o $DIR/2ki5_chains_test_all.pdb --non_interactive chains --select All &> $DIR/2ki5_chains_test_all.log
$CLI -i pdb:2ki5 -o $DIR/2ki5_chains_test_B.pdb --non_interactive chains --select B &> $DIR/2ki5_chains_test_B.log
#inscodes
#echo "Running inscodes on 104l"
$CLI -i pdb:104l -o $DIR/104l_chains_test.pdb --non_interactive inscodes &> $DIR/104l_chains_test.log
#altloc
#echo "Running altloc on 2ki5"
$CLI -i pdb:2ki5 -o $DIR/2ki5_altloc_test_occ.pdb --non_interactive altloc --select occupancy &> $DIR/2ki5_altloc_test_occ.log
$CLI -i pdb:2ki5 -o $DIR/2ki5_altloc_test_B.pdb --non_interactive altloc --select B &> $DIR/2ki5_altloc_test_B.log
#ligands
#echo "Running ligands on 2ki5"
$CLI -i pdb:2ki5 -o $DIR/2ki5_ligands_test_SO4.pdb --non_interactive ligands --remove SO4 &> $DIR/2ki5_ligands_test_SO4.log
$CLI -i pdb:2ki5 -o $DIR/2ki5_ligands_test_All.pdb --non_interactive ligands --remove All &> $DIR/2ki5_ligands_test_All.log
$CLI -i pdb:2ki5 -o $DIR/2ki5_ligands_test_None.pdb --non_interactive ligands --remove None &> $DIR/2ki5_ligands_test_None.log
#metals
#echo "Running metals on 1bqo"
$CLI -i pdb:1bqo -o $DIR/1bqo_metals_test_ZN.pdb --non_interactive metals --remove ZN &> $DIR/1bqo_metals_test_ZN.log
$CLI -i pdb:1bqo -o $DIR/1bqo_metals_test_A305.pdb --non_interactive metals --remove A305 &> $DIR/1bqo_metals_test_A305.log
$CLI -i pdb:1bqo -o $DIR/1bqo_metals_test_All.pdb --non_interactive metals --remove All &> $DIR/1bqo_metals_test_All.log
$CLI -i pdb:1bqo -o $DIR/1bqo_metals_test_None.pdb --non_interactive metals --remove None &> $DIR/1bqo_metals_test_None.log
#remwat
#echo "Running remwat on 2ki5"
$CLI -i pdb:2ki5 --non_interactive water --remove No &> $DIR/2ki5_remwat_test_None.log
$CLI -i pdb:2ki5 -o $DIR/2ki5_remwat_test_Yes.pdb --non_interactive water --remove Yes &> $DIR/2ki5_remwat_test_Yes.log
#remh
#echo "Running remh on 1ark"
$CLI -i pdb:1ark --non_interactive rem_hydrogen --remove No &> $DIR/1ark_remh_test_None.log
$CLI -i pdb:1ark -o $DIR/1ark_remh_test_Yes.pdb --non_interactive rem_hydrogen --remove Yes &> $DIR/1ark_remh_test_Yes.log
#ss bonds
#echo "Running getss on 4ku1"
$CLI -i pdb:4ku1 --check_only --non_interactive getss &> $DIR/4ku1_getss_test.log
#clashes
#echo "Running clashes on 2ki5"
$CLI -i pdb:2ki5 --check_only --non_interactive clashes &> $DIR/2ki5_clashes_test.log
#amide
#echo "Running amide on 1ubq"
$CLI -i pdb:1ubq --non_interactive amide --fix None &> $DIR/1ubq_amide_test_None.log
$CLI -i pdb:1ubq -o $DIR/1ubq_amide_test_All.pdb --non_interactive amide --fix All &> $DIR/1ubq_amide_test_All.log
#chiral
#echo "Running chiral on modified 1ubq"
$CLI -i $CHI/1ubq_chi.pdb -o $DIR/1ubq_chi_chiral_test.pdb --non_interactive chiral --fix All &> $DIR/1ubq_chi_chiral_test_All.log
$CLI -i $CHI/1ubq_chi.pdb --non_interactive chiral --fix None &> $DIR/1ubq_chi_none_test_None.log
#mutateside
#echo "Running mutateside on 2ki5"
$CLI -i pdb:2ki5  -o $DIR/1ki5_mutateside_test.pdb mutateside --mut Leu49Ile,B:arg51Lys &> $DIR/2ki5_mutateside_test.log
#fixside
#echo "Running fixside on 2ki5"
$CLI -i pdb:2ki5 fixside --fix None &> $DIR/2ki5_fixside_None_test.log
$CLI -i pdb:2ki5  -o $DIR/1ki5_fixside_All_test.pdb fixside --fix All &> $DIR/2ki5_fixside_All_test.log
#chiral_bck
#echo "Running chiral_bck on modified 1ark"
$CLI -i $CHI/1ark_m1_chica_nh.pdb -o $DIR/1ark_chiral_bck_test.pdb --non_interactive chiral_bck &> $DIR/1ark_chiral_bck_test.log
#backbone
#echo "Running backbone on 2ki5"
$CLI -i pdb:2ki5  -o $DIR/1ki5_backbone_test.pdb --non_interactive --check_only backbone &> $DIR/2ki5_ibackbone_test.log
#cistransbck
#echo "Running cistransbck on 4mdh"
$CLI -i pdb:4mdh -o $DIR/4mdh_cistransbck_test.pdb --non_interactive cistransbck &> $DIR/4mdh_cistransbck_test.log
#echo "Running all checks on 1ldn.1"
$CLI -i pdb:1ldn.1 --pdb_server mmb  --non_interactive checkall &> $DIR/1ldn1_checkall_test.log
#echo "Running all checks on 2ki5.1"
$CLI -i pdb:2ki5.1 --pdb_server mmb  --non_interactive checkall &> $DIR/2ki5_checkall_test.log
#echo "Running all checks on 1vtk.1"
$CLI -i pdb:1vtk.1 --pdb_server mmb  --non_interactive checkall &> $DIR/1vtk_checkall_test.log
#All_test
#echo -n "Running all checks on 1ark"
$CLI -i pdb:1ark -o $DIR/1ark_all_test.pdb --json $DIR/1ark_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/1ark_all_test.log
#echo -n " 2ki5"
$CLI -i pdb:2ki5 -o $DIR/2ki5_all_test.pdb --json $DIR/2ki5_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/2ki5_all_test.log
#echo -n " 1bqo"
$CLI -i pdb:1bqo -o $DIR/1bqo_all_test.pdb --json $DIR/1bqo_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/1bqo_all_test.log
#echo -n " 4ku1"
$CLI -i pdb:4ku1 -o $DIR/4ku1_all_test.pdb --json $DIR/4ku1_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/4ku1_all_test.log
#echo -n " 1ubq"
$CLI -i pdb:1ubq -o $DIR/1ubq_all_test.pdb --json $DIR/1ubq_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/1ubq_all_test.log
#echo -n " 1svc"
$CLI -i pdb:1svc -o $DIR/1svc_all_test.pdb --json $DIR/1svc_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/1svc_all_test.log
#echo " 1vtk.1"
$CLI -i pdb:1vtk.1 --pdb_server mmb -o $DIR/1vtk_1_all_test.pdb --json $DIR/1vtk_1_all_test.json --non_interactive command_list --list $DIR/all_checks &> $DIR/1vtk_1_all_test.log
#echo "Calculating diffs"

for f in $DIR/*.pdb $DIR/*.json $DIR/*.log
do
  if [ -s $f ]; then
    cp $f $DIR/file_a
    cp $REF/$(basename $f) $DIR/file_b
    sed -i '/Structure exists/d' $DIR/file_a
    sed -i '/Structure exists/d' $DIR/file_b
    sed -i '/Structure saved on/d' $DIR/file_a
    sed -i '/Structure saved on/d' $DIR/file_b
    sed -i '/Downloading PDB structure/d' $DIR/file_a
    sed -i '/Downloading PDB structure/d' $DIR/file_b
    sed -i '/Structure .* loaded/d' $DIR/file_a
    sed -i '/Structure .* loaded/d' $DIR/file_b
    sed -i '/Summary data saved on/d' $DIR/file_a
    sed -i '/Summary data saved on/d' $DIR/file_b
    sed -i '/ERROR Unknown selection SO4/d' $DIR/file_a
    sed -i '/ERROR Unknown selection SO4/d' $DIR/file_b
    diff $DIR/file_a $DIR/file_b &> $DIR/$(basename $f).diff
    if [ -s $DIR/$(basename $f).diff ]
    then
      echo $DIR/$(basename $f)
      cat $DIR/$(basename $f).diff
    fi
    rm $DIR/file_a $DIR/file_b
  fi
done

#Remove all tmp files inside script directory
rm $DIR/*.pdb $DIR/*.log $DIR/*.json $DIR/*.diff &> /dev/null
rm -rf tmpPDB &> /dev/null
