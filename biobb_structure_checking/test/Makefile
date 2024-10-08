BASEDIR = ../..
APPBIN = ${BASEDIR}/bin/check_structure
TESTDIR = ${BASEDIR}/biobb_structure_checking/test
SHELL = /bin/csh

COMMON = setenv PYTHONPATH ${BASEDIR}; cd ${TESTDIR}

clean:
	rm *pdb *log *json *diff *cmip *pqr *pdbqt *fasta

update:
	cp *pdb *log *json *cmip *pqr *pdbqt *fasta ref


all_nomodeller: models chains inscodes altloc ligands metals water rem_hydrogen getss clashes amide chiral mutateside fixside chiral_bck \
  backbone backbone_caps cistransbck biounits sequences dnamut add_hydrogen all_checks

modeller: backbone_main

all_nodebug: all_nomodeller modeller

all: debug all_nodebug

debug:
	echo "Running all checks on 1ubq with --debug"
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq_all_test_debug.pdb --json 1ubq_all_test_debug.json --debug --non_interactive command_list --list scripts/all_checks > 1ubq_all_test_debug.log

load:
	@echo "Loading structures"
	${COMMON};${APPBIN} -v -i pdb:1ark load > load_cache.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1svc load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:104l load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1bqo load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:4ku1 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1ubq load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:4mdh load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1d23 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:6m0j load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1ldn.1 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1vtk.1 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:2ki5.1 load >> load_cache.log
	${COMMON};${APPBIN} -v -i pdb:1svc.1 load >> load_cache.log

models:
	@echo "Running models on 1ark"
	${COMMON};${APPBIN} -v -i pdb:1ark -o 1ark_models_5_test.pdb --non_interactive models --select 5 > 1ark_models_test_5.log
	${COMMON};${APPBIN} -v -i pdb:1ark -o 1ark_models_1_test.pdb --non_interactive models --select 1 > 1ark_models_test_1.log
	${COMMON};${APPBIN} -v -i pdb:1ark -o 1ark_models_sup_test.pdb --non_interactive models --select 3-6 --superimpose --save_split > 1ark_models_3-6_sup_test.log
	@echo "Running models build complex on 1vtk"
	${COMMON};${APPBIN} -v -i pdb:1vtk.1 --file_format pdb -o 1vtk_models_complex.pdb --non_interactive models  --select all --build_complex > 1vtk_1_models_complex.log

chains:
	@echo "Running chains on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_chains_test_all.pdb --non_interactive chains --select All > 2ki5_chains_test_all.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_chains_test_B.pdb --non_interactive chains --select B > 2ki5_chains_test_B.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_chains_test_renum_chain.pdb --non_interactive chains --select All --renumber A:=C:,B:=D:  > 2ki5_chains_test_renum_chain.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_chains_test_fix_chains.pdb --non_interactive chains --select All --renumber A:501=W:,A:1-1=A:376,A:4-4=A:377,W:=A:501,B:504=W:,B:2-3=B:376,W:=B:504  > 2ki5_chains_test_fix_chains.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_chains_rebuild.pdb --non_interactive --coords_only chains --select all --rebuild > 2ki5_chains_test_rebuild.log
	@echo "Running chains on 1svc"
	${COMMON};${APPBIN} -v -i pdb:1svc -o 1svc_chains_test_protein.pdb --non_interactive chains --select protein > 1svc_chains_test_protein.log
	${COMMON};${APPBIN} -v -i pdb:1svc -o 1svc_chains_rebuild.pdb --non_interactive chains --select all --rebuild > 1svc_chains_test_rebuild.log
inscodes:
	@echo "Running inscodes on 104l"
	${COMMON};${APPBIN} -v -i pdb:104l -o 104l_inscodes_test.pdb --non_interactive inscodes > 104l_inscodes_test.log
	${COMMON};${APPBIN} -v -i pdb:104l -o 104l_inscodes_renum_test.pdb --non_interactive inscodes --renumber > 104l_inscodes_renum_test.log

altloc:
	@echo "Running altloc on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_altloc_test_occ.pdb --non_interactive altloc --select occupancy > 2ki5_altloc_test_occ.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_altloc_test_B.pdb --non_interactive altloc --select B > 2ki5_altloc_test_B.log

ligands:
	@echo "Running ligands on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_ligands_test_SO4.pdb --non_interactive ligands --remove SO4 > 2ki5_ligands_test_SO4.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_ligands_test_All.pdb --non_interactive ligands --remove All > 2ki5_ligands_test_All.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_ligands_test_None.pdb --non_interactive ligands --remove None > 2ki5_ligands_test_None.log

metals:
	@echo "Running metals on 1bqo"
	${COMMON};${APPBIN} -v -i pdb:1bqo -o 1bqo_metals_test_ZN.pdb --non_interactive metals --remove ZN > 1bqo_metals_test_ZN.log
	${COMMON};${APPBIN} -v -i pdb:1bqo -o 1bqo_metals_test_A305.pdb --non_interactive metals --remove A305 > 1bqo_metals_test_A305.log
	${COMMON};${APPBIN} -v -i pdb:1bqo -o 1bqo_metals_test_All.pdb --non_interactive metals --remove All > 1bqo_metals_test_All.log
	${COMMON};${APPBIN} -v -i pdb:1bqo -o 1bqo_metals_test_None.pdb --non_interactive metals --remove None > 1bqo_metals_test_None.log

water:
	@echo "Running water on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 --non_interactive water --remove No > 2ki5_remwat_test_None.log
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_remwat_test_Yes.pdb --non_interactive water --remove Yes > 2ki5_remwat_test_Yes.log

rem_hydrogen:
	@echo "Running rem_hydrogen on 1ark"
	${COMMON};${APPBIN} -v -i pdb:1ark --non_interactive rem_hydrogen --remove No > 1ark_remh_test_None.log
	${COMMON};${APPBIN} -v -i pdb:1ark -o 1ark_remh_test_Yes.pdb --non_interactive rem_hydrogen --remove Yes > 1ark_remh_test_Yes.log

getss:
	@echo "Running getss on 4ku1"
	${COMMON};${APPBIN} -v -i pdb:4ku1 --check_only --non_interactive getss > 4ku1_getss_test.log
	${COMMON};${APPBIN} -v -i pdb:4ku1 -o 4ku1_getss_mark_test.pdb --non_interactive getss --mark all > 4ku1_getss_mark_test.log

clashes:
	@echo "Running clashes on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 --check_only --non_interactive clashes > 2ki5_clashes_test.log

amide:
	@echo "Running amide on 1ubq"
	${COMMON};${APPBIN} -v -i pdb:1ubq --non_interactive amide --fix None > 1ubq_amide_test_None.log
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq_amide_test_All.pdb --non_interactive amide --fix All > 1ubq_amide_test_All.log
	@echo "Running amide auto on 6m0j"
	${COMMON};${APPBIN} -v -i pdb:6m0j -o 6m0j_amide_test_auto.pdb --non_interactive amide --fix auto > 6m0j_amide_test_auto.log

chiral:
	@echo "Running chiral on modified 1ubq"
	${COMMON};${APPBIN} -v -i chiral_pdb_test/1ubq_chi.pdb -o 1ubq_chi_chiral_test.pdb --non_interactive chiral --fix All > 1ubq_chi_chiral_test_All.log
	${COMMON};${APPBIN} -v -i chiral_pdb_test/1ubq_chi.pdb --non_interactive chiral --fix None > 1ubq_chi_none_test_None.log

mutateside:
	@echo "Running mutateside on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5  -o 2ki5_mutateside_test.pdb --json 2ki5_mutateside_test.json  mutateside --mut Leu49Ile,B:arg51Lys > 2ki5_mutateside_test.log
	${COMMON};${APPBIN} -v -i pdb:2ki5  -o 2ki5_mutateside_rb_test.pdb --json 2ki5_mutateside_rb_test.json  mutateside --mut Leu49Ile --rebuild> 2ki5_mutateside_rb_test.log

fixside:
	@echo "Running fixside on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 fixside --fix None > 2ki5_fixside_None_test.log
	${COMMON};${APPBIN} -v -i pdb:2ki5  -o 2ki5_fixside_All_test.pdb fixside --fix All > 2ki5_fixside_All_test.log
	${COMMON};${APPBIN} -v -i pdb:2ki5  -o 2ki5_fixside_All_rb_test.pdb fixside --fix All --rebuild > 2ki5_fixside_All_rb_test.log
	@echo "Running fixside on 1d23_fixtest"
	${COMMON};${APPBIN} -v -i dna_test/1d23_fixtest.pdb -o 1d23_fixtest_All_test.pdb fixside --fix All > 1d23_fixtest_All_test.log

chiral_bck:
	@echo "Running chiral_bck on modified 1ark"
	${COMMON};${APPBIN} -v -i chiral_pdb_test/1ark_m1_chica_nh.pdb -o 1ark_chiral_bck_test.pdb --non_interactive chiral_bck > 1ark_chiral_bck_test.log

backbone:
	@echo "Running backbone on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5  -o 2ki5_backbone_test.pdb --non_interactive --check_only backbone > 2ki5_backbone_test.log

backbone_main:
	@echo "Running backbone_main on 1ubq_caps"
	${COMMON};${APPBIN} -v -i backbone_test/1ubq_caps.pdb -o 1ubq_backbone_main_test.pdb --sequence backbone_test/1ubq.fasta --non_interactive  backbone --no_recheck --fix_chain all --add_caps none --fix_atoms none > 1ubq_backbone_main_test.log

backbone_caps:
	@echo "Running backbone_caps on 1ubq_caps"
	${COMMON};${APPBIN} -v -i backbone_test/1ubq_caps.pdb -o 1ubq_backbone_caps_test.pdb --non_interactive  backbone --no_recheck --fix_chain none --add_caps all --fix_atoms all > 1ubq_backbone_caps_test.log

cistransbck:
	@echo "Running cistransbck on 4mdh"
	${COMMON};${APPBIN} -v -i pdb:4mdh -o 4mdh_cistransbck_test.pdb --non_interactive cistransbck > 4mdh_cistransbck_test.log

sequences:
	@echo "Running sequences on 2ki5"
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_sequences_test.pdb --non_interactive sequences --output_fasta 2ki5_sequences.fasta > 2ki5_sequences_test.log

dnamut:
	@echo "Running dnamut on 1D23"
	${COMMON};${APPBIN} -v -i pdb:1d23 -o 1d23_dnamut_test.pdb --non_interactive mutateside --mut DA3DG > 1d23_dnamut_test.log
	${COMMON};${APPBIN} -v -i pdb:1d23 -o 1d23_poliA_test.pdb --non_interactive mutateside --na_seq AAAAAAAAAA  > 1d23_poliA_test.log

biounits:
	@echo "Running checkall on biounits"
	@echo " 1ldn.1 "
	${COMMON};${APPBIN} -v -i pdb:1ldn.1 --non_interactive checkall > 1ldn1_checkall_test.log
	@echo " 2ki5.1 "
	${COMMON};${APPBIN} -v -i pdb:2ki5.1 --non_interactive checkall > 2ki5_checkall_test.log
	@echo " 1vtk.1 "
	${COMMON};${APPBIN} -v -i pdb:1vtk.1 --non_interactive checkall > 1vtk_checkall_test.log
	@echo " 1svc.1 "
	${COMMON};${APPBIN} -v -i pdb:1svc.1 --non_interactive checkall > 1svc_checkall_test.log

add_hydrogen:
	@echo "Running add_hydrogen on 1ubq"
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq.pdbqt --non_interactive command_list --list 'water --remove yes; add_hydrogen --add_mode auto --add_charges ADT' > 1ubq_add_hydrogen_ADT_test.log
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq.pqr   --non_interactive command_list --list 'water --remove yes; add_hydrogen --add_mode auto --add_charges ADT' > 1ubq_add_hydrogen_PQR_test.log
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq.cmip  --non_interactive command_list --list 'water --remove yes; add_hydrogen --add_mode auto --add_charges CMIP' > 1ubq_add_hydrogen_CMIP_test.log
	@echo "Running add_hydrogen on 1d23"
	${COMMON};${APPBIN} -v -i pdb:1d23 -o 1d23.pdb --output_format cmip  --non_interactive command_list --list 'water --remove yes; add_hydrogen --add_mode auto --add_charges CMIP' > 1d23_add_hydrogen_CMIP_test.log
	@echo "Running add_hydrogen on 1rna"
	${COMMON};${APPBIN} -v -i pdb:1rna -o 1rna.pdbqt --output_format pdbqt  --non_interactive command_list --list 'water --remove yes; add_hydrogen --add_mode auto --add_charges adt' > 1rna_add_hydrogen_ADT_test.log

all_checks:
	@echo "Running all checks"
	@echo " 1ark "
	${COMMON};${APPBIN} -v -i pdb:1ark -o 1ark_all_test.pdb --json 1ark_all_test.json --non_interactive command_list --list scripts/all_checks > 1ark_all_test.log
	@echo " 2ki5 "
	${COMMON};${APPBIN} -v -i pdb:2ki5 -o 2ki5_all_test.pdb --json 2ki5_all_test.json --non_interactive command_list --list scripts/all_checks > 2ki5_all_test.log
	@echo " 1bqo "
	${COMMON};${APPBIN} -v -i pdb:1bqo -o 1bqo_all_test.pdb --json 1bqo_all_test.json --non_interactive command_list --list scripts/all_checks > 1bqo_all_test.log
	@echo " 4kui "
	${COMMON};${APPBIN} -v -i pdb:4ku1 -o 4ku1_all_test.pdb --json 4ku1_all_test.json --non_interactive command_list --list scripts/all_checks > 4ku1_all_test.log
	@echo " 1ubq "
	${COMMON};${APPBIN} -v -i pdb:1ubq -o 1ubq_all_test.pdb --json 1ubq_all_test.json --non_interactive command_list --list scripts/all_checks > 1ubq_all_test.log
	@echo " 1svc "
	${COMMON};${APPBIN} -v -i pdb:1svc -o 1svc_all_test.pdb --json 1svc_all_test.json --non_interactive command_list --list scripts/all_checks > 1svc_all_test.log
	@echo " 1vtk.1 "
	${COMMON};${APPBIN} -v -i pdb:1vtk.1 -o 1vtk_1_all_test.pdb --json 1vtk_1_all_test.json --non_interactive command_list --list scripts/all_checks > 1vtk_1_svc_all_test.log

diffs:
	cd ${TESTDIR};scripts/calcDiffs.csh

