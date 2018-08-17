"""
    Main Class for Structure Checking functionality

"""

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"


import argparse
import re
import sys

from   structure_checking.help_manager import HelpManager
from   structure_checking.json_writer import JSONWriter
import structure_manager.data_constants as dataCts
from   structure_manager.structure_manager import StructureManager
import structure_manager.model_utils as mu

class StructureChecking():
    def __init__(self, args):
        self.args = args
        self.summary = {}
        self.rr_dist = []

    def launch(self, sets):
        help = HelpManager(sets.help_dir_path)
        help.print_help('header')

        try:
            f = getattr(self, self.args.command)
        except AttributeError:
            print ('Error: {} command unknown or not implemented'.format(self.args.command))
            sys.exit(1)

        f(self.args.options)

        if not self.args.check_only:
            self._save_structure()
            self.stm.set_num_ats()
            print ('Structure saved on {}'.format(self.args.output_structure_path))
            self.stm.print_stats('Final')
            self.summary['final_stats'] = self.stm.get_stats()

        if self.args.json_output_path is not None:
            json_writer = JSONWriter()
            for k in self.summary:
                json_writer.set(k, self.summary[k])
            json_writer.save(self.args.json_output_path)
            print ('Summary data saved on {}'.format(self.args.json_output_path))

    def command_list(self, options):
        opts = _get_parameters(options, "command_list", "--list", "op_list", "Command List File")
        opts.op_list = _check_parameter(opts.op_list, "Command list file: ")

        try:
            fh = open(opts.op_list, "r")

        except OSError:
            print ('Error when opening file {}'.format(opts.op_list))
            sys.exit(1)

        print ('Running command_list from {}'.format(opts.op_list))

        self._load_structure()

        i = 1
        for line in fh:
            if line == "\n" or line[0:1] == '#':
                continue
            print ("\nStep {}: {}".format(i, line))
            data = line.split()
            command = data[0]
            options = data[1:]
            try:
                f = getattr(self, command)
            except AttributeError:
                print ("Error: command unknown or not implemented")
                continue
                #sys.exit(1)
            f(options)
            i += 1

        print ("Command list completed")

    def models(self, options):
        opts = _get_parameters(options, "models", "--select_model", "select_model", "Select model to keep", int)
        print ('Running models. Options: {}'.format(' '.join(options)))
        models_sum = {}

        self._load_structure()


        print ('{} Model(s) detected'.format(self.stm.nmodels))
        models_sum['detected'] = {'nmodels': self.stm.nmodels}
        if self.stm.nmodels > 1:
            models_sum['detected']['type'] =self.stm.models_type
            if self.stm.models_type['type'] == mu.ENSM:
                print ('Models superimpose, RMSd: {:8.3f} A, guessed as ensemble type (NMR / MD TRAJ)'.format(self.stm.models_type['rmsd']))
            elif self.stm.models_type['type'] == mu.BUNIT:
                print ('Models do not superimpose, RMSd: {:8.3f} A, guessed as Biounit type'.format(self.stm.models_type['rmsd']))
            else:
                print ('Models type unknown')
        if self.args.check_only:
            print ('Running with --check_only. Nothing else to do.')
        else:
            if self.stm.nmodels > 1:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.select_model = _check_parameter(opts.select_model, "Select Model Num [1-" + str(self.stm.nmodels) + "]: ")
                    opts.select_model = int(opts.select_model)
                    ok = opts.select_model > 0 or opts.select_model <= self.stm.nmodels
                    if not ok:
                        print ('Error: unknown model {}'.format(opts.select_model), file=sys.stderr)
                        opts.select_model = ''
                        if self.args.non_interactive:
                            self.summary['models'] = models_sum
                            return 1

                print ('Selecting model num. {}'.format(opts.select_model))
                self.stm.select_model(opts.select_model)
                models_sum['selected_model'] = opts.select_model

            else:
                print ("Nothing to do")

        self.summary['models'] = models_sum

    def chains(self, options):
        opts = _get_parameters(options, "chains", '--select_chains', "select_chains", "Chains (All | Chain list comma separated)")
        print ('Running chains. Options: {}'.format(' '.join(options)))
        chains_sum = {}

        self._load_structure()

        print ('{} Chains detected'.format(len(self.stm.chain_ids)))
        for ch_id in sorted(self.stm.chain_ids):
            print ('  {}: {}'.format(ch_id, mu.chain_type_labels[self.stm.chain_ids[ch_id]]))
        chains_sum['detected'] = self.stm.chain_ids

        if self.args.check_only:
            print ('Running with --check_only. Nothing else to do.')
        else:
            if len(self.stm.chain_ids) > 1:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.select_chains = _check_parameter(opts.select_chains, 'Select All or chain id(s), [{}]: '.format(','.join(sorted(self.stm.chain_ids))))
                    ok = opts.select_chains.lower() == 'all'
                    if not ok:
                        ok = True
                        for ch in opts.select_chains.split(','):
                            ok = ok and ch in self.stm.chain_ids
                    if not ok:
                        print ('Error unknown selection: {}'.format(opts.select_chains))
                        if self.args.non_interactive:
                            self.summary['chains'] = chains_sum
                            return 1
                        opts.select_chains = ''

                if opts.select_chains.lower() == 'all':
                    print ('Selecting all chains')
                    chains_sum['selected'] = self.stm.chain_ids
                else:
                    self.stm.select_chains(opts.select_chains)
                    print ('Selecting chain(s) {}'.format(','.join(opts.select_chains)))
                    chains_sum['selected'] = opts.select_chains.split(',')
                self.stm.set_chain_ids()
                chains_sum['final'] = self.stm.chain_ids
            else:
                print ("Nothing to do")

        self.summary['chains'] = chains_sum

    def altloc(self, options):
        opts = _get_parameters(options, "altloc", '--select_altloc', "select_altloc", "select altloc occupancy|alt_id")

        print ('Running altloc. Options: {}'.format(' '.join(options)))
        altloc_sum = {}

        self._load_structure()

        alt_loc_res = mu.get_altloc_residues(self._get_structure())

        if len(alt_loc_res) > 0:
            print ('Detected {} residues with alternative location labels'.format(len(alt_loc_res)))

            altloc_sum['detected'] = {}

            for r in sorted(alt_loc_res):
                if len(alt_loc_res[r]) > 1:
                    print ('{}:'.format(r))
                    altloc_sum['detected'][r] = []
                    altlocs = sorted(alt_loc_res[r][0].child_dict)
                    for at in alt_loc_res[r]:
                        s = '  {:4}'.format(at.id)
                        for alt in sorted(at.child_dict):
                            s += ' {} ({:4.2f})'.format(alt, at.child_dict[alt].occupancy)
                        print (s)
                        altloc_sum['detected'][r].append({
                                                         'atom':at.id,
                                                         'loc_label':alt,
                                                         'occupancy':at.child_dict[alt].occupancy}
                                                         )

                    if not self.args.check_only:
                        ok = opts.select_altloc in altlocs or opts.select_altloc == 'occupancy'
                        while not ok:
                            if not self.args.non_interactive:
                                opts.select_altloc = _check_parameter(opts.select_altloc, 'Select alternative (occupancy, {}): '.format(','.join(altlocs)))
                            ok = opts.select_altloc in altlocs or opts.select_altloc.lower() == 'occupancy'
                            if not ok:
                                print ('Error: Unknown selection {} '.format(opts.select_altloc), file=sys.stderr)
                                opts.select_altloc = ''
                                if self.args.non_interactive:
                                    self.summary['altloc'] = altloc_sum
                                    return 1

                        print ('Selecting location {}'.format(opts.select_altloc))
                        altloc_sum['selected'] = opts.select_altloc
                        self.stm.select_altloc_residues(r, opts.select_altloc)
        else:
            print ("No residues with alternative location labels detected")

        self.summary['altloc'] = altloc_sum

    def metals (self, options):
        opts = _get_parameters(options, "metals", '--remove', 'remove_metals', 'Remove Metal ions')

        print ('Running metals. Options: {}'.format(' '.join(options)))
        metals_sum = {}

        self._load_structure()

        met_list = mu.get_metal_atoms(self._get_structure(), dataCts.metal_ats)

        if len(met_list) > 1:
            print ('{} Metal ions found'.format(len(met_list)))
            metals_sum['detected'] = []
            met_rids = []
            at_groups = {}
            for at in sorted(met_list, key=lambda x: x.serial_number):
                print ("  ", mu.atom_id(at, self.stm.nmodels > 1))
                r = at.get_parent()
                met_rids.append(mu.residue_num(r))
                if not at.id in at_groups:
                    at_groups[at.id] = []
                at_groups[at.id].append(at)
                metals_sum['detected'].append(mu.residue_num(r))

            if self.args.check_only:
                print ('Running with --check_only. Nothing else to do.')
            else:
                ok = False
                resids = False
                atids = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_metals = _check_parameter(opts.remove_metals,
                            'Remove (All | None | any of {} | any of {}):'.format(','.join(sorted(at_groups)), ','.join(met_rids)))
                    ok = opts.remove_metals.lower() in ['all', 'none']
                    if not ok:
                        ok = True
                        for r in opts.remove_metals.split(','):
                            ok = ok and r in met_rids
                        resids = ok

                    if not ok:
                        ok = True
                        for atid in opts.remove_metals.split(','):
                            ok = ok and atid in at_groups
                            if not ok and atid.upper() in at_groups:
                                print ('Warning: Atom names are case sensitive') # TODO check all case options with metals list
                        atids = ok

                    if not ok:
                        sys.stderr.write ('Error: unknown selection {}'.format(opts.remove_metals))
                        opts.remove_metals = ''
                        if self.args.non_interactive:
                            self.summary['metals'] = metals_sum
                            return 1

                if opts.remove_metals.lower() == 'none':
                    to_remove = []

                elif opts.remove_metals.lower() == 'all':
                    to_remove = met_list

                elif resids:
                    to_remove = []
                    rid_list = opts.remove_metals.split(',')
                    for at in met_list:
                        r = at.get_parent()
                        if mu.residue_num(r) in rid_list:
                            to_remove.append(at)
                elif atids:
                    to_remove = []
                    for atid in opts.remove_metals.split(','):
                        to_remove.extend(at_groups[atid])
                metals_sum['removed'] = []

                n = 0
                for at in to_remove:
                    metals_sum['removed'].append(mu.residue_id(at.get_parent(), self.stm.nmodels > 1))
                    self.stm.remove_residue(at.get_parent())
                    n += 1

                print ('Metal Atoms removed {} ({:d})'.format(opts.remove_metals, n))
                metals_sum['n_removed'] = n

        else:
            print ("No metal ions found")

        self.summary['metals'] = metals_sum

    def remwat(self, options):
        opts = _get_parameters(options, "remwat", '--remove', 'remove_wat', 'Remove Water molecules')

        print ('Running remwat. Options: {}'.format(' '.join(options)))
        remwat_sum = {}

        self._load_structure()

        lig_list = mu.get_ligands(self._get_structure(),incl_water=True)

        wat_list = []
        for r in lig_list:
            if mu.is_wat(r):
                wat_list.append(r)

        if len(wat_list) > 0:
            print ('{} Water molecules detected'.format(len(wat_list)))
            remwat_sum['n_detected'] = len(wat_list)

            if self.args.check_only:
                print ('Running with --check_only. Nothing else to do.')
            else:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_wat = _check_parameter(opts.remove_wat, 'Remove (Yes | No): ')
                    ok = opts.remove_wat.lower() in ['yes', 'no']
                    if not ok:
                        print ('Warning: unknown option {}'.format(opts.remove_wat))
                        if self.args.non_interactive:
                            self.summary['remwat'] = remwat_sum
                            return 1

                if opts.remove_wat.lower() == 'yes':
                    n = 0
                    for r in wat_list:
                        self.stm.remove_residue(r)
                        n += 1
                    print ('{} Water molecules removed'.format(n))
                    remwat_sum['n_removed'] = n
        else:
            print ("No water molecules detected")

        self.summary['remwat'] = remwat_sum


    def ligands(self, options):
        opts = _get_parameters(options, "ligands", '--remove', 'remove_ligands', 'Remove Ligand residues')

        print ('Running ligands. Options: {}'.format(' '.join(options)))
        ligands_sum = {}

        self._load_structure()

        lig_list = mu.get_ligands(self._get_structure())

        if len(lig_list) > 0:
            print ('{} Ligands detected '.format(len(lig_list)))
            rids = set()
            rnums = []
            ligands_sum['detected'] = []

            for r in sorted(lig_list, key=lambda x: x.index):
                print (mu.residue_id(r, self.stm.nmodels > 1))
                ligands_sum['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))
                rids.add(r.get_resname())
                rnums.append(mu.residue_num(r))

            if self.args.check_only:
                print ('Running with --check_only. Nothing else to do.')
            else:
                ok = False
                byresnum = False
                byrids = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_ligands = _check_parameter(opts.remove_ligands,
                            'Remove (All | None | any of {} | any of {}): '.format(','.join(sorted(rids)), ','.join(rnums)))

                    ok = opts.remove_ligands.lower() in ['all', 'none']
                    if not ok:
                        ok = True
                        for rid in opts.remove_ligands.split(','):
                            ok = ok and (rid.upper() in rids)
                        byrids = ok
                        if ok:
                            opts.remove_ligands = opts.remove_ligands.upper()
                    if not ok:
                        ok = True
                        for rn in opts.remove_ligands.split(','):
                            ok = ok and (rn in rnums)
                        byresnum = ok
                    if not ok:
                        print ('Error: unknown selection {}'.format(opts.remove_ligands))
                        if self.args.non_interactive:
                            self.summary['ligands'] = ligands_sum
                            return 1
                        opts.remove_ligands = ''

                    ligands_sum['removed'] = {'opt':opts.remove_ligands, 'lst':[]}

                to_remove = []

                if opts.remove_ligands.lower() == 'none':
                    print ("Nothing to do")
                elif opts.remove_ligands.lower() == 'all':
                    to_remove = lig_list
                elif byrids:
                    rm = opts.remove_ligands.split(',')
                    for r in lig_list:
                        if r.get_resname() in rm:
                            to_remove.append(r)
                elif byresnum:
                    rm = opts.remove_ligands.split(',')
                    for r in lig_list:
                        if mu.residue_num(r) in rm:
                            to_remove.append(r)
                n = 0
                for r in to_remove:
                    ligands_sum['removed']['lst'].append(mu.residue_id(r, self.stm.nmodels > 1))
                    self.stm.remove_residue(r)
                    n += 1

                print ('Ligands removed {} ({})'.format(opts.remove_ligands, n))
                ligands_sum['n_removed'] = n
        else:
            print ("No ligands detected")

        self.summary['ligands'] = ligands_sum

    def remh(self, options):
        opts = _get_parameters(options, "remh", '--remove', 'remove_h', 'Remove Hydrogen atoms')

        print ('Running remh. Options: {}'.format(' '.join(options)))
        remh_sum = {}

        self._load_structure()

        remh_list = mu.get_residues_with_H(self._get_structure())

        if len(remh_list) > 0:
            print ('{} Residues containing H atoms detected'.format(len(remh_list)))
            remh_sum['n_detected'] = len(remh_list)

            if self.args.check_only:
                print ('Running with --check_only. Nothing else to do.')
            else:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_h = _check_parameter(opts.remove_h, 'Remove hydrogen atoms (Yes | No): ')
                    ok = opts.remove_h.lower() in ['yes', 'no']
                    if not ok:
                        print ('Warning: unknown option {}'.format(opts.remove_h))
                        if self.args.non_interactive:
                            self.summary['remh'] = remh_sum
                            return 1

                if opts.remove_h.lower() == 'yes':
                    n = 0
                    for r in remh_list:
                        mu.remove_H_from_r(r['r'])
                        n += 1
                    print ('Hydrogen atom removed from {} residues'.format(n))
                    #self.stm.set_num_ats()
                    #self.stm.atom_renumbering()
                    remh_sum['n_removed'] = n
        else:
            print ("No residues with hydrogen atoms detected")

        self.summary['remh'] = remh_sum

    def getss (self, options):
        print ("Running getss")
        getss_sum = {}

        self._load_structure()

        SS_bonds = mu.get_all_at2at_distances(self._get_structure(),'SG', dataCts.SS_DIST)

        if len(SS_bonds):
            print ('{} Possible SS Bonds detected'.format(len(SS_bonds)))
            getss_sum['detected'] = []
            for ssb in SS_bonds:
                print ('  {} {}{:8.3f}'.format(mu.atom_id(ssb[0], self.stm.nmodels > 1), mu.atom_id(ssb[1], self.stm.nmodels > 1), ssb[2]))
                getss_sum['detected'].append({'at1':mu.atom_id(ssb[0], self.stm.nmodels > 1), 'at2':mu.atom_id(ssb[1], self.stm.nmodels > 1), 'dist': float(ssb[2])})

        else:
            print ("No SS bonds detected")

        self.summary['getss'] = getss_sum

    def amide(self, options):
        opts = _get_parameters(options, "amide", '--fix', 'amide_fix', 'Fix Residues (All | None | List)')
        print ('Running amide. Options: {}'.format(' '.join(options)))
        amide_sum = {}

        self._load_structure()
        amide_list = []
        for r in self._get_structure().get_residues():
            if r.get_resname() in dataCts.amide_res:
                amide_list.append(r)

        amide_sum['n_amides'] = len(amide_list)

        if len(amide_list) > 0:
            if len(self.rr_dist) == 0:
                self.rr_dist = mu.get_all_r2r_distances(self._get_structure(),'all', dataCts.R_R_CUTOFF)
            res_to_fix = []
            cont_list = []
            rnums = []
            for r_pair in self.rr_dist:
                [r1, r2, d] = r_pair
                if r1 in amide_list or r2 in amide_list:
                    if r1 != r2 and mu.same_model(r1, r2) and not mu.is_wat(r1) and not mu.is_wat(r2):
                        for at_pair in mu.get_all_rr_distances(r1, r2):
                            [at1, at2, dist] = at_pair
                            r1 = at1.get_parent()
                            ch1 = r1.get_parent()
                            r2 = at2.get_parent()
                            ch2 = r2.get_parent()
                            if not at1.id in dataCts.amide_atoms and not at2.id in dataCts.amide_atoms:
                                continue
                            for cls in ['acceptor', 'donor']:
                                if dist < dataCts.CLASH_DIST[cls]:
                                    if cls == 'donor' and not \
                                        (  self.stm.is_at_in_list(at1, dataCts.polar_donor) \
                                       and self.stm.is_at_in_list(at2, dataCts.polar_donor)):
                                            continue
                                    if cls == 'acceptor' and not \
                                        (  self.stm.is_at_in_list(at1, dataCts.polar_acceptor) \
                                       and self.stm.is_at_in_list(at2, dataCts.polar_acceptor)):
                                            continue
                                    if at1.id in dataCts.amide_atoms and r1 in amide_list:
                                        res_to_fix.append(r1)
                                        rnums.append(mu.residue_num(r1))
                                    if at2.id in dataCts.amide_atoms and r2 in amide_list:
                                        res_to_fix.append(r2)
                                        rnums.append(mu.residue_num(r2))
                                    cont_list.append(at_pair)
            print (res_to_fix)
            if len(cont_list):
                print ('{} unusual contact(s) involving amide atoms found'.format(len(cont_list)))
                amide_sum['detected'] = []
                for at_pair in cont_list:
                    print (' {:12} {:12} {:8.3f}'.format(mu.atom_id(at_pair[0], self.stm.nmodels > 1), mu.atom_id(at_pair[1], self.stm.nmodels > 1), at_pair[2]))
                    amide_sum['detected'].append({'at1':mu.atom_id(at_pair[0], self.stm.nmodels > 1), 'at2':mu.atom_id(at_pair[1], self.stm.nmodels > 1), 'dist': float(at_pair[2])})

                if self.args.check_only:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    ok = False
                    while not ok:
                        if not self.args.non_interactive:
                            opts.amide_fix = _check_parameter(opts.amide_fix, 'Fix amides (All | None | {}): '.format(','.join(rnums)))
                        ok = opts.amide_fix.lower() in ['all', 'none']
                        if not ok:
                            ok = True
                            for rn in opts.amide_fix.split(','):
                                ok = ok and rn in rnums
                        if not ok:
                            print ('Warning: unknown option {}'.format(opts.amide_fix))
                            if self.args.non_interactive:
                                self.summary['amide'] = amide_sum
                                return 1
                    if opts.amide_fix.lower() == 'none':
                        to_fix = []
                    elif opts.amide_fix.lower() == 'all':
                        to_fix = res_to_fix
                    else:
                        to_fix = []
                        for r in res_to_fix:
                            if mu.residue_num(r) in opts.amide_fix.split(','):
                                to_fix.append(r)
                    n = 0
                    for r in to_fix:
                        mu.invert_side_atoms(r, dataCts.amide_res)
                        n += 1
                    print ('Amide residues fixed {} ({})'.format(opts.amide_fix, n))

        else:
            print ("No amide residues found")

        self.summary['amide'] = amide_sum

    def chiral(self, options):
        opts = _get_parameters(options, "chiral", '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
        print ('Running chiral. Options: {}'.format(' '.join(options)))
        chiral_sum = {}

        self._load_structure()
        chiral_list = []
        for r in self._get_structure().get_residues():
            if r.get_resname() in dataCts.chiral_res:
                chiral_list.append(r)

        chiral_sum['n_chirals'] = len(chiral_list)
 
        if len(chiral_list) > 0:
            res_to_fix = []
            rnums = []
            for r in chiral_list:
                if not mu.check_chiral_residue(r, dataCts.chiral_res):
                    res_to_fix.append(r)
                    rnums.append(mu.residue_num(r))
            if len(res_to_fix):
                print ('{} residues with incorrect chirality found'.format(len(res_to_fix)))
                chiral_sum['detected'] = []
                for r in res_to_fix:
                    print (' {:10}'.format(mu.residue_id(r, self.stm.nmodels > 1)))
                    chiral_sum['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))

                if self.args.check_only:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    ok = False
                    while not ok:
                        if not self.args.non_interactive:
                            opts.chiral_fix = _check_parameter(opts.chiral_fix, 'Fix chiralities (All | None | {}): '.format(','.join(rnums)))
                        ok = opts.chiral_fix.lower() in ['all', 'none']
                        if not ok:
                            ok = True
                            for rn in opts.chiral_fix.split(','):
                                ok = ok and rn in rnums
                        if not ok:
                            print ('Warning: unknown option {}'.format(opts.amide_fix))
                            if self.args.non_interactive:
                                self.summary['chiral'] = chiral_sum
                                return 1
                    if opts.chiral_fix.lower() == 'none':
                        to_fix = []
                    elif opts.chiral_fix.lower() == 'all':
                        to_fix = res_to_fix
                    else:
                        to_fix = []
                        for r in res_to_fix:
                            if mu.residue_num(r) in opts.chiral_fix.split(','):
                                to_fix.append(r)
                    n = 0
                    for r in to_fix:
                        mu.invert_side_atoms(r, dataCts.chiral_res)
                        n += 1
                    print ('Quiral residues fixed {} ({})'.format(opts.chiral_fix, n))
            else:
                print ("No residues with incorrect chirality found")
        else:
            print ("No chiral residues found")

        self.summary['chiral'] = chiral_sum

    def chiral_bck(self, options):
        opts = _get_parameters(options, "chiral_bck", '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
        print ('Running chiral. Options: {}'.format(' '.join(options)))
        chiral_bck_sum = {}
        
        self.args.check_only=True # Provisional
        
        self._load_structure()
        chiral_list = []
        for ch in self._get_structure().get_chains():
            if self.stm.chain_ids[ch.id] == mu.PROTEIN:
                for r in ch.get_residues():
                    if r.get_resname() != 'GLY' and not mu.is_hetatm(r):
                        chiral_list.append(r)
        chiral_bck_sum['n_chirals'] = len(chiral_list)
 
        if len(chiral_list) > 0:
            res_to_fix = []
            rnums = []
            for r in chiral_list:
                if not mu.check_chiral_ca(r):
                    res_to_fix.append(r)
                    rnums.append(mu.residue_num(r))
            if len(res_to_fix):
                print ('{} residues with incorrect backbone chirality found'.format(len(res_to_fix)))
                chiral_bck_sum['detected'] = []
                for r in res_to_fix:
                    print (' {:10}'.format(mu.residue_id(r, self.stm.nmodels > 1)))
                    chiral_bck_sum['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))

                if self.args.check_only:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    ok = False
                    while not ok:
                        if not self.args.non_interactive:
                            opts.chiral_fix = _check_parameter(opts.chiral_fix, 'Fix CA chiralities (All | None | {}): '.format(','.join(rnums)))
                        ok = opts.chiral_fix.lower() in ['all', 'none']
                        if not ok:
                            ok = True
                            for rn in opts.chiral_fix.split(','):
                                ok = ok and rn in rnums
                        if not ok:
                            print ('Warning: unknown option {}'.format(opts.amide_fix))
                            if self.args.non_interactive:
                                self.summary['chiral'] = chiral_sum
                                return 1
                    if opts.chiral_fix.lower() == 'none':
                        to_fix = []
                    elif opts.chiral_fix.lower() == 'all':
                        to_fix = res_to_fix
                    else:
                        to_fix = []
                        for r in res_to_fix:
                            if mu.residue_num(r) in opts.chiral_fix.split(','):
                                to_fix.append(r)
                    n = 0
                    for r in to_fix:
                        mu.stm.invert_chiral_CA(r) # TODO
                        n += 1
                    print ('Quiral residues fixed {} ({})'.format(opts.chiral_fix, n))
            else:
                print ("No residues with incorrect backbone chirality found")
        else:
            print ("No chiral residues found")

        self.summary['chiral_bck'] = chiral_bck_sum

    def clashes(self, options):
        opts = _get_parameters(options, "clashes", '--no_wat', 'discard_wat', 'Discard water molecules')
        if opts.discard_wat is None:
            opts.discard_wat = True
        print ('Running clashes: Options {}'.format(' '.join(options)))
        clashes_sum = {}
        self._load_structure()

        if len(self.rr_dist) == 0:
            self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', dataCts.R_R_CUTOFF)

        clashes = {
            'severe':{},
            'apolar':{},
            'acceptor':{},
            'donor':{},
            'positive':{},
            'negative':{}
        }

        clashes_sum = {
            'detected': {
                'severe':[],
                'apolar':[],
                'acceptor':[],
                'donor':[],
                'positive':[],
                'negative':[]
            }
        }

        for r_pair in self.rr_dist:
            [r1, r2, d] = r_pair
            if opts.discard_wat and (mu.is_wat(r1) or mu.is_wat(r2)):
                continue
            if r1 != r2 and not mu.seq_consecutive(r1, r2) and mu.same_model(r1, r2):
                rkey = mu.residue_id(r1) + '-' + mu.residue_id(r2)
                for at_pair in mu.get_all_rr_distances(r1, r2):
                    [at1, at2, dist] = at_pair
                    for cls in ['severe', 'apolar', 'acceptor', 'donor', 'positive', 'negative']:
                        if dist < dataCts.CLASH_DIST[cls]:
                            if cls == 'apolar' and (at1.element not in dataCts.apolar_elements and at2.element not in dataCts.apolar_elements):
                                continue
                            if cls == 'donor' and not \
                                (   self.stm.is_at_in_list(at1, dataCts.polar_donor) \
                                and self.stm.is_at_in_list(at2, dataCts.polar_donor)):
                                continue
                            if cls == 'acceptor' and not \
                                (   self.stm.is_at_in_list(at1, dataCts.polar_acceptor) \
                                and self.stm.is_at_in_list(at2, dataCts.polar_acceptor)):
                                continue
                            if cls == 'positive' and not (at1.id in dataCts.pos_ats and at2.id in dataCts.pos_ats):
                                continue
                            if cls == 'negative' and not (at1.id in dataCts.neg_ats and at2.id in dataCts.neg_ats):
                                continue
                            if not rkey in clashes[cls]:
                                clashes[cls][rkey] = at_pair
                            if dist < clashes[cls][rkey][2]:
                                clashes[cls][rkey] = at_pair
        for cls in ['severe', 'apolar', 'acceptor', 'donor', 'positive', 'negative']:
            if len(clashes[cls]):
                print ('{} Steric {} clashes detected'.format(len(clashes[cls]), cls))
                for rkey in sorted(clashes[cls], key=lambda x: 10000 * clashes[cls][x][0].serial_number + clashes[cls][x][1].serial_number):
                    print (' {:12} {:12} {:8.3f}'.format(mu.atom_id(clashes[cls][rkey][0], self.stm.nmodels>1), mu.atom_id(clashes[cls][rkey][1], self.stm.nmodels>1), clashes[cls][rkey][2]))
                    at_pair = clashes[cls][rkey]
                    clashes_sum['detected'][cls].append(
                                                        {
                                                        'at1':mu.atom_id(clashes[cls][rkey][0], self.stm.nmodels>1),
                                                        'at2':mu.atom_id(clashes[cls][rkey][1], self.stm.nmodels>1),
                                                        'dist': float(clashes[cls][rkey][2])
                                                        }
                                                        )
            else:
                print ('No {} clashes detected'.format(cls))

        self.summary['clashes'] = clashes_sum

#===============================================================================

    def _load_structure(self, verbose=True):
        if not hasattr(self, 'stm'):
            if not self.args.non_interactive and self.args.input_structure_path is None:
                self.args.input_structure_path = input("Enter input structure path (PDB, mmcif | pdb:pdbid): ")
            self.stm = StructureManager(self.args.input_structure_path, self.args.debug)
            if verbose:
                print ('Structure {} loaded'.format(self.args.input_structure_path))
                self.stm.print_stats()
            self.summary['stats']=self.stm.get_stats()

    def _get_structure(self):
        return self.stm.get_structure()

    def _save_structure(self):
        if not self.args.non_interactive and self.args.output_structure_path is None:
            self.args.output_structure_path = input("Enter output structure path: ")
        self.stm.save_structure(self.args.output_structure_path)

    def json(self):
        json_writer = JSONWriter()

#===============================================================================

def _get_parameters (options, this_prog, this_param, this_dest, this_help, this_type=str):
    optionsParser = argparse.ArgumentParser(prog=this_prog)
    optionsParser.add_argument(this_param, dest=this_dest, help=this_help, type=this_type)
    return optionsParser.parse_args(options)

def _check_parameter (opts_param, input_text):
    while opts_param is None or opts_param == '':
        opts_param = input (input_text)
    if opts_param is str:
        opts_param = opts_param.replace(' ', '')
    return opts_param
