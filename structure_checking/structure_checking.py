"""
    Main Class for Structure Checking functionality

"""

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"


import argparse
import re
import sys
from structure_checking.help_manager import HelpManager
from structure_checking.json_writer import JSONWriter
from structure_manager.data_lib_manager import DataLibManager
import structure_manager.model_utils as mu
from structure_manager.structure_manager import StructureManager
from structure_checking.param_input import ParamInput

class StructureChecking():
    def __init__(self, sets, args):
        self.args = args
        self.sets = sets
        self.summary = {}
        self.rr_dist = []
        self.data_library = DataLibManager(sets.data_library_path)
        if not 'Notebook' in self.args:
            self.args['Notebook']= False
        if self.args['Notebook']:
            self.args['non_interactive']=True
            self.args['check_only']=False


    def launch(self):
        help = HelpManager(self.sets.help_dir_path)
        help.print_help('header')

        try:
            f = getattr(self, self.args['command'])
        except AttributeError:
            print ('Error: {} command unknown or not implemented'.format(self.args['command']))
            sys.exit(1)

        f(self.args['options'])


        if not self.args['check_only']:
            if self.stm.modified or self.args['force_save']:
                if not self.stm.modified:
                    print ('Structure not modified, saving due to --force_save option')
                self._save_structure()
                self.stm.calc_stats()
                print ('Structure saved on {}'.format(self.args['output_structure_path']))
                self.stm.print_stats('Final')
                self.summary['final_stats'] = self.stm.get_stats()
            elif not self.stm.modified:
                print ('Structure not modified, not saving. Override with --force_save')

        if self.args['json_output_path'] is not None:
            json_writer = JSONWriter()
            for k in self.summary:
                json_writer.set(k, self.summary[k])
            json_writer.save(self.args['json_output_path'])
            print ('Summary data saved on {}'.format(self.args['json_output_path']))

           
    def print_stats(self, verbose=True):
        self._load_structure()
        self.stm.print_stats()
        
    def command_list(self, options):
        opts = _get_parameters(options, "command_list", "--list", "op_list", "Command List File")
        [input_option, opts.op_list] = ParamInput('Command List File', opts.op_list, False).run()

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
            models_sum['detected']['type'] = self.stm.models_type
            if self.stm.models_type['type'] == mu.ENSM:
                print ('Models superimpose, RMSd: {:8.3f} A, guessed as ensemble type (NMR / MD TRAJ)'.format(self.stm.models_type['rmsd']))
            elif self.stm.models_type['type'] == mu.BUNIT:
                print ('Models do not superimpose, RMSd: {:8.3f} A, guessed as Biounit type'.format(self.stm.models_type['rmsd']))
            else:
                print ('Models type unknown')
        if self.args['check_only']:
            print ('Running with --check_only. Nothing else to do.')
        else:
            if self.stm.nmodels > 1:
                input_line = ParamInput('Select Model Num', opts.select_model, self.args['non_interactive'])
                input_line.add_option_all ()
                input_line.add_option ('modelno', [], type='int', min=1, max=self.stm.nmodels)
                [input_option, opts.select_model] = input_line.run()
                
                if input_option == 'error':
                    print ('Error: unknown model {}'.format(opts.select_model), file=sys.stderr)
                    self.summary['models'] = models_sum
                    return 1

                print ('Selecting model num. {}'.format(opts.select_model))
                if input_option != 'all':
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

        if self.args['check_only']:
            print ('Running with --check_only. Nothing else to do.')
        else:
            if len(self.stm.chain_ids) > 1:
                input_line = ParamInput('Select chain', opts.select_chains, self.args['non_interactive'])
                input_line.add_option_all()
                input_line.add_option('chid',sorted(self.stm.chain_ids), multiple=True, case="sensitive")
                [input_option, opts.select_chains] = input_line.run()
                
                if input_option == 'error':
                    print ('Error unknown selection: {}'.format(opts.select_chains))
                    self.summary['chains'] = chains_sum
                    return 1

                if input_option == 'all':
                    print ('Selecting all chains')
                    chains_sum['selected'] = self.stm.chain_ids
                else:
                    self.stm.select_chains(opts.select_chains)
                    print ('Selecting chain(s) {}'.format(opts.select_chains))
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

                    if not self.args['check_only']:
                        input_line = ParamInput('Select alternative', opts.select_altloc, self.args['non_interactive'])
                        input_line.add_option_all()
                        input_line.add_option('occup',['occupancy'])
                        input_line.add_option('altids',altlocs, case='upper')
                        [input_option, opts.select_altloc] = input_line.run()
                        
                        if input_option=='error':
                            print ('Error: Unknown selection {} '.format(opts.select_altloc), file=sys.stderr)
                            self.summary['altloc'] = altloc_sum
                            return 1

                        print ('Selecting location {}'.format(opts.select_altloc))
                        altloc_sum['selected'] = opts.select_altloc
                        if input_option == 'all':
                            print ("Nothing to do")
                        else:
                            self.stm.select_altloc_residues(r, opts.select_altloc.upper())
        else:
            print ("No residues with alternative location labels detected")

        self.summary['altloc'] = altloc_sum

    def metals (self, options):
        opts = _get_parameters(options, "metals", '--remove', 'remove_metals', 'Remove Metal ions')

        print ('Running metals. Options: {}'.format(' '.join(options)))
        metals_sum = {}

        self._load_structure()

        met_list = mu.get_metal_atoms(self._get_structure(), self.data_library.get_metal_atoms())

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

            if self.args['check_only']:
                print ('Running with --check_only. Nothing else to do.')
            else:
                input_sess = ParamInput("Remove", opts.remove_metals, self.args['non_interactive'])
                input_sess.add_option_all()
                input_sess.add_option_none()
                input_sess.add_option('atids', sorted(at_groups), case='sensitive', multiple=True)
                input_sess.add_option('resids', met_rids, case='sensitive', multiple=True)
                [input_option, opts.remove_metals] = input_sess.run()

                if input_option == "error":
                    sys.stderr.write ('Error: unknown selection {}\n'.format(opts.remove_metals))
                    self.summary['metals'] = metals_sum
                    return 1
        
                if input_option == 'none':
                    print ("Nothing to do")
                else:
                    if input_option == 'all':
                        to_remove = met_list

                    elif input_option == 'resids':
                        to_remove = []
                        rid_list = opts.remove_metals.split(',')
                        for at in met_list:
                            r = at.get_parent()
                            if mu.residue_num(r) in rid_list:
                                to_remove.append(at)
                    
                    elif input_option == 'atids':
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

        lig_list = mu.get_ligands(self._get_structure(), incl_water=True)

        wat_list = []
        for r in lig_list:
            if mu.is_wat(r):
                wat_list.append(r)

        if len(wat_list) > 0:
            print ('{} Water molecules detected'.format(len(wat_list)))
            remwat_sum['n_detected'] = len(wat_list)

            if self.args['check_only']:
                print ('Running with --check_only. Nothing else to do.')
            else:
                input_line = ParamInput('Remove', opts.remove_wat, self.args['non_interactive'])
                input_line.add_option_yes_no()
                [input_option, opts.remove_wat] = input_line.run()

                if input_option == 'error':
                    print ('Warning: unknown option {}'.format(opts.remove_wat))
                    self.summary['remwat'] = remwat_sum
                    return 1

                if input_option == 'yes':
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

            if self.args['check_only']:
                print ('Running with --check_only. Nothing else to do.')
            else:
                input_line = ParamInput('Remove', opts.remove_ligands, self.args['non_interactive'])
                input_line.add_option_all()
                input_line.add_option_none()
                input_line.add_option('byrids', sorted(rids), multiple=True)
                input_line.add_option('byresnum', rnums,  case='sensitive', multiple=True)
                [input_option, opts.remove_ligands]= input_line.run()
                
                if input_option == 'error':
                    print ('Error: unknown selection {}'.format(opts.remove_ligands))
                    self.summary['ligands'] = ligands_sum
                    return 1

                ligands_sum['removed'] = {'opt':opts.remove_ligands, 'lst':[]}

                to_remove = []

                if input_option == 'none':
                    print ("Nothing to do")
                else:
                    if input_option == 'all':
                        to_remove = lig_list
                    elif input_option == 'byrids':
                        rm = opts.remove_ligands.split(',')
                        for r in lig_list:
                            if r.get_resname() in rm:
                                to_remove.append(r)
                    elif input_option == 'byresnum':
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

            if self.args['check_only']:
                print ('Running with --check_only. Nothing else to do.')
            else:
                input_line = ParamInput('Remove hydrogen atoms', opts.remove_h, self.args['non_interactive'])
                input_line.add_option_yes_no()
                [input_option, opts.remove_h] = input_line.run()
                
                if input_option == 'error':
                    print ('Warning: unknown option {}'.format(opts.remove_h))
                    self.summary['remh'] = remh_sum
                    return 1

                if input_option == 'yes':
                    n = 0
                    for r in remh_list:
                        mu.remove_H_from_r(r['r'])
                        n += 1
                    print ('Hydrogen atom removed from {} residues'.format(n))
                    self.stm.modified=True
                    remh_sum['n_removed'] = n
        else:
            print ("No residues with hydrogen atoms detected")

        self.summary['remh'] = remh_sum

    def getss (self, options):
        print ("Running getss")
        getss_sum = {}

        self._load_structure()

        SS_bonds = mu.get_all_at2at_distances(self._get_structure(), 'SG', self.data_library.get_distance('SS_DIST'))

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
        
        [amide_res, amide_atoms] = self.data_library.get_amide_data()
        CLASH_DIST = self.data_library.get_distance('CLASH_DIST')
        atoms_list = {
            'donor': self.data_library.get_atom_feature_list('polar_donor_atoms'),
            'acceptor' :self.data_library.get_atom_feature_list('polar_acceptor_atoms')
        }
        
        amide_list = []
        for r in self._get_structure().get_residues():
            if r.get_resname() in amide_res:
                amide_list.append(r)

        amide_sum['n_amides'] = len(amide_list)
        
        if len(amide_list) > 0:
            if len(self.rr_dist) == 0:
                self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distance('R_R_CUTOFF'))
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
                            if not at1.id in amide_atoms and not at2.id in amide_atoms:
                                continue
                            for cls in ['acceptor', 'donor']:
                                if dist < CLASH_DIST[cls]:
                                    if not mu.is_at_in_list(at1, atoms_list[cls]) or not mu.is_at_in_list(at2, atoms_list[cls]):
                                        continue
                                    if at1.id in amide_atoms and r1 in amide_list:
                                        res_to_fix.append(r1)
                                        rnums.append(mu.residue_num(r1))
                                    if at2.id in amide_atoms and r2 in amide_list:
                                        res_to_fix.append(r2)
                                        rnums.append(mu.residue_num(r2))
                                    cont_list.append(at_pair)
            if len(cont_list):
                print ('{} unusual contact(s) involving amide atoms found'.format(len(cont_list)))
                amide_sum['detected'] = []
                for at_pair in cont_list:
                    print (' {:12} {:12} {:8.3f} A'.format(mu.atom_id(at_pair[0], self.stm.nmodels > 1), mu.atom_id(at_pair[1], self.stm.nmodels > 1), at_pair[2]))
                    amide_sum['detected'].append({'at1':mu.atom_id(at_pair[0], self.stm.nmodels > 1), 'at2':mu.atom_id(at_pair[1], self.stm.nmodels > 1), 'dist': float(at_pair[2])})

                if self.args['check_only']:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    input_line = ParamInput('Fix amide atoms', opts.amide_fix, self.args['non_interactive'])
                    input_line.add_option_all()
                    input_line.add_option_none()
                    input_line.add_option('resnum',rnums, case='sensitive', multiple=True)
                    [input_option, opts.amide_fix] = input_line.run()
                    if input_option == 'error':
                        print ('Warning: unknown option {}'.format(opts.amide_fix))
                        self.summary['amide'] = amide_sum
                        return 1
                    
                    if input_option == 'none':
                        print ("Nothing to do")
                    else:
                        if input_option == 'all':
                            to_fix = res_to_fix
                        else:
                            to_fix = []
                            for r in res_to_fix:
                                if mu.residue_num(r) in opts.amide_fix.split(','):
                                    to_fix.append(r)
                        n = 0
                        for r in to_fix:
                            mu.invert_side_atoms(r, amide_res)
                            n += 1
                        print ('Amide residues fixed {} ({})'.format(opts.amide_fix, n))
                        self.stm.modified = True
        else:
            print ("No amide residues found")

        self.summary['amide'] = amide_sum

    def chiral(self, options):
        opts = _get_parameters(options, "chiral", '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
        print ('Running chiral. Options: {}'.format(' '.join(options)))
        chiral_sum = {}

        self._load_structure()
        
        chiral_res = self.data_library.get_chiral_data()
        chiral_list = []
        
        for r in self._get_structure().get_residues():
            if r.get_resname() in chiral_res:
                chiral_list.append(r)

        chiral_sum['n_chirals'] = len(chiral_list)
        
        if len(chiral_list) > 0:
            res_to_fix = []
            rnums = []
            for r in chiral_list:
                if not mu.check_chiral_residue(r, chiral_res):
                    res_to_fix.append(r)
                    rnums.append(mu.residue_num(r))
            if len(res_to_fix):
                print ('{} residues with incorrect chirality found'.format(len(res_to_fix)))
                chiral_sum['detected'] = []
                for r in res_to_fix:
                    print (' {:10}'.format(mu.residue_id(r, self.stm.nmodels > 1)))
                    chiral_sum['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))

                if self.args['check_only']:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    input_line = ParamInput('Fix chiralities', opts.chiral_fix, self.args['non_interactive'])
                    input_line.add_option_all()
                    input_line.add_option_none()
                    input_line.add_option('resnum', rnums, case='sensitive', multiple=True)
                    [input_option, opts.chiral_fix] = input_line.run()
                    
                    if input_option == 'error':
                        print ('Warning: unknown option {}'.format(opts.chiral_fix))
                        self.summary['chiral'] = chiral_sum
                        return1

                    if input_option == 'none':
                        print ("Nothing to do")
                    else:
                        if input_option == 'all':
                            to_fix = res_to_fix
                        else:
                            to_fix = []
                            for r in res_to_fix:
                                if mu.residue_num(r) in opts.chiral_fix.split(','):
                                    to_fix.append(r)
                        n = 0
                        for r in to_fix:
                            mu.invert_side_atoms(r, chiral_res)
                            n += 1
                        print ('Quiral residues fixed {} ({})'.format(opts.chiral_fix, n))
                        self.stm.modified= True
            else:
                print ("No residues with incorrect chirality found")
        else:
            print ("No chiral residues found")

        self.summary['chiral'] = chiral_sum

    def chiral_bck(self, options):
        opts = _get_parameters(options, "chiral_bck", '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
        print ('Running chiral. Options: {}'.format(' '.join(options)))
        chiral_bck_sum = {}

        self.args['check_only'] = True # Provisional

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

                if self.args['check_only']:
                    print ('Running with --check_only. Nothing else to do.')
                else:
                    input_line = ParamInput('Fix CA chiralities', opts.chiral_fix, elf.args['non_interactive'])
                    input_line.add_option_all()
                    input_line.add_option_none()
                    input_line.add_option('resnum', opts.chiral_fix, rnums, case='sensitive', multiple=True)
                    [input_option, opts.chiral_fix] = input_line.run()
                    if input_option == 'error':
                        print ('Warning: unknown option {}'.format(opts.amide_fix))
                        self.summary['chiral'] = chiral_sum
                        return 1

                    if input_option == 'none':
                        print ("Nothing to do")
                    else:
                        if input_option == 'all':
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
                        self.stm.modified=True
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

        CLASH_DIST = self.data_library.get_distance('CLASH_DIST')
        atom_lists = {
            'acceptor': self.data_library.get_atom_feature_list('polar_acceptor_atoms'),
            'donor': self.data_library.get_atom_feature_list('polar_donor_atoms'),
            'positive':self.data_library.get_atom_feature_list('positive_atoms'),
            'negative':self.data_library.get_atom_feature_list('negative_atoms'),
            'apolar':self.data_library.get_atom_feature_list('apolar_atoms'),
        }
        
        if len(self.rr_dist) == 0:
            self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distance('R_R_CUTOFF'))

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
                    if dist < CLASH_DIST['severe']:
                        if not rkey in clashes['severe']:
                            clashes['severe'][rkey] = at_pair
                        if dist < clashes['severe'][rkey][2]:
                            clashes['severe'][rkey] = at_pair
                    else:
                        for cls in ['apolar', 'acceptor', 'donor', 'positive', 'negative']:
                            if dist < CLASH_DIST[cls]:
                                if cls == 'apolar':
                                    #Only one of the atoms should be apolar
                                    if not mu.is_at_in_list(at1, atom_lists[cls]) and not mu.is_at_in_list(at2, atom_lists[cls]):
                                       continue
                                    #Remove n->n+2 backbone clashes. TODO Improve
                                    if abs(at1.get_parent().index - at2.get_parent().index) <= 2:
                                        continue
                                    #Remove Ca2+ looking like backbone CA's 
                                    if mu.is_hetatm(at1.get_parent()) and at1.id =='CA' or \
                                       mu.is_hetatm(at2.get_parent()) and at2.id =='CA':
                                        continue
                                else:
                                    # Both atoms should be of the same kind
                                    if not mu.is_at_in_list(at1, atom_lists[cls]) or not mu.is_at_in_list(at2, atom_lists[cls]):
                                        continue
                                if not rkey in clashes[cls]:
                                    clashes[cls][rkey] = at_pair
                                if dist < clashes[cls][rkey][2]:
                                    clashes[cls][rkey] = at_pair
        for cls in ['severe', 'apolar', 'acceptor', 'donor', 'positive', 'negative']:
            if len(clashes[cls]):
                print ('{} Steric {} clashes detected'.format(len(clashes[cls]), cls))
                for rkey in sorted(clashes[cls], key=lambda x: 10000 * clashes[cls][x][0].serial_number + clashes[cls][x][1].serial_number):
                    print (' {:12} {:12} {:8.3f} A'.format(mu.atom_id(clashes[cls][rkey][0], self.stm.nmodels > 1), mu.atom_id(clashes[cls][rkey][1], self.stm.nmodels > 1), clashes[cls][rkey][2]))
                    at_pair = clashes[cls][rkey]
                    clashes_sum['detected'][cls].append(
                        {
                            'at1':mu.atom_id(clashes[cls][rkey][0], self.stm.nmodels > 1),
                            'at2':mu.atom_id(clashes[cls][rkey][1], self.stm.nmodels > 1),
                            'dist': float(clashes[cls][rkey][2])
                        }
                    )
            else:
                print ('No {} clashes detected'.format(cls))

        self.summary['clashes'] = clashes_sum

    def fixside (self, options):
        opts = _get_parameters(options, "fixside", '--fix', 'fix_side', 'Add missing atoms to side chains')
        
        print ('Running fixside: Options {}'.format(' '.join(options)))
        fixside_sum = {}
        self._load_structure()
        
        miss_at_list = self.stm.get_missing_side_chain_atoms(self.data_library.get_valid_codes('protein'), self.data_library.get_atom_lists())
        
        if len(miss_at_list):
            rnums=[]
            fixside_sum['detected']={}
            print('{} Residues with missing side chain atoms found'.format(len(miss_at_list)))
            for r_at in miss_at_list:
                [r,at_list]= r_at
                print ('{:10} {}'.format(mu.residue_id(r), ','.join(at_list)))
                rnums.append(mu.residue_num(r))
                fixside_sum['detected'][mu.residue_id(r)] = at_list
            
            if not self.args['check_only']:
                input_line = ParamInput ('fixside', opts.fix_side, self.args['non_interactive'])
                input_line.add_option_all()
                input_line.add_option_none()
                input_line.add_option('resnum',rnums, case='sensitive', multiple=True)
                [input_option, opts.fix_side] = input_line.run()
                
                if input_option == 'error':
                    print ("Invalid option", opts.fix_side)
                    self.summary['fix_side']=fixside_sum
                    return 1
                
                fixside_sum['selected'] = opts.fix_side
                
                if input_option == 'none':
                    print ('Nothing to do')
                else:
                    if input_option == 'all':
                        to_fix = miss_at_list
                    else:
                        to_fix = []
                        for r_at in miss_at_list:
                            [r, at_list] = r_at
                            if mu.residue_num(r) in opts.fix_side.split(','):
                                to_fix.append(r_at)
                    n=0
                    fixside_sum['fixed']=[]
                    for r_at in to_fix:
                        self.stm.fix_side_chain(r_at)
                        n+=1
                        fixside_sum['fixed'].append(mu.residue_id(r_at[0]))
                    
                    print ('Fixed {} side chain(s)'.format(n))
                    
        else:
            print ("No residue with missing side chain atoms found")
        
        self.summary['fix_side']=fixside_sum
                                
                    
            
        
        
        
#===============================================================================

    def _load_structure(self, verbose=True):
        if not hasattr(self, 'stm'):
            if not self.args['non_interactive'] and self.args['input_structure_path'] is None:
                self.args['input_structure_path'] = input("Enter input structure path (PDB, mmcif | pdb:pdbid): ")
            self.stm = StructureManager(self.args['input_structure_path'], self.args['debug'])
            if verbose:
                print ('Structure {} loaded'.format(self.args['input_structure_path']))
                self.stm.print_stats()
            self.summary['stats'] = self.stm.get_stats()

    def _get_structure(self):
        return self.stm.get_structure()

    def _save_structure(self):
        if not self.args['non_interactive'] and self.args['output_structure_path'] is None:
            self.args['output_structure_path'] = input("Enter output structure path: ")
        self.stm.save_structure(self.args['output_structure_path'])

#===============================================================================

def _get_parameters (options, this_prog, this_param, this_dest, this_help, this_type=str):
    optionsParser = argparse.ArgumentParser(prog=this_prog)
    optionsParser.add_argument(this_param, dest=this_dest, help=this_help, type=this_type)
    return optionsParser.parse_args(options)
