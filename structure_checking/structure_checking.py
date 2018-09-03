"""
    Main Class for Structure Checking functionality

"""

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys

from structure_checking.help_manager import HelpManager
from structure_checking.json_writer import JSONWriter
from structure_checking.param_input import Dialog
from structure_checking.param_input import ParamInput

from structure_manager.data_lib_manager import DataLibManager
from structure_manager.structure_manager import StructureManager
from structure_manager.mutation_manager import MutationManager
from structure_manager.residue_lib_manager import ResidueLib

import structure_manager.model_utils as mu

# Interactive dialogs to complete command_line missing paramters
dialogs = Dialog()

#dialogs.add_option(command, prompt, destination, help_text, type(str))
dialogs.add_option('command_list', '--list', 'op_list', 'Command List File')
dialogs.add_option('models', '--select_model', 'select_model', 'Select model to keep', int)
dialogs.add_option('chains', '--select_chains', 'select_chains', 'Chains (All | Chain list comma separated)')
dialogs.add_option('altloc', '--select_altloc', 'select_altloc', 'Select altloc occupancy|alt_id')
dialogs.add_option('metals', '--remove', 'remove_metals', 'Remove Metal ions')
dialogs.add_option('remwat', '--remove', 'remove_wat', 'Remove Water molecules')
dialogs.add_option('ligands', '--remove', 'remove_ligands', 'Remove Ligand residues')
dialogs.add_option('remh', '--remove', 'remove_h', 'Remove Hydrogen atoms')
dialogs.add_option('amide', '--fix', 'amide_fix', 'Fix Residues (All | None | List)')
dialogs.add_option('chiral', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
dialogs.add_option('chiral_bck', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
dialogs.add_option('clashes', '--no_wat', 'discard_wat', 'Discard water molecules')
dialogs.add_option('fixside', '--fix', 'fix_side', 'Add missing atoms to side chains (All | None | List)')
dialogs.add_option('mutateside', '--mut', 'mut_list', 'Mutate side chains (Mutation List as [*:]arg234Thr)')

# Main class
class StructureChecking():
    def __init__(self, sets, args):
        self.args = args
        self.sets = sets
        self.summary = {}
        self.tmp_data = {}
        self.rr_dist = []
        self.data_library = DataLibManager(sets.data_library_path)
        if not 'Notebook' in self.args:
            self.args['Notebook'] = False
        if self.args['Notebook']:
            self.args['non_interactive'] = True
            self.args['check_only'] = False


    def launch(self):
        help_manager= HelpManager(self.sets.help_dir_path)
        help_manager.print_help('header')

        if '-h' in self.args['options'] or '--help' in self.args['options']:
            help_manager.print_help(self.args['command'])
            sys.exit()

        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        else:
            self.run_method(self.args['command'], self.args['options'])

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

    def run_method(self, command, opts):
        try:
            f_check = getattr(self, command + '_check')
        except AttributeError:
            print ('Error: {} command unknown or not implemented'.format(command))
            sys.exit(1)

        if not command in self.summary:
            self.summary[command] = {}

        msg = 'Running {}.'.format(command)

        if opts:
            self.summary[command]['opts'] = opts
            msg += ' Options: {} '.format(' '.join(opts))

        print (msg)

        self._load_structure()

    #Running checking method
        self.to_fix = f_check()

        if self.args['check_only'] or opts is None or opts == '':
            print ('Running  check_only. Nothing else to do.')
        elif self.to_fix:
            try:
                f_fix = getattr(self, command + '_fix')
            except AttributeError:
                print ('Error: {}_fix command unknown or not implemented'.format(command))
                sys.exit(1)
            if not self.args['Notebook']:
                if dialogs.exists(command):
                    opts = dialogs.get_parameter(command, opts)
                    opts = opts[dialogs.get_dialog(command)['dest']]
                else:
                    opts = ''
            f_fix(opts)

    def command_list(self, opts):
        opts = dialogs.get_parameter('command_list', opts)
        op_list = opts[dialogs.get_dialog('command_list')['dest']]

        op_list = ParamInput('Command List File', op_list, False).run()

        try:
            fh = open(op_list, "r")

        except OSError:
            print ('Error when opening file {}'.format(op_list))
            sys.exit(1)

        self._load_structure()

        i = 1
        for line in fh:
            if line == "\n" or line[0:1] == '#':
                continue
            print ("\nStep {}: {}".format(i, line))
            data = line.split()
            command = data[0]
            opts = data[1:]
            self.run_method(command, opts)
            i += 1

        print ("Command list completed")

    def print_stats(self, verbose=True):
        self._load_structure(print_stats=False)
        self.stm.print_stats()
# =============================================================================
    def models(self, opts=None):
        self.run_method('models', opts)

    def models_check(self):
        print ('{} Model(s) detected'.format(self.stm.nmodels))
        self.summary['models']['detected'] = {'nmodels': self.stm.nmodels}
        if self.stm.nmodels > 1:
            self.summary['models']['detected']['type'] = self.stm.models_type
            if self.stm.models_type['type'] == mu.ENSM:
                print ('Models superimpose, RMSd: {:8.3f} A, guessed as ensemble type (NMR / MD TRAJ)'.format(self.stm.models_type['rmsd']))
            elif self.stm.models_type['type'] == mu.BUNIT:
                print ('Models do not superimpose, RMSd: {:8.3f} A, guessed as Biounit type'.format(self.stm.models_type['rmsd']))
            else:
                print ('Models type unknown')
            return True
        else:
            print ("Single model found")
            return False

    def models_fix(self, select_model):
        input_line = ParamInput('Select Model Num', select_model, self.args['non_interactive'])
        input_line.add_option_all ()
        input_line.add_option ('modelno', [], opt_type='int', min=1, max=self.stm.nmodels)
        [input_option, select_model] = input_line.run()

        if input_option == 'error':
            print ('Error: unknown model {}'.format(select_model), file=sys.stderr)
            self.summary['models']['error'] = 'Unknown model'
            return 1

        print ('Selecting model num. {}'.format(select_model))
        if input_option != 'all':
            self.stm.select_model(select_model)

        self.summary['models']['selected_model'] = select_model

# =============================================================================
    def chains (self, opts=None):
        self.run_method('chains', opts)

    def chains_check(self):
        print ('{} Chains detected'.format(len(self.stm.chain_ids)))
        for ch_id in sorted(self.stm.chain_ids):
            print ('  {}: {}'.format(ch_id, mu.chain_type_labels[self.stm.chain_ids[ch_id]]))
        self.summary['chains']['detected'] = self.stm.chain_ids
        return len(self.stm.chains_ids) > 1

    def chains_fix(self, select_chains):
        self.summary['chains']['selected'] = {}
        input_line = ParamInput('Select chain', select_chains, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option('chid', sorted(self.stm.chain_ids), multiple=True, case="sensitive")
        [input_option, select_chains] = input_line.run()

        if input_option == 'error':
            print ('Error unknown selection: {}'.format(opts.select_chains))
            self.summary['chains']['error'] = 'Unknown selection'
            return 1

        if input_option == 'all':
            print ('Selecting all chains')
        else:
            self.stm.select_chains(select_chains)
            print ('Selecting chain(s) {}'.format(select_chains))
            self.summary['chains']['selected'] = select_chains.split(',')
            self.stm.set_chain_ids()
            self.summary['chains']['selected'] = self.stm.chain_ids

# =============================================================================
    def altloc(self, opts=None):
        self.run_method('altloc', opts)

    def altloc_check(self): #TODO improve output

        self.alt_loc_res = mu.get_altloc_residues(self._get_structure())

        if len(self.alt_loc_res) > 0:
            print ('Detected {} residues with alternative location labels'.format(len(self.alt_loc_res)))

            self.summary['altloc']['detected'] = {}
            self.alt_loc_rnums = []
            self.altlocs = {}
            for r in sorted(self.alt_loc_res, key=lambda x: x.index):
                rid = mu.residue_id(r)
                print(rid)
                self.alt_loc_rnums.append(mu.residue_num(r))
                self.summary['altloc']['detected'][rid] = []
                self.altlocs[r] = sorted(self.alt_loc_res[r][0].child_dict)
                for at in self.alt_loc_res[r]:
                    s = '  {:4}'.format(at.id)
                    for alt in sorted(at.child_dict):
                        s += ' {} ({:4.2f})'.format(alt, at.child_dict[alt].occupancy)
                    print (s)
                    self.summary['altloc']['detected'][rid].append({
                                                                   'atom':at.id,
                                                                   'loc_label':alt,
                                                                   'occupancy':at.child_dict[alt].occupancy}
                                                                   )
            return True
        else:
            print ("No residues with alternative location labels detected")
            return False

    def altloc_fix(self, select_altloc):
        #Prepare the longest possible list of alternatives
        altlocs = []
        l = 0
        for r in self.altlocs:
            if len(self.altlocs[r]) > l:
                altlocs = self.altlocs[r]
                l = len(self.altlocs[r])

        input_line = ParamInput('Select alternative', select_altloc, self.args['non_interactive'])
        input_line.add_option('occup', ['occupancy'])
        input_line.add_option('altids', altlocs, case='upper')
        input_line.add_option('resnum', self.alt_loc_rnums, opt_type='pair_list', list2=altlocs, case='sensitive', multiple=True)
        [input_option, select_altloc] = input_line.run()

        if input_option == 'error':
            print ('Error: Unknown selection {} '.format(opts.select_altloc), file=sys.stderr)
            self.summary['altloc']['error'] = "Unknown selection"
            return 1

        print ('Selecting location {}'.format(select_altloc))
        if input_option == 'occup' or input_option == 'altids':
            to_fix = {}
            for r in self.alt_loc_res.keys():
                to_fix[r] = {}
                to_fix[r]['ats'] = self.alt_loc_res[r]
                to_fix[r]['select'] = select_altloc
        elif input_option == 'resnum':
            to_fix = {}
            selected_rnums = {}
            for rsel in select_altloc.split(','):
                [rn, alt] = rsel.split(':')
                selected_rnums[rn] = alt
            for r in self.alt_loc_res.keys():
                rn0 = mu.residue_num(r)
                if rn0 in selected_rnums.keys():
                    to_fix[r] = {}
                    to_fix[r]['ats'] = self.alt_loc_res[r]
                    to_fix[r]['select'] = selected_rnums[rn0]

        self.summary['altloc']['selected'] = select_altloc
        for r in to_fix.keys():
            self.stm.select_altloc_residue(r, to_fix[r])

# =============================================================================
    def metals (self, opts=None):
        self.run_method('metals', opts)

    def metals_check (self):
        self.met_list = mu.get_metal_atoms(self._get_structure(), self.data_library.get_metal_atoms())

        if len(self.met_list) > 1:
            print ('{} Metal ions found'.format(len(self.met_list)))
            self.summary['metals']['detected'] = []
            self.met_rids = []
            self.at_groups = {}
            for at in sorted(self.met_list, key=lambda x: x.serial_number):
                print ("  ", mu.atom_id(at, self.stm.nmodels > 1))
                r = at.get_parent()
                self.met_rids.append(mu.residue_num(r))
                if not at.id in self.at_groups:
                    self.at_groups[at.id] = []
                self.at_groups[at.id].append(at)
                self.summary['metals']['detected'].append(mu.residue_num(r))
            return True
        else:
            print ("No metal ions present")
            return False

    def metals_fix(self, remove_metals):
        input_sess = ParamInput("Remove", remove_metals, self.args['non_interactive'])
        input_sess.add_option_all()
        input_sess.add_option_none()
        input_sess.add_option('atids', sorted(self.at_groups), case='sensitive', multiple=True)
        input_sess.add_option('resids', self.met_rids, case='sensitive', multiple=True)
        [input_option, remove_metals] = input_sess.run()

        if input_option == "error":
            sys.stderr.write ('Error: unknown selection {}\n'.format(opts.remove_metals))
            self.summary['metals']['error'] = 'Unknown selection'
            return 1

        if input_option == 'none':
            print ("Nothing to do")
        else:
            if input_option == 'all':
                to_remove = self.met_list
            elif input_option == 'resids':
                to_remove = []
                rid_list = remove_metals.split(',')
                for at in self.met_list:
                    r = at.get_parent()
                    if mu.residue_num(r) in rid_list:
                        to_remove.append(at)
            elif input_option == 'atids':
                to_remove = []
                for atid in remove_metals.split(','):
                    to_remove.extend(self.at_groups[atid])

            self.summary['metals']['removed'] = []

            n = 0
            for at in to_remove:
                self.summary['metals']['removed'].append(mu.residue_id(at.get_parent(), self.stm.nmodels > 1))
                self.stm.remove_residue(at.get_parent())
                n += 1

            print ('Metal Atoms removed {} ({:d})'.format(remove_metals, n))
            self.summary['metals']['n_removed'] = n

# =============================================================================
    def remwat(self, opts=None):
        self.run_method('remwat', opts)

    def remwat_check(self):
        self.lig_list = mu.get_ligands(self._get_structure(), incl_water=True)
        self.wat_list = []
        for r in self.lig_list:
            if mu.is_wat(r):
                self.wat_list.append(r)

        if len(self.wat_list) > 0:
            print ('{} Water molecules detected'.format(len(self.wat_list)))
            self.summary['remwat']['n_detected'] = len(self.wat_list)
            return True
        else:
            print ("No water molecules found")
            return False

    def remwat_fix(self, remove_wat):
        input_line = ParamInput('Remove', remove_wat, self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_wat] = input_line.run()

        if input_option == 'error':
            print ('Warning: unknown option {}'.format(opts.remove_wat))
            self.summary['remwat']['error'] = 'Unknown option'
            return 1

        if input_option == 'yes':
            n = 0
            for r in self.wat_list:
                self.stm.remove_residue(r)
                n += 1
            print ('{} Water molecules removed'.format(n))
            self.summary['remwat']['n_removed'] = n
            self.lig_list = mu.get_ligands(self._get_structure(), incl_water=True)

# =============================================================================
    def ligands(self, opts=None):
        self.run_method('ligands', opts)

    def ligands_check(self):
        self.lig_list = mu.get_ligands(self._get_structure(), incl_water=False)

        if len(self.lig_list):
            print ('{} Ligands detected '.format(len(self.lig_list)))
            self.ligand_rids = set()
            self.ligand_rnums = []
            self.summary['ligands']['detected'] = []

            for r in sorted(self.lig_list, key=lambda x: x.index):
                print (mu.residue_id(r, self.stm.nmodels > 1))
                self.summary['ligands']['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))
                self.ligand_rids.add(r.get_resname())
                self.ligand_rnums.append(mu.residue_num(r))
            return True
        else:
            print ("No ligands found")
            return False

# =============================================================================
    def ligands_fix(self, remove_ligands):
        input_line = ParamInput('Remove', remove_ligands, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('byrids', sorted(self.ligand_rids), multiple=True)
        input_line.add_option('byresnum', self.ligand_rnums, case='sensitive', multiple=True)
        [input_option, remove_ligands] = input_line.run()

        if input_option == 'error':
            print ('Error: unknown selection {}'.format(remove_ligands))
            self.summary['ligands']['error'] = 'Unknown selection'
            return 1

        self.summary['ligands']['removed'] = {'opt':remove_ligands, 'lst':[]}

        to_remove = []

        if input_option == 'none':
            print ("Nothing to do")
        else:
            if input_option == 'all':
                to_remove = self.lig_list
            elif input_option == 'byrids':
                rm = remove_ligands.split(',')
                for r in self.lig_list:
                    if r.get_resname() in rm:
                        to_remove.append(r)
            elif input_option == 'byresnum':
                rm = remove_ligands.split(',')
                for r in self.lig_list:
                    if mu.residue_num(r) in rm:
                        to_remove.append(r)
            n = 0
            for r in to_remove:
                self.summary['ligands']['removed']['lst'].append(mu.residue_id(r, self.stm.nmodels > 1))
                self.stm.remove_residue(r)
                n += 1

            print ('Ligands removed {} ({})'.format(remove_ligands, n))
            self.summary['ligands']['n_removed'] = n


# =============================================================================
    def remh(self, opts=None):
        self.run_method('remh', opts)

    def remh_check(self):
        self.remh_list = mu.get_residues_with_H(self._get_structure())

        if len(self.remh_list) > 0:
            print ('{} Residues containing H atoms detected'.format(len(self.remh_list)))
            self.summary['remh']['n_detected'] = len(self.remh_list)
            return True
        else:
            print ("No residues with Hydrogen atoms found")
            return False

    def remh_fix(self, remove_h):
        input_line = ParamInput('Remove hydrogen atoms', remove_h, self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_h] = input_line.run()

        if input_option == 'error':
            print ('Warning: unknown option {}'.format(remove_h))
            self.summary['remh']['error'] = 'Unknown option'
            return 1
        if input_option == 'yes':
            n = 0
            for r in self.remh_list:
                mu.remove_H_from_r(r['r'])
                n += 1
            print ('Hydrogen atoms removed from {} residues'.format(n))
            self.stm.modified = True
            self.summary['remh']['n_removed'] = n

# =============================================================================
    def getss (self, opts=None):
        self.run_method('getss', opts)

    def getss_check(self):
        self.SS_bonds = mu.get_all_at2at_distances(self._get_structure(), 'SG', self.data_library.get_distances('SS_DIST'))

        if len(self.SS_bonds):
            print ('{} Possible SS Bonds detected'.format(len(self.SS_bonds)))
            self.summary['getss']['detected'] = []
            for ssb in self.SS_bonds:
                print ('  {} {}{:8.3f}'.format(mu.atom_id(ssb[0], self.stm.nmodels > 1), mu.atom_id(ssb[1], self.stm.nmodels > 1), ssb[2]))
                self.summary['getss']['detected'].append({'at1':mu.atom_id(ssb[0], self.stm.nmodels > 1), 'at2':mu.atom_id(ssb[1], self.stm.nmodels > 1), 'dist': float(ssb[2])})
        else:
            print ("No SS bonds detected")

        return False

    def getss_fix(self):
        pass

# =============================================================================
    def amide(self, opts=None):
        self.run_method('amide', opts)

    def amide_check(self):
        [amide_res, amide_atoms] = self.data_library.get_amide_data()
        CLASH_DIST = self.data_library.get_distances('CLASH_DIST')
        
        atom_lists = {
            'polar_donor': self.data_library.get_atom_feature_list('polar_donor_atoms'),
            'polar_acceptor':self.data_library.get_atom_feature_list('polar_acceptor_atoms')
        }

        self.amide_list = []

        for r in self._get_structure().get_residues():
            if r.get_resname() in amide_res:
                self.amide_list.append(r)

        self.summary['amide']['n_amides'] = len(self.amide_list)

        if len(self.amide_list) > 0:
            if len(self.rr_dist) == 0:
                self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))
            self.amide_res_to_fix = []
            self.amide_cont_list = []
            self.amide_rnums = []
            for r_pair in self.rr_dist:
                [r1, r2, d] = r_pair
                if r1 != r2 and mu.same_model(r1, r2) and not mu.is_wat(r1) and not mu.is_wat(r2) \
                    and (r1 in self.amide_list or r2 in self.amide_list):
                    c_list = self._get_clashes(r1,r2,CLASH_DIST,atom_lists)
                    for cls in c_list:
                        if len(c_list[cls]):
                            [at1,at2,d]=c_list[cls]
                            add_pair=False
                            if at1.id in amide_atoms and r1 in self.amide_list:
                                self.amide_res_to_fix.append(r1)
                                self.amide_rnums.append(mu.residue_num(r1))
                                add_pair=True
                            if at2.id in amide_atoms and r2 in self.amide_list:
                                self.amide_res_to_fix.append(r2)
                                self.amide_rnums.append(mu.residue_num(r2))
                                add_pair=True
                            if add_pair:
                                self.amide_cont_list.append(c_list[cls])
            
            if len(self.amide_cont_list):
                print ('{} unusual contact(s) involving amide atoms found'.format(len(self.amide_cont_list)))
                self.summary['amide']['detected'] = []
                for at_pair in self.amide_cont_list:
                    print (' {:12} {:12} {:8.3f} A'.format(mu.atom_id(at_pair[0], self.stm.nmodels > 1), mu.atom_id(at_pair[1], self.stm.nmodels > 1), at_pair[2]))
                    self.summary['amide']['detected'].append({'at1':mu.atom_id(at_pair[0], self.stm.nmodels > 1), 'at2':mu.atom_id(at_pair[1], self.stm.nmodels > 1), 'dist': float(at_pair[2])})
                return True
            else:
                print ("No unusual contact(s) involving amide atoms found")
                return False
        else:
            print ("No amide residues found")
            return False

# =============================================================================
    def amide_fix(self, amide_fix):
        amide_res = self.data_library.get_amide_data()[0]

        input_line = ParamInput('Fix amide atoms', amide_fix, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', self.amide_rnums, case='sensitive', multiple=True)
        [input_option, amide_fix] = input_line.run()

        if input_option == 'error':
            print ('Warning: unknown option {}'.format(amide_fix))
            self.summary['amide']['error'] = 'Unknown option'
            return 1

        if input_option == 'none':
            print ("Nothing to do")
        else:
            if input_option == 'all':
                to_fix = self.amide_res_to_fix
            else:
                to_fix = []
                for r in self.amide_res_to_fix:
                    if mu.residue_num(r) in amide_fix.split(','):
                        to_fix.append(r)
            n = 0
            for r in to_fix:
                mu.invert_side_atoms(r, amide_res)
                n += 1
            print ('Amide residues fixed {} ({})'.format(amide_fix, n))
            self.stm.modified = True

# =============================================================================
    def chiral(self, opts=None):
        self.run_method('chiral', opts)

    def chiral_check(self):

        chiral_res = self.data_library.get_chiral_data()

        self.chiral_list = []

        for r in self._get_structure().get_residues():
            if r.get_resname() in chiral_res:
                self.chiral_list.append(r)

        self.summary['chiral']['n_chirals'] = len(self.chiral_list)

        if len(self.chiral_list) > 0:
            self.chiral_res_to_fix = []
            self.chiral_rnums = []
            for r in self.chiral_list:
                if not mu.check_chiral_residue(r, chiral_res):
                    self.chiral_res_to_fix.append(r)
                    self.chiral_rnums.append(mu.residue_num(r))
            if len(self.chiral_res_to_fix):
                print ('{} residues with incorrect side-chain chirality found'.format(len(self.chiral_res_to_fix)))
                self.summary['chiral']['detected'] = []
                for r in self.chiral_res_to_fix:
                    print (' {:10}'.format(mu.residue_id(r, self.stm.nmodels > 1)))
                    self.summary['chiral']['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))
                return True
            else:
                print ("No residues with incorrent side-chain chirality found")
                return False
        else:
            print ("No chiral side-chains found")
            return False

    def chiral_fix(self, chiral_fix):

        chiral_res = self.data_library.get_chiral_data()

        input_line = ParamInput('Fix chiralities', chiral_fix, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', self.chiral_rnums, case='sensitive', multiple=True)
        [input_option, chiral_fix] = input_line.run()

        if input_option == 'error':
            print ('Warning: unknown option {}'.format(chiral_fix))
            self.summary['chiral']['error'] = 'Unknown option'
            return 1

        if input_option == 'none':
            print ("Nothing to do")
        else:
            if input_option == 'all':
                to_fix = self.chiral_res_to_fix
            else:
                to_fix = []
                for r in self.chiral_res_to_fix:
                    if mu.residue_num(r) in chiral_fix.split(','):
                        to_fix.append(r)
            n = 0
            for r in to_fix:
                mu.invert_side_atoms(r, chiral_res)
                if r.get_resname() == 'ILE':
                    mu.delete_atom(r,'CD1')
                    if not hasattr(self,'residue_lib'):
                        self.residue_lib = ResidueLib(self.sets.res_library_path)
                    mu.build_atom(r, 'CD1', self.residue_lib,'ILE')
                n += 1
            print ('Quiral residues fixed {} ({})'.format(chiral_fix, n))
            self.stm.modified = True

# =============================================================================
    def chiral_bck(self, opts=None):
        self.run_method('chiral_bck', opts)

    def chiral_bck_check(self):

        self.chiral_bck_list = []
        for ch in self._get_structure().get_chains():
            if self.stm.chain_ids[ch.id] == mu.PROTEIN:
                for r in ch.get_residues():
                    if r.get_resname() != 'GLY' and not mu.is_hetatm(r):
                        self.chiral_bck_list.append(r)
        self.summary['chiral_bck']['n_chirals'] = len(self.chiral_bck_list)

        if len(self.chiral_bck_list) > 0:
            self.chiral_bck_res_to_fix = []
            self.chiral_bck_rnums = []
            for r in self.chiral_bck_list:
                if not mu.check_chiral_ca(r):
                    self.chiral_bck_res_to_fix.append(r)
                    self.chiral_bck_rnums.append(mu.residue_num(r))
            if len(self.chiral_bck_res_to_fix):
                print ('{} residues with incorrect backbone chirality found'.format(len(self.chiral_bck_res_to_fix)))
                self.summary['chiral_bck']['detected'] = []
                for r in self.chiral_bck_res_to_fix:
                    print (' {:10}'.format(mu.residue_id(r, self.stm.nmodels > 1)))
                    self.summary['chiral_bck']['detected'].append(mu.residue_id(r, self.stm.nmodels > 1))
                return False
            else:
                print ("No residues with incorrect backbone chirality found")
                return False
        else:
            print ("No residues with chiral CA found (weird!!)")
            return False


    def chiral_bck_fix(self, chiral_fix):
        return 0 #TODO
        input_line = ParamInput('Fix CA chiralities', chiral_fix, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', self.chiral_bck_rnums, case='sensitive', multiple=True)
        [input_option, chiral_fix] = input_line.run()
        if input_option == 'error':
            print ('Warning: unknown option {}'.format(amide_fix))
            self.summary['chiral']['error'] = 'Unknown option'
            return 1

        if input_option == 'none':
            print ("Nothing to do")
        else:
            if input_option == 'all':
                to_fix = self.chiral_bck_res_to_fix
            else:
                to_fix = []
                for r in self.chiral_bck_res_to_fix:
                    if mu.residue_num(r) in chiral_fix.split(','):
                        to_fix.append(r)
            n = 0
            for r in to_fix:
                mu.stm.invert_chiral_CA(r) # TODO
                n += 1
            print ('Quiral residues fixed {} ({})'.format(chiral_fix, n))
            self.stm.modified = True

# =============================================================================
    def clashes(self, opts=None):
        self.run_method('clashes', opts)

    def clashes_check(self):
        contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']
        
        CLASH_DIST = self.data_library.get_distances('CLASH_DIST')
        
        atom_lists={}
        self.clash_list={}
        self.summary['clashes']['detected']={}        
        for cls in contact_types:
            self.clash_list[cls]={}
            self.summary['clashes']['detected'][cls]=[]
            if cls == 'severe':
                continue
            atom_lists[cls]= self.data_library.get_atom_feature_list(cls + '_atoms')
        
        if len(self.rr_dist) == 0:
            self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))

        for r_pair in self.rr_dist:
            [r1, r2, d] = r_pair
            if mu.is_wat(r1) or mu.is_wat(r2):
                continue
            c_list = self._get_clashes(r1,r2,CLASH_DIST, atom_lists)
            rkey = mu.residue_id(r1)+'-'+mu.residue_id(r2)
            for cls in c_list:
                if len(c_list[cls]):
                    self.clash_list[cls][rkey] = c_list[cls]
        
        
        for cls in contact_types:
            if len(self.clash_list[cls]):
                print ('{} Steric {} clashes detected'.format(len(self.clash_list[cls]), cls))
                for rkey in sorted(self.clash_list[cls], key=lambda x: 10000 * self.clash_list[cls][x][0].serial_number + self.clash_list[cls][x][1].serial_number):
                    print (' {:12} {:12} {:8.3f} A'.format(mu.atom_id(self.clash_list[cls][rkey][0], self.stm.nmodels > 1), mu.atom_id(self.clash_list[cls][rkey][1], self.stm.nmodels > 1), self.clash_list[cls][rkey][2]))
                    at_pair = self.clash_list[cls][rkey]
                    self.summary['clashes']['detected'][cls].append(
                                                                    {
                                                                    'at1':mu.atom_id(self.clash_list[cls][rkey][0], self.stm.nmodels > 1),
                                                                    'at2':mu.atom_id(self.clash_list[cls][rkey][1], self.stm.nmodels > 1),
                                                                    'dist': float(self.clash_list[cls][rkey][2])
                                                                    }
                                                                    )
            else:
                print ('No {} clashes detected'.format(cls))

        return False

    def clashes_fix(self, res_to_fix):
        pass
#===============================================================================
    def _get_clashes(self, r1, r2, CLASH_DIST, atom_lists):

        clash_list={}
        min_dist={}
        for cls in atom_lists:
            clash_list[cls]=[]
            min_dist[cls]=999.
        if r1 != r2 and not mu.seq_consecutive(r1, r2) and mu.same_model(r1, r2):
            for at_pair in mu.get_all_rr_distances(r1, r2):
                [at1, at2, dist] = at_pair
                if 'severe' in atom_lists and dist < CLASH_DIST['severe']:
                    if dist < min_dist:
                        clash_list['severe'] = at_pair
                        min_dist['severe'] = dist
                else:
                    for cls in atom_lists:
                        if cls == 'apolar':
                            #Only one of the atoms should be apolar
                            if not mu.is_at_in_list(at1, atom_lists[cls]) and not mu.is_at_in_list(at2, atom_lists[cls]):
                                continue
                            #Remove n->n+2 backbone clashes. TODO Improve
                            if abs(at1.get_parent().index - at2.get_parent().index) <= 2:
                                continue
                            #Remove Ca2+ looking like backbone CA's
                            if mu.is_hetatm(at1.get_parent()) and at1.id == 'CA' or \
                                mu.is_hetatm(at2.get_parent()) and at2.id == 'CA':
                                    continue
                        else:
                            # Both atoms should be of the same kind
                            if not mu.is_at_in_list(at1, atom_lists[cls]) or not mu.is_at_in_list(at2, atom_lists[cls]):
                                continue
                        if dist < CLASH_DIST[cls]:
                            if dist < min_dist[cls]:
                                clash_list[cls] = at_pair
                                min_dist[cls] = dist

        return clash_list

# =============================================================================
    def fixside (self, opts=None):
        self.run_method('fixside', opts)

    def fixside_check (self):
        self.miss_at_list = self.stm.get_missing_side_chain_atoms(self.data_library.get_valid_codes('protein'), self.data_library.get_atom_lists())

        if len(self.miss_at_list):
            self.fixside_rnums = []
            self.summary['fixside']['detected'] = {}
            print('{} Residues with missing side chain atoms found'.format(len(self.miss_at_list)))
            for r_at in self.miss_at_list:
                [r, at_list] = r_at
                print ('{:10} {}'.format(mu.residue_id(r), ','.join(at_list)))
                self.fixside_rnums.append(mu.residue_num(r))
                self.summary['fixside']['detected'][mu.residue_id(r)] = at_list
            return True
        else:
            print ("No residues with missing side chain atoms found")
            return False

    def fixside_fix(self, fix_side):
        input_line = ParamInput ('fixside', fix_side, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', self.fixside_rnums, case='sensitive', multiple=True)
        [input_option, fix_side] = input_line.run()

        if input_option == 'error':
            print ("Invalid option", fix_side)
            self.summary['fix_side']['error'] = "Unknown option"
            return 1

        self.summary['fixside']['selected'] = fix_side

        if input_option == 'none':
            print ('Nothing to do')
        else:
            if input_option == 'all':
                to_fix = self.miss_at_list
            else:
                to_fix = []
                for r_at in self.miss_at_list:
                    if mu.residue_num(r_at[0]) in fix_side.split(','):
                        to_fix.append(r_at)
            print ("Fixing side chains")
            n = 0
            self.summary['fixside']['fixed'] = []
            if not hasattr(self,'residue_lib'):
                self.residue_lib = ResidueLib(self.sets.res_library_path)
            for r_at in to_fix:
                print (mu.residue_id(r_at[0]))
                mu.remove_H_from_r(r_at[0], verbose= True)
                self.stm.fix_side_chain(r_at, self.residue_lib)
                n += 1
                self.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))

            print ('Fixed {} side chain(s)'.format(n))
# =============================================================================
    def mutateside(self, mut_list):
        self.run_method('mutateside', mut_list)

    def mutateside_check(self):
        return True

    def mutateside_fix(self, mut_list):
        input_line = ParamInput ('Mutation list', mut_list, self.args['non_interactive'])
        mut_list = input_line.run()
        self.mutations = MutationManager(mut_list)

        self.mutations.prepare_mutations(self.stm.st)
        print ('Mutations to perform')
        print (self.mutations)
        if not hasattr(self,'mutation_rules'):
            self.mutation_rules = self.data_library.get_mutation_map()

        if not hasattr(self,'residue_lib'):
            self.residue_lib = ResidueLib(self.sets.res_library_path)

        self.mutations.apply_mutations (self.stm.st, self.mutation_rules, self.residue_lib)

        self.stm.modified=True


#===============================================================================

    def _load_structure(self, verbose=True, print_stats=True):
        if not hasattr(self, 'stm'):
            if not self.args['non_interactive'] and self.args['input_structure_path'] is None:
                self.args['input_structure_path'] = input("Enter input structure path (PDB, mmcif | pdb:pdbid): ")
            self.stm = StructureManager(self.args['input_structure_path'], self.args['debug'])
            if verbose:
                print ('Structure {} loaded'.format(self.args['input_structure_path']))
            if print_stats:
                self.stm.print_stats()
            self.summary['stats'] = self.stm.get_stats()

    def _get_structure(self):
        return self.stm.get_structure()

    def _save_structure(self):
        if not self.args['non_interactive'] and self.args['output_structure_path'] is None:
            self.args['output_structure_path'] = input("Enter output structure path: ")
        self.stm.save_structure(self.args['output_structure_path'])

#===============================================================================

