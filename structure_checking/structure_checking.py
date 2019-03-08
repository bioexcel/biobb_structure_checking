"""
    Main Class for Structure Checking functionality

"""

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys
import numpy as np

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

#dialogs.add_option(command, prompt, destinmore ation, help_text, type(str))
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
dialogs.add_option('backbone', '--fix', 'fix_back', 'Add missing O atoms to backbone (All | None | List)')
dialogs.add_option('mutateside', '--mut', 'mut_list', 'Mutate side chains (Mutation List as [*:]arg234Thr)')
dialogs.add_option('addH', '--mode', 'mode', 'Selection mode (None | auto | interactive | interactive_his | ph )')

AVAILABLE_METHODS=[
    'models','chains','inscodes','altloc','remh','remwat', 'metals','ligands',
    'getss','amide','chiral','chiral_bck','fixside','backbone','cistransbck',
    'clashes']

# Main class
class StructureChecking():
    def __init__(self, sets, args):
        self.args = args
        self.sets = sets
        self.summary = {}
        self.tmp_data = {}
        self.rr_dist = []
        self.pdb_server=self.args['pdb_server']
        self.cache_dir=self.args['cache_dir']
        self.data_library = DataLibManager(sets.data_library_path)
        if not 'Notebook' in self.args:
            self.args['Notebook'] = False
        if self.args['Notebook']:
            self.args['non_interactive'] = True
            self.args['check_only'] = False

        self.CLASH_DISTS = self.data_library.get_distances('CLASH_DIST')

    def launch(self):
        help_manager= HelpManager(self.sets.help_dir_path)
        help_manager.print_help('header')

        if '-h' in self.args['options'] or '--help' in self.args['options']:
            help_manager.print_help(self.args['command'])
            sys.exit()

        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        elif self.args['command'] == 'checkall':
            self.checkall(self.args['options'])
        elif self.args['command'] == 'load':
            self._load_structure()
            sys.exit(0)
        elif self.args['command'] == 'stats':
            self._load_structure()
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

        if not self.args['quiet']:
            print (msg)
        self._load_structure()

    #Running checking method
        self.to_fix = f_check()

        if self.args['check_only'] or opts is None or opts == '':
            if not self.args['quiet']:
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
            if not self.args['quiet']:
                print ("\nStep {}: {}".format(i, line))
            data = line.split()
            command = data[0]
            opts = data[1:]
            self.run_method(command, opts)
            i += 1

        print ("Command list completed")

    def checkall(self, opts):
        self.args['check_only'] = True
        self._load_structure()
        for meth in AVAILABLE_METHODS:
            self.run_method(meth, opts)

    def print_stats(self, verbose=True):
        self._load_structure(print_stats=False)
        self.stm.print_stats()
# =============================================================================
    def models(self, opts=None):
        self.run_method('models', opts)

    def models_check(self):
        print ('{} Model(s) detected'.format(self.stm.nmodels))
        self.summary['models'] = {'nmodels': self.stm.nmodels}
        if self.stm.nmodels > 1:
            self.summary['models']['type'] = self.stm.models_type
            if self.stm.models_type['type'] == mu.ENSM:
                print ('Models superimpose, RMSd: {:8.3f} A, guessed as ensemble type (NMR / MD TRAJ)'.format(self.stm.models_type['rmsd']))
            elif self.stm.models_type['type'] == mu.BUNIT:
                print ('Models do not superimpose, RMSd: {:8.3f} A, guessed as Biounit type'.format(self.stm.models_type['rmsd']))
            else:
                print ('Models type unknown')
            return True
        else:
            if not self.args['quiet']:
                print ("Single model found")
            return False

    def models_fix(self, select_model):
        input_line = ParamInput('Select Model Num', select_model, self.args['non_interactive'])
        input_line.add_option_all ()
        input_line.add_option ('modelno', [], opt_type='int', min_val=1, max_val=self.stm.nmodels)
        [input_option, select_model] = input_line.run()

        if input_option == 'error':
            print ('Error: unknown model {}'.format(select_model), file=sys.stderr)
            self.summary['models']['error'] = 'Unknown model '+select_model
            return 1

        print ('Selecting model num. {}'.format(select_model))
        if input_option != 'all':
            self.stm.select_model(select_model)

        self.summary['models']['selected'] = select_model

# =============================================================================
    def chains (self, opts=None):
        self.run_method('chains', opts)

    def chains_check(self):
        print ('{} Chain(s) detected'.format(len(self.stm.chain_ids)))
        for ch_id in sorted(self.stm.chain_ids):
            if isinstance(self.stm.chain_ids[ch_id],list):
                print (' {}: Unknown (PROTEIN: {s[0]:4.2f} DNA: {s[1]:4.2f} RNA: {s[2]:4.2f} Other: {s[3]:4.2f})'.format(ch_id, s=self.stm.chain_ids[ch_id]))
            else:
                print (' {}: {}'.format(ch_id, mu.chain_type_labels[self.stm.chain_ids[ch_id]]))
        self.summary['chains'] = {'ids':self.stm.chain_ids}
        return len(self.stm.chains_ids) > 1

    def chains_fix(self, select_chains):
        self.summary['chains']['selected'] = {}
        input_line = ParamInput('Select chain', select_chains, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option('chid', sorted(self.stm.chain_ids), multiple=True, case="sensitive")
        [input_option, select_chains] = input_line.run()

        if input_option == 'error':
            print ('Error unknown selection: {}'.format(select_chains))
            self.summary['chains']['error'] = 'Unknown selection ' + select_chains
            return 1

        if input_option == 'all':
            print ('Selecting all chains')
        else:
            self.stm.select_chains(select_chains)
            print ('Selecting chain(s) {}'.format(select_chains))
            self.stm.set_chain_ids()
            self.summary['chains']['selected'] = self.stm.chain_ids

# =============================================================================
    def inscodes(self, opts=None):
        self.run_method('inscodes', opts)

    def inscodes_check(self):
        self.stm.get_ins_codes()
        if self.stm.ins_codes_list:
            print ('{} Residues with insertion codes found'.format(len(self.stm.ins_codes_list)))
            self.summary['inscodes']=[]
            for r in self.stm.ins_codes_list:
                print (mu.residue_id(r))
                self.summary['inscodes'].append(mu.residue_id(r))
        else:
            if not self.args['quiet']:
                print ("No residues with insertion codes found")

    def inscodes_fix(self,select_codes):
        pass
# =============================================================================
    def altloc(self, opts=None):
        self.run_method('altloc', opts)

    def altloc_check(self): #TODO improve output
        self.alt_loc_res = mu.get_altloc_residues(self._get_structure())
        if len(self.alt_loc_res) > 0:
            print ('Detected {} residues with alternative location labels'.format(len(self.alt_loc_res)))

            self.summary['altloc'] = {}
            self.alt_loc_rnums = []
            self.altlocs = {}
            for r in sorted(self.alt_loc_res, key=lambda x: x.index):
                rid = mu.residue_id(r)
                print(rid)
                self.alt_loc_rnums.append(mu.residue_num(r))
                self.summary['altloc'][rid] = {}
                self.altlocs[r] = sorted(self.alt_loc_res[r][0].child_dict)
                for at in self.alt_loc_res[r]:
                    self.summary['altloc'][rid][at.id]=[]
                    s = '  {:4}'.format(at.id)
                    for alt in sorted(at.child_dict):
                        s += ' {} ({:4.2f})'.format(alt, at.child_dict[alt].occupancy)
                        self.summary['altloc'][rid][at.id].append({
                                'loc_label':alt,
                                'occupancy':at.child_dict[alt].occupancy
                            }
                            )
                    print (s)
            return True
        else:
            if not self.args['quiet']:
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
            print ('Error: Unknown selection {} '.format(select_altloc), file=sys.stderr)
            self.summary['altloc']['error'] = "Unknown selection " + select_altloc
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
            self.summary['metals'] = {'detected':[]}
            self.met_rids = []
            self.at_groups = {}
            for at in sorted(self.met_list, key=lambda x: x.serial_number):
                print (" {:12}".format(mu.atom_id(at)))
                r = at.get_parent()
                self.met_rids.append(mu.residue_num(r))
                if not at.id in self.at_groups:
                    self.at_groups[at.id] = []
                self.at_groups[at.id].append(at)
                self.summary['metals']['detected'].append(mu.residue_num(r))
            return True
        else:
            if not self.args['quiet']:
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
            sys.stderr.write ('Error: unknown selection {}\n'.format(remove_metals))
            self.summary['metals']['error'] = 'Unknown selection ' + remove_metals
            return 1

        if input_option == 'none':
            if not self.args['quiet']:
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
            if not self.args['quiet']:
                print ("No water molecules found")
            return False

    def remwat_fix(self, remove_wat):
        input_line = ParamInput('Remove', remove_wat, self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_wat] = input_line.run()

        if input_option == 'error':
            print ('Warning: unknown option {}'.format(remove_wat))
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
            self.summary['ligands'] = {'detected': []}

            for r in sorted(self.lig_list, key=lambda x: x.index):
                if self.stm.has_models():
                    if r.get_parent().get_parent().id > 0:
                        continue
                    print (' '+mu.residue_id(r,False)+'/*')
                else:
                    print (' '+mu.residue_id(r,False))
                self.summary['ligands']['detected'].append(mu.residue_id(r, self.stm.has_models()))
                self.ligand_rids.add(r.get_resname())
                self.ligand_rnums.append(mu.residue_num(r))
            return True
        else:
            if not self.args['quiet']:
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
            if not self.args['quiet']:
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
            if not self.args['quiet']:
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
        self.SS_bonds = mu.get_all_at2at_distances(
            self._get_structure(),
            'SG',
            self.data_library.get_distances('SS_DIST'),
            not self.stm.has_superimp_models()
        )

        if len(self.SS_bonds):
            print ('{} Possible SS Bonds detected'.format(len(self.SS_bonds)))
            self.summary['getss'] = []
            for ssb in self.SS_bonds:
                print (' {:12} {:12} {:8.3f}'.format(
                    mu.atom_id(ssb[0]),
                    mu.atom_id(ssb[1]),
                    ssb[2]))
                self.summary['getss'].append({
                    'at1':mu.atom_id(ssb[0]),
                    'at2':mu.atom_id(ssb[1]),
                    'dist': round(float(ssb[2]),4)}
                )
        else:
            if not self.args['quiet']:
                print ("No SS bonds detected")

        return False

    def getss_fix(self):
        pass

# =============================================================================
    def amide(self, opts=None):
        self.run_method('amide', opts)

    def amide_check(self):

        [amide_res, amide_atoms] = self.data_library.get_amide_data()
        contact_types = ['polar_donor','polar_acceptor']
        atom_lists = self.data_library.get_atom_lists(contact_types)

        self.amide_list = set()

        for r in self._get_structure().get_residues():
            if r.get_resname() in amide_res:
                self.amide_list.add(r)

        self.summary['amide']['n_amides'] = len(self.amide_list)

        if len(self.amide_list) > 0:
            if len(self.rr_dist) == 0:
                self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))
            self.amide_res_to_fix = []
            self.amide_cont_list = []
            self.amide_rnums = []

            c_list = mu.check_r_list_clashes(
                self.amide_list,
                self.rr_dist,
                self.CLASH_DISTS,
                atom_lists,
                not self.stm.has_superimp_models()
            )
            for cls in c_list:
                for rkey in c_list[cls]:
                    [at1,at2,d]=c_list[cls][rkey]
                    r1 = at1.get_parent()
                    r2 = at2.get_parent()
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
                        self.amide_cont_list.append(c_list[cls][rkey])

            if len(self.amide_cont_list):
                print ('{} unusual contact(s) involving amide atoms found'.format(len(self.amide_cont_list)))
                self.summary['amide']['detected'] = []
                for at_pair in sorted(self.amide_cont_list, key=lambda x:  _key_sort_atom_pairs(x)):
                    print (' {:12} {:12} {:8.3f} A'.format(
                        mu.atom_id(at_pair[0]),
                        mu.atom_id(at_pair[1]),
                        np.sqrt(at_pair[2]))
                    )
                    self.summary['amide']['detected'].append({
                        'at1':mu.atom_id(at_pair[0]),
                        'at2':mu.atom_id(at_pair[1]),
                        'dist': round(float(np.sqrt(at_pair[2])),4)}
                    )
                return True
            else:
                if not self.args['quiet']:
                    print ("No unusual contact(s) involving amide atoms found")
                return False
        else:
            if not self.args['quiet']:
                print ("No amide residues found")
            return False

# =============================================================================
    def amide_fix(self, amide_fix, recheck=True):
        amide_res = self.data_library.get_amide_data()[0]

        no_int_recheck = amide_fix is not None or self.args['non_interactive']

        fix=True
        while fix:
            input_line = ParamInput('Fix amide atoms', amide_fix, self.args['non_interactive'])
            input_line.add_option_all()
            input_line.add_option_none()
            input_line.add_option('resnum', sorted(self.amide_rnums), case='sensitive', multiple=True)
            [input_option, amide_fix] = input_line.run()

            if input_option == 'error':
                print ('Warning: unknown option {}'.format(amide_fix))
                self.summary['amide']['error'] = 'Unknown option'
                return 1

            if input_option == 'none':
                if not self.args['quiet']:
                    print ("Nothing to do")
                fix=False
            else:
                if input_option == 'all':
                    to_fix = self.amide_res_to_fix
                else:
                    to_fix = []
                    for r in self.amide_res_to_fix:
                        if mu.residue_num(r) in amide_fix.split(','):
                            to_fix.append(r)
                n = 0
                done=[]
                for r in to_fix:
                    if r not in done:
                        mu.invert_side_atoms(r, amide_res)
                        done.append(r)
                    n += 1
                print ('Amide residues fixed {} ({})'.format(amide_fix, n))
                fix = False
                if recheck:
                    if not self.args['quiet']:
                        print ("Rechecking")
                    fix=self.amide_check() #TODO reduce check to fixed residues if necessary
                    amide_fix=''
                    if no_int_recheck:
                        fix=False
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
                    print (' {:10}'.format(mu.residue_id(r)))
                    self.summary['chiral']['detected'].append(mu.residue_id(r))
                return True
            else:
                if not self.args['quiet']:
                       print ("No residues with incorrect side-chain chirality found")
                return False
        else:
            if not self.args['quiet']:
                print ("No chiral side-chains found")
            return False

    def chiral_fix(self, chiral_fix, check_clashes=True):

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
            if not self.args['quiet']:
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

            if check_clashes:
                if not self.args['quiet']:
                    print ("Checking for steric clashes")

                contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']
                atom_lists = self.data_library.get_atom_lists(contact_types)

                if len(self.rr_dist) == 0:
                    self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))

                self.summary['chiral_clashes'] = self._clash_report(
                    contact_types,
                    mu.check_r_list_clashes(to_fix, self.rr_dist, self.CLASH_DISTS, atom_lists)
                )

            self.stm.modified = True

# =============================================================================
    def chiral_bck(self, opts=None):
        self.run_method('chiral_bck', opts)

    def chiral_bck_check(self):

        self.chiral_bck_list = []
        prot_chains=0
        for ch in self._get_structure().get_chains():
            if self.stm.chain_ids[ch.id] == mu.PROTEIN:
                prot_chains += 1
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
                    print (' {:10}'.format(mu.residue_id(r)))
                    self.summary['chiral_bck']['detected'].append(mu.residue_id(r))
                return False
            else:
                if not self.args['quiet']:
                    print ("No residues with incorrect backbone chirality found")
                return False
        else:
            if not prot_chains:
                print ("No protein chains detected, skipping")
            else:
                print ("Protein chains detected but no residues with chiral CA found (weird!!)")
                if not self.args['quiet']:
                    print ("No residues with chiral CA found (weird!!)")
            return False


    def chiral_bck_fix(self, chiral_fix, check_clashes=True):
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
            if not self.args['quiet']:
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
            if check_clashes:
                if not self.args['quiet']:
                    print ("Checking new steric clashes")
                #TODO

            self.stm.modified = True

# =============================================================================
    def clashes(self, opts=None):
        self.run_method('clashes', opts)

    def clashes_check(self):
        contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']

        self.clash_list={}
        atom_lists = self.data_library.get_atom_lists(contact_types)
        for cls in contact_types:
            self.clash_list[cls]={}

        if len(self.rr_dist) == 0:
            self.rr_dist = mu.get_all_r2r_distances(
                self._get_structure(),
                'all',
                self.data_library.get_distances('R_R_CUTOFF'),
                not self.stm.has_superimp_models()
            )

        for r_pair in self.rr_dist:
            [r1, r2, d] = r_pair
            if mu.is_wat(r1) or mu.is_wat(r2):
                continue
            c_list = mu.check_rr_clashes(r1,r2,self.CLASH_DISTS, atom_lists)
            rkey = mu.residue_id(r1)+'-'+mu.residue_id(r2)
            for cls in c_list:
                if len(c_list[cls]):
                    self.clash_list[cls][rkey] = c_list[cls]

        self.summary['clashes']['detected'] = self._clash_report(contact_types, self.clash_list)

        return False

    def clashes_fix(self, res_to_fix):
        pass
# =============================================================================
    def fixside (self, opts=None):
        self.run_method('fixside', opts)

    def fixside_check (self):
        self.miss_at_list = self.stm.get_missing_side_chain_atoms(self.data_library.get_valid_codes('protein'), self.data_library.get_all_atom_lists())

        if len(self.miss_at_list):
            self.fixside_rnums = []
            self.summary['fixside']['detected'] = {}
            print('{} Residues with missing side chain atoms found'.format(len(self.miss_at_list)))
            for r_at in self.miss_at_list:
                [r, at_list] = r_at
                print (' {:10} {}'.format(mu.residue_id(r), ','.join(at_list)))
                self.fixside_rnums.append(mu.residue_num(r))
                self.summary['fixside']['detected'][mu.residue_id(r)] = at_list
            return True
        else:
            if not self.args['quiet']:
                print ("No residues with missing side chain atoms found")
            return False

    def fixside_fix(self, fix_side, check_clashes=True):
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
            if not self.args['quiet']:
                print ('Nothing to do')
        else:
            if input_option == 'all':
                to_fix = self.miss_at_list
            else:
                to_fix = []
                for r_at in self.miss_at_list:
                    if mu.residue_num(r_at[0]) in fix_side.split(','):
                        to_fix.append(r_at)
            if not self.args['quiet']:
                print ("Fixing side chains")
            n = 0
            self.summary['fixside']['fixed'] = []
            if not hasattr(self,'residue_lib'):
                self.residue_lib = ResidueLib(self.sets.res_library_path)
            fixed_res=[]
            for r_at in to_fix:
                mu.remove_H_from_r(r_at[0], verbose= True)
                self.stm.fix_side_chain(r_at, self.residue_lib)
                n += 1
                self.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])

            print ('Fixed {} side chain(s)'.format(n))
            # Checking new clashes
            if check_clashes:
                print ("Checking possible new clashes")
                contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']
                atom_lists= self.data_library.get_atom_lists(contact_types)

                if len(self.rr_dist) == 0:
                    self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))

                self.summary['fixside_clashes'] = self._clash_report(
                    contact_types,
                    mu.check_r_list_clashes(fixed_res, self.rr_dist, self.CLASH_DISTS, atom_lists)
                )
# =============================================================================
    def addH(self, opts=None):
        self.run_method('addH', opts)

    def addH_check(self):

        self.ion_res_list = self.stm.get_ion_res_list(self.data_library.ion_res, self.data_library.get_hydrogen_atoms())

        if len(self.ion_res_list):
            self.addH_rnums = []
            self.summary['addH']['detected'] = {}
            print('{} Residues requiring selection on adding H'.format(len(self.ion_res_list)))
            for r_at in self.ion_res_list:
                [r, at_list] = r_at
                #print (' {:10} {}'.format(mu.residue_id(r), ','.join(at_list)))
                self.addH_rnums.append(mu.residue_num(r))
                self.summary['addH']['detected'][mu.residue_id(r)] = at_list
            return True
        else:
            if not self.args['quiet']:
                print ("No residues requiring selection on adding H")
            return True

    def addH_fix(self,mode):
        input_line = ParamInput ('Mode', mode, self.args['non_interactive'])
        input_line.add_option_none()
        input_line.add_option('auto', ['auto'], case="lower")
        input_line.add_option('inter', ['interactive'], case="lower")
        input_line.add_option('inter_His', ['interactive_his'], case="lower")
        input_line.add_option('ph', ['pH'])
        #input_line.add_option('resnum', sorted(self.addH_rnums), case='sensitive', multiple=True)
        [input_option, addH_mode] = input_line.run()

        if input_option == 'error':
            print ("Invalid option", addH_mode)
            self.summary['addH']['error'] = "Unknown option"
            return 1
        if input_option== "none":
            if not self.args['quiet']:
                print ('Nothing to do')
        else:
            if not hasattr(self,'residue_lib'):
                self.residue_lib = ResidueLib(self.sets.res_library_path)
            to_fix=[]
            std_ion= self.data_library.get_ion_data()
            if input_option == 'auto':
                if not self.args['quiet']:
                    print ('Selection: auto')
                for r_at in self.ion_res_list:
                    [r, at_list] = r_at
                    rcode=r.get_resname()
                    to_fix.append([r,std_ion[rcode]['std']])
            elif input_option == 'ph':
                ph_value=None
                input_line = ParamInput("pH Value", ph_value,self.args['non_interactive'])
                input_line.add_option("pH", [], opt_type="float", min_val=0., max_val=14.)
                [input_option, ph_value] = input_line.run()
                if not self.args['quiet']:
                    print ('Selection: pH',ph_value)
                for r_at in self.ion_res_list:
                    [r, at_list] = r_at
                    rcode=r.get_resname()
                    if ph_value <= std_ion[rcode]['pK']:
                        to_fix.append([r,std_ion[rcode]['lowpH']])
                    else:
                        to_fix.append([r,std_ion[rcode]['highpH']])
            elif input_option == 'interactive':
                if not self.args['quiet']:
                    print ('Selection: interactive')
                pass
            elif input_option == 'interactive_His':
                if not self.args['quiet']:
                    print ('Selection: interactive-his')
                pass
            else:
                print ("Not Valid")
                return 1
            backbone_atoms = self.data_library.get_all_atom_lists()['GLY']['backbone']
            addH_rules = self.data_library.get_addH_rules()
            self.stm.add_hydrogens(to_fix,self.data_library.get_hydrogen_atoms(), self.residue_lib, backbone_atoms, self.data_library.get_distances('COVLNK'), addH_rules)

        return ""



# =============================================================================
    def mutateside(self, mut_list):
        self.run_method('mutateside', mut_list)

    def mutateside_check(self):
        return True

    def mutateside_fix(self, mut_list, check_clashes=True):
        input_line = ParamInput ('Mutation list', mut_list, self.args['non_interactive'])
        mut_list = input_line.run()
        self.mutations = MutationManager(mut_list)

        self.mutations.prepare_mutations(self.stm.st)
        print ('Mutations to perform')
        for m in self.mutations.mutation_list:
            print (m)
        if not hasattr(self,'mutation_rules'):
            self.mutation_rules = self.data_library.get_mutation_map()

        if not hasattr(self,'residue_lib'):
            self.residue_lib = ResidueLib(self.sets.res_library_path)

        mutated_res = self.mutations.apply_mutations (self.stm.st, self.mutation_rules, self.residue_lib)
        self.stm.atom_renumbering()
        if check_clashes:
            if not self.args['quiet']:
                print ("Checking new clashes")

            contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']
            atom_lists= self.data_library.get_atom_lists(contact_types)

            if len(self.rr_dist) == 0:
                self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))

            self.summary['mutateside_clashes'] = self._clash_report(
                contact_types,
                mu.check_r_list_clashes(mutated_res, self.rr_dist, self.CLASH_DISTS, atom_lists)
            )

        self.stm.modified=True
#===============================================================================
    def backbone(self, opts):
        self.run_method('backbone',opts)

    def backbone_check(self):
        backbone_atoms = self.data_library.get_all_atom_lists()['GLY']['backbone']
        # Residues with missing backbone
        self.miss_bck_at_list = self.stm.get_missing_backbone_atoms(
            self.data_library.get_valid_codes('protein'),
            self.data_library.get_all_atom_lists()
        )
        if len(self.miss_bck_at_list):
            self.summary['backbone']['missing_atoms'] = {}
            self.fixbck_rnums=[]
            print('{} Residues with missing backbone atoms found'.format(len(self.miss_bck_at_list)))
            for r_at in self.miss_bck_at_list:
                [r, at_list] = r_at
                print (' {:10} {}'.format(mu.residue_id(r), ','.join(at_list)))
                self.fixbck_rnums.append(mu.residue_num(r))
                self.summary['backbone']['missing_atoms'][mu.residue_id(r)] = at_list

        self.stm.check_backbone_connect(backbone_atoms, self.data_library.get_distances('COVLNK'))

        #Not bound consecutive residues
        self.stm.get_backbone_breaks()
        if self.stm.bck_breaks_list:
            print ("{} Backbone breaks found".format(len(self.stm.bck_breaks_list)))
            self.summary['backbone']['breaks']=[]
            for br in self.stm.bck_breaks_list:
                print (" {:10} - {:10} ".format(mu.residue_id(br[0]),mu.residue_id(br[1])))
                self.summary['backbone']['breaks'].append([mu.residue_id(br[0]),mu.residue_id(br[1])])

        if self.stm.wrong_link_list:
            print ("Unexpected backbone links found")
            self.summary['backbone']['wrong_links']=[]
            for br in self.stm.wrong_link_list:
                print (" {:10} linked to {:10}, expected {:10} ".format(
                    mu.residue_id(br[0]),
                    mu.residue_id(br[1]),
                    mu.residue_id(br[2]))
                )
                self.summary['backbone']['wrong_links'].append([mu.residue_id(br[0]),mu.residue_id(br[1]),mu.residue_id(br[2])])

        if self.stm.not_link_seq_list:
            print ("Consecutive residues too far away to be covalently linked")
            self.summary['backbone']['long_links']=[]
            for br in self.stm.not_link_seq_list:
                print (" {:10} - {:10}, bond distance {:8.3f} ".format(
                    mu.residue_id(br[0]),
                    mu.residue_id(br[1]),
                    br[2])
                )
                self.summary['backbone']['long_links'].append(
                    [mu.residue_id(br[0]),
                    mu.residue_id(br[1]),
                    round(float(br[2]),5)]
                )
        #TODO move this section to ligands
        if self.stm.modified_residue_list:
            print ("Modified residues found")
            self.summary['backbone']['mod_residues']=[]
            for br in self.stm.modified_residue_list:
                print (" {:10} ".format(mu.residue_id(br)))
                self.summary['backbone']['mod_residues'].append(mu.residue_id(br))
        #Provisional only missing atoms can be fixed
        if len(self.miss_bck_at_list):
            return True
        else:
            return False

    def backbone_fix(self, fix_back, check_clashes=True):
        input_line = ParamInput ('fixbck', fix_back, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', self.fixbck_rnums, case='sensitive', multiple=True)
        [input_option, fix_back] = input_line.run()

        if input_option == 'error':
            print ("Invalid option", fix_back)
            self.summary['backbone']['missing_atoms']['error'] = "Unknown option"
            return 1

        self.summary['backbone']['missing_atoms']['selected'] = fix_back

        if input_option == 'none':
            if not self.args['quiet']:
                print ('Nothing to do')
        else:
            if input_option == 'all':
                to_fix = self.miss_bck_at_list
            else:
                to_fix = []
                for r_at in self.miss_bck_at_list:
                    if mu.residue_num(r_at[0]) in fix_back.split(','):
                        to_fix.append(r_at)
            if not self.args['quiet']:
                print ("Adding missing backbone atoms")
            n = 0
            self.summary['backbone']['missing_atoms']['fixed'] = []
            fixed_res=[]
            for r_at in to_fix:
                if self.stm.fix_backbone_atoms(r_at):
                    n += 1
                self.summary['backbone']['missing_atoms']['fixed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])

            print ('Fixed {} backbone atom(s)'.format(n))
            # Checking new clashes
            if check_clashes:
                print ("Checking possible new clashes")
                contact_types = ['severe','apolar','polar_acceptor','polar_donor', 'positive', 'negative']
                atom_lists= self.data_library.get_atom_lists(contact_types)

                if len(self.rr_dist) == 0:
                    self.rr_dist = mu.get_all_r2r_distances(self._get_structure(), 'all', self.data_library.get_distances('R_R_CUTOFF'))

                self.summary['backbone']['missing_atoms']['clashes'] = self._clash_report(
                    contact_types,
                    mu.check_r_list_clashes(fixed_res, self.rr_dist, self.CLASH_DISTS, atom_lists)
                )

        #TODO Chain fix

        pass
#===============================================================================
    def cistransbck(self, opts):
        self.run_method('cistransbck',opts)

    def cistransbck_check(self):
        self.stm.check_cis_backbone(self.data_library.get_distances('COVLNK'), not self.stm.biounit)
        if self.stm.cis_backbone_list:
            self.summary['cistransbck']['cis']=[]
            print ('{} cis peptide bonds'.format(len(self.stm.cis_backbone_list)))
            for lnk in self.stm.cis_backbone_list:
                [r1,r2,dih] = lnk
                print ('{:10} {:10} Dihedral: {:8.3f}'.format(mu.residue_id(r1), mu.residue_id(r2), dih))
                self.summary['cistransbck']['cis'].append([mu.residue_id(r1), mu.residue_id(r2), round(float(dih),3)])
        else:
            if not self.args['quiet']:
                print ("No cis peptide bonds found")

        if self.stm.lowtrans_backbone_list:
            self.summary['cistransbck']['unusual_trans']=[]
            print ('{} trans peptide bonds with unusual omega dihedrals'.format(len(self.stm.lowtrans_backbone_list)))
            for lnk in self.stm.lowtrans_backbone_list:
                [r1,r2,dih] = lnk
                print ('{:10} {:10} Dihedral: {:8.3f}'.format(mu.residue_id(r1), mu.residue_id(r2), dih))
                self.summary['cistransbck']['unusual_trans'].append([mu.residue_id(r1), mu.residue_id(r2), round(float(dih),3)])
        else:
            if not self.args['quiet']:
                print ("No trans peptide bonds with unusual omega dihedrals found")

    def cistransbck_fix(self,option):
        pass
#===============================================================================

    def _load_structure(self, verbose=True, print_stats=True):
        if not hasattr(self, 'stm'):
            if not self.args['non_interactive'] and self.args['input_structure_path'] is None:
                self.args['input_structure_path'] = input("Enter input structure path (PDB, mmcif | pdb:pdbid): ")
            self.stm = StructureManager(
                self.args['input_structure_path'],
                pdb_server=self.pdb_server,
                cache_dir=self.cache_dir,
                file_format=self.args['file_format']
            )
            if verbose:
                print ('Structure {} loaded'.format(self.args['input_structure_path']))
                self.stm.print_headers()
                print()
                self.summary['headers'] = self.stm.meta
            if print_stats:
                self.stm.print_stats()
                print()
            self.summary['stats'] = self.stm.get_stats()

            backbone_atoms = self.data_library.get_all_atom_lists()['GLY']['backbone']
            self.stm.check_backbone_connect(backbone_atoms, self.data_library.get_distances('COVLNK'))

    def _get_structure(self):
        return self.stm.get_structure()

    def _save_structure(self):
        if not self.args['non_interactive'] and self.args['output_structure_path'] is None:
            self.args['output_structure_path'] = input("Enter output structure path: ")
        self.stm.save_structure(self.args['output_structure_path'])

    def _clash_report(self, contact_types, clash_list):
        summary={}
        for cls in contact_types:
            summary[cls]=[]
            if len(clash_list[cls]):
                print ('{} Steric {} clashes detected'.format(len(clash_list[cls]), cls))
                for rkey in sorted(clash_list[cls], key=lambda x: _key_sort_atom_pairs(clash_list[cls][x])):
                    print (' {:12} {:12} {:8.3f} A'.format(
                        mu.atom_id(clash_list[cls][rkey][0]),
                        mu.atom_id(clash_list[cls][rkey][1]),
                        np.sqrt(clash_list[cls][rkey][2])
                        )
                    )
                    summary[cls].append(
                        {
                            'at1':mu.atom_id(clash_list[cls][rkey][0]),
                            'at2':mu.atom_id(clash_list[cls][rkey][1]),
                            'dist': round(float(np.sqrt(clash_list[cls][rkey][2])),4)
                        }
                    )
            else:
                if not self.args['quiet']:
                    print ('No {} clashes detected'.format(cls))
        return summary

#===============================================================================
def _key_sort_atom_pairs(at_pair):
    return at_pair[0].serial_number * 10000 + at_pair[1].serial_number
