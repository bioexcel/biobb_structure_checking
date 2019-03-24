"""
    Main Class for Structure Checking functionality

"""

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys
import numpy as np

from structure_checking.json_writer import JSONWriter
from structure_checking.param_input import Dialog, ParamInput
from structure_manager.structure_manager import StructureManager
from structure_manager.data_lib_manager import DataLibManager
from structure_manager.residue_lib_manager import ResidueLib
import structure_manager.model_utils as mu

# Interactive DIALOGS to complete command_line missing parameters
DIALOGS = Dialog()

#DIALOGS.add_option(command, prompt, destinmore ation, help_text, type(str))
DIALOGS.add_option('command_list', '--list', 'op_list', 'Command List File')
DIALOGS.add_option('models', '--select_model', 'select_model', \
    'Select model to keep', int)
DIALOGS.add_option('chains', '--select_chains', 'select_chains',\
    'Chains (All | Chain list comma separated)')
DIALOGS.add_option('altloc', '--select_altloc', 'select_altloc', \
    'Select altloc occupancy|alt_id')
DIALOGS.add_option('metals', '--remove', 'remove_metals', 'Remove Metal ions')
DIALOGS.add_option('remwat', '--remove', 'remove_wat', 'Remove Water molecules')
DIALOGS.add_option('ligands', '--remove', 'remove_ligands', 'Remove Ligand residues')
DIALOGS.add_option('remh', '--remove', 'remove_h', 'Remove Hydrogen atoms')
DIALOGS.add_option('amide', '--fix', 'amide_fix', 'Fix Residues (All | None | List)')
DIALOGS.add_option('chiral', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
DIALOGS.add_option('chiral_bck', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')
DIALOGS.add_option('clashes', '--no_wat', 'discard_wat', 'Discard water molecules')
DIALOGS.add_option('fixside', '--fix', 'fix_side',\
    'Add missing atoms to side chains (All | None | List)')
DIALOGS.add_option('backbone', '--fix', 'fix_back',\
    'Add missing O atoms to backbone (All | None | List)')
DIALOGS.add_option('mutateside', '--mut', 'mut_list',\
    'Mutate side chains (Mutation List as [*:]arg234Thr)')
DIALOGS.add_option('addH', '--mode', 'mode',\
    'Selection mode (None | auto | interactive | interactive_his | ph )')

AVAILABLE_METHODS = [
    'models', 'chains', 'inscodes', 'altloc', 'remh', 'addH', 'remwat', 'metals', 'ligands',
    'getss', 'amide', 'chiral', 'chiral_bck', 'fixside', 'backbone', 'cistransbck',
    'clashes']

ALL_CONTACT_TYPES = [
    'severe',
    'apolar',
    'polar_acceptor',
    'polar_donor',
    'positive',
    'negative'
]
AMIDE_CONTACT_TYPES = [
    'polar_acceptor',
    'polar_donor',
]

# Main class
class StructureChecking():
    """Main class to control structure checking front end"""
    def __init__(self, sets, args):
        self.args = args
        self.sets = sets
        self.summary = {}

        self.pdb_server = self.args['pdb_server']
        self.cache_dir = self.args['cache_dir']

        self.residue_lib = ResidueLib(self.sets.res_library_path)

        self.stm = self._load_structure()

        if not 'Notebook' in self.args:
            self.args['Notebook'] = False

        if self.args['Notebook']:
            self.args['non_interactive'] = True
            self.args['check_only'] = False


    def launch(self):
        """main method to run checking"""
        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        elif self.args['command'] == 'checkall':
            self.checkall(self.args['options'])
        elif self.args['command'] == 'load':
            sys.exit(0)
        else:
            self._run_method(self.args['command'], self.args['options'])

        if not self.args['check_only']:
            if self.stm.modified or self.args['force_save']:
                if not self.stm.modified:
                    print('Structure not modified, saving due to --force_save option')
                self._save_structure()
                self.stm.calc_stats()
                print('Structure saved on {}'.format(self.args['output_structure_path']))
                self.stm.print_stats('Final')
                self.summary['final_stats'] = self.stm.get_stats()
            elif not self.stm.modified:
                print('Structure not modified, not saving. Override with --force_save')

        if self.args['json_output_path'] is not None:
            json_writer = JSONWriter()
            for k in self.summary:
                json_writer.set(k, self.summary[k])
            json_writer.save(self.args['json_output_path'])
            print('Summary data saved on {}'.format(self.args['json_output_path']))

    def command_list(self, opts):
        """manages command_list worlflows"""
        opts = DIALOGS.get_parameter('command_list', opts)
        op_list = opts[DIALOGS.get_dialog('command_list')['dest']]

        op_list = ParamInput('Command List File', op_list, False).run()

        try:
            list_file_h = open(op_list, "r")

        except OSError:
            print('Error when opening file {}'.format(op_list))
            sys.exit(1)

        #self._load_structure()

        i = 1
        for line in list_file_h:
            if line == "\n" or line[0:1] == '#':
                continue
            if not self.args['quiet']:
                print("\nStep {}: {}".format(i, line))
            data = line.split()
            command = data[0]
            opts = data[1:]
            self._run_method(command, opts)
            i += 1

        print("Command list completed")

    def checkall(self, opts):
        """Predefined workflow for complete checking"""
        self.args['check_only'] = True
        #self._load_structure()
        for meth in AVAILABLE_METHODS:
            self._run_method(meth, opts)

    def _run_method(self, command, opts):
        """run check and fix methods for specific command"""
        try:
            f_check = getattr(self, '_' + command + '_check')
        except AttributeError:
            print('Error: {} command unknown or not implemented'.format(command))
            sys.exit(1)

        if not command in self.summary:
            self.summary[command] = {}

        msg = 'Running {}.'.format(command)

        if opts:
            self.summary[command]['opts'] = opts
            msg += ' Options: {} '.format(' '.join(opts))

        if not self.args['quiet']:
            print(msg)
        # self._load_structure()

    #Running checking method
        fix_data = f_check()
    #Running fix method if needed
        if self.args['check_only'] or opts is None or opts == '':
            if not self.args['quiet']:
                print('Running  check_only. Nothing else to do.')
        elif fix_data:
            try:
                f_fix = getattr(self, '_' + command + '_fix')
            except AttributeError:
                print('Error: {}_fix command unknown or not implemented'.format(command))
                sys.exit(1)
            if not self.args['Notebook']:
                if DIALOGS.exists(command):
                    opts = DIALOGS.get_parameter(command, opts)
                    opts = opts[DIALOGS.get_dialog(command)['dest']]
                else:
                    opts = ''
            f_fix(opts, fix_data)

# ==============================================================================q
    def models(self, opts=None):
        """ run models command """
        self._run_method('models', opts)

    def _models_check(self):
        fix_data = False
        print('{} Model(s) detected'.format(self.stm.nmodels))
        self.summary['models'] = {'nmodels': self.stm.nmodels}
        if self.stm.nmodels > 1:
            self.summary['models']['type'] = self.stm.models_type
            supimp = ''
            if self.stm.models_type['type'] == mu.BUNIT:
                supimp = 'do not'
            print(
                'Models {} superimpose, RMSd: {:8.3f} A, guessed as {} '.format(
                    supimp,
                    self.stm.models_type['rmsd'],
                    mu.MODEL_TYPE_LABELS[self.stm.models_type['type']]
                )
            )

            fix_data = True
        else:
            if not self.args['quiet']:
                print("Single model found")
        return fix_data

    def _models_fix(self, select_model, fix_data=None):
        input_line = ParamInput('Select Model Num', select_model, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option('modelno', [], opt_type='int', min_val=1, max_val=self.stm.nmodels)
        [input_option, select_model] = input_line.run()

        if input_option == 'error':
            print('Error: unknown model {}'.format(select_model), file=sys.stderr)
            self.summary['models']['error'] = 'Unknown model ' + select_model
            return

        print('Selecting model num. {}'.format(select_model))
        if input_option != 'all':
            self.stm.select_model(select_model)

        self.summary['models']['selected'] = select_model

# =============================================================================
    def chains(self, opts=None):
        """ run chains command """
        self._run_method('chains', opts)

    def _chains_check(self):
        print('{} Chain(s) detected'.format(len(self.stm.chain_ids)))
        for ch_id in sorted(self.stm.chain_ids):
            if isinstance(self.stm.chain_ids[ch_id], list):
                print(
                    ' {}: Unknown (PROTEIN: {s[0]:4.2f} DNA: {s[1]:4.2f} RNA: {s[2]:4.2f} Other: {s[3]:4.2f})'.format(
                        ch_id, s=self.stm.chain_ids[ch_id]
                    )
                )
            else:
                print(' {}: {}'.format(ch_id, mu.CHAIN_TYPE_LABELS[self.stm.chain_ids[ch_id]]))
        self.summary['chains'] = {'ids':self.stm.chain_ids}
        return len(self.stm.chains_ids) > 1

    def _chains_fix(self, select_chains, fix_data=None):
        self.summary['chains']['selected'] = {}
        input_line = ParamInput('Select chain', select_chains, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option('chid', sorted(self.stm.chain_ids), multiple=True, case="sensitive")
        [input_option, select_chains] = input_line.run()

        if input_option == 'error':
            print('Error unknown selection: {}'.format(select_chains))
            self.summary['chains']['error'] = 'Unknown selection ' + select_chains
            return

        if input_option == 'all':
            print('Selecting all chains')
        else:
            self.stm.select_chains(select_chains)
            print('Selecting chain(s) {}'.format(select_chains))
            self.stm.set_chain_ids()
            self.summary['chains']['selected'] = self.stm.chain_ids

# =============================================================================
    def inscodes(self, opts=None):
        """ run inscodes command """
        self._run_method('inscodes', opts)

    def _inscodes_check(self):
        self.stm.get_ins_codes()
        if self.stm.ins_codes_list:
            print('{} Residues with insertion codes found'.format(len(self.stm.ins_codes_list)))
            self.summary['inscodes'] = []
            for res in self.stm.ins_codes_list:
                print(mu.residue_id(res))
                self.summary['inscodes'].append(mu.residue_id(res))
        else:
            if not self.args['quiet']:
                print("No residues with insertion codes found")
        return {}

    def _inscodes_fix(self, select_codes, fix_data=None):
        pass
# =============================================================================
    def altloc(self, opts=None):
        """ run altloc command """
        self._run_method('altloc', opts)

    def _altloc_check(self): #TODO improve output
        fix_data = {}
        alt_loc_res = mu.get_altloc_residues(self._get_structure())
        if alt_loc_res:
            fix_data['alt_loc_res'] = alt_loc_res
            print(
                'Detected {} residues with alternative location labels'.format(
                    len(alt_loc_res)
                )
            )

            self.summary['altloc'] = {}
            fix_data['alt_loc_rnums'] = []
            fix_data['altlocs'] = {}
            for res in sorted(alt_loc_res, key=lambda x: x.index):
                rid = mu.residue_id(res)
                print(rid)
                fix_data['alt_loc_rnums'].append(mu.residue_num(res))
                self.summary['altloc'][rid] = {}
                fix_data['altlocs'][res] = sorted(alt_loc_res[res][0].child_dict)
                for atm in alt_loc_res[res]:
                    self.summary['altloc'][rid][atm.id] = []
                    alt_str = '  {:4}'.format(atm.id)
                    for alt in sorted(atm.child_dict):
                        alt_str += ' {} ({:4.2f})'.format(alt, atm.child_dict[alt].occupancy)
                        self.summary['altloc'][rid][atm.id].append({
                            'loc_label':alt,
                            'occupancy':atm.child_dict[alt].occupancy
                        })
                    print(alt_str)
        else:
            if not self.args['quiet']:
                print("No residues with alternative location labels detected")
        return fix_data

    def _altloc_fix(self, select_altloc, fix_data=None):
        #Prepare the longest possible list of alternatives
        altlocs = []
        max_al_len = 0
        for res in fix_data['altlocs']:
            if len(fix_data['altlocs'][res]) > max_al_len:
                altlocs = fix_data['altlocs'][res]
                max_al_len = len(fix_data['altlocs'][res])

        input_line = ParamInput('Select alternative', select_altloc, self.args['non_interactive'])
        input_line.add_option('occup', ['occupancy'])
        input_line.add_option('altids', altlocs, case='upper')
        input_line.add_option(
            'resnum',
            fix_data['alt_loc_rnums'],
            opt_type='pair_list',
            list2=altlocs,
            case='sensitive',
            multiple=True
        )
        [input_option, select_altloc] = input_line.run()

        if input_option == 'error':
            print('Error: Unknown selection {} '.format(select_altloc), file=sys.stderr)
            self.summary['altloc']['error'] = "Unknown selection " + select_altloc
            return

        print('Selecting location {}'.format(select_altloc))
        if input_option in ('occup', 'altids'):
            to_fix = {}
            for res in fix_data['alt_loc_res']:
                to_fix[res] = {}
                to_fix[res]['ats'] = fix_data['alt_loc_res'][res]
                to_fix[res]['select'] = select_altloc
        elif input_option == 'resnum':
            to_fix = {}
            selected_rnums = {}
            for rsel in select_altloc.split(','):
                [rnum, alt] = rsel.split(':')
                selected_rnums[rnum] = alt
            for res in fix_data['alt_loc_res']:
                rnum0 = mu.residue_num(res)
                if rnum0 in selected_rnums:
                    to_fix[res] = {}
                    to_fix[res]['ats'] = fix_data['alt_loc_res'][res]
                    to_fix[res]['select'] = selected_rnums[rnum0]

        self.summary['altloc']['selected'] = select_altloc
        for res in to_fix:
            self.stm.select_altloc_residue(res, to_fix[res])

# =============================================================================
    def metals(self, opts=None):
        """ run metals command """
        self._run_method('metals', opts)

    def _metals_check(self):
        fix_data = {}

        met_list = self.stm.get_metal_atoms()

        if len(met_list) > 1:
            fix_data['met_list'] = met_list
            print('{} Metal ions found'.format(len(met_list)))
            self.summary['metals'] = {'detected':[]}
            fix_data['met_rids'] = []
            fix_data['at_groups'] = {}
            for atm in sorted(met_list, key=lambda x: x.serial_number):
                print(" {:12}".format(mu.atom_id(atm)))
                res = atm.get_parent()
                fix_data['met_rids'].append(mu.residue_num(res))
                if not atm.id in fix_data['at_groups']:
                    fix_data['at_groups'][atm.id] = []
                fix_data['at_groups'][atm.id].append(atm)
                self.summary['metals']['detected'].append(mu.residue_num(res))
        else:
            if not self.args['quiet']:
                print("No metal ions present")
        return fix_data

    def _metals_fix(self, remove_metals, fix_data=None):
        input_sess = ParamInput("Remove", remove_metals, self.args['non_interactive'])
        input_sess.add_option_all()
        input_sess.add_option_none()
        input_sess.add_option(
            'atids',
            sorted(fix_data['at_groups']),
            case='sensitive',
            multiple=True
        )
        input_sess.add_option('resids', fix_data['met_rids'], case='sensitive', multiple=True)
        [input_option, remove_metals] = input_sess.run()

        if input_option == "error":
            print('Error: unknown selection {}\n'.format(remove_metals), file=sys.stderr)
            self.summary['metals']['error'] = 'Unknown selection ' + remove_metals
            return

        if input_option == 'none':
            if not self.args['quiet']:
                print("Nothing to do")
        else:
            if input_option == 'all':
                to_remove = fix_data['met_list']
            elif input_option == 'resids':
                to_remove = []
                rid_list = remove_metals.split(',')
                for atm in fix_data['met_list']:
                    res = atm.get_parent()
                    if mu.residue_num(res) in rid_list:
                        to_remove.append(atm)
            elif input_option == 'atids':
                to_remove = []
                for atid in remove_metals.split(','):
                    to_remove.extend(fix_data['at_groups'][atid])

            self.summary['metals']['removed'] = []

            rmm_num = 0
            for atm in to_remove:
                self.summary['metals']['removed'].append(
                    mu.residue_id(
                        atm.get_parent(), self.stm.has_models()
                    )
                )
                self.stm.remove_residue(atm.get_parent())
                rmm_num += 1

            print('Metal Atoms removed {} ({:d})'.format(remove_metals, rmm_num))
            self.summary['metals']['n_removed'] = rmm_num

# =============================================================================
    def remwat(self, opts=None):
        """ run remwat command """
        self._run_method('remwat', opts)

    def _remwat_check(self):
        fix_data = {}
        lig_list = mu.get_ligands(self._get_structure(), incl_water=True)
        wat_list = []
        for res in lig_list:
            if mu.is_wat(res):
                wat_list.append(res)

        if wat_list:
            fix_data['wat_list'] = wat_list
            print('{} Water molecules detected'.format(len(wat_list)))
            self.summary['remwat']['n_detected'] = len(wat_list)
        else:
            if not self.args['quiet']:
                print("No water molecules found")
        return fix_data

    def _remwat_fix(self, remove_wat, fix_data=None):
        input_line = ParamInput('Remove', remove_wat, self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_wat] = input_line.run()

        if input_option == 'error':
            print(' unknown option {}'.format(remove_wat), file=sys.stderr)
            self.summary['remwat']['error'] = 'Unknown option'
            return

        if input_option == 'yes':
            rmw_num = 0
            for res in fix_data['wat_list']:
                self.stm.remove_residue(res)
                rmw_num += 1
            print('{} Water molecules removed'.format(rmw_num))
            self.summary['remwat']['n_removed'] = rmw_num

# =============================================================================
    def ligands(self, opts=None):
        """ run ligands command """
        self._run_method('ligands', opts)

    def _ligands_check(self):
        fix_data = {}
        lig_list = mu.get_ligands(self._get_structure(), incl_water=False)

        if lig_list:
            fix_data['lig_list'] = lig_list
            print('{} Ligands detected '.format(len(lig_list)))
            fix_data['ligand_rids'] = set()
            fix_data['ligand_rnums'] = []
            self.summary['ligands'] = {'detected': []}

            for res in sorted(lig_list, key=lambda x: x.index):
                if self.stm.has_models():
                    if res.get_parent().get_parent().id > 0:
                        continue
                    print(' ' + mu.residue_id(res, False) + '/*')
                else:
                    print(' ' + mu.residue_id(res, False))
                self.summary['ligands']['detected'].append(
                    mu.residue_id(res, self.stm.has_models())
                )
                fix_data['ligand_rids'].add(res.get_resname())
                fix_data['ligand_rnums'].append(mu.residue_num(res))
        else:
            if not self.args['quiet']:
                print("No ligands found")
        return fix_data

# =============================================================================
    def _ligands_fix(self, remove_ligands, fix_data=None):
        input_line = ParamInput('Remove', remove_ligands, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('byrids', sorted(fix_data['ligand_rids']), multiple=True)
        input_line.add_option('byresnum', fix_data['ligand_rnums'], case='sensitive', multiple=True)
        [input_option, remove_ligands] = input_line.run()

        if input_option == 'error':
            print('Error: unknown selection {}'.format(remove_ligands))
            self.summary['ligands']['error'] = 'Unknown selection'
            return

        self.summary['ligands']['removed'] = {'opt':remove_ligands, 'lst':[]}

        to_remove = []

        if input_option == 'none':
            if not self.args['quiet']:
                print("Nothing to do")
        else:
            if input_option == 'all':
                to_remove = fix_data['lig_list']
            elif input_option == 'byrids':
                rm_lig = remove_ligands.split(',')
                for res in fix_data['lig_list']:
                    if res.get_resname() in rm_lig:
                        to_remove.append(res)
            elif input_option == 'byresnum':
                rm_lig = remove_ligands.split(',')
                for res in fix_data['lig_list']:
                    if mu.residue_num(res) in rm_lig:
                        to_remove.append(res)
            rl_num = 0
            for res in to_remove:
                self.summary['ligands']['removed']['lst'].append(
                    mu.residue_id(res, self.stm.has_models())
                )
                self.stm.remove_residue(res)
                rl_num += 1

            print('Ligands removed {} ({})'.format(remove_ligands, rl_num))
            self.summary['ligands']['n_removed'] = rl_num

# =============================================================================
    def remh(self, opts=None):
        """ run remh command """
        self._run_method('remh', opts)

    def _remh_check(self):
        fix_data = {}
        remh_list = mu.get_residues_with_H(self._get_structure())

        if remh_list:
            fix_data['remh_list'] = remh_list
            print('{} Residues containing H atoms detected'.format(len(remh_list)))
            self.summary['remh']['n_detected'] = len(remh_list)
        else:
            if not self.args['quiet']:
                print("No residues with Hydrogen atoms found")
        return fix_data

    def _remh_fix(self, remove_h, fix_data=None):
        input_line = ParamInput('Remove hydrogen atoms', remove_h, self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_h] = input_line.run()

        if input_option == 'error':
            print('Warning: unknown option {}'.format(remove_h))
            self.summary['remh']['error'] = 'Unknown option'
            return
        if input_option == 'yes':
            rmh_num = 0
            for resh in fix_data['remh_list']:
                mu.remove_H_from_r(resh['r'])
                rmh_num += 1
            print('Hydrogen atoms removed from {} residues'.format(rmh_num))
            self.stm.modified = True
            self.summary['remh']['n_removed'] = rmh_num

# =============================================================================
    def getss(self, opts=None):
        """ run getss command """
        self._run_method('getss', opts)

    def _getss_check(self):

        SS_bonds = self.stm.get_SS_bonds()

        if SS_bonds:
            print('{} Possible SS Bonds detected'.format(len(SS_bonds)))
            self.summary['getss'] = []
            for ssb in SS_bonds:
                print(
                    ' {:12} {:12} {:8.3f}'.format(
                        mu.atom_id(ssb[0]), mu.atom_id(ssb[1]), ssb[2]
                    )
                )
                self.summary['getss'].append(
                    {
                        'at1':mu.atom_id(ssb[0]),
                        'at2':mu.atom_id(ssb[1]),
                        'dist': round(float(ssb[2]), 4)
                    }
                )
        else:
            if not self.args['quiet']:
                print("No SS bonds detected")

        return {}

    def _getss_fix(self):
        pass

# =============================================================================
    def amide(self, opts=None):
        """ run amide command """
        self._run_method('amide', opts)

    def _amide_check(self):
        fix_data = {}

        amide_atoms = self.stm.data_library.get_amide_data()[1]

        amide_list = self.stm.get_amide_list()

        self.summary['amide']['n_amides'] = len(amide_list)

        if amide_list:
            c_list = self.stm.check_r_list_clashes(
                amide_list,
                AMIDE_CONTACT_TYPES
            )
            amide_res_to_fix = []
            amide_rnums = []
            amide_cont_list = []
            for cls in c_list:
                for rkey in c_list[cls]:
                    [at1, at2] = c_list[cls][rkey][0:2]
                    res1 = at1.get_parent()
                    res2 = at2.get_parent()
                    add_pair = False
                    if at1.id in amide_atoms and res1 in amide_list:
                        amide_res_to_fix.append(res1)
                        amide_rnums.append(mu.residue_num(res1))
                        add_pair = True
                    if at2.id in amide_atoms and res2 in amide_list:
                        amide_res_to_fix.append(res2)
                        amide_rnums.append(mu.residue_num(res2))
                        add_pair = True
                    if add_pair:
                        amide_cont_list.append(c_list[cls][rkey])

            if amide_cont_list:
                print(
                    '{} unusual contact(s) involving amide atoms found'.format(
                        len(amide_cont_list)
                    )
                )
                self.summary['amide']['detected'] = []
                for at_pair in sorted(amide_cont_list, key=_key_sort_atom_pairs):
                    print(
                        ' {:12} {:12} {:8.3f} A'.format(
                            mu.atom_id(at_pair[0]),
                            mu.atom_id(at_pair[1]),
                            np.sqrt(at_pair[2])
                        )
                    )
                    self.summary['amide']['detected'].append(
                        {
                            'at1': mu.atom_id(at_pair[0]),
                            'at2': mu.atom_id(at_pair[1]),
                            'dist': round(float(np.sqrt(at_pair[2])), 4)
                        }
                    )
                fix_data['amide_list'] = amide_list
                fix_data['amide_res_to_fix'] = amide_res_to_fix
                fix_data['amide_cont_list'] = amide_cont_list
                fix_data['amide_rnums'] = amide_rnums
            else:
                if not self.args['quiet']:
                    print("No unusual contact(s) involving amide atoms found")
        else:
            if not self.args['quiet']:
                print("No amide residues found")
        return fix_data

# =============================================================================
    def _amide_fix(self, amide_fix, fix_data=None, recheck=True, ):
        amide_res = self.stm.data_library.get_amide_data()[0]

        no_int_recheck = amide_fix is not None or self.args['non_interactive']

        fix = True
        while fix:
            input_line = ParamInput(
                'Fix amide atoms', amide_fix, self.args['non_interactive']
            )
            input_line.add_option_all()
            input_line.add_option_none()
            input_line.add_option(
                'resnum', sorted(fix_data['amide_rnums']), case='sensitive', multiple=True
            )
            [input_option, amide_fix] = input_line.run()

            if input_option == 'error':
                print('Warning: unknown option {}'.format(amide_fix))
                self.summary['amide']['error'] = 'Unknown option'
                return

            if input_option == 'none':
                if not self.args['quiet']:
                    print("Nothing to do")
                fix = False
            else:
                if input_option == 'all':
                    to_fix = fix_data['amide_res_to_fix']
                else:
                    to_fix = []
                    for res in fix_data['amide_res_to_fix']:
                        if mu.residue_num(res) in amide_fix.split(','):
                            to_fix.append(res)
                fix_num = 0
                done = []
                for res in to_fix:
                    if res not in done:
                        mu.invert_side_atoms(res, amide_res)
                        done.append(res)
                    fix_num += 1
                print('Amide residues fixed {} ({})'.format(amide_fix, fix_num))
                fix = {}
                if recheck:
                    if not self.args['quiet']:
                        print("Rechecking")
                    fix = self._amide_check() #TODO reduce check to fixed residues if necessary
                    amide_fix = ''
                    if no_int_recheck:
                        fix = {}
                self.stm.modified = True

# =============================================================================
    def chiral(self, opts=None):
        """ run chiral command """
        self._run_method('chiral', opts)

    def _chiral_check(self):
        fix_data = {}

        chiral_res = self.stm.data_library.get_chiral_data()

        chiral_list = self.stm.get_chiral_list()

        self.summary['chiral']['n_chirals'] = len(chiral_list)

        if chiral_list:
            chiral_res_to_fix = []
            chiral_rnums = []
            for res in chiral_list:
                if not mu.check_chiral_residue(res, chiral_res):
                    chiral_res_to_fix.append(res)
                    chiral_rnums.append(mu.residue_num(res))
            if chiral_res_to_fix:
                print(
                    '{} residues with incorrect side-chain chirality found'.format(
                        len(chiral_res_to_fix)
                    )
                )
                self.summary['chiral']['detected'] = []
                for res in chiral_res_to_fix:
                    print(' {:10}'.format(mu.residue_id(res)))
                    self.summary['chiral']['detected'].append(mu.residue_id(res))
                fix_data['chiral_res_to_fix'] = chiral_res_to_fix
                fix_data['chiral_rnums'] = chiral_rnums
            else:
                if not self.args['quiet']:
                    print("No residues with incorrect side-chain chirality found")
        else:
            if not self.args['quiet']:
                print("No chiral side-chains found")
        return fix_data

    def _chiral_fix(self, chiral_fix, fix_data=None, check_clashes=True):

        chiral_res = self.stm.data_library.get_chiral_data()

        input_line = ParamInput('Fix chiralities', chiral_fix, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', fix_data['chiral_rnums'], case='sensitive', multiple=True)
        [input_option, chiral_fix] = input_line.run()

        if input_option == 'error':
            print('Warning: unknown option {}'.format(chiral_fix))
            self.summary['chiral']['error'] = 'Unknown option'
            return

        if input_option == 'none':
            if not self.args['quiet']:
                print("Nothing to do")
        else:
            if input_option == 'all':
                to_fix = fix_data['chiral_res_to_fix']
            else:
                to_fix = []
                for res in fix_data['chiral_res_to_fix']:
                    if mu.residue_num(res) in chiral_fix.split(','):
                        to_fix.append(res)
            fix_num = 0
            for res in to_fix:
                mu.invert_side_atoms(res, chiral_res)
                if res.get_resname() == 'ILE':
                    mu.delete_atom(res, 'CD1')
                    mu.build_atom(res, 'CD1', self.residue_lib, 'ILE')
                fix_num += 1
            print('Quiral residues fixed {} ({})'.format(chiral_fix, fix_num))

            if check_clashes:
                if not self.args['quiet']:
                    print("Checking for steric clashes")

                self.summary['chiral_clashes'] = self._clash_report(
                    ALL_CONTACT_TYPES,
                    self.stm.check_r_list_clashes(
                        to_fix,
                        ALL_CONTACT_TYPES
                    )
                )

            self.stm.modified = True

# =============================================================================
    def chiral_bck(self, opts=None):
        """ run chiral_bck command """
        self._run_method('chiral_bck', opts)

    def _chiral_bck_check(self):
        fix_data = {}

        chiral_bck_list = self.stm.get_chiral_bck_list()

        self.summary['chiral_bck']['n_chirals'] = len(chiral_bck_list)

        if chiral_bck_list:
            chiral_bck_res_to_fix = []
            chiral_bck_rnums = []
            for res in chiral_bck_list:
                if not mu.check_chiral_ca(res):
                    chiral_bck_res_to_fix.append(res)
                    chiral_bck_rnums.append(mu.residue_num(res))
            if chiral_bck_res_to_fix:
                print(
                    '{} residues with incorrect backbone chirality found'.format(
                        len(chiral_bck_res_to_fix)
                    )
                )
                self.summary['chiral_bck']['detected'] = []
                for res in chiral_bck_res_to_fix:
                    print(' {:10}'.format(mu.residue_id(res)))
                    self.summary['chiral_bck']['detected'].append(mu.residue_id(res))
                fix_data['chiral_bck_res_to_fix'] = chiral_bck_res_to_fix
                fix_data['chiral_bck_rnums'] = chiral_bck_rnums
            else:
                if not self.args['quiet']:
                    print("No residues with incorrect backbone chirality found")

        return fix_data


    def _chiral_bck_fix(self, chiral_fix, fix_data=None, check_clashes=True):
        return
#TODO chiral_bck_fix
#        input_line = ParamInput('Fix CA chiralities', chiral_fix, self.args['non_interactive'])
#        input_line.add_option_all()
#        input_line.add_option_none()
#        input_line.add_option('resnum', self.chiral_bck_rnums, case='sensitive', multiple=True)
#        [input_option, chiral_fix] = input_line.run()
#        if input_option == 'error':
#            print('Warning: unknown option {}'.format(amide_fix))
#            self.summary['chiral']['error'] = 'Unknown option'
#            return
#
#        if input_option == 'none':
#            if not self.args['quiet']:
#                print("Nothing to do")
#        else:
#            if input_option == 'all':
#                to_fix = self.chiral_bck_res_to_fix
#            else:
#                to_fix = []
#                for res in self.chiral_bck_res_to_fix:
#                    if mu.residue_num(res) in chiral_fix.split(','):
#                        to_fix.append(res)
#            n = 0
#            for res in to_fix:
#                mu.stm.invert_chiral_CA(res)
#                n += 1
#            print('Quiral residues fixed {} ({})'.format(chiral_fix, n))
#            if check_clashes:
#                if not self.args['quiet']:
#                    print("Checking new steric clashes")
#
#            self.stm.modified = True

# =============================================================================
    def clashes(self, opts=None):
        """ run clashes command """
        self._run_method('clashes', opts)

    def _clashes_check(self):

        clash_list = {}
        for cls in ALL_CONTACT_TYPES:
            clash_list[cls] = {}

        #Recalc when models are separated molecules
        if not self.stm.has_superimp_models():
            rr_dist = mu.get_all_r2r_distances(
                self._get_structure(),
                'all',
                self.stm.data_library.get_distances('R_R_CUTOFF'),
                join_models=True #not self.stm.has_superimp_models()
            )
        else:
            rr_dist = self.stm.rr_dist

        for r_pair in rr_dist:
            [res1, res2] = r_pair[0:2]

            if mu.is_wat(res1) or mu.is_wat(res2):
                continue

            c_list = self.stm.check_rr_clashes(res1, res2, ALL_CONTACT_TYPES)

            rkey = mu.residue_id(res1) + '-' + mu.residue_id(res2)
            for cls in c_list:
                if c_list[cls]:
                    clash_list[cls][rkey] = c_list[cls]

        self.summary['clashes']['detected'] = self._clash_report(
            ALL_CONTACT_TYPES, clash_list
        )

        return False

    def _clashes_fix(self, res_to_fix, fix_data=None):
        pass
# =============================================================================
    def fixside(self, opts=None):
        """ run fixside command """
        self._run_method('fixside', opts)

    def _fixside_check(self):
        fix_data = {}

        miss_at_list = self.stm.get_missing_atoms('side')

        if miss_at_list:
            fix_data['miss_at_list'] = miss_at_list
            fix_data['fixside_rnums'] = []
            self.summary['fixside']['detected'] = {}
            print(
                '{} Residues with missing side chain atoms found'.format(
                    len(miss_at_list)
                )
            )
            for r_at in miss_at_list:
                [res, at_list] = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                fix_data['fixside_rnums'].append(mu.residue_num(res))
                self.summary['fixside']['detected'][mu.residue_id(res)] = at_list
        else:
            if not self.args['quiet']:
                print("No residues with missing side chain atoms found")
        return fix_data

    def _fixside_fix(self, fix_side, fix_data=None, check_clashes=True):
        input_line = ParamInput('fixside', fix_side, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', fix_data['fixside_rnums'], case='sensitive', multiple=True)
        [input_option, fix_side] = input_line.run()

        if input_option == 'error':
            print("Invalid option", fix_side)
            self.summary['fix_side']['error'] = "Unknown option"
            return

        self.summary['fixside']['selected'] = fix_side

        if input_option == 'none':
            if not self.args['quiet']:
                print('Nothing to do')
        else:
            if input_option == 'all':
                to_fix = fix_data['miss_at_list']
            else:
                to_fix = []
                for r_at in fix_data['miss_at_list']:
                    if mu.residue_num(r_at[0]) in fix_side.split(','):
                        to_fix.append(r_at)
            if not self.args['quiet']:
                print("Fixing side chains")
            fix_num = 0
            self.summary['fixside']['fixed'] = []
            fixed_res = []
            for r_at in to_fix:
                mu.remove_H_from_r(r_at[0], verbose=True)
                self.stm.fix_side_chain(r_at, self.residue_lib)
                fix_num += 1
                self.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])

            print('Fixed {} side chain(s)'.format(fix_num))
            # Checking new clashes
            if check_clashes:
                print("Checking possible new clashes")

                self.summary['fixside_clashes'] = self._clash_report(
                    ALL_CONTACT_TYPES,
                    self.stm.check_r_list_clashes(fixed_res, ALL_CONTACT_TYPES)
                )
# =============================================================================
    def addH(self, opts=None):
        """ run addH command """
        self._run_method('addH', opts)

    def _addH_check(self):
        fix_data = {}

        ion_res_list = self.stm.get_ion_res_list()

        if ion_res_list:
            fix_data['ion_res_list'] = ion_res_list
            fix_data['add_h_rnums'] = []
            self.summary['addH']['detected'] = {}
            print(
                '{} Residues requiring selection on adding H atoms'.format(
                    len(ion_res_list)
                )
            )
            for r_at in ion_res_list:
                [res, at_list] = r_at
                #print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                fix_data['add_h_rnums'].append(mu.residue_num(res))
                self.summary['addH']['detected'][mu.residue_id(res)] = at_list
        else:
            if not self.args['quiet']:
                print("No residues requiring selection on adding H atoms")
        return fix_data

    def _addH_fix(self, mode, fix_data=None):
        input_line = ParamInput('Mode', mode, self.args['non_interactive'])
        input_line.add_option_none()
        input_line.add_option('auto', ['auto'], case="lower")
        input_line.add_option('inter', ['interactive'], case="lower")
        input_line.add_option('inter_His', ['interactive_his'], case="lower")
        input_line.add_option('ph', ['pH'])
        #input_line.add_option('resnum', sorted(self.add_h_rnums), case='sensitive', multiple=True)
        [input_option, add_h_mode] = input_line.run()

        if input_option == 'error':
            print("Invalid option", add_h_mode)
            self.summary['addH']['error'] = "Unknown option"
            return
        if input_option == "none":
            if not self.args['quiet']:
                print('Nothing to do')
        else:
            to_fix = []
            std_ion = self.stm.data_library.get_ion_data()
            if input_option == 'auto':
                if not self.args['quiet']:
                    print('Selection: auto')
                for r_at in fix_data['ion_res_list']:
                    to_fix.append([r_at[0], std_ion[r_at[0].get_resname()]['std']])
            elif input_option == 'ph':
                ph_value = None
                input_line = ParamInput("pH Value", ph_value, self.args['non_interactive'])
                input_line.add_option("pH", [], opt_type="float", min_val=0., max_val=14.)
                [input_option, ph_value] = input_line.run()
                if not self.args['quiet']:
                    print('Selection: pH', ph_value)
                for r_at in fix_data['ion_res_list']:
                    res = r_at[0]
                    rcode = res.get_resname()
                    if ph_value <= std_ion[rcode]['pK']:
                        to_fix.append([res, std_ion[rcode]['lowpH']])
                    else:
                        to_fix.append([res, std_ion[rcode]['highpH']])
            elif input_option == 'interactive':
                if not self.args['quiet']:
                    print('Selection: interactive')
                #TODO
            elif input_option == 'interactive_His':
                if not self.args['quiet']:
                    print('Selection: interactive-his')
                #TODO
            else:
                print("Not Valid")
                return
            self.stm.add_hydrogens(to_fix, self.residue_lib)

# =============================================================================
    def mutateside(self, mut_list):
        """ run mutateside command """
        self._run_method('mutateside', mut_list)

    def _mutateside_check(self):
        return True

    def _mutateside_fix(self, mut_list, fix_data=None, check_clashes=True):
        input_line = ParamInput('Mutation list', mut_list, self.args['non_interactive'])
        mut_list = input_line.run()

        mutations = self.stm.prepare_mutations(mut_list)
#        mutations = MutationManager(mut_list)
#        mutations.prepare_mutations(self.stm.st)

        print('Mutations to perform')
        for mut in mutations.mutation_list:
            print(mut)

        mutated_res = self.stm.apply_mutations(mutations, self.residue_lib)

        if check_clashes:
            if not self.args['quiet']:
                print("Checking new clashes")

            self.summary['mutateside_clashes'] = self._clash_report(
                ALL_CONTACT_TYPES,
                self.stm.check_r_list_clashes(mutated_res, ALL_CONTACT_TYPES)
            )

        self.stm.modified = True
#===============================================================================
    def backbone(self, opts):
        """ run backbone command """
        self._run_method('backbone', opts)

    def _backbone_check(self):
        fix_data = {}
        #backbone_atoms = self.stm.data_library.get_all_atom_lists()['GLY']['backbone']
        # Residues with missing backbone
        miss_bck_at_list = self.stm.get_missing_atoms('backbone')
        if miss_bck_at_list:
            fix_data['miss_bck_at_list'] = miss_bck_at_list
            self.summary['backbone']['missing_atoms'] = {}
            fix_data['fixbck_rnums'] = []
            print(
                '{} Residues with missing backbone atoms found'.format(
                    len(miss_bck_at_list)
                )
            )
            for r_at in miss_bck_at_list:
                [res, at_list] = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                fix_data['fixbck_rnums'].append(mu.residue_num(res))
                self.summary['backbone']['missing_atoms'][mu.residue_id(res)] = at_list

        #Not bound consecutive residues
        self.stm.get_backbone_breaks()
        if self.stm.bck_breaks_list:
            print("{} Backbone breaks found".format(len(self.stm.bck_breaks_list)))
            self.summary['backbone']['breaks'] = []
            for brk in self.stm.bck_breaks_list:
                print(
                    " {:10} - {:10} ".format(
                        mu.residue_id(brk[0]),
                        mu.residue_id(brk[1])
                    )
                )
                self.summary['backbone']['breaks'].append([
                    mu.residue_id(brk[0]), mu.residue_id(brk[1])
                ])
            fix_data['bck_breaks_list'] = True
        if self.stm.wrong_link_list:
            print("Unexpected backbone links found")
            self.summary['backbone']['wrong_links'] = []
            for brk in self.stm.wrong_link_list:
                print(
                    " {:10} linked to {:10}, expected {:10} ".format(
                        mu.residue_id(brk[0]),
                        mu.residue_id(brk[1]),
                        mu.residue_id(brk[2])
                    )
                )
                self.summary['backbone']['wrong_links'].append([
                    mu.residue_id(brk[0]),
                    mu.residue_id(brk[1]),
                    mu.residue_id(brk[2])
                ])
            fix_data['wrong_link_list'] = True

        if self.stm.not_link_seq_list:
            print("Consecutive residues too far away to be covalently linked")
            self.summary['backbone']['long_links'] = []
            for brk in self.stm.not_link_seq_list:
                print(
                    " {:10} - {:10}, bond distance {:8.3f} ".format(
                        mu.residue_id(brk[0]),
                        mu.residue_id(brk[1]),
                        brk[2]
                    )
                )
                self.summary['backbone']['long_links'].append([
                    mu.residue_id(brk[0]),
                    mu.residue_id(brk[1]),
                    round(float(brk[2]), 5)
                ])
            fix_data['not_link_seq_list'] = True
        #TODO move this section to ligands
        if self.stm.modified_residue_list:
            print("Modified residues found")
            self.summary['backbone']['mod_residues'] = []
            for brk in self.stm.modified_residue_list:
                print(" {:10} ".format(mu.residue_id(brk)))
                self.summary['backbone']['mod_residues'].append(mu.residue_id(brk))
        #Provisional only missing atoms can be fixed
            fix_data['modified_residue_list'] = True
        return fix_data

    def _backbone_fix(self, fix_back, fix_data=None, check_clashes=True):
        if 'miss_bck_at_list' not in fix_data:
            return
        input_line = ParamInput('fixbck', fix_back, self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option('resnum', fix_data['fixbck_rnums'], case='sensitive', multiple=True)
        [input_option, fix_back] = input_line.run()

        if input_option == 'error':
            print("Invalid option", fix_back)
            self.summary['backbone']['missing_atoms']['error'] = "Unknown option"
            return

        self.summary['backbone']['missing_atoms']['selected'] = fix_back

        if input_option == 'none':
            if not self.args['quiet']:
                print('Nothing to do')
        else:
            if input_option == 'all':
                to_fix = fix_data['miss_bck_at_list']
            else:
                to_fix = []
                for r_at in fix_data['miss_bck_at_list']:
                    if mu.residue_num(r_at[0]) in fix_back.split(','):
                        to_fix.append(r_at)
            if not self.args['quiet']:
                print("Adding missing backbone atoms")
            fix_num = 0
            self.summary['backbone']['missing_atoms']['fixed'] = []
            fixed_res = []
            for r_at in to_fix:
                if self.stm.fix_backbone_atoms(r_at):
                    fix_num += 1
                self.summary['backbone']['missing_atoms']['fixed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])

            print('Fixed {} backbone atom(s)'.format(fix_num))
            # Checking new clashes
            if check_clashes:
                print("Checking possible new clashes")
                self.summary['backbone']['missing_atoms']['clashes'] = self._clash_report(
                    ALL_CONTACT_TYPES,
                    self.stm.check_r_list_clashes(fixed_res, ALL_CONTACT_TYPES)
                )

        #TODO Chain fix
#===============================================================================
    def cistransbck(self, opts):
        """ run cistransbck command """
        self._run_method('cistransbck', opts)

    def _cistransbck_check(self):
        (cis_backbone_list, lowtrans_backbone_list) = self.stm.check_cis_backbone()
        if cis_backbone_list:
            self.summary['cistransbck']['cis'] = []
            print('{} cis peptide bonds'.format(len(cis_backbone_list)))
            for lnk in cis_backbone_list:
                [res1, res2, dih] = lnk
                print(
                    '{:10} {:10} Dihedral: {:8.3f}'.format(
                        mu.residue_id(res1),
                        mu.residue_id(res2),
                        dih
                    )
                )
                self.summary['cistransbck']['cis'].append([
                    mu.residue_id(res1),
                    mu.residue_id(res2),
                    round(float(dih), 3)
                ])
        else:
            if not self.args['quiet']:
                print("No cis peptide bonds found")

        if lowtrans_backbone_list:
            self.summary['cistransbck']['unusual_trans'] = []
            print(
                '{} trans peptide bonds with unusual omega dihedrals'.format(
                    len(lowtrans_backbone_list)
                )
            )
            for lnk in lowtrans_backbone_list:
                [res1, res2, dih] = lnk
                print(
                    '{:10} {:10} Dihedral: {:8.3f}'.format(
                        mu.residue_id(res1), mu.residue_id(res2), dih
                    )
                )
                self.summary['cistransbck']['unusual_trans'].append([
                    mu.residue_id(res1), mu.residue_id(res2), round(float(dih), 3)
                ])
        else:
            if not self.args['quiet']:
                print("No trans peptide bonds with unusual omega dihedrals found")
        return {}

    def _cistransbck_fix(self, option):
        pass
#===============================================================================
    def _load_structure(self, verbose=True, print_stats=True):
        if not self.args['non_interactive'] and self.args['input_structure_path'] is None:
            self.args['input_structure_path'] =\
                input("Enter input structure path (PDB, mmcif | pdb:pdbid): ")

        stm = StructureManager(
            self.args['input_structure_path'],
            DataLibManager(self.sets.data_library_path),
            pdb_server=self.pdb_server,
            cache_dir=self.cache_dir,
            file_format=self.args['file_format']
        )

        if verbose:
            print('Structure {} loaded'.format(self.args['input_structure_path']))
            stm.print_headers()
            print()
            self.summary['headers'] = stm.meta

        if print_stats:
            stm.print_stats()
            print()

        self.summary['stats'] = stm.get_stats()

        return stm

    def _get_structure(self):
        return self.stm.get_structure()

    def _save_structure(self):
        if not self.args['non_interactive'] and self.args['output_structure_path'] is None:
            self.args['output_structure_path'] = input("Enter output structure path: ")
        self.stm.save_structure(self.args['output_structure_path'])

    def _clash_report(self, contact_types, clash_list):
        summary = {}
        for cls in contact_types:
            summary[cls] = []
            if clash_list[cls]:
                print('{} Steric {} clashes detected'.format(len(clash_list[cls]), cls))
                for rkey in sorted(
                        clash_list[cls],
                        key=lambda x: _key_sort_atom_pairs(clash_list[cls][x])):
                    print(
                        ' {:12} {:12} {:8.3f} A'.format(
                            mu.atom_id(clash_list[cls][rkey][0]),
                            mu.atom_id(clash_list[cls][rkey][1]),
                            np.sqrt(clash_list[cls][rkey][2])
                        )
                    )
                    summary[cls].append({
                        'at1':mu.atom_id(clash_list[cls][rkey][0]),
                        'at2':mu.atom_id(clash_list[cls][rkey][1]),
                        'dist': round(float(np.sqrt(clash_list[cls][rkey][2])), 4)
                    })
            else:
                if not self.args['quiet']:
                    print('No {} clashes detected'.format(cls))
        return summary

#===============================================================================
def _key_sort_atom_pairs(at_pair):
    return at_pair[0].serial_number * 10000 + at_pair[1].serial_number
