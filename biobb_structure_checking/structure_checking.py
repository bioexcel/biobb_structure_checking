"""
    Main Class for Structure Checking functionality

"""
__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys
import numpy as np

import biobb_structure_checking.constants as cts

from biobb_structure_checking.json_writer import JSONWriter
from biobb_structure_checking.param_input import ParamInput, NoDialogAvailableError

import biobb_structure_manager.structure_manager as stm
import biobb_structure_manager.model_utils as mu

# Main class
class StructureChecking():
    """Main class to control structure checking front end"""
    def __init__(self, sets, args):
        self.args = args
        self.sets = sets
        self.summary = {}

        self.pdb_server = self.args['pdb_server']
        self.cache_dir = self.args['cache_dir']

        try:
            self.strucm = self._load_structure(self.args['input_structure_path'])
        except IOError:
            print(
                'ERROR: fetching/parsing structure from {}'.format(
                    self.args['input_structure_path']
                ), file=sys.stderr
            )
            sys.exit(2)
#        except OSError:
#            print("ERROR: parsing PDB", file=sys.stderr)
#            sys.exit(2)
        except (stm.WrongServerError, stm.UnknownFileTypeError) as err:
            print(err.message, file=sys.stderr)
            sys.exit(1)

        if 'Notebook' not in self.args:
            self.args['Notebook'] = False

        if self.args['Notebook']:
            self.args['non_interactive'] = True
            self.args['check_only'] = False

    def launch(self):
        """main method to run checking"""
        if self.args['command'] == 'load':
            sys.exit()

        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        elif self.args['command'] == 'checkall':
            self.checkall(self.args['options'])
        else:
            self._run_method(self.args['command'], self.args['options'])

        if not self.args['check_only']:
            if self.strucm.modified or self.args['force_save']:
                if not self.strucm.modified:
                    print(cts.MSGS['FORCE_SAVE_STRUCTURE'])
                self.strucm.calc_stats()
                self.strucm.print_stats('Final')
                self.summary['final_stats'] = self.strucm.get_stats()
                try:
                    output_structure_path = self._save_structure(
                        self.args['output_structure_path']
                    )
                    print(cts.MSGS['STRUCTURE_SAVED'], output_structure_path)

                except OSError:
                    print(
                        'ERROR: unable to save PDB data on ',
                        output_structure_path,
                        file=sys.stderr
                    )
                except stm.OutputPathNotProvidedError as err:
                    print(err.message, file=sys.stderr)
            elif not self.strucm.modified:
                print(cts.MSGS['NON_MODIFIED_STRUCTURE'])

        if self.args['json_output_path'] is not None:
            json_writer = JSONWriter()
            for k in self.summary:
                json_writer.set(k, self.summary[k])
            try:
                json_writer.save(self.args['json_output_path'])
                print(cts.MSGS['JSON_SAVED'], self.args['json_output_path'])
            except IOError:
                print(cts.MSGS['JSON_NOT_SAVED'], self.args['json_output_path'])

    def command_list(self, opts):
        """ Manages command_list workflows"""
        try:
            opts = cts.DIALOGS.get_parameter('command_list', opts)
            op_list = opts[cts.DIALOGS.get_dialog('command_list')['dest']]
        except NoDialogAvailableError as err:
            print(err.message)

        op_list = ParamInput('Command List File', False).run(op_list)

        try:
            list_file_h = open(op_list, "r")
        except OSError:
            print(cts.MSGS['ERROR_OPEN_FILE'], op_list)
            sys.exit(1)
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

        print(cts.MSGS['COMMAND_LIST_COMPLETED'])

    def checkall(self, opts):
        """ Predefined workflow for complete checking"""
        self.args['check_only'] = True

        for meth in cts.AVAILABLE_METHODS:
            self._run_method(meth, opts)

    def _run_method(self, command, opts):
        """ Run check and fix methods for specific command"""
        try:
            f_check = getattr(self, '_' + command + '_check')
        except AttributeError:
            print(cts.MSGS['COMMAND_NOT_FOUND'].format(command))
            sys.exit(1)

        if command not in self.summary:
            self.summary[command] = {}

        msg = 'Running {}.'.format(command)

        if opts:
            self.summary[command]['opts'] = opts
            msg += ' Options: {} '.format(' '.join(opts))

        if not self.args['quiet']:
            print(msg)
    #Running checking method
        data_to_fix = f_check()
    #Running fix method if needed
        if self.args['check_only'] or opts in (None, ''):
            if not self.args['quiet']:
                print(cts.MSGS['CHECK_ONLY_DONE'])

        elif data_to_fix:
            try:
                f_fix = getattr(self, '_' + command + '_fix')
            except AttributeError:
                print(cts.MSGS['FIX_COMMAND_NOT_FOUND'].format(command))
                sys.exit(1)

            if not self.args['Notebook']:
                if cts.DIALOGS.exists(command):
                    opts = cts.DIALOGS.get_parameter(command, opts)
                    opts = opts[cts.DIALOGS.get_dialog(command)['dest']]
                else:
                    opts = ''
            error_status = f_fix(opts, data_to_fix)
            if error_status:
                print('ERROR', ' '.join(error_status), file=sys.stderr)
                self.summary[command]['error'] = ' '.join(error_status)


# ==============================================================================q
    def models(self, opts=None):
        """ direct entry to run models command """
        self._run_method('models', opts)

    def _models_check(self):
        print(cts.MSGS['MODELS_FOUND'].format(self.strucm.nmodels))
        self.summary['models'] = {'nmodels': self.strucm.nmodels}
        if self.strucm.nmodels == 1:
            if not self.args['quiet']:
                print(cts.MSGS['SINGLE_MODEL'])
            return {}

        self.summary['models']['type'] = self.strucm.models_type
        supimp = ''
        if self.strucm.models_type['type'] == mu.BUNIT:
            supimp = 'do not'
        print(
            cts.MSGS['MODELS_GUESS'].format(
                supimp,
                self.strucm.models_type['rmsd'],
                mu.MODEL_TYPE_LABELS[self.strucm.models_type['type']]
            )
        )

        return True

    def _models_fix(self, select_model, fix_data=None):
        input_line = ParamInput('Select Model Num', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_numeric(
            'modelno', [], opt_type='int', min_val=1, max_val=self.strucm.nmodels
        )
        [input_option, select_model] = input_line.run(select_model)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], select_model]

        print(cts.MSGS['SELECT_MODEL'], select_model)

        if input_option != 'all':
            self.strucm.select_model(select_model)

        self.summary['models']['selected'] = select_model

        return False

# =============================================================================
    def chains(self, opts=None):
        """ Run chains command """
        self._run_method('chains', opts)

    def _chains_check(self):
        print(cts.MSGS['CHAINS_DETECTED'].format(len(self.strucm.chain_ids)))
        for ch_id in sorted(self.strucm.chain_ids):
            if isinstance(self.strucm.chain_ids[ch_id], list):
                print(
                    cts.MSGS['UNKNOWN_CHAINS'].format(
                        ch_id, s=self.strucm.chain_ids[ch_id]
                    )
                )
            else:
                print(
                    ' {}: {}'.format(
                        ch_id, mu.CHAIN_TYPE_LABELS[self.strucm.chain_ids[ch_id]]
                    )
                )
        self.summary['chains'] = {'ids':self.strucm.chain_ids}
        return len(self.strucm.chain_ids) > 1

    def _chains_fix(self, select_chains, fix_data=None):
        self.summary['chains']['selected'] = {}
        input_line = ParamInput('Select chain', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_list(
            'chid', sorted(self.strucm.chain_ids), multiple=True, case="sensitive"
        )
        [input_option, select_chains] = input_line.run(select_chains)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], select_chains]

        if input_option == 'all':
            print(cts.MSGS['SELECT_ALL_CHAINS'])
            return False

        self.strucm.select_chains(select_chains)
        print(cts.MSGS['SELECT_CHAINS'], select_chains)
        self.summary['chains']['selected'] = self.strucm.chain_ids

        return False

# =============================================================================
    def inscodes(self, opts=None):
        """ Run inscodes command """
        self._run_method('inscodes', opts)

    def _inscodes_check(self):
        ins_codes_list = self.strucm.get_ins_codes()
        if not ins_codes_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_INSCODES_FOUND'])
            return {}

        print(cts.MSGS['INSCODES_FOUND'].format(len(ins_codes_list)))
        self.summary['inscodes'] = []
        for res in ins_codes_list:
            print(mu.residue_id(res))
            self.summary['inscodes'].append(mu.residue_id(res))
        return {}

#    def _inscodes_fix(self, select_codes, fix_data=None):
#        Renumber residues??
#        pass
# =============================================================================
    def altloc(self, opts=None):
        """ run altloc command """
        self._run_method('altloc', opts)

    def _altloc_check(self): #TODO improve output
        alt_loc_res = mu.get_altloc_residues(self._get_structure())
        if not alt_loc_res:
            if not self.args['quiet']:
                print(cts.MSGS['NO_ALTLOC_FOUND'])
            return {}

        print(cts.MSGS['ALTLOC_FOUND'].format(len(alt_loc_res)))

        self.summary['altloc'] = {}

        fix_data = {
            'alt_loc_res' : alt_loc_res,
            'alt_loc_rnums' : [],
            'altlocs' : {}
        }

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

        return fix_data

    def _altloc_fix(self, select_altloc, fix_data=None):
        #Prepare the longest possible list of alternatives
        altlocs = []
        max_al_len = 0
        for res in fix_data['altlocs']:
            if len(fix_data['altlocs'][res]) > max_al_len:
                altlocs = fix_data['altlocs'][res]
                max_al_len = len(fix_data['altlocs'][res])

        input_line = ParamInput('Select alternative', self.args['non_interactive'])
        input_line.add_option_list('occup', ['occupancy'])
        input_line.add_option_list('altids', altlocs, case='upper')
        input_line.add_option_list(
            'resnum',
            fix_data['alt_loc_rnums'],
            opt_type='pair_list',
            list2=altlocs,
            case='sensitive',
            multiple=True
        )
        [input_option, select_altloc] = input_line.run(select_altloc)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], select_altloc]

        print('Selecting location {}'.format(select_altloc))
        if input_option in ('occup', 'altids'):
            to_fix = {}
            select_altloc = select_altloc.upper()
            for res in fix_data['alt_loc_res']:
                to_fix[res] = {
                    'ats': fix_data['alt_loc_res'][res],
                    'select' : select_altloc
                }
        elif input_option == 'resnum':
            to_fix = {}
            selected_rnums = {}
            for rsel in select_altloc.split(','):
                [rnum, alt] = rsel.split(':')
                selected_rnums[rnum] = alt
            for res in fix_data['alt_loc_res']:
                rnum0 = mu.residue_num(res)
                if rnum0 in selected_rnums:
                    to_fix[res] = {
                        'ats' : fix_data['alt_loc_res'][res],
                        'select' : selected_rnums[rnum0]
                    }

        self.summary['altloc']['selected'] = select_altloc
        for res in to_fix:
                self.strucm.select_altloc_residue(res, to_fix[res])
        return False
# =============================================================================
    def metals(self, opts=None):
        """ Run metals command """
        self._run_method('metals', opts)

    def _metals_check(self):
        met_list = self.strucm.get_metal_atoms()

        if not met_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_METALS_FOUND'])
            return {}

        print(cts.MSGS['METALS_FOUND'].format(len(met_list)))
        self.summary['metals'] = {'detected':[]}
        fix_data = {
            'met_list': met_list,
            'met_rids': [],
            'at_groups': {}
        }
        for atm in sorted(met_list, key=lambda x: x.serial_number):
            print(" {:12}".format(mu.atom_id(atm)))
            res = atm.get_parent()
            fix_data['met_rids'].append(mu.residue_num(res))
            if not atm.id in fix_data['at_groups']:
                fix_data['at_groups'][atm.id] = []
            fix_data['at_groups'][atm.id].append(atm)
            self.summary['metals']['detected'].append(mu.residue_num(res))

        return fix_data

    def _metals_fix(self, remove_metals, fix_data=None):
        input_sess = ParamInput("Remove", self.args['non_interactive'])
        input_sess.add_option_all()
        input_sess.add_option_none()
        input_sess.add_option_list(
            'atids',
            sorted(fix_data['at_groups']),
            case='sensitive',
            multiple=True
        )
        input_sess.add_option_list('resids', fix_data['met_rids'], case='sensitive', multiple=True)
        [input_option, remove_metals] = input_sess.run(remove_metals)

        if input_option == "error":
            return [cts.MSGS['UNKNOWN_SELECTION'], remove_metals]

        if input_option == 'none':
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False

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
            self.summary['metals']['removed'].append(mu.residue_id(atm.get_parent()))
            self.strucm.remove_residue(atm.get_parent(), False)
            rmm_num += 1
        self.strucm.update_internals()
        print(cts.MSGS['METALS_REMOVED'].format(remove_metals, rmm_num))
        self.summary['metals']['n_removed'] = rmm_num
        return False
# =============================================================================
    def water(self, opts=None):
        """ Run water command """
        self._run_method('water', opts)

    def _water_check(self):
        lig_list = mu.get_ligands(self._get_structure(), incl_water=True)
        wat_list = []
        for res in lig_list:
            if mu.is_wat(res):
                wat_list.append(res)

        if not wat_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_WATERS'])
            return {}

        print(cts.MSGS['WATERS_FOUND'].format(len(wat_list)))
        self.summary['water']['n_detected'] = len(wat_list)

        return {'wat_list': wat_list}

    def _water_fix(self, remove_wat, fix_data=None):
        input_line = ParamInput('Remove', self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_wat] = input_line.run(remove_wat)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], remove_wat]

        if input_option == 'yes':
            rmw_num = 0
            for res in fix_data['wat_list']:
                self.strucm.remove_residue(res, False)
                rmw_num += 1
            self.strucm.update_internals()
            print(cts.MSGS['WATER_REMOVED'].format(rmw_num))
            self.summary['water']['n_removed'] = rmw_num
        return False
# =============================================================================
    def ligands(self, opts=None):
        """ Run ligands command """
        self._run_method('ligands', opts)

    def _ligands_check(self):
        lig_list = mu.get_ligands(self._get_structure(), incl_water=False)

        if not lig_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_LIGANDS_FOUND'])
            return {}

        print(cts.MSGS['LIGANDS_DETECTED'].format(len(lig_list)))
        fix_data = {
            'lig_list': lig_list,
            'ligand_rids': set(),
            'ligand_rnums': []
        }
        self.summary['ligands'] = {'detected': []}

        for res in sorted(lig_list, key=lambda x: x.index):
            if self.strucm.has_models():
                if res.get_parent().get_parent().id > 0:
                    continue
                print(' {}/*'.format(mu.residue_id(res, False)))
            else:
                print(' {}'.format(mu.residue_id(res, False)))

            self.summary['ligands']['detected'].append(mu.residue_id(res))
            fix_data['ligand_rids'].add(res.get_resname())
            fix_data['ligand_rnums'].append(mu.residue_num(res))

        return fix_data

# =============================================================================
    def _ligands_fix(self, remove_ligands, fix_data=None):
        input_line = ParamInput('Remove', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list('byrids', sorted(fix_data['ligand_rids']), multiple=True)
        input_line.add_option_list(
            'byresnum', fix_data['ligand_rnums'], case='sensitive', multiple=True
        )
        [input_option, remove_ligands] = input_line.run(remove_ligands)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], remove_ligands]

        self.summary['ligands']['removed'] = {'opt':remove_ligands, 'lst':[]}

        to_remove = []

        if input_option == 'none':
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False

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
                mu.residue_id(res, self.strucm.has_models())
            )
            self.strucm.remove_residue(res, False)
            rl_num += 1
        self.strucm.update_internals()
        print(cts.MSGS['LIGANDS_REMOVED'].format(remove_ligands, rl_num))
        self.summary['ligands']['n_removed'] = rl_num
        return False

# =============================================================================
    def rem_hydrogen(self, opts=None):
        """ Run rem_hydrogen command """
        self._run_method('rem_hydrogen', opts)

    def _rem_hydrogen_check(self):

        remh_list = mu.get_residues_with_H(self._get_structure())

        if not remh_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_RESIDUES_H_FOUND'])
            return {}

        print(cts.MSGS['RESIDUES_H_FOUND'].format(len(remh_list)))
        self.summary['rem_hydrogen']['n_detected'] = len(remh_list)
        return {'remh_list': remh_list}

    def _rem_hydrogen_fix(self, remove_h, fix_data=None):
        input_line = ParamInput('Remove hydrogen atoms', self.args['non_interactive'])
        input_line.add_option_yes_no()
        [input_option, remove_h] = input_line.run(remove_h)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], remove_h]

        if input_option == 'yes':
            rmh_num = 0
            for resh in fix_data['remh_list']:
                mu.remove_H_from_r(resh['r'])
                rmh_num += 1
            print(cts.MSGS['REMOVED_H'].format(rmh_num))
            self.strucm.modified = True
            self.summary['rem_hydrogen']['n_removed'] = rmh_num
        return False

# =============================================================================
    def getss(self, opts=None):
        """ run getss command """
        self._run_method('getss', opts)

    def _getss_check(self):

        SS_bonds = self.strucm.get_SS_bonds()

        if not SS_bonds:
            if not self.args['quiet']:
                print(cts.MSGS['NO_SS'])
            return {}
        print(cts.MSGS['POSSIBLE_SS'].format(len(SS_bonds)))
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
        return {}

#    def _getss_fix(self):
#        pass

# =============================================================================
    def amide(self, opts=None):
        """ run amide command """
        self._run_method('amide', opts)

    def _amide_check(self):

        amide_atoms = self.strucm.data_library.get_amide_data()[1]

        amide_list = self.strucm.get_amide_list()

        self.summary['amide']['n_amides'] = len(amide_list)

        if not amide_list:
            if not self.args['quiet']:
                print('NO_AMIDES')
            return {}

        c_list = self.strucm.check_r_list_clashes(
            amide_list,
            stm.AMIDE_CONTACT_TYPES
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

        if not amide_cont_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_UNUSUAL_AMIDES'])
            return {}

        print(cts.MSGS['UNUSUAL_AMIDES'].format(len(amide_cont_list)))
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
        return {
            'amide_list':  amide_list,
            'amide_res_to_fix': amide_res_to_fix,
            'amide_cont_list': amide_cont_list,
            'amide_rnums': amide_rnums
        }

# =============================================================================
    def _amide_fix(self, amide_fix, fix_data=None, recheck=True):
        no_int_recheck = amide_fix is not None or self.args['non_interactive']

        while fix_data:
            input_line = ParamInput(
                'Fix amide atoms', self.args['non_interactive']
            )
            input_line.add_option_all()
            input_line.add_option_none()
            input_line.add_option_list(
                'resnum', sorted(fix_data['amide_rnums']), case='sensitive', multiple=True
            )
            [input_option, amide_fix] = input_line.run(amide_fix)

            if input_option == 'error':
                return [cts.MSGS['UNKNOWN_SELECTION'], amide_fix]

            if input_option == 'none':
                if not self.args['quiet']:
                    print(cts.MSGS['DO_NOTHING'])
                return False

            to_fix = []
            for res in fix_data['amide_res_to_fix']:
                if mu.residue_num(res) in amide_fix.split(',')\
                        or input_option == 'all':
                    to_fix.append(res)
            fix_num = 0
            done = set()
            for res in to_fix:
                if res not in done:
                    try:
                        self.strucm.invert_amide_atoms(res)
                    except stm.NotAValidResidueError as err:
                        return [err.message]
                    done.add(res) # To avoid double change
                    fix_num += 1
                    
            print(cts.MSGS['AMIDES_FIXED'].format(amide_fix, fix_num))
            self.strucm.modified = True
            fix_data = {}
            if recheck:
                if not self.args['quiet']:
                    print(cts.MSGS['AMIDES_RECHECK'])
                fix_data = self._amide_check()
                amide_fix = ''
                if no_int_recheck:
                    fix_data = {}
        return False

# =============================================================================
    def chiral(self, opts=None):
        """ run chiral command """
        self._run_method('chiral', opts)

    def _chiral_check(self):

        chiral_list = self.strucm.get_chiral_list()

        self.summary['chiral']['n_chirals'] = len(chiral_list)

        if not chiral_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_CHIRALS'])
            return {}

        chiral_res = self.strucm.data_library.get_chiral_data()

        chiral_res_to_fix = []
        chiral_rnums = []

        for res in chiral_list:
            if not mu.check_chiral_residue(res, chiral_res):
                chiral_res_to_fix.append(res)
                chiral_rnums.append(mu.residue_num(res))

        if not chiral_res_to_fix:
            if not self.args['quiet']:
                print(cts.MSGS['NO_WRONG_CHIRAL_SIDE'])
            return {}

        print(cts.MSGS['WRONG_CHIRAL_SIDE'].format(len(chiral_res_to_fix)))
        self.summary['chiral']['detected'] = []
        for res in chiral_res_to_fix:
            print(' {:10}'.format(mu.residue_id(res)))
            self.summary['chiral']['detected'].append(mu.residue_id(res))

        return {'chiral_res_to_fix': chiral_res_to_fix, 'chiral_rnums' : chiral_rnums}

    def _chiral_fix(self, chiral_fix, fix_data=None, check_clashes=True):

        input_line = ParamInput('Fix chiralities', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum', fix_data['chiral_rnums'], case='sensitive', multiple=True
        )
        [input_option, chiral_fix] = input_line.run(chiral_fix)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], chiral_fix]

        if input_option == 'none':
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'all':
            to_fix = fix_data['chiral_res_to_fix']
        else:
            to_fix = []
            for res in fix_data['chiral_res_to_fix']:
                if mu.residue_num(res) in chiral_fix.split(','):
                    to_fix.append(res)
        fix_num = 0
        for res in to_fix:
            self.strucm.fix_chiral_chains(res)
            fix_num += 1
        print(cts.MSGS['CHIRAL_SIDE_FIXED'].format(chiral_fix, fix_num))

        if check_clashes:
            if not self.args['quiet']:
                print(cts.MSGS['CHECKING_CLASHES'])

            self.summary['chiral_clashes'] = self._check_report_clashes(to_fix)

        self.strucm.modified = True
        return False

# =============================================================================
    def chiral_bck(self, opts=None):
        """ run chiral_bck command """
        self._run_method('chiral_bck', opts)

    def _chiral_bck_check(self):

        chiral_bck_list = self.strucm.get_chiral_bck_list()

        self.summary['chiral_bck']['n_chirals'] = len(chiral_bck_list)

        if chiral_bck_list:
            chiral_bck_res_to_fix = []
            chiral_bck_rnums = []
            for res in chiral_bck_list:
                if not mu.check_chiral_ca(res):
                    chiral_bck_res_to_fix.append(res)
                    chiral_bck_rnums.append(mu.residue_num(res))

            if not chiral_bck_res_to_fix:
                if not self.args['quiet']:
                    print(cts.MSGS['NO_CHIRAL_BCK_RESIDUES'])
                return {}

            print(cts.MSGS['CHIRAL_BCK_RESIDUES'].format(len(chiral_bck_res_to_fix)))
            self.summary['chiral_bck']['detected'] = []
            for res in chiral_bck_res_to_fix:
                print(' {:10}'.format(mu.residue_id(res)))
                self.summary['chiral_bck']['detected'].append(mu.residue_id(res))

#            return {
#                'chiral_bck_res_to_fix': chiral_bck_res_to_fix,
#                'chiral_bck_rnums': chiral_bck_rnums
#            }

        return {}

#    def _chiral_bck_fix(self, chiral_fix, fix_data=None, check_clashes=True):
#        return False
#TODO chiral_bck_fix
#        input_line = ParamInput('Fix CA chiralities', self.args['non_interactive'])
#        input_line.add_option_all()
#        input_line.add_option_none()
#        input_line.add_option_list('resnum', self.chiral_bck_rnums,
#            case='sensitive', multiple=True)
#        [input_option, chiral_fix] = input_line.run(chiral_fix)
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
#                mu.strucm.invert_chiral_CA(res)
#                n += 1
#            print('Quiral residues fixed {} ({})'.format(chiral_fix, n))
#            if check_clashes:
#                if not self.args['quiet']:
#                    print("Checking new steric clashes")
#
#            self.strucm.modified = True

# =============================================================================
    def clashes(self, opts=None):
        """ run clashes command """
        self._run_method('clashes', opts)

    def _clashes_check(self):

        clash_list = {}
        for cls in stm.ALL_CONTACT_TYPES:
            clash_list[cls] = {}

        #Recalc when models are separated molecules
        if not self.strucm.has_superimp_models():
            rr_dist = mu.get_all_r2r_distances(
                self._get_structure(),
                'all',
                self.strucm.data_library.get_distances('R_R_CUTOFF'),
                join_models=True #not self.strucm.has_superimp_models()
            )
        else:
            rr_dist = self.strucm.rr_dist

        for r_pair in rr_dist:
            [res1, res2] = r_pair[0:2]

            if mu.is_wat(res1) or mu.is_wat(res2):
                continue

            c_list = self.strucm.check_rr_clashes(res1, res2, stm.ALL_CONTACT_TYPES)

            rkey = mu.residue_id(res1) + '-' + mu.residue_id(res2)
            for cls in c_list:
                if c_list[cls]:
                    clash_list[cls][rkey] = c_list[cls]

        self.summary['clashes']['detected'] = self._clash_report(
            stm.ALL_CONTACT_TYPES, clash_list
        )

        return False

#    def _clashes_fix(self, res_to_fix, fix_data=None):
#        pass
# =============================================================================
    def fixside(self, opts=None):
        """ run fixside command """
        self._run_method('fixside', opts)

    def _fixside_check(self):
        miss_at_list = self.strucm.get_missing_atoms('side')

        if not miss_at_list:
            if not self.args['quiet']:
                print("No residues with missing side chain atoms found")
            return {}

        fix_data = {
            'miss_at_list': miss_at_list,
            'fixside_rnums': []
        }

        self.summary['fixside']['detected'] = {}
        print(cts.MSGS['MISSING_SIDE_ATOMS'].format(len(miss_at_list)))
        for r_at in miss_at_list:
            [res, at_list] = r_at
            print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
            fix_data['fixside_rnums'].append(mu.residue_num(res))
            self.summary['fixside']['detected'][mu.residue_id(res)] = at_list

        return fix_data

    def _fixside_fix(self, fix_side, fix_data=None, check_clashes=True):
        input_line = ParamInput('fixside', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum', fix_data['fixside_rnums'], case='sensitive', multiple=True
        )
        [input_option, fix_side] = input_line.run(fix_side)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], fix_side]

        self.summary['fixside']['selected'] = fix_side

        if input_option == 'none':
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'all':
            to_fix = fix_data['miss_at_list']
        else:
            to_fix = []
            for r_at in fix_data['miss_at_list']:
                if mu.residue_num(r_at[0]) in fix_side.split(','):
                    to_fix.append(r_at)

        if not self.args['quiet']:
            print(cts.MSGS['FIXING_SIDE_CHAINS'])

        fix_num = 0
        self.summary['fixside']['fixed'] = []
        fixed_res = []
        for r_at in to_fix:
            mu.remove_H_from_r(r_at[0], verbose=True)
            self.strucm.fix_side_chain(r_at)
            fix_num += 1
            self.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))
            fixed_res.append(r_at[0])

        print(cts.MSGS['SIDE_CHAIN_FIXED'].format(fix_num))
        # Checking new clashes
        if check_clashes:
            print(cts.MSGS['CHECKING_CLASHES'])

            self.summary['fixside_clashes'] = self._check_report_clashes(fixed_res)
        return False
# =============================================================================
    def add_hydrogen(self, opts=None):
        """ Run add_hydrogen command """
        self._run_method('add_hydrogen', opts)

    def _add_hydrogen_check(self):

        ion_res_list = self.strucm.get_ion_res_list()

        if not ion_res_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_SELECT_ADDH'])
            return {'ion_res_list':[]}

        fix_data = {
            'ion_res_list': ion_res_list,
            'add_h_rnums' : []
        }
        self.summary['add_hydrogen']['detected'] = {}
        print(cts.MSGS['SELECT_ADDH_RESIDUES'].format(len(ion_res_list)))
        for r_at in ion_res_list:
            [res, at_list] = r_at
            #print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
            fix_data['add_h_rnums'].append(mu.residue_num(res))
            self.summary['add_hydrogen']['detected'][mu.residue_id(res)] = at_list
        return fix_data

    def _add_hydrogen_fix(self, mode, fix_data=None):
        input_line = ParamInput('Mode', self.args['non_interactive'])
        input_line.add_option_none()
        input_line.add_option_list('auto', ['auto'], case="lower")
        input_line.add_option_list('ph', ['pH'])
        if fix_data['ion_res_list']:
            input_line.add_option_list('inter', ['interactive'], case="lower")
            input_line.add_option_list('inter_His', ['interactive_his'], case="lower")
        #input_line.add_option('resnum', sorted(self.add_h_rnums), case='sensitive', multiple=True)
        [input_option, add_h_mode] = input_line.run(mode)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], add_h_mode]

        if input_option == "none":
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        to_fix = []
        std_ion = self.strucm.data_library.get_ion_data()
        if input_option == 'auto':
            if not self.args['quiet']:
                print('Selection: auto')
            for r_at in fix_data['ion_res_list']:
                to_fix.append([r_at[0], std_ion[r_at[0].get_resname()]['std']])
        elif input_option == 'ph':
            ph_value = None
            input_line = ParamInput("pH Value", self.args['non_interactive'])
            input_line.add_option_numeric("pH", [], opt_type="float", min_val=0., max_val=14.)
            [input_option, ph_value] = input_line.run(ph_value)
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
        self.strucm.add_hydrogens(to_fix)
        return False

# =============================================================================
    def mutateside(self, mut_list):
        """ Run mutateside command """
        self._run_method('mutateside', mut_list)

    def _mutateside_check(self):
        return True

    def _mutateside_fix(self, mut_list, fix_data=None, check_clashes=True):
        input_line = ParamInput('Mutation list', self.args['non_interactive'])
        mut_list = input_line.run(mut_list)

        mutations = self.strucm.prepare_mutations(mut_list)

        print(cts.MSGS['MUTATIONS_TO_DO'])
        for mut in mutations.mutation_list:
            print(mut)

        mutated_res = self.strucm.apply_mutations(mutations)

        if check_clashes:
            if not self.args['quiet']:
                print(cts.MSGS['CHECKING_CLASHES'])

            self.summary['mutateside_clashes'] = self._check_report_clashes(mutated_res)

        self.strucm.modified = True
        return False
#===============================================================================
    def backbone(self, opts):
        """ Run backbone command """
        self._run_method('backbone', opts)

    def _backbone_check(self):
        fix_data = {}
        # Residues with missing backbone
        miss_bck_at_list = self.strucm.get_missing_atoms('backbone')
        if miss_bck_at_list:
            fix_data['miss_bck_at_list'] = miss_bck_at_list
            self.summary['backbone']['missing_atoms'] = {}
            fix_data['fixbck_rnums'] = []
            print(cts.MSGS['BCK_MISSING_RESIDUES'].format(len(miss_bck_at_list)))
            for r_at in miss_bck_at_list:
                [res, at_list] = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                fix_data['fixbck_rnums'].append(mu.residue_num(res))
                self.summary['backbone']['missing_atoms'][mu.residue_id(res)] = at_list

        #Not bound consecutive residues
        bck_check  = self.strucm.get_backbone_breaks()
        if bck_check['bck_breaks_list']:
            print(cts.MSGS['BACKBONE_BREAKS'].format(len(bck_check['bck_breaks_list'])))
            self.summary['backbone']['breaks'] = []
            for brk in bck_check['bck_breaks_list']:
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
        if bck_check['wrong_link_list']:
            print(cts.MSGS['UNEXPECTED_BCK_LINKS'])
            self.summary['backbone']['wrong_links'] = []
            for brk in bck_check['wrong_link_list']:
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

        if bck_check['not_link_seq_list']:
            print(cts.MSGS['CONSEC_RES_FAR'])
            self.summary['backbone']['long_links'] = []
            for brk in bck_check['not_link_seq_list']:
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
        if self.strucm.modified_residue_list:
            print(cts.MSGS['MODIF_RESIDUES'])
            self.summary['backbone']['mod_residues'] = []
            for brk in self.strucm.modified_residue_list:
                print(" {:10} ".format(mu.residue_id(brk)))
                self.summary['backbone']['mod_residues'].append(mu.residue_id(brk))
        #Provisional only missing atoms can be fixed
            fix_data['modified_residue_list'] = True
        return fix_data

    def _backbone_fix(self, fix_back, fix_data=None, check_clashes=True):
        if 'miss_bck_at_list' not in fix_data:
            return False
        input_line = ParamInput('fixbck', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum', fix_data['fixbck_rnums'], case='sensitive', multiple=True
        )
        [input_option, fix_back] = input_line.run(fix_back)

        if input_option == 'error':
            return [cts.MSGS['UNKNOWN_SELECTION'], fix_back]

        self.summary['backbone']['missing_atoms']['selected'] = fix_back

        if input_option == 'none':
            if not self.args['quiet']:
                print(cts.MSGS['DO_NOTHING'])
            return False
        
        to_fix = []
        for r_at in fix_data['miss_bck_at_list']:
            if mu.residue_num(r_at[0]) in fix_back.split(',')\
                    or input_option == 'all':
                to_fix.append(r_at)
        
        if not self.args['quiet']:
            print(cts.MSGS['ADDING_BCK_ATOMS'])
        
        fix_num = 0
        self.summary['backbone']['missing_atoms']['fixed'] = []
        fixed_res = []
        for r_at in to_fix:
            try:
                if self.strucm.fix_backbone_O_atoms(r_at):
                    fix_num += 1
            except stm.NotEnoughAtomsError as err:
                print(err.message)

            self.summary['backbone']['missing_atoms']['fixed'].append(
                mu.residue_id(r_at[0])
            )
            fixed_res.append(r_at[0])

        print(cts.MSGS['BCK_ATOMS_FIXED'].format(fix_num))
        # Checking new clashes
        if fix_num and check_clashes:
            print(cts.MSGS['CHECKING_CLASHES'])
            self.summary['backbone']['missing_atoms']['clashes'] =\
                self._check_report_clashes(fixed_res)

        #TODO Chain fix
        return False
#===============================================================================
    def cistransbck(self, opts):
        """ Run cistransbck command """
        self._run_method('cistransbck', opts)

    def _cistransbck_check(self):
        (cis_backbone_list, lowtrans_backbone_list) = self.strucm.check_cis_backbone()
        if cis_backbone_list:
            self.summary['cistransbck']['cis'] = []
            print(cts.MSGS['CIS_BONDS'].format(len(cis_backbone_list)))
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
                print(cts.MSGS['NO_CIS_BONDS'])

        if lowtrans_backbone_list:
            self.summary['cistransbck']['unusual_trans'] = []
            print(cts.MSGS['LOWTRANS_BONDS'].format(len(lowtrans_backbone_list)))
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
                print(cts.MSGS['NO_LOWTRANS_BONDS'])

        return {}

#    def _cistransbck_fix(self, option):
#        pass
#===============================================================================
    def _load_structure(self, input_structure_path, verbose=True, print_stats=True):

        input_line = ParamInput(
            "Enter input structure path (PDB, mmcif | pdb:pdbid)",
            self.args['non_interactive']
        )
        input_structure_path = input_line.run(input_structure_path)

        strucm = stm.StructureManager(
            input_structure_path,
            self.sets.data_library_path,
            self.sets.res_library_path,
            pdb_server=self.pdb_server,
            cache_dir=self.cache_dir,
            file_format=self.args['file_format']
        )

        if verbose:
            print(cts.MSGS['STRUCTURE_LOADED'].format(input_structure_path))
            strucm.print_headers()
            print()
            self.summary['headers'] = strucm.meta

        if print_stats:
            strucm.print_stats()
            print()

        self.summary['stats'] = strucm.get_stats()

        return strucm

    def _get_structure(self):
        return self.strucm.get_structure()

    def _save_structure(self, output_structure_path):
        input_line = ParamInput(
            "Enter output structure path",
            self.args['non_interactive']
        )
        output_structure_path = input_line.run(output_structure_path)
        self.strucm.save_structure(output_structure_path)
        return output_structure_path

    def _check_report_clashes(self, residue_list, contact_types=None):
        if contact_types is None:
            contact_types = stm.ALL_CONTACT_TYPES
        return self._clash_report(
            contact_types,
            self.strucm.check_r_list_clashes(residue_list, contact_types)
        )

    def _clash_report(self, contact_types, clash_list):
        summary = {}
        for cls in contact_types:
            summary[cls] = []
            if clash_list[cls]:
                print(cts.MSGS['CLASHES_DETECTED'].format(len(clash_list[cls]), cls))
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
                    print(cts.MSGS['NO_CLASHES_DETECTED'].format(cls))
        return summary

#===============================================================================
def _key_sort_atom_pairs(at_pair):
    return at_pair[0].serial_number * 10000 + at_pair[1].serial_number
