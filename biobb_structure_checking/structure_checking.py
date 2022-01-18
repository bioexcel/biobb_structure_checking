"""
    Class for Structure Checking functionality
"""
__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys
import os
import time
import numpy as np

import biobb_structure_checking.constants as cts

from biobb_structure_checking.json_writer import JSONWriter
from biobb_structure_checking.param_input import ParamInput, NoDialogAvailableError

import biobb_structure_checking.structure_manager as stm
import biobb_structure_checking.model_utils as mu


# Main class
class StructureChecking():
    """
    | biobb_structure_checking.StructureChecking
    | Main class to control structure checking functionality
    | Provides support for to check_structure command line 
    | Load directly for Jupyter Notebook or python scripts. 
    
    Args:
        base_dir_path (str): Base directory path where application resides. 
        args (dict): Arguments dictionary see https://biobb-structure-checking.readthedocs.io/en/latest/command_line_usage.html. Recommended 'non-interactive':True for Notebook use.
    """
    def __init__(self, base_dir_path, args):

        if args is None:
            args = {}

        self.args = cts.set_defaults(base_dir_path, args)

        self.summary = {}

        if self.args['debug']:

            import psutil

            self.start_time = time.time()
            self.timings = []
            self.summary['elapsed_times'] = {}
            self.summary['memsize'] = []

        try:
            self.strucm = self._load_structure(self.args['input_structure_path'], self.args['fasta_seq_path'])
        except IOError:
            sys.exit(
                'ERROR: fetching/parsing structure from {}'.format(
                    self.args['input_structure_path']
                )
            )
        except (stm.WrongServerError, stm.UnknownFileTypeError, stm.ParseError) as err:
            sys.exit(err.message)

        if self.args['debug']:
            self.timings.append(['load', time.time() - self.start_time])
            process = psutil.Process(os.getpid())
            memsize = process.memory_info().rss/1024/1024
            self.summary['memsize'].append(['load', memsize])
            print(
                "#DEBUG Memory used after structure load: {:f} MB ".format(memsize)
            )

        if self.args['atom_limit'] and self.strucm.num_ats > self.args['atom_limit']:
            sys.exit(cts.MSGS['ATOM_LIMIT'].format(self.strucm.num_ats, self.args['atom_limit']))

    def help(self, command=None):
        """ StructureChecking.help
        Provides help on StructureChecking commands
        
        Args:
            command (str) : (None) Requested command. If empty returns all commands help
        """
        return cts.help(command)

    def launch(self):
        """ StructureChecking.launch
        Method run from the command line invocation
        """
        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        elif self.args['command'] == 'checkall':
            self.checkall()
        elif self.args['command'] == 'fixall':
            self.fixall(self.args['options'])
        elif self.args['command'] != 'load':
            self._run_method(self.args['command'], self.args['options'])
        if not self.args['check_only'] or self.args['force_save']:
            if self.strucm.modified or self.args['force_save']:
                if not self.strucm.modified:
                    print(cts.MSGS['FORCE_SAVE_STRUCTURE'])
                if not self.args['check_only']:
                    self.print_stats('Final')
                    self.summary['final_stats'] = self.strucm.get_stats()
                try:
                    output_structure_path = self._save_structure(
                        self.args['output_structure_path'],
                        self.args['rename_terms'],
                        split_models = '--save_split' in self.args['options']
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

        if self.args['debug']:
            total = time.time() - self.start_time
            self.summary['elapsed_times']['total'] = total
            print("#DEBUG TIMINGS")
            print("#DEBUG =======")
            ant = 0.
            for op, timing in self.timings:
                elapsed = timing - ant
                self.summary['elapsed_times'][op] = elapsed
                print(
                    '#DEBUG {:15s}: {:10.4f} s ({:6.2f}%)'.format(
                        op,
                        elapsed,
                        elapsed / total * 100.
                    )
                )
                ant = timing
            print('#DEBUG {:15s}: {:10.4F} s'.format('TOTAL', total))
            print()
            print("#DEBUG MEMORY USAGE EVOLUTION")
            print("#DEBUG ======================")
            for op, memsize in self.summary['memsize']:
                print('#DEBUG {:15s}: {:.2f} MB'.format(op, memsize))

        if self.args['json_output_path'] is not None:
            json_writer = JSONWriter()
            json_writer.data = self.summary
            try:
                json_writer.save(self.args['json_output_path'])
                print(cts.MSGS['JSON_SAVED'], self.args['json_output_path'])
            except IOError:
                print(cts.MSGS['JSON_NOT_SAVED'], self.args['json_output_path'])
    
    def print_stats(self, prefix=None):
        """ StructureChecking.print_stats
        Print statistics on the loaded structure
        
        Args:
            prefix (str): (None) Prefix to add to the output lines for identification.
        """
        self.strucm.calc_stats()
        if prefix is None: 
            prefix = ''
        self.strucm.print_stats(prefix)

    def command_list(self, opts):
        """ StructureChecking.command_list
        Manages command_list workflows
        
        Args:
            opts (str | list(str)): Command options as str or str list. 
        """
        try:
            opts = cts.DIALOGS.get_parameter('command_list', opts)
            op_list = opts['op_list']
        except NoDialogAvailableError as err:
            print(err.message)
        
        if not op_list:
            if not self.args['non_interactive']:
                op_list = ParamInput('Command List File', False).run(op_list)
            else:
                sys.exit(f'ERROR: command list not provided and non_interactive')
            

        if os.path.isfile(op_list):
            command_list = []
            try:
                with open(op_list, "r") as list_file_h:
                    for line in list_file_h:
                        if line == "\n" or line[0:1] == '#':
                            continue
                        command_list.append(line)
            except OSError:
                sys.exit('{} {}'.format(cts.MSGS['ERROR_OPEN_FILE'], op_list))
        else:
            command_list = op_list.split(';')
            
        i = 1
        for line in command_list:
            if not self.args['quiet']:
                print("\nStep {}: {}".format(i, line))
            data = line.split()
            command = data[0]
            opts = data[1:]
            self._run_method(command, opts)
            i += 1

        print(cts.MSGS['COMMAND_LIST_COMPLETED'])

    def checkall(self, opts=None):
        """ StructureChecking.checkall
        Predefined workflow for complete checking
        """
        #Required for interactive run in Notebooks
        old_check_only = self.args['check_only']
        self.args['check_only'] = True

        for meth in cts.AVAILABLE_METHODS:
            self._run_method(meth, None)

        self.args['check_only'] = old_check_only

    # def fixall(self):
    #     """ StructureChecking.fixall
    #     Fix all using defaults. Not implemented (yet)
    #     """
    #     # TODO Implement method fixall
    #     print("Fixall not implemented (yet)")

    def revert_changes(self):
        """ StructureChecking.revert_changes
        Reload original structure. Used in Pipelines or Notebooks to revert changes.
        """
        self.strucm = self._load_structure(self.args['input_structure_path'], self.args['fasta_seq_path'])
        self.summary = {}
        print(cts.MSGS['ALL_UNDO'])

    def _run_method(self, command, opts):
        """ Private. StructureChecking._run_method
        Run check and fix methods for specific command
        
        Args:
            command (str): Command to run
            opts (str | list(str) | dict): Command options, passed from callers
        """
        try:
            f_check = getattr(self, '_' + command + '_check')
        except AttributeError:
            sys.exit(cts.MSGS['COMMAND_NOT_FOUND'].format(command))

        if command not in self.summary:
            self.summary[command] = {}

        msg = 'Running {}.'.format(command)
        if opts:
            if isinstance(opts, list):
                opts_str = ' '.join(opts)
            elif isinstance(opts, dict):
                opts_str = str(opts)
            else:
                opts_str = opts
            msg += ' Options: ' + opts_str
            self.summary[command]['opts'] = opts_str

        if not self.args['quiet'] or self.args['verbose']:
            print(msg.strip())

    # Running checking method
        data_to_fix = f_check()
    # Running fix method if needed
        if self.args['check_only'] or opts in (None, ''):
            if self.args['verbose']:
                print(cts.MSGS['CHECK_ONLY_DONE'])
        elif data_to_fix:
            try:
                f_fix = getattr(self, '_' + command + '_fix')
            except AttributeError:
                sys.exit(cts.MSGS['FIX_COMMAND_NOT_FOUND'].format(command))
            if isinstance(opts, str) or isinstance(opts, list):
                if cts.DIALOGS.exists(command):
                    opts = cts.DIALOGS.get_parameter(command, opts)
                else:
                    opts = {}
            else:
                #Adding default parameters
                if cts.DIALOGS.exists(command):
                    defs = cts.DIALOGS.get_parameter(command, '')
                    for k in defs:
                        if k not in opts:
                            opts[k] = defs[k]
            error_status = f_fix(opts, data_to_fix)
            
            if error_status:
                if isinstance(error_status, tuple):
                    if error_status[1] is None:
                        error_status = [error_status[0]]
                    
                print('ERROR', ' '.join(error_status), file=sys.stderr)
                self.summary[command]['error'] = ' '.join(error_status)

        if self.args['debug']:
            import psutil
            self.timings.append([command, time.time() - self.start_time])
            process = psutil.Process(os.getpid())
            memsize = process.memory_info().rss/1024/1024
            self.summary['memsize'].append([command, memsize])
            print(
                "#DEBUG Memory used after {}: {:f} MB ".format(
                    command, memsize
                )
            )
# ==============================================================================
    def sequences(self):
        """ StructureChecking.sequences
        Print canonical and structure sequences in FASTA format
        """
        self._run_method('sequences', None)

    def _sequences_check(self):
        if self.strucm.sequence_data.has_canonical:
            print('Canonical sequence')
            can_seq = self.strucm.sequence_data.get_canonical()
            print(can_seq)
            pdb_seq = self.strucm.sequence_data.get_pdbseq()
            print('Structure sequence')
            print(pdb_seq)
            self.summary['FASTA'] = {
                'canonical': can_seq,
                'structure': pdb_seq
            }
        else:
            print(cts.MSGS['NO_CANONICAL'])

        return {}

    def models(self, opts=None):
        """ StructureChecking.models
        Detect/Select Models. Check only with no options. Options accepted as command-line string, or python dictionary.
        
        Args:
            opts (str | dict - Options dictionary):
                * select (int) - model(s) to select
                * superimpose (bool) - superimpose models
        """
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
        if not self.strucm.has_superimp_models():
            supimp = 'do not'
        print(
            cts.MSGS['MODELS_GUESS'].format(
                supimp,
                self.strucm.models_type['rmsd'],
                mu.MODEL_TYPE_LABELS[self.strucm.models_type['type']]
            )
        )
        return True

    def _models_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            select_model = opts
        else:
            select_model = opts['select']
        
        input_line = ParamInput('Select Model Num', self.args['non_interactive'], set_none='All')
        input_line.add_option_all()
        input_line.add_option_numeric(
            'modelno', [], opt_type='int', min_val=1, max_val=self.strucm.nmodels, multiple=True
        )
        input_option, select_model = input_line.run(select_model)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], select_model

        print(cts.MSGS['SELECT_MODEL'], select_model)

        if input_option != 'all':
            self.strucm.select_model(select_model)

        if (opts['superimpose']):
            self.strucm.superimpose_models()
            print(cts.MSGS['SUPIMP_MODELS'].format(self.strucm.models_type["rmsd"]))
            self.summary['models']['superimp_rmsd'] = self.strucm.models_type["rmsd"]

        if (opts['save_split']): # tag as modified to force save
            self.strucm.modified = True

        self.summary['models']['selected'] = select_model

        return False
# =============================================================================

    def chains(self, opts=None):
        """ StructureChecking.chains
        Detect/Select Chains. Check only with no options. Options accepted as command-line string, or python dictionary.
        
        Args:
            opts (str | dict - Options dictionary): 
                * select: 
                    * **chain_id_list** - List of chains to retain (comma separated, case sensitive), 
                    * **protein** - Select all protein chains,
                    * **na** - Select all NA chains,
                    * **rna** - Select all RNA chains,
                    * **dna** - Select all DNA chains.
        """
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

    def _chains_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            select_chains = opts
        else:
            select_chains = opts['select']
        
        self.summary['chains']['selected'] = {}
        input_line = ParamInput('Select chain', self.args['non_interactive'], set_none='All')
        input_line.add_option_all()
        input_line.add_option_list(
            'type', ['protein'], multiple=False
        )
        input_line.add_option_list(
            'type', ['na'], multiple=False
        )
        input_line.add_option_list(
            'type', ['dna'], multiple=False
        )
        input_line.add_option_list(
            'type', ['rna'], multiple=False
        )
        input_line.add_option_list(
            'chid', sorted(self.strucm.chain_ids), multiple=True, case="sensitive"
        )
        input_line.set_default('All')

        input_option, select_chains = input_line.run(select_chains)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], select_chains

        if input_option == 'all':
            print(cts.MSGS['SELECT_ALL_CHAINS'])
            return False

        self.strucm.select_chains(select_chains)
        print(cts.MSGS['SELECT_CHAINS'], select_chains)
        self.summary['chains']['selected'] = self.strucm.chain_ids

        return False

# =============================================================================
    def inscodes(self):
        """ StructureChecking.inscodes
        Detects residues with insertion codes. No fix provided (yet)
        """
        self._run_method('inscodes', None)

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
        return {'ins_codes_list': ins_codes_list}

    def _inscodes_fix(self, opts, fix_data=None):
        # TODO implement method _inscodes_fix
        if opts['renum']:
            print("--renum option not implemented (yet)")
        return False
# =============================================================================

    def altloc(self, opts=None):
        """ StructureChecking.altloc
        Detect/Select Alternative Locations. Check only with no options. Options accepted as command-line string, or python dictionary.
        
        Args:
            opts (str | dict - Options dictionary):
                * select: 
                    * **occupancy** - select higher occupancy,
                    * **alt_id** - All atoms of the indicated alternative
                    * **list** of res_id:alt_id - Indicate selection per atom                
        """
        self._run_method('altloc', opts)

    def _altloc_check(self): #TODO improve output
        alt_loc_res = self.strucm.get_altloc_residues()
        if not alt_loc_res:
            if not self.args['quiet']:
                print(cts.MSGS['NO_ALTLOC_FOUND'])
            return {}

        print(cts.MSGS['ALTLOC_FOUND'].format(len(alt_loc_res)))

        self.summary['altloc'] = {}

        fix_data = {
            'alt_loc_res' : alt_loc_res,
            'altlocs' : {}
        }

        for res in sorted(alt_loc_res, key=lambda x: x.index):
            rid = mu.residue_id(res)
            print(rid)
            self.summary['altloc'][rid] = {}
            fix_data['altlocs'][res] = sorted(alt_loc_res[res][0].child_dict)
            for atm in alt_loc_res[res]:
                self.summary['altloc'][rid][atm.id] = []
                alt_str = '  {:4}'.format(atm.id)
                for alt in sorted(atm.child_dict):
                    alt_str += ' {} ({:4.2f})'.format(alt, atm.child_dict[alt].occupancy)
                    self.summary['altloc'][rid][atm.id].append({
                        'loc_label': alt,
                        'occupancy': atm.child_dict[alt].occupancy
                    })
                print(alt_str)

        return fix_data

    def _altloc_fix(self, opts, fix_data=None):

        if isinstance(opts, str):
            select_altloc = opts
        else:
            select_altloc = opts['select']
            
        # Prepare the longest possible list of alternatives
        altlocs = []
        max_al_len = 0
        for res in fix_data['altlocs']:
            if len(fix_data['altlocs'][res]) > max_al_len:
                altlocs = fix_data['altlocs'][res]
                max_al_len = len(fix_data['altlocs'][res])

        input_line = ParamInput('Select alternative', self.args['non_interactive'], set_none='All')
        input_line.add_option_all()
        input_line.add_option_list('occup', ['occupancy'])
        input_line.add_option_list('altids', altlocs, case='upper')
        input_line.add_option_list(
            'resnum',
            mu.prep_rnums_list(fix_data['altlocs']),
            opt_type='pair_list',
            list2=altlocs,
            case='sensitive',
            multiple=True
        )
        input_line.set_default('occupancy')

        input_option, select_altloc = input_line.run(select_altloc)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], select_altloc
        
        if input_option != 'all':
            print('Selecting location {}'.format(select_altloc))
            if input_option in ('occup', 'altids'):
                select_altloc = select_altloc.upper()
                to_fix = {
                    res: {
                        'ats': value,
                        'select' : select_altloc
                    } for res, value in fix_data['alt_loc_res'].items()
                }

            elif input_option == 'resnum':
                to_fix = {}
                selected_rnums = {}
                for rsel in select_altloc.split(','):
                    rnum, alt = rsel.split(':')
                    selected_rnums[rnum] = alt
                to_fix = {
                    res : {
                        'ats' : value,
                        'select' : selected_rnums[mu.residue_num(res)]
                    }
                    for res, value in fix_data['alt_loc_res'].items()
                    if mu.residue_num(res) in selected_rnums
                }
            for res in to_fix:
                self.strucm.select_altloc_residue(res, to_fix[res])

        self.summary['altloc']['selected'] = select_altloc

        return False
# =============================================================================

    def metals(self, opts=None):
        """ StructureChecking.metals
        Detect/Remove Metals. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * remove:
                    * **all** - Remove all metal atoms,
                    * **atom_type_list**: Remove all Metals of listed types
                    * **residue_list**: Remove indicated residues
        """
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
            if atm.id not in fix_data['at_groups']:
                fix_data['at_groups'][atm.id] = []
            fix_data['at_groups'][atm.id].append(atm)
            self.summary['metals']['detected'].append(mu.residue_num(res))

        return fix_data

    def _metals_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            remove_metals = opts
        else:
            remove_metals = opts['remove']
        
        input_line = ParamInput("Remove", self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'atids',
            sorted(fix_data['at_groups']),
            case='sensitive',
            multiple=True
        )
        input_line.add_option_list('resids', fix_data['met_rids'], case='sensitive', multiple=True)
        input_line.set_default('All')
        input_option, remove_metals = input_line.run(remove_metals)

        if input_option == "error":
            return cts.MSGS['UNKNOWN_SELECTION'], remove_metals

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'all':
            to_remove = fix_data['met_list']
        elif input_option == 'resids':
            rid_list = remove_metals.split(',')
            to_remove = [
                atm for atm in fix_data['met_list']
                if mu.residue_num(atm.get_parent()) in rid_list
            ]
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
        """ StructureChecking.water
        Detect/Select Remove Water molecules. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * remove: Yes - Remove All Water molecules
        """
        self._run_method('water', opts)

    def _water_check(self):
        wat_list = [
            res
            for res in mu.get_ligands(self.strucm.st, incl_water=True)
            if mu.is_wat(res)
        ]

        if not wat_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_WATERS'])
            return {}

        print(cts.MSGS['WATERS_FOUND'].format(len(wat_list)))
        self.summary['water']['n_detected'] = len(wat_list)

        return {'wat_list': wat_list}

    def _water_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            remove_wat = opts
        else:
            remove_wat = opts['remove']
        
        input_line = ParamInput('Remove', self.args['non_interactive'], set_none='no')
        input_line.add_option_yes_no()
        input_line.set_default('yes')
        input_option, remove_wat = input_line.run(remove_wat)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], remove_wat

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

    def hetatm(self, opts=None):
        """  StructureChecking.hetatm
        Manages hetero atoms. Not implemented yet. See Ligands
        """
        print("Warning: hatatm function not implemented yet, running ligands instead")
        self._run_method('ligands', opts)

    def ligands(self, opts=None):
        """ StructureChecking.ligands
        Detect/Remove Ligands. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * remove:
                    * **all** - Remove all hetatm,
                    * **res_type_list** - Remove Hetatm of given types,
                    * **residue_list** - Remove indicated residues
        """
        self._run_method('ligands', opts)

    def _ligands_check(self):
        lig_list = mu.get_ligands(self.strucm.st, incl_water=False)

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

    def _ligands_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            remove_ligands = opts
        else:
            remove_ligands = opts['remove']
        input_line = ParamInput('Remove', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'byrids', sorted(fix_data['ligand_rids']), multiple=True
        )
        input_line.add_option_list(
            'byresnum', fix_data['ligand_rnums'], case='sensitive', multiple=True
        )
        input_line.set_default('All')
        input_option, remove_ligands = input_line.run(remove_ligands)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], remove_ligands

        self.summary['ligands']['removed'] = {'opt':remove_ligands, 'lst':[]}

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'all':
            to_remove = fix_data['lig_list']
        elif input_option == 'byrids':
            to_remove = [
                res
                for res in fix_data['lig_list']
                if res.get_resname() in remove_ligands.split(',')
            ]
        elif input_option == 'byresnum':
            to_remove = [
                res
                for res in fix_data['lig_list']
                if mu.residue_num(res) in remove_ligands.split(',')
            ]
        rl_num = 0
        for res in to_remove:
            self.summary['ligands']['removed']['lst'].append(mu.residue_id(res))
            self.strucm.remove_residue(res, False)
            rl_num += 1
        self.strucm.update_internals()
        print(cts.MSGS['LIGANDS_REMOVED'].format(remove_ligands, rl_num))
        self.summary['ligands']['n_removed'] = rl_num
        return False
# =============================================================================

    def rem_hydrogen(self, opts=None):
        """ StructureChecking.add_hydrogen
        Remove Hydrogen atoms from structure. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * remove (str): Yes - remove all hydrogen atoms
        """
        self._run_method('rem_hydrogen', opts)

    def _rem_hydrogen_check(self):
        remh_list = mu.get_residues_with_H(self.strucm.st)
        if not remh_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_RESIDUES_H_FOUND'])
            return {}

        print(cts.MSGS['RESIDUES_H_FOUND'].format(len(remh_list)))
        self.summary['rem_hydrogen']['n_detected'] = len(remh_list)
        return {'remh_list': remh_list}

    def _rem_hydrogen_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            remove_h = opts
        else:
            remove_h = opts['remove']

        input_line = ParamInput('Remove hydrogen atoms', self.args['non_interactive'], set_none='no')
        input_line.add_option_yes_no()
        input_line.set_default('yes')

        input_option, remove_h = input_line.run(remove_h)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], remove_h

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
        """ StructureChecking.getss
        Detect SS Bonds. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * mark: 
                    * **all** - Rename all reported cys residues as CYX,
                    * **residue_list** - Rename indicated residues
        """
        self._run_method('getss', opts)

    def _getss_check(self):
        SS_bonds = self.strucm.get_SS_bonds()
        if not SS_bonds:
            if not self.args['quiet']:
                print(cts.MSGS['NO_SS'])
            return {}
        print(cts.MSGS['POSSIBLE_SS'].format(len(SS_bonds)))
        self.summary['getss'] = {'found':[]}
        for ssb in SS_bonds:
            print(
                ' {:12} {:12} {:8.3f}'.format(
                    mu.atom_id(ssb[0]), mu.atom_id(ssb[1]), ssb[2]
                )
            )
            self.summary['getss']['found'].append(
                {
                    'at1':mu.atom_id(ssb[0]),
                    'at2':mu.atom_id(ssb[1]),
                    'dist': round(float(ssb[2]), 4)
                }
            )
        return SS_bonds

    def _getss_fix(self, opts, fix_data=None):
        if not fix_data:
            return False
        if isinstance(opts, str):
            getss_mark = opts
        else:
            getss_mark = opts['mark']

        pairs_list = [
            mu.residue_num(a[0].get_parent()) + "-" + mu.residue_num(a[1].get_parent())
            for a in fix_data
        ]
        input_line = ParamInput('Mark SS', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'bypair', pairs_list, multiple=True
        )
        input_line.set_default('All')
        input_option, getss_mark = input_line.run(getss_mark)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], getss_mark

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        cys_to_mark = []

        for pair in fix_data:
            if (input_option == 'all') or (mu.residue_num(pair[0].get_parent()) + "-" + mu.residue_num(pair[1].get_parent()) in getss_mark.split(',')):
                cys_to_mark.append(pair[0].get_parent())
                cys_to_mark.append(pair[1].get_parent())
        self.summary['getss']['marked'] = [mu.residue_id(a) for a in cys_to_mark]
        self.strucm.mark_ssbonds(cys_to_mark)
        self.strucm.update_internals()
        return False

# =============================================================================
    def amide(self, opts=None):
        """  StructureChecking.amide
        Detect/Fix Amide atoms Assignment. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * fix:
                    * **all** - Fix all residues,
                    * **residue_list** - Fix indicated residues
                * no_recheck (bool) - (False) Do not recheck amide residues after modification.
        """
        self._run_method('amide', opts)

    def _amide_check(self):
        amide_check = self.strucm.check_amide_contacts()
        if 'list' not in amide_check:
            if not self.args['quiet']:
                print(cts.MSGS['NO_AMIDES'])
            return {}
        self.summary['amide']['n_amides'] = len(amide_check['list'])

        if not amide_check['cont_list']:
            if not self.args['quiet']:
                print(cts.MSGS['NO_UNUSUAL_AMIDES'])
            return {}

        print(cts.MSGS['UNUSUAL_AMIDES'].format(len(amide_check['cont_list'])))

        self.summary['amide']['detected'] = []
        for at_pair in sorted(amide_check['cont_list'], key=_key_sort_atom_pairs):
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
        return amide_check
# =============================================================================

    def _amide_fix(self, opts, fix_data=None):
        if not fix_data:
            return False
        if isinstance(opts, str):
            amide_fix = opts
        else:
            amide_fix = opts['fix']
        no_int_recheck = amide_fix is not None or self.args['non_interactive']
        while fix_data:
            input_line = ParamInput('Fix amide atoms', self.args['non_interactive'])
            input_line.add_option_all()
            input_line.add_option_none()
            input_line.add_option_list(
                'resnum',
                sorted(mu.prep_rnums_list(fix_data['res_to_fix'])),
                case='sensitive',
                multiple=True
            )
            input_line.set_default('All')
            input_option, amide_fix = input_line.run(amide_fix)

            if input_option == 'error':
                return cts.MSGS['UNKNOWN_SELECTION'], amide_fix

            if input_option == 'none':
                if self.args['verbose']:
                    print(cts.MSGS['DO_NOTHING'])
                return False
            to_fix = [
                res
                for res in fix_data['res_to_fix']
                if mu.residue_num(res) in amide_fix.split(',') or input_option == 'all'
            ]
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
            if not opts['no_recheck']:
                if not self.args['quiet']:
                    print(cts.MSGS['AMIDES_RECHECK'])
                fix_data = self._amide_check()
                amide_fix = ''
                if no_int_recheck:
                    fix_data = {}
        return False
# =============================================================================

    def chiral(self, opts=None):
        """ StructureChecking.chiral
        Detect/Fix Improper side chain chirality. Check only with no options.  Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * fix:
                    * **All** - Fix all residues
                    * **residue_list** - Fix indicates residues
                * no_check_clashes (bool) - (False) Do not check generated clashes
        """
        self._run_method('chiral', opts)

    def _chiral_check(self):

        chiral_check = self.strucm.check_chiral_sides()

        if 'list' not in chiral_check:
            if not self.args['quiet']:
                print(cts.MSGS['NO_CHIRALS'])
            return {}
        self.summary['chiral']['n_chirals'] = len(chiral_check['list'])

        if not chiral_check['res_to_fix']:
            if not self.args['quiet']:
                print(cts.MSGS['NO_WRONG_CHIRAL_SIDE'])
            return {}

        print(cts.MSGS['WRONG_CHIRAL_SIDE'].format(len(chiral_check['res_to_fix'])))
        self.summary['chiral']['detected'] = []
        for res in chiral_check['res_to_fix']:
            print(' {:10}'.format(mu.residue_id(res)))
            self.summary['chiral']['detected'].append(mu.residue_id(res))

        return chiral_check

    def _chiral_fix(self, opts, fix_data=None):
        if not fix_data:
            return False
        if isinstance(opts, str):
            chiral_fix = opts
        else:
            chiral_fix = opts['fix']
        input_line = ParamInput('Fix chiralities', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum',
            mu.prep_rnums_list(fix_data['res_to_fix']),
            case='sensitive',
            multiple=True
        )
        input_line.set_default('All')
        input_option, chiral_fix = input_line.run(chiral_fix)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], chiral_fix

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option == 'all':
            to_fix = fix_data['res_to_fix']
        else:
            to_fix = [
                res
                for res in fix_data['res_to_fix']
                if mu.residue_num(res) in chiral_fix.split(',')
            ]
        fix_num = 0
        for res in to_fix:
            self.strucm.fix_chiral_chains(res)
            fix_num += 1
        print(cts.MSGS['CHIRAL_SIDE_FIXED'].format(chiral_fix, fix_num))

        if not opts['no_check_clashes']:
            if not self.args['quiet']:
                print(cts.MSGS['CHECKING_CLASHES'])

            self.summary['chiral_clashes'] = self._check_report_clashes(to_fix)

        self.strucm.modified = True
        return False
# =============================================================================

    def chiral_bck(self):
        """  StructureChecking.chiral_bck
        Detect/Fix Improper CA chirality. No fix.
        """
        self._run_method('chiral_bck', None)

    def _chiral_bck_check(self):
        check = self.strucm.get_chiral_bck_list()
        if 'list' not in check:
            if not self.args['quiet']:
                print(cts.MSGS['NO_BCK_CHIRALS'])
            self.summary['chiral_bck']['n_chirals'] = 0
            return {}

        self.summary['chiral_bck']['n_chirals'] = len(check['list'])

        if not check['res_to_fix']:
            if not self.args['quiet']:
                print(cts.MSGS['NO_CHIRAL_BCK_RESIDUES'])
            return {}

        print(cts.MSGS['CHIRAL_BCK_RESIDUES'].format(len(check['res_to_fix'])))
        self.summary['chiral_bck']['detected'] = []
        for res in check['res_to_fix']:
            print(' {:10}'.format(mu.residue_id(res)))
            self.summary['chiral_bck']['detected'].append(mu.residue_id(res))

#       return check
        return {}

#    def _chiral_bck_fix(self, chiral_fix, fix_data=None, check_clashes=True):
#        return False
# TODO chiral_bck_fix
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
#            if self.args['verbose']:
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
    def clashes(self):
        """ StructureChecking.clashes
        Detect steric clashes in groups: Severe, Apolar, Polar Donors, Polar Acceptors, Ionic Positive, Ionic Negative
        """
        self._run_method('clashes', None)

    def _clashes_check(self):
        self.summary['clashes']['detected'] = self._check_report_clashes()
        return False

# =============================================================================
    def fixside(self, opts=None):
        """  StructureChecking.fixside
        Complete side chains (heavy atoms, protein only). Check only with no options.  Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * fix: 
                    * **all** - Fix all residues
                    * **residue_list** - Fix indicated residues
                * no_check_clashes (bool) - (False) Do not check for generated clashes
                * rebuild (bool) - (False) Rebuild side chains using Modeller
        """
        self._run_method('fixside', opts)

    def _fixside_check(self):
        miss_at_list = self.strucm.get_missing_atoms('side')
        extra_at_list = self.strucm.check_extra_atoms()

        if not miss_at_list and not extra_at_list:
            if not self.args['quiet']:
                print("No residues with missing or unknown side chain atoms found")
            return {}

        fix_data = {
            'miss_at_list': miss_at_list,
            'extra_at_list': extra_at_list
        }
        if miss_at_list:
            self.summary['fixside']['detected_missing'] = {}
            print(cts.MSGS['MISSING_SIDE_ATOMS'].format(len(miss_at_list)))
            for r_at in miss_at_list:
                res, at_list = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                self.summary['fixside']['detected_missing'][mu.residue_id(res)] = at_list
        if extra_at_list:
            self.summary['fixside']['detected_unknown'] = {}
            print(cts.MSGS['UNKNOWN_SIDE_ATOMS'].format(len(extra_at_list)))
            for r_at in extra_at_list:
                res, at_list = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                self.summary['fixside']['detected_unknown'][mu.residue_id(res)] = at_list

        return fix_data

    def _fixside_fix(self, opts, fix_data=None):

        if not fix_data:
            return False
        if isinstance(opts, str):
            fix_side = opts
        else:
            fix_side = opts['fix']

        fixside_rnums = [mu.residue_num(r_at[0]) for r_at in fix_data['miss_at_list']]

        input_line = ParamInput('fixside', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum', fixside_rnums, case='sensitive', multiple=True
        )
        input_line.set_default('All')
        input_option_fix, fix_side = input_line.run(fix_side)

        if input_option_fix == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], fix_side

        self.summary['fixside']['fix'] = fix_side

        if input_option_fix == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        if input_option_fix == 'all':
            to_add = fix_data['miss_at_list']
            to_remove = fix_data['extra_at_list']
        else:
            to_add = [
                r_at
                for r_at in fix_data['miss_at_list']
                if mu.residue_num(r_at[0]) in fix_side.split(',')
            ]

            to_remove = [
                r_at
                for r_at in fix_data['extra_at_list']
                if mu.residue_num(r_at[0]) in fix_side.split(',')
            ]

        if not self.args['quiet']:
            print(cts.MSGS['FIXING_SIDE_CHAINS'])

        fix_num = 0
        rem_num = 0
        self.summary['fixside']['fixed'] = []
        self.summary['fixside']['removed'] = []
        fixed_res = []
        if not opts['no_rem_extra']:
            for r_at in to_remove:
                print(mu.residue_id(r_at[0]))
                for at_id in r_at[1]:
                    print("  Removing", at_id)
                    mu.remove_atom_from_res(r_at[0], at_id)
                    rem_num += 1
                self.summary['fixside']['removed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])
        if not opts['rebuild']:
            for r_at in to_add:
                mu.remove_H_from_r(r_at[0], verbose=True)
                self.strucm.fix_side_chain(r_at)
                fix_num += 1
                self.summary['fixside']['fixed'].append(mu.residue_id(r_at[0]))
                fixed_res.append(r_at[0])
        else:
            self.strucm.rebuild_side_chains(to_add)
            fixed_res = [r_at[0] for r_at in to_add]

        print(cts.MSGS['SIDE_CHAIN_FIXED'].format(fix_num))
        self.strucm.fixed_side = True
        self.strucm.modified = True
        # Checking new clashes
        if not opts['no_check_clashes']:
            print(cts.MSGS['CHECKING_CLASHES'])
            self.summary['fixside_clashes'] = self._check_report_clashes(fixed_res)
        return False
# =============================================================================

    def add_hydrogen(self, opts=None):
        """ StructureChecking.add_hydrogen
        Add Hydrogen Atoms to the structure. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * add_mode: 
                    * **auto** - Add hydrogen atom considering pH 7.0.
                    * **pH** (float) - Set explicit pH value.
                    * **list** (str) - Explicit residue list as [*:]HisXXHid.
                * no_fix_side (bool) - (False) Do not fix side chains.
                * keep_h (bool) - (False) Keep original Hydrigen atoms.
                * add_charges FF (str) - Add charges and atom types for the selected FF.
        """
        self._run_method('add_hydrogen', opts)

    def _add_hydrogen_check(self):

        ion_res_list = self.strucm.get_ion_res_list()

        if not ion_res_list:
            if not self.args['quiet']:
                print(cts.MSGS['NO_SELECT_ADDH'])
            return {'ion_res_list':[]}

        fix_data = {
            'ion_res_list': ion_res_list,
        }
        self.summary['add_hydrogen']['detected'] = {}
        print(cts.MSGS['SELECT_ADDH_RESIDUES'].format(len(ion_res_list)))
        for res_type in self.strucm.data_library.residue_codes['protein']:
            residue_list = [
                mu.residue_num(r_at[0])
                for r_at in ion_res_list
                if r_at[0].get_resname() == res_type
            ]
            if residue_list:
                print(' {} {}'.format(res_type, ','.join(residue_list)))

        self.summary['add_hydrogen']['ionic_detected'] = [
            mu.residue_id(r_at[0]) for r_at in ion_res_list
        ]
        # print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
        return fix_data

    def _add_hydrogen_fix(self, opts, fix_data=None):
        if not fix_data:
            return False
        
        if not self.strucm.fixed_side and not opts['no_fix_side']:
            print("WARNING: fixing side chains, override with --no_fix_side")
            self.fixside(['--fix', 'all'])
        
        # Fixing previously marked N and C terms
        self.strucm.revert_terms()

        if isinstance(opts, str):
            add_h_mode = opts
        else:
            add_h_mode = opts['add_mode']

        input_line = ParamInput('Mode', self.args['non_interactive'])
        input_line.add_option_none()
        sel_options = ['auto']
        if fix_data['ion_res_list']:
            sel_options += ['pH', 'list', 'int', 'int_his']
        input_line.add_option_list('selection', sel_options)
        input_line.set_default('auto')
        input_option, add_h_mode = input_line.run(add_h_mode)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], add_h_mode

        if input_option == "none":
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False

        std_ion = self.strucm.data_library.std_ion

        add_h_mode = add_h_mode.lower()

        ion_to_fix = {
            r_at[0]: std_ion[r_at[0].get_resname()]['std']
            for r_at in fix_data['ion_res_list']
        }

        if add_h_mode == 'auto':
            if not self.args['quiet']:
                print('Selection: auto')
        else:
            if add_h_mode == 'ph':
                ph_value = opts['pH']
                input_line = ParamInput("pH Value", self.args['non_interactive'], set_none=7.0)
                input_line.add_option_numeric("pH", [], opt_type="float", min_val=0., max_val=14.)
                input_line.set_default(7.0)
                input_option, ph_value = input_line.run(ph_value)
                if not self.args['quiet']:
                    print('Selection: pH', ph_value)
                for r_at in fix_data['ion_res_list']:
                    res = r_at[0]
                    rcode = res.get_resname()
                    if ph_value <= std_ion[rcode]['pK']:
                        ion_to_fix[res] = std_ion[rcode]['lowpH']
                    else:
                        ion_to_fix[res] = std_ion[rcode]['highpH']
            else:
                if add_h_mode == 'list':
                    ions_list = opts['list']
                    if not self.args['quiet']:
                        print('Selection: list')
                    ions_list = ParamInput(
                        "Enter Forms list as [*:]his22hip",
                        self.args['non_interactive']
                    ).run(ions_list)

                    # Changes in tautomeric forms made as mutationts eg.HisXXHip
                    mutations = self.strucm.prepare_mutations(ions_list)
                    for mut in mutations.mutation_list:
                        for mut_res in mut.mutations:
                            ion_to_fix[mut_res['resobj']] = mut_res['new_id']
                else:
                    if add_h_mode == 'int':
                        if not self.args['quiet']:
                            print('Selection: interactive')
                        res_list = fix_data['ion_res_list']
                    elif add_h_mode == 'int_his':
                        if not self.args['quiet']:
                            print('Selection: int_his')
                        res_list = [
                            r_at for r_at in fix_data['ion_res_list']
                            if r_at[0].get_resname() == 'HIS'
                        ]

                    for r_at in res_list:
                        rcode = r_at[0].get_resname()
                        input_line = ParamInput(
                            "Select residue form for " + mu.residue_id(r_at[0]),
                            self.args['non_interactive']
                        )
                        input_line.add_option_list('list', r_at[1].keys())
                        input_line.default = std_ion[rcode]['std']

                        form = None
                        input_option, form = input_line.run(form)
                        ion_to_fix[r_at[0]] = form.upper()
        self.strucm.add_hydrogens(ion_to_fix, add_charges=opts['add_charges'])
        self.strucm.modified = True
        return False
# =============================================================================

    def mutateside(self, mut_list):
        """ StructureChecking.mutateside
        Mutate side chain with minimal atom replacement. Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * mut (str) - List of mutations
                * no_check_clashes (bool) - (False) Do not check for generated clashes
                * rebuild (bool) - (False) - Optimize new side chains using Modeller
        """
        self._run_method('mutateside', mut_list)

    def _mutateside_check(self):
        # TODO Check _mutateside_check function
        return True

    def _mutateside_fix(self, opts, fix_data=None):
        if isinstance(opts, str):
            mut_list = opts
        else:
            mut_list = opts['mut']

        if opts['na_seq']:
            mut_list = self.strucm.prepare_mutations_from_na_seq(opts['na_seq'])

        input_line = ParamInput('Mutation list', self.args['non_interactive'], set_none='')

        mut_list = input_line.run(mut_list)
        
        if mut_list == '':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return False   
        
        mutations = self.strucm.prepare_mutations(mut_list)

        print(cts.MSGS['MUTATIONS_TO_DO'])
        for mut in mutations.mutation_list:
            print(mut)
        if opts['rebuild']:
            if self.strucm.has_NA:
                print(cts.MSGS['WARN_NOBUILD_NA'])
            mutated_res = self.strucm.rebuild_mutations(mutations)
        else:
            mutated_res = self.strucm.apply_mutations(mutations)
        if not opts['no_check_clashes']:
            if not self.args['quiet']:
                print(cts.MSGS['CHECKING_CLASHES'])

            self.summary['mutateside_clashes'] = self._check_report_clashes(mutated_res)

        self.strucm.modified = True
        return False
#===============================================================================

    def backbone(self, opts=None):
        """ StructureChecking.backbone
        Analyze/Fix main chain missing atoms and fragments (protein only). Check only with no options. Options accepted as command-line string, or python dictionary.
                
        Args:
            opts (str | dict - Options dictionary):
                * fix_atoms (str - Fix missing O, OXT backbone atoms):
                    * **all** - Fix all residues
                    * **residue List** - Fix indicated residues
                * fix_main (str - Fix backbone main chain):
                    * **all** - All detected breaks
                    * **break list** - Indicated breaks
                * add_caps (str - Add ACE and NME residues):
                    * **all** - All detected terminals
                    * **residue_list** - Indicated terminals
                * extra_gap (int) - ('0') Recover addiciontal residues from the model to improve match (experimental)
                * no_recheck (bool) - (False) Do not recheck backbone after fixing 
                * no_check_clashes (bool) - (False) Do not check for generated clashes
        """
        self._run_method('backbone', opts)

    def _backbone_check(self):
        fix_data = {}
        # Residues with missing backbone
        miss_bck_at_list = self.strucm.get_missing_atoms('backbone')
        fix_data['miss_bck_at_list'] = miss_bck_at_list
        if miss_bck_at_list:
            self.summary['backbone']['missing_atoms'] = {}
            print(cts.MSGS['BCK_MISSING_RESIDUES'].format(len(miss_bck_at_list)))
            for r_at in miss_bck_at_list:
                [res, at_list] = r_at
                print(' {:10} {}'.format(mu.residue_id(res), ','.join(at_list)))
                self.summary['backbone']['missing_atoms'][mu.residue_id(res)] = at_list
        else:
            if not self.args['quiet']:
                print(cts.MSGS['NO_BCK_MISSING'])

        #Not bound consecutive residues
        bck_check = self.strucm.get_backbone_breaks()
        if bck_check['bck_breaks_list']:
            print(cts.MSGS['BACKBONE_BREAKS'].format(len(bck_check['bck_breaks_list'])))
            self.summary['backbone']['breaks'] = {'detected': []}
            for brk in bck_check['bck_breaks_list']:
                print(
                    " {:10} - {:10} ".format(
                        mu.residue_id(brk[0]),
                        mu.residue_id(brk[1])
                    )
                )
                self.summary['backbone']['breaks']['detected'].append([
                    mu.residue_id(brk[0]), mu.residue_id(brk[1])
                ])
            fix_data['bck_breaks_list'] = bck_check['bck_breaks_list']
        else:
            fix_data['bck_breaks_list'] = []
            if not self.args['quiet']:
                print(cts.MSGS['NO_BCK_BREAKS'])


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
        else:
            if not self.args['quiet']:
                print(cts.MSGS['NO_BCK_LINKS'])

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
        # TODO move this section to ligands
        if self.strucm.modified_residue_list:
            print(cts.MSGS['MODIF_RESIDUES'])
            self.summary['backbone']['mod_residues'] = []
            for brk in self.strucm.modified_residue_list:
                print(" {:10} ".format(mu.residue_id(brk)))
                self.summary['backbone']['mod_residues'].append(mu.residue_id(brk))
        # Provisional only missing atoms can be fixed
            fix_data['modified_residue_list'] = True
        return fix_data

    def _backbone_fix(self, opts, fix_data=None):

        no_int_recheck = opts['fix_atoms'] is not None or self.args['non_interactive']

        fix_done = not fix_data['bck_breaks_list']
        fixed_main_res = []
        while not fix_done:
            if opts['extra_gap'] is None:
                opts['extra_gap'] = 0
            fixed_main = self._backbone_fix_main_chain(
                opts['fix_chain'],
                fix_data['bck_breaks_list'],
                self.args['modeller_key'],
                opts['extra_gap']
            )
            if not fixed_main:
                fix_done = True
                continue
            else:
                fixed_main_res += fixed_main
            
            self.summary['backbone']['main_chain_fix'] = [mu.residue_id(r) for r in fixed_main_res]
            if fixed_main:
                self.strucm.modified = True

            if no_int_recheck or not fixed_main or opts['no_recheck']:
                fix_done = True
                # Force silent re-checking to update modified residue pointers
                fix_data = self._backbone_check()
            else:
                if not self.args['quiet']:
                    print(cts.MSGS['BACKBONE_RECHECK'])
                fix_data = self._backbone_check()
                fix_done = not fix_data['bck_breaks_list']

        # Add CAPS
        fixed_caps = self._backbone_add_caps(
            opts['add_caps'],
            fix_data['bck_breaks_list']
        )

        self.summary['backbone']['added_caps'] = [mu.residue_id(r) for r in fixed_caps]

        if fixed_caps:
            print('Added caps:', ', '.join(map(mu.residue_num, fixed_caps)))
            self.strucm.modified = True
            fix_data = self._backbone_check()
        else:
            print('No caps added')

        # Add missing atoms
        fixed_res = []
        if fix_data['miss_bck_at_list']:
            fixed_res = self._backbone_fix_missing(
                opts['fix_atoms'],
                fix_data['miss_bck_at_list']
            )
        if not isinstance(fixed_res, list):
            fixed_res = []
        if not isinstance(fixed_caps, list):
            fixed_caps = []
        if not isinstance(fixed_main_res, list):
            fixed_main_res = []
        res_to_check = fixed_res + fixed_caps + fixed_main_res
        if res_to_check and not opts['no_check_clashes']:
            if not self.args['quiet']:
                print(cts.MSGS['CHECKING_CLASHES'])
            self.summary['backbone']['clashes'] = self._check_report_clashes(res_to_check)

        return False

    def _backbone_fix_main_chain(self, fix_main_bck, breaks_list, modeller_key, extra_gap):
        print("Main chain fixes")

        brk_rnums = [
            (mu.residue_num(resp[0]) + '-' + mu.residue_num(resp[1])).replace(' ', '')
            for resp in breaks_list
        ]
        input_line = ParamInput('Fix Main Breaks', self.args['non_interactive'])
        input_line.add_option_none()
        input_line.add_option_all()
        input_line.add_option_list(
            'brk', brk_rnums, case='sensitive', multiple=True
        )
        input_line.set_default('All')
        input_option, fix_main_bck = input_line.run(fix_main_bck)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], fix_main_bck

        self.summary['backbone']['breaks']['option_selected'] = fix_main_bck

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return []

        # Checking for canonical sequence
        if not self.strucm.sequence_data.has_canonical:
            read_ok = False
            while not read_ok:
                input_line = ParamInput(
                    "Enter canonical sequence path (FASTA)",
                    self.args['non_interactive']
                )
                self.args['fasta_seq_path'] = input_line.run(self.args['fasta_seq_path'])
                if not self.args['fasta_seq_path']:
                    print(cts.MSGS['FASTA_MISSING'])

                read_ok = self.strucm.sequence_data.load_sequence_from_fasta(self.args['fasta_seq_path'])
                if not read_ok:
                    self.args['fasta_seq_path'] = None
                if self.args['non_interactive'] and not read_ok:
                    print(cts.MSGS['FASTA_MISSING'])
                    return []
            self.strucm.sequence_data.read_canonical_seqs(self.strucm, False)
            self.strucm.sequence_data.match_sequence_numbering()

        to_fix = [
            rpair
            for rpair in breaks_list
            if (mu.residue_num(rpair[0]) + '-' + mu.residue_num(rpair[1])).replace(' ', '')\
                in fix_main_bck.split(',') or input_option == 'all'
        ]
        return self.strucm.fix_backbone_chain(to_fix, modeller_key, extra_gap)

    def _backbone_add_caps(self, add_caps, bck_breaks_list):
        print("Capping terminal ends")
        term_res = self.strucm.get_term_res()
        term_rnums = [mu.residue_num(p[1]) for p in term_res]
        brk_res = set()

        for r0, r1 in bck_breaks_list:
            brk_res.add(r0)
            brk_res.add(r1)
        true_term = []
        for r in term_res:
            if r[1] not in brk_res:
                true_term.append(r)

        print("True terminal residues: ", ','.join([mu.residue_num(r[1]) for r in true_term]))
        if bck_breaks_list:
            print(
                "Terminal residues from backbone breaks: ",
                ','.join(
                    [
                        mu.residue_num(r0) + '-' + mu.residue_num(r1)
                        for r0, r1 in bck_breaks_list
                    ]
                )
            )

        input_line = ParamInput('Cap residues', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list('terms', ['Terms'])
        if bck_breaks_list:
            input_line.add_option_list('brks', ['Breaks'])
        input_line.add_option_list(
            'resnum', term_rnums, case='sensitive', multiple=True
        )
        input_line.set_default('none')
        input_option, add_caps = input_line.run(add_caps)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], add_caps

        self.summary['backbone']['add_caps'] = add_caps

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return []
        elif input_option == 'terms':
            to_fix = true_term
        elif input_option == 'brks':
            to_fix = [pair for pair in term_res if pair[1] in brk_res]
        else:
            to_fix = [
                pair
                for pair in term_res
                if mu.residue_num(pair[1]) in add_caps.split(',') or input_option == 'all'
            ]

        return self.strucm.add_main_chain_caps(to_fix)

    def _backbone_fix_missing(self, fix_back, fix_at_list):
        print("Fixing missing backbone atoms")
        fixbck_rnums = [mu.residue_num(r_at[0]) for r_at in fix_at_list]
        input_line = ParamInput('Fix bck missing', self.args['non_interactive'])
        input_line.add_option_all()
        input_line.add_option_none()
        input_line.add_option_list(
            'resnum', fixbck_rnums, case='sensitive', multiple=True
        )
        input_line.set_default('none')
        input_option, fix_back = input_line.run(fix_back)

        if input_option == 'error':
            return cts.MSGS['UNKNOWN_SELECTION'], fix_back

        self.summary['backbone']['missing_atoms']['selected'] = fix_back

        if input_option == 'none':
            if self.args['verbose']:
                print(cts.MSGS['DO_NOTHING'])
            return []

        to_fix = [
            r_at
            for r_at in fix_at_list
            if mu.residue_num(r_at[0]) in fix_back.split(',') or input_option == 'all'
        ]

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

        return fixed_res
#===============================================================================

    def cistransbck(self):
        """ StructureChecking.cistransbck
        Analyzes cis-trans dihedrals on backbone atoms
        """
        self._run_method('cistransbck', None)

    def _cistransbck_check(self):
        (cis_backbone_list, lowtrans_backbone_list) = self.strucm.check_cis_backbone()
        if cis_backbone_list:
            self.summary['cistransbck']['cis'] = []
            print(cts.MSGS['CIS_BONDS'].format(len(cis_backbone_list)))
            for lnk in cis_backbone_list:
                res1, res2, dih = lnk
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
                res1, res2, dih = lnk
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
    def _load_structure(self, input_structure_path, fasta_seq_path=None, verbose=True, print_stats=True):
        """ Private. StructureChecking._load_structure
        Prepares Structure Manager and load structure

        Args:
            input_structure_path (str): Path to structure file or pdb:{pdbid]
            fasta_seq_path (str): (None) Path to sequence FASTA file
            verbose (bool): (True) Output progress information.
            print_stats (bool): (True) Print structure statistics
        """

        input_line = ParamInput(
            "Enter input structure path (PDB, mmcif | pdb:pdbid)",
            self.args['non_interactive']
        )
        input_structure_path = input_line.run(input_structure_path)

        strucm = stm.StructureManager(
            input_structure_path,
            self.args['data_library_path'],
            self.args['res_library_path'],
            pdb_server=self.args['pdb_server'],
            cache_dir=self.args['cache_dir_path'],
            file_format=self.args['file_format'],
            fasta_sequence_path=fasta_seq_path
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
    
    def save_structure(self, output_structure_path, rename_terms=False, split_models=False):
        """ StuctureChecking.save_structure
        Saving the current structure in a the output file
        
        Args:
            output_structure_path (str): Path to saved File
            rename_terms (bool): (False) Rename terminal residues as NXXX, CXXX
            split_models (bool): (False) Save models in separated output files
        """
        return self._save_structure(output_structure_path, rename_terms=rename_terms, split_models=split_models)
    
    #Kept for back compatibility
    def _save_structure(self, output_structure_path, rename_terms=False, split_models=False):
        """ Private. StuctureChecking._save_structure
        Saving the current structure in a the output file
        
        Args:
            output_structure_path (str): Path to saved File
            rename_terms (bool): (False) Rename terminal residues as NXXX, CXXX
            split_models (bool): (False) Save models in separated output files
        """
        input_line = ParamInput(
            "Enter output structure path",
            self.args['non_interactive'],
            set_none='fixed_structure.pdb'
        )
        
        output_structure_path = input_line.run(output_structure_path)
        if self.args['output_format']:
            output_format = self.args['output_format']
        else:
            output_format = os.path.splitext(output_structure_path)[1][1:]
        
        if not split_models:
            self.strucm.save_structure(output_structure_path, rename_terms=rename_terms, output_format=output_format)
        else:
            for mod in self.strucm.st.get_models():
                output_path = f'{output_structure_path}_{mod.serial_num}.{output_format}'
                self.strucm.save_structure(output_path, mod_id = mod.id, rename_terms=rename_terms, output_format=output_format)
            print(cts.MSGS["SPLIT_MODELS"])
        
        return output_structure_path

    def _check_report_clashes(self, residue_list=None, contact_types=None):
        if contact_types is None:
            contact_types = stm.ALL_CONTACT_TYPES
        if not residue_list:
            residue_list = self.strucm.all_residues
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
