"""
    Class for Structure Checking functionality
"""
__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import importlib
import sys
import os
import time
from numpy import sqrt

import biobb_structure_checking.constants as cts

from biobb_structure_checking.io.json_writer import JSONWriter
from biobb_structure_checking.io.param_input import ParamInput, NoDialogAvailableError

import biobb_structure_checking.structure_manager as stm
import biobb_structure_checking.modelling.utils as mu


# Main class
class StructureChecking():
    """
    | biobb_structure_checking.StructureChecking
    | Main class to control structure checking functionality
    | Provides support for to check_structure command line
    | Load directly for Jupyter Notebook or python scripts.

    Args:
        base_dir_path (str): Base directory path where application resides.
        args (dict): Arguments dictionary see https://biobb-structure-checking.readthedocs.io/en/latest/command_line_usage.html.
    """
    def __init__(self, base_dir_path, args):

        if args is None:
            args = {}

        self.args = cts.set_defaults(base_dir_path, args)

        self.summary = {}
        if self.args['debug'] or args['time_limit']:
            import psutil
            self.start_time = time.time()

        if self.args['debug']:
            self.timings = []
            self.summary['elapsed_times'] = {}
            self.summary['memsize'] = []

        try:
            self.strucm = self._load_structure(
                self.args['input_structure_path'],
                self.args['fasta_seq_path']
            )
        except IOError:
            sys.exit(
                "ERROR: fetching/parsing structure from "
                f"{self.args['input_structure_path']}"
            )
        except (
            stm.WrongServerError,
            stm.UnknownFileTypeError,
            stm.ParseError
        ) as err:
            sys.exit(err.message)


        if self.args['debug']:
            self.timings.append(['load', time.time() - self.start_time])
            process = psutil.Process(os.getpid())
            memsize = process.memory_info().rss/1024/1024
            self.summary['memsize'].append(['load', memsize])
            print(f"#DEBUG Memory used after structure load: {memsize:f} MB ")

        if self.args['time_limit'] and self._check_time_limit():
            sys.exit(1)

        # if self.args['atom_limit'] and \
        #         self.strucm.st_data.stats['num_ats'] > self.args['atom_limit']:
        #     sys.exit(
        #         cts.MSGS['ATOM_LIMIT'].format(
        #             self.strucm.st_data.stats['num_ats'],
        #             self.args['atom_limit']
        #         )
        #     )

    def launch(self):
        """ StructureChecking.launch
        Method run from the command line invocation
        """
        if self.args['command'] == 'command_list':
            self.command_list(self.args['options'])
        elif self.args['command'] == 'checkall':
            self.checkall(self.args['options'])
        elif self.args['command'] == 'fixall':
            self.fixall(self.args['options'])
        elif self.args['command'] == 'load':
            if self.args['nocache'] and \
                    not self.args['force_save'] and \
                    not self.args['copy_input']:
                print(
                    "WARNING: load with --nocache will not "
                    "have any effect unless --copy_input is set"
                )
        else:
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
                        split_models='--save_split' in self.args['options']
                    )
                    print(cts.MSGS['STRUCTURE_SAVED'], output_structure_path)
                    self.summary['saved_structure'] = output_structure_path
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
            self.summary['modified_structure'] = self.strucm.modified

        if self.args['debug']:
            total = time.time() - self.start_time
            self.summary['elapsed_times']['total'] = total
            print("#DEBUG TIMINGS")
            print("#DEBUG =======")
            ant = 0.
            for operation, timing in self.timings:
                elapsed = timing - ant
                self.summary['elapsed_times'][operation] = elapsed
                print(
                    f"#DEBUG {operation:15s}: "
                    f"{elapsed:10.4f} s "
                    f"({elapsed / total * 100:6.2f}%)"
                )
                ant = timing
            print(f"#DEBUG TOTAL          : {total:10.4F} s")
            print()
            print("#DEBUG MEMORY USAGE EVOLUTION")
            print("#DEBUG ======================")
            for operation, memsize in self.summary['memsize']:
                print(f"#DEBUG {operation:15s}: {memsize:.2f} MB")

        if self.args['json_output_path'] is not None:
            json_writer = JSONWriter()
            json_writer.data = self.summary
            try:
                json_writer.save(self.args['json_output_path'])
                print(
                    cts.MSGS['JSON_SAVED'],
                    self.args['json_output_path']
                )
            except IOError:
                print(
                    cts.MSGS['JSON_NOT_SAVED'],
                    self.args['json_output_path']
                )

    def print_stats(self, prefix=None):
        """ StructureChecking.print_stats
        Print statistics on the loaded structure

        Args:
            prefix (str): (None) Prefix to add to the output lines
                for identification.
        """
        self.strucm.st_data.calc_stats()
        if prefix is None:
            prefix = ''
        self.strucm.print_stats(prefix)

    def command_list(self, opts):
        """ StructureChecking.command_list
        Manages command_list workflows

        Args:
            opts (str | list(str)): Command options as
                str, file or str list (';' separated).
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
                sys.exit('ERROR: command list not provided and non_interactive')

        if os.path.isfile(op_list):
            command_list = []
            try:
                with open(op_list, "r") as list_file_h:
                    for line in list_file_h:
                        if line == "\n" or line[0:1] == '#':
                            continue
                        command_list.append(line)
            except OSError:
                sys.exit(f"{cts.MSGS['ERROR_OPEN_FILE']} {op_list}")
        else:
            command_list = op_list.split(';')

        i = 1
        for line in command_list:
            if not self.args['quiet']:
                print(f"\nStep {i}: {line}")
            data = line.split()
            command = data[0]
            opts = data[1:]
            if self._run_method(command, opts):
                break
            i += 1

        print(cts.MSGS['COMMAND_LIST_COMPLETED'])

    def checkall(self, opts=None):
        """ StructureChecking.checkall
        Predefined workflow for complete checking
        """
        # Required to allow for interactive run in Notebooks
        old_check_only = self.args['check_only']
        self.args['check_only'] = True

        for meth in cts.AVAILABLE_METHODS:
            if self._run_method(meth, None):
                break

        self.args['check_only'] = old_check_only

    def fixall(self, opts=None):
        """ StructureChecking.fixall
         Fix all using defaults. Not implemented (yet)
        """
        # TODO Implement method fixall
        print("Fixall not implemented (yet)")

    def revert_changes(self):
        """ StructureChecking.revert_changes
            Reload original structure. Used in Pipelines or Notebooks
            to revert changes.
        """
        self.strucm = self._load_structure(
            self.args['input_structure_path'],
            self.args['fasta_seq_path']
        )
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
            importlib.import_module(
                f'biobb_structure_checking.commands.{command}'
            )
            f_check = sys.modules[f'biobb_structure_checking.commands.{command}'].check
            f_fix = sys.modules[f'biobb_structure_checking.commands.{command}'].fix
        except ImportError as e:
            print(command, e)
            sys.exit(cts.MSGS['COMMAND_NOT_FOUND'].format(command))

        if command not in self.summary:
            self.summary[command] = {}

        msg = f"Running {command}."
        if opts:
            if isinstance(opts, list):
                opts_str = ' '.join(opts)
            elif isinstance(opts, dict):
                opts_str = str(opts)
            else:
                opts_str = opts
            msg += f' Options: {opts_str}'
            self.summary[command]['opts'] = opts_str

        if not self.args['quiet'] or self.args['verbose']:
            print(msg.strip())

    # Running checking method
        data_to_fix = f_check(self)
    # Running fix method if needed
        if self.args['check_only'] or opts in (None, ''):
            if self.args['verbose']:
                print(cts.MSGS['CHECK_ONLY_DONE'])
        elif data_to_fix:
            if isinstance(opts, (str, list)):
                if cts.DIALOGS.exists(command):
                    opts = cts.DIALOGS.get_parameter(command, opts)
                else:
                    opts = {}
            else:
                # Adding default parameters
                if cts.DIALOGS.exists(command):
                    defaults = cts.DIALOGS.get_parameter(command, '')
                    for k in defaults:
                        if k not in opts:
                            opts[k] = defaults[k]
            error_status = f_fix(self, opts, data_to_fix)

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
            print(f"#DEBUG Memory used after {command}: {memsize:f} MB ")

        return self.args['time_limit'] and self._check_time_limit()

# ==============================================================================
    def _load_structure(
        self,
        input_structure_path,
        fasta_seq_path=None,
        verbose=True,
        print_stats=True
    ):
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
            nocache=self.args['nocache'],
            copy_dir=self.args['copy_input'],
            fasta_sequence_path=fasta_seq_path,
            nowarn=not self.args['build_warnings'],
            coords_only=self.args['coords_only'],
            overwrite=self.args['overwrite'],
            atom_limit=self.args['atom_limit']
        )

        self.summary['loaded_structure'] = input_structure_path

        if verbose:
            print(cts.MSGS['STRUCTURE_LOADED'].format(input_structure_path))
            strucm.st_data.print_headers()
            print()

        self.summary['headers'] = strucm.st_data.meta

        if print_stats:
            strucm.print_stats()

        self.summary['stats'] = strucm.get_stats()

        return strucm

    def save_structure(
        self,
        output_structure_path,
        rename_terms=False,
        split_models=False
    ):
        """ StuctureChecking.save_structure
        Saving the current structure in a the output file

        Args:
            output_structure_path (str): Path to saved File
            rename_terms (bool): (False) Rename terminal residues as NXXX, CXXX
            split_models (bool): (False) Save models in separated output files
        """
        return self._save_structure(
            output_structure_path,
            rename_terms=rename_terms,
            split_models=split_models
        )

    # Kept for back compatibility
    def _save_structure(
        self,
        output_structure_path,
        rename_terms=False,
        split_models=False
    ):
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
        if output_format == 'mmCif':
            output_format = 'cif'
        if not split_models:
            self.strucm.save_structure(
                output_structure_path,
                rename_terms=rename_terms,
                output_format=output_format,
                keep_resnames=self.args['keep_canonical']
            )
        else:
            for mod in self.strucm.st.get_models():
                output_path = f'{output_structure_path}_{mod.serial_num}.{output_format}'
                self.strucm.save_structure(
                    output_path,
                    mod_id=mod.id,
                    rename_terms=rename_terms,
                    output_format=output_format,
                    keep_resnames=self.args['keep_canonical']
                )
            print(cts.MSGS["SPLIT_MODELS"])

        return output_structure_path

    def check_report_clashes(self, residue_list=None, contact_types=None):
        """ StructureChecking.check_report_clashes
            Check and reports clashes

            Args:
                residue_list (res (list)) : Residues to check
                contact_types (int (list)): Types of contacts to consider
        """
        if contact_types is None:
            contact_types = mu.ALL_CONTACT_TYPES
        if not residue_list:
            residue_list = self.strucm.st_data.all_residues
        return self._clash_report(
            contact_types,
            self.strucm.check_r_list_clashes(residue_list, contact_types)
        )

    def _clash_report(self, contact_types, clash_list):
        summary = {}
        for cls in contact_types:
            summary[cls] = []
            if clash_list[cls]:
                print(
                    cts.MSGS['CLASHES_DETECTED'].format(
                        len(clash_list[cls]),
                        cls
                    )
                )
                for rkey in sorted(
                        clash_list[cls],
                        key=lambda x: mu.key_sort_atom_pairs(clash_list[cls][x])
                    ):
                    print(
                        f" {mu.atom_id(clash_list[cls][rkey][0]):12}"
                        f" {mu.atom_id(clash_list[cls][rkey][1]):12}"
                        f" {sqrt(clash_list[cls][rkey][2]):8.3f} A"
                    )
                    summary[cls].append({
                        'at1':mu.atom_id(clash_list[cls][rkey][0]),
                        'at2':mu.atom_id(clash_list[cls][rkey][1]),
                        'dist': round(float(sqrt(clash_list[cls][rkey][2])), 4)
                    })
            else:
                if not self.args['quiet']:
                    print(cts.MSGS['NO_CLASHES_DETECTED'].format(cls))
        return summary

    def _check_time_limit(self):
        if time.time() - self.start_time > self.args['time_limit']:
            print(cts.MSGS['TIME_LIMIT'].format(self.args['time_limit']), file=sys.stderr)
            return True
        return False
# ==============================================================================
# Entry points for direct notebook calls and Docstrings for commands' help

    def help(self, command=None):
        """ StructureChecking.help
        Provides help on StructureChecking commands

        Args:
            command (str) : (None) Requested command. If empty returns all commands help
        """
        return cts.help(command)

    def sequences(self, opts=None):
        """ StructureChecking.sequences
        Print canonical and structure sequences in FASTA format

        Args:
            opts (str | dict - Options dictionary):
                * output_fasta (str) - File name to output (FASTA format)
        """
        self._run_method('sequences', opts)

    def models(self, opts=None):
        """ StructureChecking.models
        Detect/Select Models. Check only with no options. Options accepted as command-line string, or python dictionary.

        Args:
            opts (str | dict - Options dictionary):
                * select (int) - model(s) to select
                * superimpose (bool) - superimpose models
                * build_complex (bool) - Build a complex from selected models
        """
        self._run_method('models', opts)

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
                * rename:
                    * **auto** - Add first possible label staring on A to unlabeled chains
                    * **label** - Use indicated label
                * renumber:
                    * **auto** - Renumbers all residues from 1 without repeating residue numbers. Chains are preserved but relabelled from A
                    * **str** - Specific renumbering recipe indicated as a list of tasks: [OldChain:]i0[-j0]=[NewChain:]i1. No j0 implies to the end of chain. No chain implies do the transformation in all chains.
                * rebuild: - Creates chain labels and renumbers residues based on backbone connectivity
        """
        self._run_method('chains', opts)

    def inscodes(self, opts=None):
        """ StructureChecking.inscodes
        Detects residues with insertion codes.

        Args:
            opts (str | dict - Options dictionary):
                * renumber (bool): Renumber residues to remove insertion codes.
        """
        self._run_method('inscodes', opts)

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

    def water(self, opts=None):
        """ StructureChecking.water
        Detect/Select Remove Water molecules. Check only with no options. Options accepted as command-line string, or python dictionary.

        Args:
            opts (str | dict - Options dictionary):
                * remove: Yes - Remove All Water molecules
        """
        self._run_method('water', opts)

    def hetatm(self, opts=None):
        """  StructureChecking.hetatm
        Manages hetero atoms. Not implemented yet. See Ligands
        """
        print("Warning: hetatm function not implemented yet, running ligands instead")
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

    def rem_hydrogen(self, opts=None):
        """ StructureChecking.add_hydrogen
        Remove Hydrogen atoms from structure. Check only with no options. Options accepted as command-line string, or python dictionary.

        Args:
            opts (str | dict - Options dictionary):
                * remove (str): Yes - remove all hydrogen atoms
        """
        self._run_method('rem_hydrogen', opts)

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

    def amide(self, opts=None):
        """  StructureChecking.amide
        Detect/Fix Amide atoms Assignment. Check only with no options. Options accepted as command-line string, or python dictionary.

        Args:
            opts (str | dict - Options dictionary):
                * fix:
                    * **all** - Fix all residues,
                    * **residue_list** - Fix indicated residues
                    * **auto** - Find the best combination to minimize amide contacts .
                * no_recheck (bool) - (False) Do not recheck amide residues after modification.
        """
        self._run_method('amide', opts)

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

    def chiral_bck(self):
        """  StructureChecking.chiral_bck
        Detect/Fix Improper CA chirality. No fix.
        """
        self._run_method('chiral_bck', None)

    def clashes(self):
        """ StructureChecking.clashes
        Detect steric clashes in groups: Severe, Apolar, Polar Donors, Polar Acceptors, Ionic Positive, Ionic Negative
        """
        self._run_method('clashes', None)

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

    def backbone(self, opts=None):
        """ StructureChecking.backbone
        Analyze/Fix main chain missing atoms and fragments (protein only). Check only with no options. Options accepted as command-line string, or python dictionary.

        Args:
            opts (str | dict - Options dictionary):
                * fix_atoms (str - Fix missing O, OXT backbone atoms):
                    * **all** - Fix all residues
                    * **residue List** - Fix indicated residues
                * fix_chain (str - Fix backbone main chain):
                    * **all** - All detected breaks
                    * **break list** - Indicated breaks
                * add_caps (str - Add ACE and NME residues):
                    * **all** - All detected terminals
                    * **residue_list** - Indicated terminals
                    * **breaks** - Add caps to backbone breaks
                    * **terms** - Add caps to true terminals
                * extra_gap (int) - ('0') Recover addiciontal residues from the model to improve match (experimental)
                * no_recheck (bool) - (False) Do not recheck backbone after fixing
                * no_check_clashes (bool) - (False) Do not check for generated clashes
        """
        self._run_method('backbone', opts)

    def cistransbck(self):
        """ StructureChecking.cistransbck
        Analyzes cis-trans dihedrals on backbone atoms
        """
        self._run_method('cistransbck', None)
# =======================================================================
