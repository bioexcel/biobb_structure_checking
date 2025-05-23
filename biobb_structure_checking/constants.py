"""
    Global constants for structure_checking module
"""
import argparse
import os
from os.path import dirname
from os.path import join as opj


from biobb_structure_checking.pdbio.param_input import Dialog

VERSION = '3.15.6'

# Default locations and settings
DATA_DIR_DEFAULT_PATH = 'dat'
RES_LIBRARY_DEFAULT_PATH = 'all_residues.in'
DATA_LIBRARY_DEFAULT_PATH = 'data_lib.json'
CACHE_DIR_DEFAULT_PATH = 'tmpPDB'
COMMANDS_HELP_PATH = 'commands.hlp'
PDB_DOWNLOAD_PREFIX = 'https://files.wwpdb.org'
FASTA_DOWNLOAD_PREFIX = 'https://www.rcsb.org/fasta/entry'
ATOM_LIMIT = 1000000
TIME_LIMIT = 3600

ALT_SERVERS = {
    'mmb': 'https://mmb.irbbarcelona.org/api/pdb',
    'bsc': 'http://mdb-login.bsc.es/api/pdb'
}

DEFAULTS = {
    'file_format': 'mmCif',
    'pdb_server': PDB_DOWNLOAD_PREFIX,
    'quiet': False,
    'force_save': False,
    'check_only': False,
    'non_interactive': False,
    'atom_limit': ATOM_LIMIT,
    'mem_check': False,
    'rename_terms': False,
    'keep_canonical': False,
    'verbose': False,
    'debug': False,
    'options': '',
    'modeller_key': None,
    'time_limit': ATOM_LIMIT,
    'nocache': False,
    'copy_input': None,
    'build_warnings': False,
    'coords_only': False,
    'overwrite_cache': False,
    'overwrite': False
}


def set_defaults(base_dir_path, args=None):
    """
    | Constants set_defaults
    | Checks input args and complete with defaults if necessary

    Args:
        base_dir_path (str) : Directory where application resides
        args (dict) : Arguments as passed from the command line
    """
    if args is None:
        args = {}

    data_dir_path = opj(base_dir_path, DATA_DIR_DEFAULT_PATH)

    if 'res_library_path' not in args or args['res_library_path'] is None:
        args['res_library_path'] = opj(
            data_dir_path, RES_LIBRARY_DEFAULT_PATH
        )

    if 'data_library_path' not in args or args['data_library_path'] is None:
        args['data_library_path'] = opj(
            data_dir_path, DATA_LIBRARY_DEFAULT_PATH
        )

    if 'cache_dir_path' not in args or args['cache_dir_path'] is None:
        args['cache_dir_path'] = CACHE_DIR_DEFAULT_PATH

    if 'commands_help_path' not in args or args['commands_help_path'] is None:
        args['commands_help_path'] = opj(data_dir_path, COMMANDS_HELP_PATH)

    if 'fasta_seq_path' not in args:
        args['fasta_seq_path'] = None

    for param, value in DEFAULTS.items():
        if param not in args or args[param] is None:
            args[param] = value

    return args


# Main Command Line Management
CMD_LINE = argparse.ArgumentParser(
    description='Basic Structure checking based on MDWeb'
)

# Input
CMD_LINE.add_argument(
    '-i', '--input',
    dest='input_structure_path',
    help='Input structure. '
         '1) Fetch from remote pdb:{pdbid}[.format] '
         '2) Direct fetch from custom URL '
         '3) Local file: Formats pdb(qt)|pqr|cif. '
         'For 2) and 3) Format assumed from extension.'
         '.gz files are automatically decompressed'
)

CMD_LINE.add_argument(
    '--file_format',
    dest='file_format',
    help='Format for retrieving remote structures '
         '(mmCif(=cif, default)|pdb|xml)',
    choices=['mmCif', 'pdb', 'cif', 'xml']
)

CMD_LINE.add_argument(
    '--coords_only',
    action='store_true',
    help='Ignores chain labels and residue ids. '
         'Used to recover faulty PDB/cif files',
)

CMD_LINE.add_argument(
    '--sequence',
    dest='fasta_seq_path',
    help='Canonical sequence in FASTA format. Accepted pdb:{pdbid} '
         'for downloading. Note that when downloading PDB files, '
         'sequences are also automatically downloaded'
         '.gz files are automatically decompressed'
)

CMD_LINE.add_argument(
    '--pdb_server',
    dest='pdb_server',
    help='Server for retrieving structures (default(ftp at wwPDB)|MMB|BSC), '
        'use alternate server for http(s) download'
)

CMD_LINE.add_argument(
    '--cache_dir',
    dest='cache_dir_path',
    help='Path for structure\'s cache directory (default: ./tmpPDB)'
)

CMD_LINE.add_argument(
    '--nocache',
    dest='nocache',
    action='store_true',
    help='Do not cache downloaded structures'
)

CMD_LINE.add_argument(
    '--overwrite_cache',
    dest='overwrite',
    action='store_true',
    help='Overwrite cached file if any'
)

CMD_LINE.add_argument(
    '--copy_input',
    dest='copy_input',
    help='Copy the downloaded structure in the indicated folder'
)

CMD_LINE.add_argument(
    '--modeller_key',
    dest='modeller_key',
    help='User key for modeller, required for backbone fix, '
         'register at https://salilab.org/modeller/registration.html'
)

# Settings, reference data
CMD_LINE.add_argument(
    '--res_lib',
    dest='res_library_path',
    help="Override settings default residue library (AMBER prep format)",
)

CMD_LINE.add_argument(
    '--data_lib',
    dest='data_library_path',
    help="Override settings default data library"
)

# Output
CMD_LINE.add_argument(
    '-o', '--output',
    dest='output_structure_path',
    help='Output structure. '
         'Formats available cif|pdb|pdbqt|pqr|cmip '
         '(use file extension or --output_format to set format). '
         'Note than cif output format includes atom records only'
)

CMD_LINE.add_argument(
    '--output_format',
    dest='output_format',
    help='Format for the output. When empty, output file extension is used.',
    choices=['pdb', 'pdbqt', 'pqr', 'cmip', 'cif', 'mmCif']
)

CMD_LINE.add_argument(
    '--keep_canonical_resnames',
    action='store_true',
    dest='keep_canonical',
    help='Keep canonical names for ionized residues in output files'
)

CMD_LINE.add_argument(
    '--rename_terms',
    action='store_true',
    dest='rename_terms',
    help='Show protein terminal residues as NXXX, CXXX in output files'
)

CMD_LINE.add_argument(
    '--json',
    dest='json_output_path',
    help='Summary checking results on a json file'
)

CMD_LINE.add_argument(
    '-nv', '--quiet',
    action='store_true',
    dest='quiet',
    help='Minimal output, removing labels and progress info'
)

CMD_LINE.add_argument(
    '-v', '--verbose',
    action='store_true',
    dest='verbose',
    help='Add extra progress info'
)

CMD_LINE.add_argument(
    '--limit',
    dest='atom_limit',
    type=int,
    default=ATOM_LIMIT,
    help=f'Set limit on the number of atoms, 0: nolimit. Default: {ATOM_LIMIT}'
)

CMD_LINE.add_argument(
    '--time_limit',
    dest='time_limit',
    type=int,
    default=TIME_LIMIT,
    help=f'Set limit on the execution time (sec), 0: nolimit. Default: {TIME_LIMIT}'
)

CMD_LINE.add_argument(
    '--debug',
    dest='debug',
    action='store_true',
    help='Add debug information (timings, resources)'
)

CMD_LINE.add_argument(
    '--build_warnings',
    dest='build_warnings',
    action='store_true',
    help='Show structure building warnings (may indicate PDB errors)'
)

CMD_LINE.add_argument(
    '--force_save',
    action='store_true',
    dest='force_save',
    help='Force saving an output file even if no modification'
)

# Operations
CMD_LINE.add_argument(
    '--check_only',
    action='store_true',
    dest='check_only',
    help='Perform checks only, structure is not modified'
)

CMD_LINE.add_argument(
    '--non_interactive',
    action='store_true',
    dest='non_interactive',
    help='Do not prompt for missing parameters'
)

CMD_LINE.add_argument(
    'command',
    help="Command to execute (help: check_structure commands) ",
)

CMD_LINE.add_argument(
    'options',
    nargs=argparse.REMAINDER,
    help="Specific command options"
)

CMD_LINE.add_argument(
    '--version',
    action='version',
    version="%(prog)s " + VERSION
)

# Interactive DIALOGS to complete command_line missing parameters
# ===============================================================
DIALOGS = Dialog()

# DIALOGS.add_option(command, prompt, destination, help_text, type(str))
# Multiple parameters should come as separate lines with a unique "command"

DIALOGS.add_entry('command_list', 'Runs a list of commands')
DIALOGS.add_option(
    'command_list',
    '--list',
    'op_list',
    'Command List File or string (; separated)'
)

DIALOGS.add_entry('models', 'Checks and selects models')
DIALOGS.add_option(
    'models',
    '--select',
    'select',
    'Select model(s) to keep'
)
DIALOGS.add_option(
    'models',
    '--save_split',
    'save_split',
    'Save each model in a separated PDB file',
    'bool'
)
DIALOGS.add_option(
    'models',
    '--superimpose',
    'superimpose',
    'Superimpose models',
    'bool'
)
DIALOGS.add_option(
    'models',
    '--build_complex',
    'build_complex',
    'Build a complex from selected models (Assemblies as models)',
    'bool'
)

DIALOGS.add_entry('chains', 'Checks and selects chains')
DIALOGS.add_option(
    'chains',
    '--select',
    'select',
    'Chains (All | protein | na | dna | rna | Chain list comma separated)'
)
DIALOGS.add_option(
    'chains',
    '--rename',
    'rename',
    'Rename unlabelled chains (auto | label)'
)
DIALOGS.add_option(
    'chains',
    '--renumber',
    'renumber',
    'Renumber residues (auto | [A:]ini0[-fin0]=[B:]ini1)'
)
DIALOGS.add_option(
    'chains',
    '--rebuild',
    'rebuild',
    'Rebuild chains from coordinates',
    'bool'
)
DIALOGS.add_option(
    'chains',
    '--rem_inscodes',
    'rem_inscodes',
    'Remove insertion codes',
    'bool'
)

DIALOGS.add_entry('altloc', 'Checks and selects alternative locations')
DIALOGS.add_option(
    'altloc',
    '--select',
    'select',
    'Select altloc occupancy|alt_id'
)

DIALOGS.add_entry('inscodes', 'Checks residues with insertion codes')
DIALOGS.add_option(
    'inscodes',
    '--renumber',
    'renumber',
    'Renumber chain residues',
    'bool'
)

DIALOGS.add_entry('metals', 'Checks and optionally removes metal atoms (will be deprecated)')
DIALOGS.add_option('metals', '--remove', 'remove', 'Remove Metal ions')

DIALOGS.add_entry('water', 'Checks and optionally removes water molecules')
DIALOGS.add_option('water', '--remove', 'remove', 'Remove All Water molecules')

DIALOGS.add_entry(
    'ligands',
    'Checks and optionally removes ligand residues (will be deprecated)'
)
DIALOGS.add_option('ligands', '--remove', 'remove', 'Remove Ligand residues')

DIALOGS.add_entry('rem_hydrogen', 'Checks and optionally removes hydrogen atoms')
DIALOGS.add_option(
    'rem_hydrogen',
    '--remove',
    'remove',
    'Remove Hydrogen atoms'
)

DIALOGS.add_entry('amide', 'Checks and optionally fixes wrong amide contacts')
DIALOGS.add_option(
    'amide',
    '--fix',
    'fix',
    'Fix Residues (All | None | auto | List)'
)
DIALOGS.add_option(
    'amide',
    '--no_recheck',
    'no_recheck',
    'Re-check after modification',
    'bool'
)

DIALOGS.add_entry('chiral', 'Checks and optionally fixes side chains with wrong chirality')
DIALOGS.add_option('chiral', '--fix', 'fix', 'Fix Residues (All | None | List)')
DIALOGS.add_option(
    'chiral',
    '--no_check_clashes',
    'no_check_clashes',
    'Do not check for new clashes',
    'bool'
)

DIALOGS.add_entry('chiral_bck', 'Checks residues with wrong CA quiral')
# DIALOGS.add_option('chiral_bck', '--fix', 'fix', 'Fix Residues (All | None | List)')

DIALOGS.add_entry('fixside', 'Checks and fixes missing side chain atoms')
DIALOGS.add_option(
    'fixside',
    '--fix',
    'fix',
    'Add missing atoms to side chains (All | None | List)'
)
DIALOGS.add_option(
    'fixside',
    '--no_rem_extra',
    'no_rem_extra',
    'Do not remove unknown atoms',
    'bool'
)
DIALOGS.add_option(
    'fixside',
    '--no_check_clashes',
    'no_check_clashes',
    'Do not check for new clashes',
    'bool'
)
DIALOGS.add_option(
    'fixside',
    '--rebuild',
    'rebuild',
    'Rebuild complete side chain',
    'bool'
)

DIALOGS.add_entry('backbone', 'Checks and fixes several backbone issues')
DIALOGS.add_option(
    'backbone',
    '--fix_atoms',
    'fix_atoms',
    'Add missing Oxygen atoms to backbone (All | None | List)'
)
DIALOGS.add_option(
    'backbone',
    '--fix_chain',
    'fix_chain',
    'Fixes missing main chain segments (All | None | List)'
)
DIALOGS.add_option(
    'backbone',
    '--add_caps',
    'add_caps',
    'Adds ACE and NME caps to missing main chain segments '
    '(All | None | Breaks | Terms)'
)
DIALOGS.add_option(
    'backbone',
    '--no_check_clashes',
    'no_check_clashes',
    'Do not check for new clashes',
    'bool'
)
DIALOGS.add_option(
    'backbone',
    '--extra_gap',
    'extra_gap',
    'Replace extra_gap additional flanking residues',
    int
)
DIALOGS.add_option(
    'backbone',
    '--no_recheck',
    'no_recheck',
    'Do not re-check after modification',
    'bool'
)

DIALOGS.add_entry('mutateside', 'Performs side chain mutations')
DIALOGS.add_option(
    'mutateside',
    '--mut',
    'mut',
    'Mutate side chains (Mutation List as [*:]arg234Thr)'
)
DIALOGS.add_option(
    'mutateside',
    '--no_check_clashes',
    'no_check_clashes',
    'Do not check for generated clashes',
    'bool'
)
DIALOGS.add_option(
    'mutateside',
    '--rebuild',
    'rebuild',
    'Rebuild complete side chain',
    'bool'
)
DIALOGS.add_option(
    'mutateside',
    '--na_seq',
    'na_seq',
    'Mutate DNA duplex to generate sequence'
)

DIALOGS.add_entry(
    'add_hydrogen',
    'Add hydrogen atoms with tautomer/ion selection'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--add_mode',
    'add_mode',
    'Selection mode (None | auto | list | ph | int | int_his )'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--pH',
    'pH',
    'pH (0-14)'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--list',
    'list',
    'Ionic form selection (list as [*:]his234hip)'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--no_fix_side',
    'no_fix_side',
    'Do not fix side chains',
    'bool'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--keep_h',
    'keep_h',
    'Keep original hydrogen atoms',
    'bool'
)
DIALOGS.add_option(
    'add_hydrogen',
    '--add_charges',
    'add_charges',
    'Update atom partial charges and add atom types '
    'from given forcefield (ADT|CMIP)',
    default=''
)

DIALOGS.add_entry('clashes', 'Checks atom clashes')
# DIALOGS.add_option('clashes', '--no_wat', 'discard_wat',
# 'Discard water molecules', 'bool')

DIALOGS.add_entry('getss', 'Checks SS bonds by distance')
DIALOGS.add_option(
    'getss',
    '--mark',
    'mark',
    'Mark Cys pairs as SS bond (All | None | List)'
)

DIALOGS.add_entry('cistransbck', 'Checks or cis peptide bonds')

DIALOGS.add_entry('checkall', 'Runs all checks, no modification')

DIALOGS.add_entry('fixall', 'Fix all found issues with default options')

DIALOGS.add_entry(
    'sequences',
    'Print Canonical and Structure sequences on FASTA format'
)
DIALOGS.add_option(
    'sequences',
    '--output_fasta',
    'output_fasta',
    'Store sequences as FASTA file'
)

# All methods to perform checkall
AVAILABLE_METHODS = [
    'models', 'chains', 'inscodes', 'altloc', 'rem_hydrogen', 'add_hydrogen',
    'water', 'metals', 'ligands', 'getss', 'amide', 'chiral', 'chiral_bck',
    'fixside', 'backbone', 'cistransbck', 'clashes', 'sequences'
]

MSGS = {
    # management
    'NON_MODIFIED_STRUCTURE': 'Structure not modified, not saving. '
                              'Override with --force_save',
    'STRUCTURE_SAVED': 'Structure saved on',
    'FORCE_SAVE_STRUCTURE': 'Structure not modified, saving due to '
                            '--force_save option',
    'JSON_SAVED': 'Summary data saved on',
    'JSON_NOT_SAVED': 'Unable to save JSON data on ',
    'UNKNOWN_SELECTION': 'Unknown selection',
    'DO_NOTHING': 'Nothing to do',
    'ALL_UNDO': 'All changes reverted to original structure',
    'ATOM_LIMIT': 'Number of atoms limit exceeded ({} > {}), '
                  'use --limit to adjust',
    'TIME_LIMIT': 'Execution time limit ({}s) exceeded, aborting, use --time_limit to adjust',
    'CA_ONLY_STRUCTURE': 'CA-Only structure, skipping',
    # command line
    'ERROR_OPEN_FILE': 'Error when opening file',
    'COMMAND_LIST_COMPLETED': 'Command list completed',
    # run method
    'COMMAND_NOT_FOUND': 'Error: {} command unknown or not implemented',
    'FIX_COMMAND_NOT_FOUND': 'Error: {} command fix not implemented',
    'CHECK_ONLY_DONE': 'Running  check_only. Nothing else to do.',
    # sequences
    'NO_CANONICAL': 'Canonical Sequences requires either mmCIF input or '
                    '--sequence',
    'NO_WRITE_FASTA': 'Error while writing FASTA on {}',
    # models
    'MODELS_FOUND': 'Detected {} Model(s)',
    'MODELS_GUESS': 'Models {} superimpose, RMSd: {:8.3f} A, guessed as {} ',
    'SINGLE_MODEL': 'Found Single model',
    'SELECT_MODEL': 'Selecting model num.',
    'SPLIT_MODELS': 'Splitting models for output',
    'SUPIMP_MODELS': 'Models superimposed: final RMSd {:8.3f} A',
    'ADDED_CHAINS_TO_COMPLEX': '{} chains added to new complex',
    # chains
    'CHAINS_DETECTED': 'Detected {} Chain(s)',
    'UNKNOWN_CHAINS': ' {}: Unknown (PROTEIN: {s[0]:4.2f} DNA: {s[1]:4.2f} '
                      'RNA: {s[2]:4.2f} Other: {s[3]:4.2f})',
    'SELECT_ALL_CHAINS': 'Selecting all chains',
    'SELECT_CHAINS': 'Selecting chain(s)',
    # inscodes
    'INSCODES_FOUND': 'Found {} Residues with insertion codes',
    'NO_INSCODES_FOUND': 'Found no residues with insertion codes',
    # altloc
    'ALTLOC_FOUND': 'Detected {} residues with alternative location labels',
    'NO_ALTLOC_FOUND': 'Detected no residues with alternative location labels',
    # metals
    'METALS_FOUND': 'Found {} Metal ions',
    'NO_METALS_FOUND': 'No metal ions found',
    'METALS_REMOVED': 'Metal Atoms removed {} ({:d})',
    # Waters
    'WATERS_FOUND': 'Detected {} Water molecules',
    'WATER_REMOVED': 'Removed {} Water molecules',
    'NO_WATERS': 'No water molecules found',
    # ligands
    'LIGANDS_DETECTED': 'Detected {} Ligands',
    'NO_LIGANDS_FOUND': 'No ligands found',
    'LIGANDS_REMOVED': 'Ligands removed {} ({})',
    # Hydrogens
    'RESIDUES_H_FOUND': 'Detected {} Residues containing H atoms',
    'NO_RESIDUES_H_FOUND': 'No residues with Hydrogen atoms found',
    'REMOVED_H': 'Hydrogen atoms removed from {} residues',
    # SS bonds
    'POSSIBLE_SS': 'Detected {} Possible SS Bonds',
    'NO_SS': 'No SS bonds detected',
    # amide
    'UNUSUAL_AMIDES': 'Found {} unusual contact(s) involving amide atoms',
    'NO_UNUSUAL_AMIDES': 'Found no unusual contact(s) involving amide atoms',
    'NO_AMIDES': 'No amide residues found',
    'AMIDES_FIXED': 'Amide residues fixed {} ({})',
    'AMIDES_RECHECK': 'Rechecking',
    # chiral
    'WRONG_CHIRAL_SIDE': 'Found {} residues with incorrect side-chain chirality',
    'NO_WRONG_CHIRAL_SIDE': 'Found no residues with incorrect side-chain chirality',
    'NO_CHIRALS': 'No chiral side-chains found',
    'CHIRAL_SIDE_FIXED': 'Chiral side chains fixed {} ({})',
    # chiral Backbone
    'CHIRAL_BCK_RESIDUES': 'Found {} residues with incorrect backbone chirality',
    'NO_CHIRAL_BCK_RESIDUES': 'Found no residues with incorrect backbone chirality',
    'NO_BCK_CHIRALS': 'No residues with chiral backbone found',
    # fixside
    'MISSING_SIDE_ATOMS': 'Found {} Residues with missing side chain atoms',
    'NO_MISSING_SIDE_ATOMS': 'Found no residues with missing or unknown side chain atoms',
    'UNKNOWN_SIDE_ATOMS': 'Found {} Residues with unknown atoms',
    'FIXING_SIDE_CHAINS': 'Fixing side chains',
    'SIDE_CHAIN_FIXED': 'Fixed {} side chain(s)',
    # add hydrogens
    'NO_SELECT_ADDH': 'No residues requiring selection on adding H atoms',
    'SELECT_ADDH_RESIDUES': 'Found {} Residues requiring selection '
                            'on adding H atoms',
    # Mutations
    'MUTATIONS_TO_DO': 'Mutations to perform',
    # Backbone
    'BCK_MISSING_RESIDUES': 'Found {} Residues with missing backbone atoms',
    'BACKBONE_BREAKS': 'Found {} Backbone breaks',
    'UNEXPECTED_BCK_LINKS': 'Found Unexpected backbone links',
    'CONSEC_RES_FAR': 'Consecutive residues too far away to be covalently linked',
    'MODIF_RESIDUES': 'Modified residues found',
    'ADDING_BCK_ATOMS': 'Adding missing backbone atoms',
    'BCK_ATOMS_FIXED': 'Fixed {} backbone atom(s)',
    'FASTA_MISSING': 'Canonical sequence unavailable, '
                     'use --sequence seq.fasta, '
                     'skipping backbone rebuilding',
    'NO_BCK_MISSING': 'Found No residues with missing backbone atoms',
    'NO_BCK_LINKS': 'No unexpected backbone links',
    'NO_BCK_BREAKS': 'No backbone breaks',
    'BACKBONE_RECHECK': 'Updating backbone check',
    # cistrasnbck
    'CIS_BONDS': 'Found {} cis peptide bonds',
    'NO_CIS_BONDS': 'No cis peptide bonds found',
    'LOWTRANS_BONDS': 'Found {} trans peptide bonds with unusual omega dihedrals',
    'NO_LOWTRANS_BONDS': 'No trans peptide bonds with unusual omega dihedrals found',
    # clashes
    'CHECKING_CLASHES': 'Checking for steric clashes',
    'CLASHES_DETECTED': '{} {} detected',
    'NO_CLASHES_DETECTED': 'No {} detected',
    # load
    'STRUCTURE_LOADED': 'Structure {} loaded',
    # NA related
    'NO_NA': 'No NA chains found, skipping',
    'WARN_NOBUILD_NA': 'Warning: --rebuild only available for protein chains',
    # Model utils
    'ATOM_NOT_FOUND': 'Warning: atom {:3} not found in {}',
    'NO_BACKBONE_ATOMS': 'Warning: No backbone atoms defined',
    'RESIDUE_NOT_VALID': "Warning: Residue not valid in this context ",
    'NOT_ENOUGH_ATOMS': "Warning: not enough atoms to build {} hydrogen atoms on"
}
# Help handler


def help(command=None):
    """
    | constants help
    | Handler for getting help on commands

    Args:
        command (str) : (None) Command requested, if empty help on all
        commands is provided.
    """
    if not command:
        help_path = opj(
            dirname(__file__),
            DATA_DIR_DEFAULT_PATH,
            COMMANDS_HELP_PATH
        )
        with open(help_path) as help_file:
            print(help_file.read())
    else:
        DIALOGS.get_parameter(command, '', print_help=True)
