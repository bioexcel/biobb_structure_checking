""" 
    Global constants for structure_checking module
"""
from biobb_structure_checking.param_input import Dialog

# Interactive DIALOGS to complete command_line missing parameters
DIALOGS = Dialog()

#DIALOGS.add_option(command, prompt, destinmore ation, help_text, type(str))
#Multiple parameters should come as separate lines with a unique "command"
DIALOGS.add_entry('command_list','Runs a list of commands')
DIALOGS.add_option('command_list', '--list', 'op_list', 'Command List File')

DIALOGS.add_entry('models','Checks and selects models')
DIALOGS.add_option('models', '--select', 'select_model', \
    'Select model to keep', int)

DIALOGS.add_entry('chains','Checks and selects chains')
DIALOGS.add_option('chains', '--select', 'select_chains',\
    'Chains (All | Chain list comma separated)')

DIALOGS.add_entry('altloc','Checks and selects alternative locations')
DIALOGS.add_option('altloc', '--select', 'select_altloc', \
    'Select altloc occupancy|alt_id')

DIALOGS.add_entry('inscodes','Checks residues with insertion codes')
DIALOGS.add_option('inscodes', '--renum', 'renum', 'Renumber residues', 'bool')

DIALOGS.add_entry('metals','Checks and optionally removes metal atoms (will be deprecated in v1.1')
DIALOGS.add_option('metals', '--remove', 'remove_metals', 'Remove Metal ions')

DIALOGS.add_entry('water','Checks and optionally removes water molecules')
DIALOGS.add_option('water', '--remove', 'remove_wat', 'Remove All Water molecules')

DIALOGS.add_entry('ligands','Checks and optionally removes ligand residues (will be deprecated in v1.1')
DIALOGS.add_option('ligands', '--remove', 'remove_ligands', 'Remove Ligand residues')

DIALOGS.add_entry('rem_hydrogen','Checks and optionally removes hydrogen atoms')
DIALOGS.add_option('rem_hydrogen', '--remove', 'remove_h', 'Remove Hydrogen atoms')

DIALOGS.add_entry('amide','Checks and optionally fixes wrong amide contacts')
DIALOGS.add_option('amide', '--fix', 'amide_fix', 'Fix Residues (All | None | List)')
DIALOGS.add_option('amide', '--recheck', 'recheck',\
    'Re-check after modification', 'bool')
    
DIALOGS.add_entry('chiral','Checks and optionally fixes side chains with wrong chirality')
DIALOGS.add_option('chiral', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')

DIALOGS.add_entry('chiral_bck','Checks residues with wrong CA quiral')
#DIALOGS.add_option('chiral_bck', '--fix', 'chiral_fix', 'Fix Residues (All | None | List)')

DIALOGS.add_entry('fixside','Checks and fixes missing side chain atoms')
DIALOGS.add_option('fixside', '--fix', 'fix_side',\
    'Add missing atoms to side chains (All | None | List)')
DIALOGS.add_option('fixside', '--no_rem', 'fix_extra',\
    'Do not remove unknown atoms', 'bool')

DIALOGS.add_entry('backbone','Checks and fixes several backbone issues')
DIALOGS.add_option('backbone', '--fix_atoms', 'fix_back',\
    'Add missing Oxygen atoms to backbone (All | None | List)')
DIALOGS.add_option('backbone', '--fill_main', 'fix_main',\
    'Fixes missing main chain segments (All | None | List)')
DIALOGS.add_option('backbone', '--add_caps', 'add_caps',\
    'Adds ACE and NME caps to missing main chain segments (All | None | List)')
DIALOGS.add_option('backbone', '--no_check_clashes', 'no_check_clashes',\
    'Do not check for new clashes', 'bool')
DIALOGS.add_option('backbone', '--no_recheck', 'recheck',\
    'Do not re-check after modification', 'bool')

DIALOGS.add_entry('mutateside','Performs side chain mutations')        
DIALOGS.add_option('mutateside', '--mut', 'mut_list',\
    'Mutate side chains (Mutation List as [*:]arg234Thr)')
DIALOGS.add_option('mutateside', '--no_check_clashes', 'no_check_clashes',\
    'Do not check for generated clashes', 'bool')
    
DIALOGS.add_entry('add_hydrogen','Add hydrogen atoms with tautomer/ion selection')
DIALOGS.add_option('add_hydrogen', '--add_mode', 'mode',\
    'Selection mode (None | auto | list | ph | int | int_his )')
DIALOGS.add_option('add_hydrogen', '--pH', 'pH_value',\
    'pH (0-14)')
DIALOGS.add_option('add_hydrogen', '--list', 'ions_list',\
    'Ionic form selection (list as [*:]his234hip)')
DIALOGS.add_option('add_hydrogen', '--no_fix_side', 'no_fix_side',\
    'Do not fix side chains', 'bool')
DIALOGS.add_option('add_hydrogen', '--keep_h', 'keep_h',\
    'Keep original hydrogen atoms', 'bool')

DIALOGS.add_entry('clashes', 'Checks atom clashes')
#DIALOGS.add_option('clashes', '--no_wat', 'discard_wat', 'Discard water molecules', 'bool')

DIALOGS.add_entry('getss', 'Checks SS bonds by distance')
DIALOGS.add_entry('cistransbck', 'Checks or cis peptide bonds')
DIALOGS.add_entry('checkall', 'Runs all checks, no modification')
DIALOGS.add_entry('fixall', 'Fix all found issues with default options')

AVAILABLE_METHODS = [
    'models', 'chains', 'inscodes', 'altloc', 'rem_hydrogen', 'add_hydrogen',
    'water', 'metals', 'ligands', 'getss', 'amide', 'chiral', 'chiral_bck',
    'fixside', 'backbone', 'cistransbck', 'clashes']

MSGS = {
    #management
    'NON_MODIFIED_STRUCTURE': 'Structure not modified, not saving. Override with --force_save',
    'STRUCTURE_SAVED': 'Structure saved on',
    'FORCE_SAVE_STRUCTURE': 'Structure not modified, saving due to --force_save option',
    'JSON_SAVED': 'Summary data saved on',
    'JSON_NOT_SAVED': 'Unable to save JSON data on ',
    'UNKNOWN_SELECTION': 'Unknown selection',
    'DO_NOTHING': 'Nothing to do',
    #command line
    'ERROR_OPEN_FILE': 'Error when opening file',
    'COMMAND_LIST_COMPLETED': 'Command list completed',
    #run method
    'COMMAND_NOT_FOUND': 'Error: {} command unknown or not implemented',
    'FIX_COMMAND_NOT_FOUND': 'Error: {} command fix not implemented',
    'CHECK_ONLY_DONE': 'Running  check_only. Nothing else to do.',
    #models
    'MODELS_FOUND': '{} Model(s) detected',
    'MODELS_GUESS': 'Models {} superimpose, RMSd: {:8.3f} A, guessed as {} ',
    'SINGLE_MODEL': 'Single model found',
    'SELECT_MODEL': 'Selecting model num.',
    #chains
    'CHAINS_DETECTED': '{} Chain(s) detected',
    'UNKNOWN_CHAINS':   ' {}: Unknown (PROTEIN: {s[0]:4.2f} DNA: {s[1]:4.2f} ' +\
                        'RNA: {s[2]:4.2f} Other: {s[3]:4.2f})',
    'SELECT_ALL_CHAINS': 'Selecting all chains',
    'SELECT_CHAINS': 'Selecting chain(s)',
    #inscodes
    'INSCODES_FOUND': '{} Residues with insertion codes found',
    'NO_INSCODES_FOUND': 'No residues with insertion codes found',
    #altloc
    'ALTLOC_FOUND': 'Detected {} residues with alternative location labels',
    'NO_ALTLOC_FOUND': 'No residues with alternative location labels detected',
    #metals
    'METALS_FOUND': '{} Metal ions found',
    'NO_METALS_FOUND' : 'No metal ions found',
    'METALS_REMOVED': 'Metal Atoms removed {} ({:d})',
    # Waters
    'WATERS_FOUND': '{} Water molecules detected',
    'WATER_REMOVED': '{} Water molecules removed',
    'NO_WATERS': 'No water molecules found',
    # ligands
    'LIGANDS_DETECTED' : '{} Ligands detected',
    'NO_LIGANDS_FOUND': 'No ligands found',
    'LIGANDS_REMOVED' : 'Ligands removed {} ({})',
    # Hydrogens
    'RESIDUES_H_FOUND' : '{} Residues containing H atoms detected',
    'NO_RESIDUES_H_FOUND': 'No residues with Hydrogen atoms found',
    'REMOVED_H': 'Hydrogen atoms removed from {} residues',
    # SS bonds
    'POSSIBLE_SS': '{} Possible SS Bonds detected',
    'NO_SS': 'No SS bonds detected',
    #amide
    'UNUSUAL_AMIDES': '{} unusual contact(s) involving amide atoms found',
    'NO_UNUSUAL_AMIDES': 'No unusual contact(s) involving amide atoms found',
    'NO_AMIDES': 'No amide residues found',
    'AMIDES_FIXED': 'Amide residues fixed {} ({})',
    'AMIDES_RECHECK': 'Rechecking',
    #chiral
    'WRONG_CHIRAL_SIDE': '{} residues with incorrect side-chain chirality found',
    'NO_WRONG_CHIRAL_SIDE': 'No residues with incorrect side-chain chirality found',
    'NO_CHIRALS' : 'No chiral side-chains found',
    'CHIRAL_SIDE_FIXED': 'Chiral side chains fixed {} ({})',
    #chiral Backbone
    'CHIRAL_BCK_RESIDUES': '{} residues with incorrect backbone chirality found',
    'NO_CHIRAL_BCK_RESIDUES': 'No residues with incorrect backbone chirality found',
    'NO_BCK_CHIRALS': 'No residues with chiral backbone found',
    #fixside
    'MISSING_SIDE_ATOMS': '{} Residues with missing side chain atoms found',
    'UNKNOWN_SIDE_ATOMS': '{} Residues with unknown atoms found',
    'FIXING_SIDE_CHAINS': '"Fixing side chains"',
    'SIDE_CHAIN_FIXED': 'Fixed {} side chain(s)',
    #add hydrogens
    'NO_SELECT_ADDH' : 'No residues requiring selection on adding H atoms',
    'SELECT_ADDH_RESIDUES' : '{} Residues requiring selection on adding H atoms',
    #Mutations
    'MUTATIONS_TO_DO': 'Mutations to perform',
    #Backbone
    'BCK_MISSING_RESIDUES': '{} Residues with missing backbone atoms found',
    'BACKBONE_BREAKS': '{} Backbone breaks found',
    'UNEXPECTED_BCK_LINKS': 'Unexpected backbone links found',
    'CONSEC_RES_FAR': 'Consecutive residues too far away to be covalently linked',
    'MODIF_RESIDUES': 'Modified residues found',
    'ADDING_BCK_ATOMS': 'Adding missing backbone atoms',
    'BCK_ATOMS_FIXED': 'Fixed {} backbone atom(s)',
    #cistrasnbck
    'CIS_BONDS': '{} cis peptide bonds',
    'NO_CIS_BONDS': 'No cis peptide bonds found',
    'LOWTRANS_BONDS': '{} trans peptide bonds with unusual omega dihedrals',
    'NO_LOWTRANS_BONDS': 'No trans peptide bonds with unusual omega dihedrals found',
    #clashes
    'CHECKING_CLASHES' : 'Checking for steric clashes',
    'CLASHES_DETECTED': '{} Steric {} clashes detected',
    'NO_CLASHES_DETECTED': 'No {} clashes detected',
    #load
    'STRUCTURE_LOADED': 'Structure {} loaded'
}
