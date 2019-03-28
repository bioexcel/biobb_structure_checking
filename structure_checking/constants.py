""" Global constants for structure_checking module
"""
from structure_checking.param_input import Dialog

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

MSGS = {
    #management
    'NON_MODIFIED_STRUCTURE': 'Structure not modified, not saving. Override with --force_save',
    'STRUCTURE_SAVED': 'Structure saved on',
    'FORCE_SAVE_STRUCTURE': 'Structure not modified, saving due to --force_save option',
    'JSON_SAVED': 'Summary data saved on',
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
    'UNKNOWN_MODEL': 'Unknown model',
    'SELECT_MODEL': 'Selecting model num.',
    #chains
    'CHAINS_DETECTED': '{} Chain(s) detected',
    'UNKNOWN_CHAINS':   ' {}: Unknown (PROTEIN: {s[0]:4.2f} DNA: {s[1]:4.2f} ' +\
                        'RNA: {s[2]:4.2f} Other: {s[3]:4.2f})',
    'SELECTION NOT VALID': 'Selection not valid',
    'SELECT_ALL_CHAINS': 'Selecting all chains',
    'SELECT_CHAINS': 'Selecting chain(s)',
    #inscodes
    'INSCODES_FOUND': '{} Residues with insertion codes found',
    'NO_INSCODES_FOUND': 'No residues with insertion codes found',
    #altloc
    'ALTLOC_FOUND': 'Detected {} residues with alternative location labels',
    'NO_ALTLOC_FOUND': 'No residues with alternative location labels detected',
    

}