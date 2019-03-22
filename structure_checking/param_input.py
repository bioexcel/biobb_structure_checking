"""
  Module to manage command line for simple parameter input
"""
import argparse

def _get_input(value, prompt_text):
    """ get input data """
    while value is None or value == '':
        value = input(prompt_text)
    if isinstance(value, str):
        value = value.replace(' ', '')
    return value

class ParamInput():
    """Class to generate and manage interactive parameter dialogs"""
    def __init__(self, prefix, opt_value, non_interactive=True):
        self.options = []
        self.opt_value = opt_value
        self.non_interactive = non_interactive
        self.prefix = prefix

    def add_option(self, label, opt_list, case=False, opt_type='list', \
            multiple=False, min_val=0, max_val=0, list2=None):
        """Add option to dialog"""
        self.options.append({
            'label':label,
            'opt_list':opt_list,
            'case':case,
            'type':opt_type,
            'min':min_val,
            'max':max_val,
            'multiple':multiple,
            'list2':list2
        })

    def add_option_all(self):
        """add All option to dialog"""
        self.add_option('all', ['All'])

    def add_option_none(self):
        """add None option to dialog"""
        self.add_option('none', ['None'])

    def add_option_yes_no(self):
        """add Yes/No option to dialog"""
        self.add_option('yes', ['Yes'])
        self.add_option('nos', ['No'])

    def run(self):
        """build and execute dialog"""
        prompt_str = self.prefix
        if self.options:
            opt_strs = []
            for opt in self.options:
                if opt['type'] == 'list':
                    opt_strs.append(','.join(opt['opt_list']))
                elif opt['type'] == 'int' or opt['type'] == 'float':
                    if opt['min'] != 0 or opt['max'] != 0:
                        opt_strs.append('{} - {}'.format(opt['min'], opt['max']))
                elif opt['type'] == 'input':
                    opt_strs.append('')
                elif opt['type'] == 'pair_list':
                    str_opt = []
                    for oper in opt['opt_list']:
                        str_opt.append('{}:[{}]'.format(oper, '|'.join(opt['list2'])))
                    opt_strs.append(','.join(str_opt))
                else:
                    opt_strs.append('?')
            prompt_str += ' (' + ' | '.join(opt_strs) + ')'
        prompt_str += ': '
# No options, simple input
        if not self.options:
            if not self.non_interactive:
                self.opt_value = _get_input(self.opt_value, prompt_str)
            return self.opt_value
# Check alternative options
        input_ok = False
        while not input_ok:
            if not self.non_interactive:
                self.opt_value = _get_input(self.opt_value, prompt_str)
            iopt = 0
            if self.opt_value is None:
                return ['error', '']
            else:
                while not input_ok and iopt < len(self.options):
                    opt = self.options[iopt]
                    if (opt['type'] == 'list' or opt['type'] == 'pair_list') and\
                            isinstance(self.opt_value, str):
                        if not opt['multiple']:
                            values = [self.opt_value]
                        else:
                            values = self.opt_value.split(',')
                        group_ok = True
                        for val in values:
                            if opt['type'] == 'pair_list':
                                val_sp = val.split(':')
                                val = val_sp[0]
                            if opt['case'] == 'upper':
                                val = val.upper()
                                group_ok = group_ok and val in opt['opt_list']
                            elif opt['case'] == 'lower':
                                val = val.lower()
                                group_ok = group_ok and val in opt['opt_list']
                            elif opt['case'] == 'sensitive':
                                group_ok = group_ok and val in opt['opt_list']
                            else:
                                group_ok = group_ok and val.lower()\
                                   in list(map(lambda x: x.lower(), opt['opt_list']))
                            if opt['type'] == 'pair_list':
                                group_ok = group_ok and val_sp[1] in opt['list2']
                        input_ok = group_ok
                    elif opt['type'] == 'int' or opt['type'] == 'float':
                        if opt['type'] == "int":
                            self.opt_value = int(self.opt_value)
                        if opt['type'] == "float":
                            self.opt_value = float(self.opt_value)
                        input_ok = self.opt_value >= opt['min'] and self.opt_value <= opt['max']
                    if not input_ok:
                        iopt += 1
                if not input_ok:
                    print("Input not valid")
                    if self.non_interactive:
                        self.options.append({'label':'error'})
                        input_ok = True
                    else:
                        self.opt_value = ''
        return [self.options[iopt]['label'], self.opt_value]

class Dialog():
    """Dialog to complete command options"""
    def __init__(self):
        self.options = {}

    def add_option(self, command, opt_prompt, opt_dest, opt_help, opt_type=str):
        """add Dialog option"""
        self.options[command] = {
            'command':command,
            'prompt': opt_prompt,
            'dest':   opt_dest,
            'help':   opt_help,
            'type':   opt_type
        }

    def get_parameter(self, command, opts):
        """Generate Dialog"""
        dialog = self.options[command]
        if not dialog:
            print("Error: no dialog available for {}".format(command))
            return 1
        opts_parser = argparse.ArgumentParser(prog=dialog['command'])
        opts_parser.add_argument(
            dialog['prompt'],
            dest=dialog['dest'],
            help=dialog['help'],
            type=dialog['type']
        )
        return vars(opts_parser.parse_args(opts))

    def get_dialog(self, command):
        """Get Dialog data"""
        return self.options[command]

    def exists(self, command):
        """Checks whether a Dialog for **command** is defined"""
        return command in self.options
