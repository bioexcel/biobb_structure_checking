"""
  Module to manage command line for simple parameter input
"""

class ParamInput():
    def __init__(self, prefix, opt_value, non_interactive=True):
        self.options=[]        
        self.opt_value=opt_value
        self.non_interactive = non_interactive
        self.prefix = prefix
        
    def add_option(self, label, opt_list, case=False, type='list', multiple=False, min=0, max=0):
        self.options.append({'label':label,'opt_list':opt_list,'case':case,'type':type, 'min':min, 'max':max, 'multiple':multiple})
    
    def add_option_all(self):
        self.add_option('all', ['All'])

    def add_option_none(self):
        self.add_option('none', ['None'])
        
    def add_option_yes_no(self):
        self.add_option('yes', ['Yes'])
        self.add_option('nos', ['No'])
        
    def run(self):
        help_str = self.prefix
        opt_strs=[]
        for opt in self.options:
            if opt['type'] == 'list':
                opt_strs.append(','.join(opt['opt_list']))
            elif opt['type'] == 'int' or opt['type'] == 'float':
                if opt['min'] != 0 or opt['max'] != 0:
                    opt_strs.append('{} - {}'.format(opt['min'], opt['max']))
            elif opt['type'] == 'input':
                opt_strs.append('')
            else: 
                opt_strs.append('?')
                
        help_str += ' ('+' | '.join(opt_strs)+'): '
        input_ok = False
        while not input_ok:
            if not self.non_interactive:
                self.opt_value = _check_parameter(self.opt_value, help_str)
            iopt=0
            if self.opt_value is None:
                return ['error','']
            else:
                while not input_ok and iopt < len(self.options):
                    opt = self.options[iopt]
                    if opt['type'] == 'list':
                        if not opt['multiple']:
                            values = [self.opt_value]
                        else:
                            values = self.opt_value.split(',')
                        for val in values:
                            group_ok = False
                            if opt['case'] == 'upper':
                                val=val.upper()
                                group_ok = group_ok and val in opt['opt_list']
                            elif opt['case'] == 'lower':
                                val=val.lower()
                                group_ok = group_ok and val in opt['opt_list']
                            elif opt['case'] == 'sensitive':
                                group_ok = val in opt['opt_list']
                            else:
                                group_ok = val.lower() in list(map(lambda x: x.lower(),opt['opt_list']))
                        input_ok = group_ok
                    elif opt['type'] == 'int': 
                        self.opt_value = int(self.opt_value)
                        input_ok = self.opt_value >= opt['min'] and self.opt_value <= opt['max']
                    if not input_ok:
                        iopt +=1
                if not input_ok:
                    print ("Input not valid")
                    self.opt_value = ''
                    if self.non_interactive:
                        self.options.append({'label':'error'})
                        input_ok=True
        return [self.options[iopt]['label'],self.opt_value]
    
def _check_parameter (opts_param, input_text):
    while opts_param is None or opts_param == '':
        opts_param = input (input_text)
    if opts_param is str:
        opts_param = opts_param.replace(' ', '')
    return opts_param
