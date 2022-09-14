""" Module supporting clashes command"""

import biobb_structure_checking.constants as cts
import biobb_structure_checking.model_utils as mu

def _check(strcheck):
    strcheck.summary['clashes']['detected'] = strcheck._check_report_clashes()
    return False

def _fix(strcheck, opts, fix_data):
    pass
