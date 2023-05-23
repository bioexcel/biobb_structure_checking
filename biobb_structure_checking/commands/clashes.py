""" Module supporting clashes command"""

import biobb_structure_checking.constants as cts

def check(strcheck):
    if strcheck.strucm.st_data.ca_only:
        print(cts.MSGS['CA_ONLY_STRUCTURE'])
        return None

    strcheck.summary['clashes']['detected'] = strcheck.check_report_clashes()
    return False


def fix(strcheck, opts, fix_data):
    pass
