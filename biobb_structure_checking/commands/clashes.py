""" Module supporting clashes command"""

def _check(strcheck):
    strcheck.summary['clashes']['detected'] = strcheck._check_report_clashes()
    return False

def _fix(strcheck, opts, fix_data):
    pass
