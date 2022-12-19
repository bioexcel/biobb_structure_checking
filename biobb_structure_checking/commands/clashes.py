""" Module supporting clashes command"""

def check(strcheck):
    strcheck.summary['clashes']['detected'] = strcheck.check_report_clashes()
    return False

def fix(strcheck, opts, fix_data):
    pass
