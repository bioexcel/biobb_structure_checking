#!/usr/bin/env python3

"""
 check_structure: Command line replacement for MDWeb structure check
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"


import sys
import pydoc

import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking

BANNER = "===============================================================================\n"\
    "=                   BioBB structure checking utility v{}                   =\n"\
    "=                 A. Hospital, P. Andrio, J.L. Gelpi 2018-20                  =\n"\
    "===============================================================================\n"
def header():
    """ Prints general application headers"""
    return BANNER.format(cts.VERSION)

def main():
    """ Command-line version of MDWeb's structure checking facility (BioBB suite)"""

    base_dir_path = biobb_structure_checking.__path__[0]

    args = cts.CMD_LINE.parse_args()

    if args.command == 'commands':
        help_str = header()
        with open(base_dir_path + "/" + cts.COMMANDS_HELP_PATH) as help_file:
            help_str += help_file.read()
        pydoc.pager(help_str)
        sys.exit(0)

    print(header())

    if '-h' in args.options or '--help' in args.options:
        cts.DIALOGS.get_parameter(args.command, args.options)

    StructureChecking(base_dir_path, vars(args)).launch()

if __name__ == "__main__":
    main()
