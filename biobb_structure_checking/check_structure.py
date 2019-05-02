#!/usr/bin/env python3

"""
 check_structure: Command line replacement for MDWeb structure check
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"


import sys
from os.path import join as opj
import biobb_structure_checking
from biobb_structure_checking.cmd_line import CmdLine
from biobb_structure_checking.help_manager import HelpManager
from biobb_structure_checking.structure_checking import StructureChecking
import biobb_structure_checking.constants as cts

VERSION = 'v1.0.7a'


def print_header():
    print(
        "===============================================================================\n"
        "=                   MDWeb structure checking utility {}                  =\n"
        "=                 A. Hospital, P. Andrio, J.L. Gelpi 2018-19                  =\n"
        "===============================================================================\n".format(VERSION)
    )

def main():
    """ Command-line version of MDWeb's structure checking facility (BioBB suite)"""

    sets = {'base_dir_path': biobb_structure_checking.__path__[0]}

    #Default locations
    sets['help_dir_path'] = sets['base_dir_path']
    sets['data_dir_path'] = opj(sets['base_dir_path'], cts.DATA_DIR_DEFAULT)
    sets['res_library_path'] = opj(sets['data_dir_path'], cts.RES_LIBRARY_DEFAULT)
    sets['data_library_path'] = opj(sets['data_dir_path'], cts.DATA_LIBRARY_DEFAULT)

    help_manager = HelpManager(sets['help_dir_path'])
    cmd_line = CmdLine(defaults={'version': VERSION})
    args = cmd_line.parse_args()

    if args.command == 'commands':
        help_manager.print_help("commands", header=True, pager=True)
        sys.exit(0)

    print_header()

    if '-h' in args.options or '--help' in args.options:
        cts.DIALOGS.get_parameter(args.command, args.options)

    StructureChecking(sets, vars(args)).launch()

if __name__ == "__main__":
    main()



