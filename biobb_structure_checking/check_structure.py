#!/usr/bin/env python3

"""
 check_structure: Command line replacement for MDWeb structure check
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"


import os
import sys
import biobb_structure_checking
from biobb_structure_checking.cmd_line import CmdLine
from biobb_structure_checking.help_manager import HelpManager
from biobb_structure_checking.structure_checking import StructureChecking
from biobb_structure_checking.default_settings import DefaultSettings

VERSION = 'v1.0.6'

def main():
    """ Command-line version of MDWeb's structure checking facility (BioBB suite)"""

    sets = DefaultSettings(biobb_structure_checking.__path__[0])
#
    help_manager = HelpManager(sets.help_dir_path)
    cmd_line = CmdLine(defaults={'version': VERSION})
    args = cmd_line.parse_args()

    if args.command == 'commands':
        help_manager.print_help("commands", header=True, pager=True)
        sys.exit(0)

    if 'help' in args.options:
        help_manager.print_help(args.command, header=True)
        sys.exit(0)

    if '-h' in args.options or '--help' in args.options:
        help_manager.print_help('header')
        help_manager.print_help(args.command)
        sys.exit(0)

    help_manager.print_help('header')
    print()

    StructureChecking(sets, vars(args)).launch()

if __name__ == "__main__":
    main()
