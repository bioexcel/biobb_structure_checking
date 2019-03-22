#! /usr/bin/python3

"""
 structureCheck: Command line replacement for MDWeb structure check
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"


import os
import sys

from structure_checking.cmd_line import CmdLine
from structure_checking.help_manager import HelpManager
from structure_checking.structure_checking import StructureChecking
from structure_checking.default_settings import DefaultSettings

VERSION = 'v1.0.5'

def main():
    """Command-line version of MDWeb's structure checking facility"""

    sets = DefaultSettings(os.path.dirname(__file__))
#
    help_manager = HelpManager(sets.help_dir_path)
    cmd_line = CmdLine(defaults={'version': VERSION})
    args = cmd_line.parse_args()

    if args.command == 'commands':
        help_manager.print_help("general", header=True, pager=True)
        sys.exit(0)

    if 'help' in args.options:
        help_manager.print_help(args.command, header=True)
        sys.exit(0)

    StructureChecking(sets, vars(args)).launch()

if __name__ == "__main__":
    main()
