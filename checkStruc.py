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

def main():

    sets = DefaultSettings(os.path.dirname(__file__))
#
    help = HelpManager(sets.help_dir_path)
    cmd_line = CmdLine(defaults={})
    args = cmd_line.parse_args()

    if args.command == 'commands':
        help.print_help("general", header=True, pager=True)
        sys.exit(0)

    if 'help' in args.options:
        help.print_help(args.command, header=True)
        sys.exit(0)

    
    StructureChecking(vars(args)).launch(sets)

if __name__ == "__main__":
    main()
