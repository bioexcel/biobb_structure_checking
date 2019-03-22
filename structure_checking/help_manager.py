"""
   Simple module for managing help texts
"""
import os
import pydoc
import sys

class HelpManager():
    """ Manager for interactive help"""
    def __init__(self, help_dir_path):
        self.help_dir_path = help_dir_path
        self.files = os.listdir(help_dir_path)

    def get(self, what, header=False):
        """get help text"""
        help_str = ''
        if header:
            if 'header.hlp' in self.files:
                with open(self.help_dir_path + "/header.hlp") as header_file:
                    help_str += header_file.read()
        if what + '.hlp' in self.files:
            with open(self.help_dir_path + "/" + what + ".hlp") as help_file:
                help_str += help_file.read()
        else:
            print("Error: no help available for " + what)
            sys.exit()

        return help_str

    def print_help(self, what, header=False, pager=False):
        """print help text"""
        if pager:
            pydoc.pager(self.get(what, header))
        else:
            print(self.get(what, header))
