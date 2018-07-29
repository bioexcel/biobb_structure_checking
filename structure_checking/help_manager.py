"""
   Simple module for managing help texts
"""
import os
import pydoc
import sys

class HelpManager():

    def __init__(self, help_dir_path):
        self.help_dir_path = help_dir_path
        self.files = os.listdir(help_dir_path)

    def get(self, what, header=False):
        s = ''
        if header:
            if 'header.hlp' in self.files:
                with open(self.help_dir_path + "/header.hlp") as f:
                    s += f.read()
        if what + '.hlp' in self.files:
            with open(self.help_dir_path + "/" + what + ".hlp") as f:
                s += f.read()
        else:
            print ("Error: no help available for " + what)
            sys.exit()

        return s

    def print_help(self, what, header=False, pager=False):
        if pager:
            pydoc.pager(self.get(what, header))
        else:
            print (self.get(what, header))
