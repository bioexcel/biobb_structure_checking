"""
  Module to manage command line arguments
"""

import argparse

class CmdLine():
    def __init__(self, defaults):
        self.argparser = argparse.ArgumentParser (
            description='Basic Structure checking based on MDWeb'
        )

        self.argparser.add_argument (
            'command',
            help="Command to execute (help: checkStruc commands) ",
            default="interactive"
        )

        self.argparser.add_argument(
            '-i', '--input',
            dest='input_structure_path',
            help='Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}'
        )

        self.argparser.add_argument(
            '-o', '--output',
            dest='output_structure_path',
            help='Output structure. Format PDB'
        )

        self.argparser.add_argument(
            '-d', '--debug',
            action="store_true",
            dest="debug",
            help='Print DEBUG information'
        )

        self.argparser.add_argument(
            '--version',
            action='version',
            version="%(prog)s v0.1"
        )

        self.argparser.add_argument(
            '--data_dir',
            dest='data_dir',
            help="Override settings default data dir"
        )

        self.argparser.add_argument(
            '--res_lib',
            dest='res_lib_path',
            help="Override settings default residue library (AMBER prep format)"
        )
        
        self.argparser.add_argument(
            '--json',
            dest='json_output_path',
            help='Cummulated checking results on a json file'
            )
        
        self.argparser.add_argument(
            '--check_only',
            action="store_true",
            dest='check_only',
            help='Perform checks, structures is not modified'
        )
        
        self.argparser.add_argument(
            '--no_interactive',
            action='store_true',
            dest='no_interactive',
            help='Do not prompt for missing parameters'
            )

        self.argparser.add_argument(
            'options',
            nargs=argparse.REMAINDER,
            help="Specific command options"
        )

    def parse_args(self):
        args = self.argparser.parse_args()
        return args

    @classmethod
    def printArgs(cls, args):
        print ("Arguments")
        print ("=========")
        print ("Command:", args.command)
        if args.debug:
            print (' DEBUG mode on')
        print()

