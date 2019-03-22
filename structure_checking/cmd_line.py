"""
  Module to manage command line arguments, adapted to the application
"""

import argparse

class CmdLine():
    """Class to manage command line input parameters"""
    def __init__(self, defaults):
        self.argparser = argparse.ArgumentParser(
            description='Basic Structure checking based on MDWeb'
        )
## Input
        self.argparser.add_argument(
            '-i', '--input',
            dest='input_structure_path',
            help='Input structure. Formats PDB|mmCIF. Remote pdb:{pdbid}'
        )

        self.argparser.add_argument(
            '--file_format',
            action='store',
            dest='file_format',
            default='mmCif',
            help='Format for retrieving structures (default=mmCif|pdb|xml)'
        )

        self.argparser.add_argument(
            '--cache_dir',
            action='store',
            dest='cache_dir',
            default='tmpPDB',
            help='Path for structure\'s cache directory (default: ./tmpPDB)'
        )

        self.argparser.add_argument(
            '--pdb_server',
            action='store',
            dest='pdb_server',
            default='ftp://ftp.wwpdb.org',
            help='Server for retrieving structures (default|MMB)'
        )


#output
        self.argparser.add_argument(
            '-o', '--output',
            dest='output_structure_path',
            help='Output structure. Format PDB'
        )
        self.argparser.add_argument(
            '--json',
            dest='json_output_path',
            help='Summary checking results on a json file'
        )

        self.argparser.add_argument(
            '--quiet',
            action="store_true",
            dest='quiet',
            help='Reduces output, removing labels and progress info'
        )

        self.argparser.add_argument(
            '--force_save',
            action='store_true',
            dest='force_save',
            help='Force saving an output file even if no modification'
        )

#Settings, reference data
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
            '--data_lib',
            dest='data_library_path',
            help="Override settings default data library"
        )

#Operations

        self.argparser.add_argument(
            '--check_only',
            action="store_true",
            dest='check_only',
            help='Perform checks only, structure is not modified'
        )

        self.argparser.add_argument(
            '--non_interactive',
            action='store_true',
            dest='non_interactive',
            help='Do not prompt for missing parameters'
        )

        self.argparser.add_argument(
            'command',
            help="Command to execute (help: checkStruc commands) ",
            default="interactive"
        )

        self.argparser.add_argument(
            'options',
            nargs=argparse.REMAINDER,
            help="Specific command options"
        )

        self.argparser.add_argument(
            '--version',
            action='version',
            version="%(prog)s " + defaults['version']
        )


    def parse_args(self):
        """Parse command line arguments"""
        args = self.argparser.parse_args()
        return args

    @classmethod
    def print_args(cls, args):
        """Print list of arguments"""
        print("Arguments")
        print("=========")
        print("Command:", args.command)
        print()
