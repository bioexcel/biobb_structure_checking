
import argparse

class cmdLine():
    def __init__(self, defaults=[]):
        self.argparser = argparse.ArgumentParser(
            prog='bnsTopology',
        description='Basic topology builder for BNS'
        )

        self.argparser.add_argument(
            '--id',
            action='store',
            dest='id',
            help='PDB Id (for output)'
        )
        self.argparser.add_argument(
            '--debug', '-d',
            action='store_true',
            dest='debug',
            help='Produce DEBUG output'
        )

        self.argparser.add_argument(
            '--useChains',
            action='store_true',
            dest='useChains',
            help='Use PDB file chain ids'
        )

        self.argparser.add_argument(
            '--graphml',
            dest='graphml',
            help='Produce GraphML output file'
        )

        self.argparser.add_argument(
            '--json',
            dest='json',
            help='Produce Json output file'
        )

        self.argparser.add_argument(
            '--bpthres',
            type = float,
            action='store',
            help='BP Score min value ('+str(defaults['BPTHRESDEF'])+')',
            dest='bpthres',
            default = defaults['BPTHRESDEF'],
        )
    
        self.argparser.add_argument(
            '--contacts',
            action='store_true',
            dest='contacts',
            help='Calculate polar contacts between chains'
        )
        
        self.argparser.add_argument(
            '--interface',
            action='store_true',
            dest='interface',
            help='Determine interface residues within a distance cutoffCalculate polar contacts between chains'
        )
        
        self.argparser.add_argument(
            '--intdist',
            type = float,
            action='store',
            help='Interface distance cutoff ('+str(defaults['INTDIST'])+')',
            dest='intdist',
            default = defaults['INTDIST'],
        )

        self.argparser.add_argument(
	    '--limit',
            action='store',
	    help='Max. number of atoms, Default 10000',
            dest='limit',
	    default = defaults['LIMIT']
        )
        
        self.argparser.add_argument(
            '--useModels',
            action='store',
            dest='useModels',
            help='Use all Models [no, auto, force]',
            default='auto'
            )
    
        
        self.argparser.add_argument('pdb_path')

    def parse_args(self):    
        args = self.argparser.parse_args()
        if not args.pdb_path:
            argparser.print_help()
            sys.exit(2)
        return args

    def printArgs(self,args):
        print ('#ARGS: Arguments list')
        print ('#ARGS: ==============')
        print ('#ARGS: pdb_path:', args.pdb_path)
        print ('#ARGS: pdb_id:', args.id)
        print ('#ARGS: input format:', args.format)
        print ('#ARGS: Use PDB Chain ids:', args.useChains)
        print ('#ARGS: Use PDB Models:', args.useModels)
        print ('#ARGS: Get polar contacts:', args.contacts)
        print ('#ARGS: Get Interface residues:', args.interface)
        print ('#ARGS: Interface Distance:', args.intdist)
        print ('#ARGS: bp min score:', args.bpthres)
        print ('#ARGS: Max num atoms:', args.limit)
        print ('#ARGS: JSON output:', args.json)
        print ('#ARGS: GRAPHML output:', args.graphml)
        print()
