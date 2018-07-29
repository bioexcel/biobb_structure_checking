# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import argparse
import re
import settings as sets
from structure_checking.help_manager import HelpManager
from structure_checking.structure_manager import StructureManager
import structure_checking.util as util
import sys



class StructureChecking():
    def __init__(self, args):
        self.args = args

    def launch(self):
        help = HelpManager(sets.help_dir_path)
        help.print_help('header')

        try:
            f = getattr(self, self.args.command)
        except AttributeError:
            print ("Error: command unknown or not implemented")
            sys.exit(1)

        f(self.args.options)

        self._save_structure()
        print ("Structure saved on", self.args.output_structure_path)


    def command_list(self, options):

        opts = _get_parameter(options, "command_list", "--list", "op_list", "Command List File")

        opts.op_list = _check_parameter(opts.op_list, "Command list file: ")

        try:
            fh = open(opts.op_list, "r")

        except OSError:
            print ("Error opening", opts.op_list)
            sys.exit(1)

        print ("Running command_list from", opts.op_list)

        self._load_structure(self.args.input_structure_path)

        i = 1
        for line in fh:
            print ("Step " + str(i) + ":", line)
            data = line.split()
            command = data[0]
            options = data[1:]
            try:
                f = getattr(self, command)
            except AttributeError:
                print ("Error: command unknown or not implemented")
                sys.exit(1)
            f(options)
            i += 1

        print ("Command list completed")


    def models(self, options):

        opts = _get_parameters(options, "models", "--select_model", "select_model", "Select model to keep", int)
        print ("Running models", ' '.join(options))
        self._load_structure(self.args.input_structure_path)
        self.nmodels = self.struc_man.get_nmodels()
        print ("Models detected: ", self.nmodels)

        if self.nmodels > 1:
            opts.select_model = int (_check_parameter(opts.select_model, "Select Model Num [1-" + str(self.nmodels) + "]: "))
            if opts.select_model < 1 or opts.select_model > self.nmodels:
                print ("Error: unknown model", opts.select_model, file=sys.stderr)
                opts.select_model = ''
            else:
                print ("Selecting model", opts.select_model)
                self.struc_man.select_model(opts.select_model)
                self.nmodels = self.struc_man.get_nmodels()
        else:
            print ("Nothing to do")

        print()

    def chains(self, options):

        opts = _get_parameter(options, "chains", '--select_chains', "select_chains", "Chains (*| list comma separated)")
        print ("Running chains", ' '.join(options))
        self._load_structure(self.args.input_structure_path)
        self.chain_ids = self.struc_man.get_chain_ids()
        print ("Chains detected: ", len(self.chain_ids), '(', ', '.join(self.chain_ids), ')')

        if len(self.chain_ids) > 1:
            opts.select_chains = _check_parameter(opts.select_chains, "Select chain id(s) or * for all [*," + ",".join(self.chain_ids) + "]: ")
            if opts.select_chains != '*':
                self.struc_man.select_chains(opts.select_chains)
                print ("Selecting chain(s) ", ','.join(ch_ok))
                self.chain_ids = self.struc_man.get_chain_ids()
            else:
                print ("Selecting all chains")
        else:
            print ("Nothing to do")
        print ()


    def altloc(self, options):

        opts = _get_parameters(options, "altloc", '--select_altloc', "select_altloc", "select altloc occupancy|alt_id")

        print ("Running altlocs", ' '.join(options))

        self._load_structure(self.args.input_structure_path)

        alt_loc_res = self.struc_man.get_altloc_residues()

        if len(alt_loc_res) > 0:
            print ("Detected alternative locations")

            for r in alt_loc_res.keys():
                if len(alt_loc_res[r]) > 1:
                    print (r, ":")
                    altlocs = sorted(alt_loc_res[r][0].child_dict.keys())
                    for at in alt_loc_res[r]:
                        s = "  " + at.id
                        for alt in sorted(at.child_dict.keys()):
                            s += " " + alt + "(" + str(at.child_dict[alt].occupancy) + ")"
                        print (s)
                    ok = False
                    while not ok:
                        opts.select_altloc = _check_parameter('', "Select alternative (occupancy, " + ','.join(altlocs) + '): ')
                        ok = opts.select_altloc in altlocs or opts.select_altloc == 'occupancy'
                        if not ok:
                            print ("Error: Unknown ", opts.select_altloc, file=sys.stderr)
                            opts.select_altloc = ''

                    print ("Selecting location", opts.select_altloc)
                    self.struc_man.select_altloc_residues(r, opts.select_altloc)
        else:
            print ("No alternative locations detected")

    def metals (self, options):

        opts = _get_parameters(options, "metals", '--remove', 'remove_metals', 'Remove Metal ions')

        print ("Running metals", ' '.join(options))

        self._load_structure(self.args.input_structure_path)

        [met_list, met_rid, at_groups] = self.struc_man.get_metals(sets.metal_ats)

        if len(met_list) > 1:
            print ("Metal ions found")
            for at in met_list:
                print ("  ", util.atomid(at))

            ok = False
            resids = False
            atids = False
            while not ok:
                opts.remove_metals = _check_parameter(opts.remove_metals,
                                                      'Remove (All | None | any of  ' + ','.join(sorted(at_groups.keys())) + ' | any of ' + ','.join(met_rid) + ': ')
                ok = opts.remove_metals in ['All', 'None']
                if not ok:
                    ok = True
                    for r in opts.remove_metals.split(','):
                        ok = ok and r in met_rid
                    resids = ok

                if not ok:
                    ok = True
                    for atid in opts.remove_metals.split(','):
                        ok = ok and atid in at_groups.keys()
                    atids = ok

                if not ok:
                    print ("Error unknown", opts.remove_metals, file=sys.stderr)
                    opts.remove_metals = ''

            if opts.remove_metals == 'None':
                to_remove = []

            elif opts.remove_metals == 'All':
                to_remove = met_rid

            elif resids:
                to_remove = opts.remove_metals.split(',')

            elif atids:
                to_remove = []
                for atid in opts.remove_metals.split(','):
                    to_remove.extend(at_groups[atid])

            selt.struc_man.remove_metals(met_list, to_remove)
            print ("Metals removed", opts.remove_metals)
        else:
            print ("No metal ions found")


#===============================================================================

    def _load_structure(self, input_structure_path):
        if not hasattr(self, 'struc_man'):
            if input_structure_path is None:
                self.args.input_structure_path = input("Enter input structure path (PDB, mmcif): ")
            self.struc_man = StructureManager()
            self.struc_man.loadStructure(self.args.input_structure_path, 'force', False, self.args.debug)
            print ("Structure", self.args.input_structure_path, "loaded")
            print()

    def _get_structure(self):
        return self.struc_man.get_structure()

    def _save_structure(self):
        if self.args.output_structure_path is None:
            self.args.output_structure_path = input("Enter output structure path: ")
        self.struc_man.saveStructure(self.args.output_structure_path)

#===============================================================================

def _get_parameters (options, this_prog, this_param, this_dest, this_help, this_type=str):
    optionsParser = argparse.ArgumentParser(prog=this_prog)
    optionsParser.add_argument(this_param, dest=this_dest, help=this_help, type=this_type)
    return optionsParser.parse_args(options)

def _check_parameter (opts_param, input_text):
    while opts_param is None or opts_param == '':
        opts_param = input (input_text)
    if opts_param is str:
        opts_param = opts_param.replace(' ', '')
    return opts_param
