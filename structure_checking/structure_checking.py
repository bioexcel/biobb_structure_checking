# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import sys
import argparse
import re

import settings as sets
from structure_checking.structure_manager import StructureManager
from structure_checking.help_manager import HelpManager
import structure_checking.util as util



class StructureChecking():
    def __init__(self, args):
        self.args=args
 
    def launch(self):
        help = HelpManager(sets.help_dir_path)
        help.print_help('header')
        
        try:
            f=getattr(self, self.args.command)
        except AttributeError:
            print ("Error: command unknown or not implemented")
            sys.exit(1)

        f(self.args.options)
        
        self._save_structure()
        print ("Structure saved on", self.args.output_structure_path)

        
    def command_list(self,options):
        optionsParser = argparse.ArgumentParser(prog="command_list")
        optionsParser.add_argument('--list', dest="op_list", )
        opts = optionsParser.parse_args(options)

        opts.op_list = _check_parameter(opts.op_list, "Command list file: ")
        
        
        try:
            fh = open(opts.op_list, "r")
        
        except OSError:
            print ("Error opening",opts.op_list)
            sys.exit(1)
        
        print ("Running command_list from", opts.op_list)

        self._load_structure(self.args.input_structure_path)
        
        i=1
        for line in fh:
            print ("Step "+ str(i) + ":",line)
            data = line.split()
            command = data[0]
            options = data[1:]
            try:
                f=getattr(self, command)
            except AttributeError:
                print ("Error: command unknown or not implemented")
                sys.exit(1)
            f(options)
            i +=1
            
        print ("Command list completed")
            
    
    def models(self, options):
        optionsParser = argparse.ArgumentParser(prog="models")
        optionsParser.add_argument('--select_model', dest="select_model", type=int)
        opts = optionsParser.parse_args(options)
        
        print ("Running models", ' '.join(options))

        self._load_structure(self.args.input_structure_path)

        self.nmodels = len(self.struc_man.get_structure())
        
        print ("Models detected: ", self.nmodels)
        
        if self.nmodels > 1:
            opts.select_model = _check_parameter(opts.select_model,"Select Model Num [1-"+str(self.nmodels)+"]: ")
            opts.select_model = int(opts.select_model)
            if opts.select_model < 1 or opts.select_model > self.nmodels:
                print ("Error: unknown model", opts.select_model, file=sys.stderr)
                opts.select_model=''
            else:
                print ("Selecting model", opts.select_model)
                self.struc_man.st = self.struc_man.st[int(opts.select_model)-1]
                
        else:
            print ("Nothing to do")
        
        print()
                
    def chains(self, options):
        optionsParser = argparse.ArgumentParser(prog="chains")
        optionsParser.add_argument('--select_chains', dest="select_chains", help="Chains (*| list comma separated)", )
        opts = optionsParser.parse_args(options)

        print ("Running chains", ' '.join(options))
        
        self._load_structure(self.args.input_structure_path)

        self.chain_ids=[]
        
        for ch in self._get_structure().get_chains():
            self.chain_ids.append(ch.id)

        print ("Chains detected: ", len(self.chain_ids), '(', ', '.join(self.chain_ids),')')
        
        if len(self.chain_ids)>1:
            opts.select_chains = _check_parameter(opts.select_chains,"Select chain id(s) or * for all [*," +",".join(self.chain_ids)+"]: ")
            if opts.select_chains != '*':
                ch_ok = opts.select_chains.split(',')
                for ch in ch_ok:
                    if not ch in self.chain_ids:
                        print ("Error: request chain not present", ch, file=sys.stderr)
                        opts.select_chains=''
                print ("Selecting chain(s) ", ','.join(ch_ok))
                for ch in self.chain_ids:
                    if ch not in ch_ok:
                        self._get_structure()[0].detach_child(ch)
            else:
                print ("Selecting all chains")
        else:
            print ("Nothing to do")
        print ()
        
    
    def altloc(self, options):
        optionsParser = argparse.ArgumentParser(prog="altloc")
        optionsParser.add_argument('--select_altloc', dest="select_altloc", help="select altloc occupancy|alt_id")
        opts = optionsParser.parse_args(options)
        
        print ("Running altlocs", ' '.join(options))
        
        self._load_structure(self.args.input_structure_path)
    
        alt_loc_res=self.struc_man.get_altloc_residues()
        
        if len(alt_loc_res) > 0:
            print ("Detected alternative locations")
            for r in alt_loc_res.keys():
                if len(alt_loc_res[r])>1:
                    print (r, ":")
                    altlocs=sorted(alt_loc_res[r][0].child_dict.keys())
                    for at in alt_loc_res[r]:
                        s = "  " + at.id
                        for alt in sorted(at.child_dict.keys()):
                            s += " "+ alt + "("+ str(at.child_dict[alt].occupancy) +")"
                        print (s)
                    ok=False
                    while not ok:
                        opts.select_altloc=_check_parameter('', "Select alternative (occupancy, "+ ','.join(altlocs) +'): ')
                        ok = opts.select_altloc in altlocs or opts.select_altloc=='occupancy'
                        if not ok:
                            print ("Error: Unknown ",opts.select_altloc, file=sys.stderr)
                            opts.select_altloc=''
                    
                        print ("Selecting location", opts.select_altloc)
                    for at in alt_loc_res[r]:
                        if opts.select_altloc == 'occupancy':
                            newat = at.selected_child
                        else:
                            newat = at.child_dict[opts.select_altloc]
                        newat.disordered_flag=0
                        newat.altloc=' '
                        res = at.get_parent()
                        res.detach_child(at.id)
                        res.add (newat)
                    res.disordered=0
        else:
            print ("No alternative location detected")
    
    def metals (self,options):
        optionsParser = argparse.ArgumentParser(prog="metals")
        optionsParser.add_argument(
            '--remove', 
            action='store',
            dest='remove_metals', 
            help='Remove Metal ions'
        )
        opts = optionsParser.parse_args(options)
        
        print ("Running metals", ' '.join(options))
        
        self._load_structure(self.args.input_structure_path)
        
        met_list=[]
        for at in self.struc_man.get_structure().get_atoms():
            if not re.match('H_', at.get_parent().id[0]):
                continue
            if at.id in sets.metal_ats:
                met_list.append(at)
        if len(met_list) > 1:
            print ("Metal ions found")
            met_rid=[]
            at_groups={}
            for at in met_list:
                print ("  ", util.atomid(at))
                r = at.get_parent()
                rid = r.get_parent().id + str(r.id[1])
                met_rid.append(rid)
                if at.id not in at_groups.keys():
                    at_groups[at.id]=[]
                at_groups[at.id].append(rid)
            print (at_groups)
            
            ok = False
            resids=False
            atids=False
            while not ok:
                opts.remove_metals = _check_parameter(opts.remove_metals,
                'Remove (All | None | any of  ' + ','.join(sorted(at_groups.keys())) + ' | any of ' + ','.join(met_rid)+': ')
                ok = opts.remove_metals in ['All', 'None']
                if not ok:
                    ok= True
                    for r in opts.remove_metals.split(','):
                        ok = ok and r in met_rid
                    resids=ok
                
                if not ok:
                    ok= True
                    for atid in opts.remove_metals.split(','):
                        ok = ok and atid in at_groups.keys()
                    atids=ok
            
                if not ok:
                    print ("Error unknown", opts.remove_metals, file=sys.stderr)
                    opts.remove_metals=''
                
                if opts.remove_metals == 'None':
                    continue
                
                elif opts.remove_metals == 'All':
                    to_remove = met_rid

                elif resids:
                    to_remove = opts.remove_metals.split(',')
                
                elif atids:
                    to_remove=[]
                    for atid in opts.remove_metals.split(','):
                        to_remove.extend(at_groups[atid])

                print (to_remove)
                
                for at in met_list:
                    r = at.get_parent()
                    rid = r.get_parent().id + str(r.id[1])
                    if rid in to_remove:
                        r.detach_child(at.id)
                print ("Metals removed", opts.remove_metals)
        else:
            print ("No metal ions found")
        
        
#===============================================================================
    
    def _load_structure(self,input_structure_path):
        if not hasattr(self,'struc_man'):
            if input_structure_path is None:
                self.args.input_structure_path = input("Enter input structure path (PDB, mmcif): ")
            self.struc_man=StructureManager()
            self.struc_man.loadStructure(self.args.input_structure_path, 'force', False, self.args.debug)
            print ("Structure",self.args.input_structure_path,"loaded")
            print()
        
    def _get_structure(self):
        return self.struc_man.get_structure()
    
    def _save_structure(self):
        if self.args.output_structure_path is None:
            self.args.output_structure_path = input("Enter output structure path: ")
        self.struc_man.saveStructure(self.args.output_structure_path)
    
def _check_parameter (opts_param, input_text):
    while opts_param is None or opts_param=='':
        opts_param = input (input_text)
    opts_param = opts_param.replace(' ','')
    return opts_param
