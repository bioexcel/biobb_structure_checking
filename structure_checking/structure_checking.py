# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "gelpi"
__date__ = "$26-jul-2018 14:34:51$"

import argparse
import re
import sys

from structure_checking.help_manager import HelpManager
from structure_manager.structure_manager import StructureManager
from structure_checking.json_writer import JSONWriter
import structure_manager.util as util
import structure_manager.data_constants as dataCts

class StructureChecking():
    def __init__(self, args):
        self.args = args
        self.summary={}

    def launch(self, sets):
        help = HelpManager(sets.help_dir_path)
        help.print_help('header')

        try:
            f = getattr(self, self.args.command)
        except AttributeError:
            print ("Error: command unknown or not implemented")
            sys.exit(1)

        f(self.args.options)

        if not self.args.check_only: 
            self._save_structure()
            print ("Structure saved on", self.args.output_structure_path)
        
        if self.args.json_output_path is not None:
            json_writer=JSONWriter()
            for k in self.summary.keys():
                json_writer.set(k,self.summary[k])
            json_writer.save(self.args.json_output_path)
        
    def command_list(self, options):

        opts = _get_parameters(options, "command_list", "--list", "op_list", "Command List File")

        opts.op_list = _check_parameter(opts.op_list, "Command list file: ")

        try:
            fh = open(opts.op_list, "r")

        except OSError:
            print ("Error opening", opts.op_list)
            sys.exit(1)

        print ("Running command_list from", opts.op_list)

        self._load_structure()

        i = 1
        for line in fh:
            if line == "\n" or line[0:1] == '#':
                continue
            print ("Step " + str(i) + ":", line)
            data = line.split()
            command = data[0]
            options = data[1:]
            try:
                f = getattr(self, command)
            except AttributeError:
                print ("Error: command unknown or not implemented")
                continue
                #sys.exit(1)
            f(options)
            i += 1

        print ("Command list completed")


    def models(self, options):

        opts = _get_parameters(options, "models", "--select_model", "select_model", "Select model to keep", int)
        print ("Running models", ' '.join(options))
        models_sum = {}

        self._load_structure()
        self.nmodels = self.struc_man.get_nmodels()

        print ("Models detected: ", self.nmodels)
        models_sum['detected'] = self.nmodels
        
        if not self.args.check_only:
            if self.nmodels > 1:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.select_model = _check_parameter(opts.select_model, "Select Model Num [1-" + str(self.nmodels) + "]: ")
                    opts.select_model = int(opts.select_model)
                    print (opts.select_model)
                    ok = opts.select_model > 0 or opts.select_model <= self.nmodels
                    if not ok:
                        print ("Error: unknown model", opts.select_model, file=sys.stderr)
                        opts.select_model = ''
                        if self.args.non_interactive:
                            self.summary['models'] = models_sum
                            return 1
        
                print ("Selecting model", opts.select_model)
                self.struc_man.select_model(opts.select_model)
                self.nmodels = self.struc_man.get_nmodels()
                models_sum['selected'] = opts.select_model
        
            else:
                print ("Nothing to do")

        self.summary['models'] = models_sum

        print()
        
    def chains(self, options):

        opts = _get_parameters(options, "chains", '--select_chains', "select_chains", "Chains (*| list comma separated)")
        print ("Running chains", ' '.join(options))
        chains_sum = {}
        
        self._load_structure()
        self.chain_ids = self.struc_man.get_chain_ids()
        print ("Chains detected: ", len(self.chain_ids), '(', ', '.join(self.chain_ids), ')')
        
        chains_sum['detected'] = self.chain_ids

        if not self.args.check_only:
            if len(self.chain_ids) > 1:
                ok= False
                while not ok:
                    if not self.args.non_interactive:
                        opts.select_chains = _check_parameter(opts.select_chains, "Select chain id(s) orAll," + ",".join(self.chain_ids) + "]: ")
                    ok = opts.select_chains == 'All'
                    if not ok:
                        ok = True
                        for ch in opts.select_chains.split(','):
                            ok = ok and ch in self.chain_ids
                    if not ok:
                        print ("Error unknown", opts.select_chains)
                        if self.args.non_interactive:
                            self.summary['chains'] = chains_sum
                            return 1
                if opts.select_chains != 'All':
                    self.struc_man.select_chains(opts.select_chains)
                    print ("Selecting chain(s) ", ','.join(opts.select_chains))
                    chains_sum['selected']=opts.select_chains.split(',')
                else:
                    print ("Selecting all chains")
                    chains_sum['selected']=self.chain_ids
                self.chain_ids = self.struc_man.get_chain_ids()
                chains_sum['final']=self.chain_ids            
            else:
                print ("Nothing to do")

        print ()
        
        self.summary['chains'] = chains_sum

    def altloc(self, options):

        opts = _get_parameters(options, "altloc", '--select_altloc', "select_altloc", "select altloc occupancy|alt_id")

        print ("Running altlocs", ' '.join(options))
        altloc_sum={}
        
        self._load_structure()

        alt_loc_res = self.struc_man.get_altloc_residues()
        if len(alt_loc_res) > 0:
            print ("Detected alternative locations")

            altloc_sum['detected']={}

            for r in alt_loc_res.keys():
                if len(alt_loc_res[r]) > 1:
                    print (r, ":")
                    altloc_sum['detected'][r]=[]
                    altlocs = sorted(alt_loc_res[r][0].child_dict.keys())
                    for at in alt_loc_res[r]:
                        s = "  " + at.id
                        for alt in sorted(at.child_dict.keys()):
                            s += " " + alt + "(" + str(at.child_dict[alt].occupancy) + ")"
                        print (s)
                        altloc_sum['detected'][r].append({'atom':at.id, 'loc_label':alt,'occupancy':at.child_dict[alt].occupancy})

                    if not self.args.check_only:
                        ok = opts.select_altloc in altlocs or opts.select_altloc == 'occupancy'
                        while not ok:
                            if not self.args.non_interactive:
                                opts.select_altloc = _check_parameter(opts.select_altloc, "Select alternative (occupancy, " + ','.join(altlocs) + '): ')
                            ok = opts.select_altloc in altlocs or opts.select_altloc == 'occupancy'
                            if not ok:
                                print ("Error: Unknown ", opts.select_altloc, file=sys.stderr)
                                opts.select_altloc = ''
                                if self.args.non_interactive:
                                    self.summary['altloc']=altloc_sum
                                    return 1

                        print ("Selecting location", opts.select_altloc)
                        altloc_sum['selected']=opts.select_altloc
                        self.struc_man.select_altloc_residues(r, opts.select_altloc)
        else:
            print ("No alternative locations detected")
        
        self.summary['altloc']=altloc_sum
    
    
    def metals (self, options):

        opts = _get_parameters(options, "metals", '--remove', 'remove_metals', 'Remove Metal ions')

        print ("Running metals", ' '.join(options))
        metals_sum={}
        
        self._load_structure()

        met_list = self.struc_man.get_metals(dataCts.metal_ats)

        if len(met_list) > 1:
            print ("Metal ions found")
            metals_sum['detected']=[]
            met_rids=[]
            at_groups={}
            for at in met_list:
                print ("  ", util.atomid(at))
                r=at.get_parent()
                met_rids.append(r.get_parent().id + str(r.id[1]))
                if not at.id in at_groups:
                    at_groups[at.id]=[]
                at_groups[at.id].append(at)
                metals_sum['detected'].append(r.get_parent().id + str(r.id[1]))

            if not self.args.check_only:
                ok = False
                resids = False
                atids = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_metals = _check_parameter(opts.remove_metals,
                             'Remove (All | None | any of  ' 
                             + ','.join(sorted(at_groups.keys())) 
                             + ' | any of ' + ','.join(met_rids) + ': ')
                    ok = opts.remove_metals in ['All', 'None']
                    if not ok:
                        ok = True
                        for r in opts.remove_metals.split(','):
                            ok = ok and r in met_rids
                        resids = ok

                    if not ok:
                        ok = True
                        for atid in opts.remove_metals.split(','):
                            ok = ok and atid in at_groups.keys()
                        atids = ok

                    if not ok:
                        print ("Error unknown", opts.remove_metals, file=sys.stderr)
                        opts.remove_metals = ''
                        if self.args.non_interactive:
                            self.summary['metals']=metals_sum
                            return 1

                if opts.remove_metals == 'None':
                    to_remove = []

                elif opts.remove_metals == 'All':
                    to_remove = met_list

                elif resids:
                    to_remove = []
                    rid_list = opts.remove_metals.split(',')
                    for at in met_list:
                        r=at.get_parent()
                        if r.get_parent().id + str(r.id[1]) in rid_list:
                            to_remove.append(at)
                elif atids:
                    to_remove = []
                    for atid in opts.remove_metals.split(','):
                        to_remove.extend(at_groups[atid])
                metals_sum['removed']=[]
                n=0
                for at in to_remove:
                    metals_sum['removed'].append(util.residueid(at.get_parent()))
                    self.struc_man.remove_residue(at.get_parent())
                    n+= 1
        
                print ("Metal Atoms removed", opts.remove_metals,"("+str(n)+")")
                metals_sum['n_removed'] = n
        
        else:
            print ("No metal ions found")
        
        self.summary['metals']=metals_sum
        
    def remwat(self, options):
        
        opts = _get_parameters(options, "remwat", '--remove', 'remove_wat', 'Remove Water molecules')

        print ("Running remwat", ' '.join(options))
        remwat_sum={}
        
        self._load_structure()
        
        lig_list = self.struc_man.get_ligands(incl_water=True)
        
        wat_list=[]
        for r in lig_list:
            if r.id[0] == 'W':
                wat_list.append(r)
        if len(wat_list) > 0:
            print ("Water molecules detected")
            remwat_sum['n_detected']= len(wat_list)
            
            if not self.args.check_only:
                ok = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_wat = _check_parameter(opts.remove_wat,'Remove (Yes | No): ')
                    ok = opts.remove_wat in ['Yes', 'No']                
                    if not ok:
                        print ("Warning: unknown", opts.remove_wat)
                        if self.args.non_interactive:
                            self.summary['remwat']=remwat_sum
                            return 1                        
            
                if opts.remove_wat == 'Yes':
                    n=0
                    for r in wat_list:
                        self.struc_man.remove_residue(r)
                        n+= 1
                    print ("Water molecules removed", '('+str(n)+')')
                    remwat_sum['n_removed'] = n
        else:
            print ("No water molecules detected")
        
        self.summary['remwat']=remwat_sum
        print ()

    def ligands(self, options):

        opts = _get_parameters(options, "ligands", '--remove', 'remove_ligands', 'Remove Ligand residues')

        print ("Running ligands", ' '.join(options))
        ligands_sum={}
        
        self._load_structure()
        
        lig_list = self.struc_man.get_ligands()
        if len(lig_list) > 0:
            print ("Ligands detected")
            rids = set()
            rnums = []
            ligands_sum['detected']=[]
            for r in lig_list:
                print (util.residueid(r))
                ligands_sum['detected'].append(util.residueid(r))
                rids.add(r.resname)
                rnums.append(r.get_parent().id + str(r.id[1]))
            
            if not self.args.check_only:
                ok = False
                byresnum = False
                byrids = False
                while not ok:
                    if not self.args.non_interactive:
                        opts.remove_ligands = _check_parameter(opts.remove_ligands,
                            'Remove (All | None | any of  ' 
                            + ','.join(sorted(rids)) 
                            + ' | any of ' + ','.join(rnums) + '): ')
                    ok = opts.remove_ligands in ['All', 'None']
                    if not ok:
                        ok = True
                        for rid in opts.remove_ligands.split(','):
                            ok = ok and (rid in rids)
                        byrids = ok
                    if not ok:
                        ok = True
                        for rn in opts.remove_ligands.split(','):
                            ok = ok and (rn in rnums)
                        byresnum = ok
                    if not ok:
                        print ("Error: unknown", opts.remove_ligands)
                        if self.args.non_interactive:
                            self.summary['ligands']=ligands_sum
                            return 1
                
                    ligands_sum['removed']={'opt':opts.remove_ligands,'lst':[]}
            
                to_remove=[]
            
                if opts.remove_ligands == 'None':
                    print ("Nothing to do")
                elif opts.remove_ligands == 'All':
                    to_remove = lig_list
                elif byrids:
                    rm = opts.remove_ligands.split(',')
                    for r in lig_list:
                        if r.resname in rm:
                            to_remove.append(r)
                elif byresnum:
                    rm = opts.remove_ligands.split(',')
                    for r in lig_list:
                        if (r.get_parent().id + str(r.id[1])) in rm:
                            to_remove.append(r)
                n=0
                for r in to_remove:
                    ligands_sum['removed']['lst'].append(util.residueid(r))
                    self.struc_man.remove_residue(r)
                    n+= 1
                    
                print ("Ligands removed", opts.remove_ligands, '('+str(n)+')')
                ligands_sum['n_removed']=n
        else:
            print ("No ligands detected")

        self.summary['ligands']=ligands_sum
        print()
    
    def getss (self, options):
        print ("Running getss")
        getss_sum={}
        
        self._load_structure()
        
        SS_bonds = self.struc_man.get_all_at2at_distances('SG', dataCts.SS_DIST)
        
        if len(SS_bonds):            
            print ("Possible SS Bonds detected")
            getss_sum['detected']=[]
            for ssb in SS_bonds:
                print ("  ", util.atomid(ssb[0]), util.atomid(ssb[1]), ssb[2])
                getss_sum['detected'].append({'at1':util.atomid(ssb[0]), 'at2':util.atomid(ssb[1]), 'dist': float(ssb[2])})
        
        else:
            print ("No SS bonds detected")

        self.summary['getss']=getss_sum
        print ()
        
    def clashes(self, options):
        opts = _get_parameters(options, "clashes", '--no_wat', 'discard_wat', 'Discard water molecules') 
        if opts.discard_wat is None: 
            opts.discard_wat=True
        print ("Running clashes")
        clashes_sum={}
        self._load_structure()
        
        rr_dist = self.struc_man.get_all_r2r_distances('All', dataCts.R_R_CUTOFF)
        
        clashes={
            'severe':{},
            'apolar':{},
            'acceptor':{},
            'donor':{},
            'positive':{},
            'negative':{}
        }
        
        clashes_sum = {
            'detected': {
                'severe':[],
                'apolar':[],
                'acceptor':[],
                'donor':[],
                'positive':[],
                'negative':[]
            }
        }
        
        for r_pair in rr_dist:
            [r1,r2,d] = r_pair
            if opts.discard_wat and (r1.id[0] == 'W' or r2.id[0] == 'W'):
                continue
            if r1 != r2 and not util.seq_consecutive(r1,r2) and util.same_model(r1,r2):
                rkey = util.residueid(r1)+'-'+util.residueid(r2)
                for at_pair in util.get_all_rr_distances(r1,r2):
                    [at1,at2,dist] = at_pair
                    for cls in ['severe', 'apolar', 'acceptor', 'donor', 'positive', 'negative']:
                        if dist < dataCts.CLASH_DIST[cls]:
                            if cls == 'apolar' and (at1.element not in dataCts.apolar_elements and at2.element not in dataCts.apolar_elements):
                                continue
                            if cls == 'donor' and not (at1.id in dataCts.polar_donor and at2.id in dataCts.polar_donor):
                                continue
                            if cls == 'acceptor' and not (at1.id in dataCts.polar_acceptor and at2.id in dataCts.polar_acceptor):
                                continue
                            if cls == 'positive' and not (at1.id in dataCts.pos_ats and at2.id in dataCts.pos_ats):
                                continue
                            if cls == 'negative' and not (at1.id in dataCts.neg_ats and at2.id in dataCts.neg_ats):
                                continue
                            if not rkey in clashes[cls].keys():
                                clashes[cls][rkey]=at_pair
                            if dist < clashes[cls][rkey][2]:
                                clashes[cls][rkey]=at_pair
        for cls in ['severe', 'apolar', 'acceptor', 'donor', 'positive', 'negative']:
            if len(clashes[cls]):
                print ("Steric", cls, "clashes detected")
                for rkey in clashes[cls].keys():
                    print (" ",util.atomid(clashes[cls][rkey][0]), util.atomid(clashes[cls][rkey][1]), clashes[cls][rkey][2])
                    at_pair = clashes[cls][rkey]
                    clashes_sum['detected'][cls].append({'at1':util.atomid(clashes[cls][rkey][0]), 'at2':util.atomid(clashes[cls][rkey][1]), 'dist': float(clashes[cls][rkey][2])})
            else:
                print ("No", cls, "clashes detected")
        
        self.summary['clashes'] = clashes_sum
        print ()
#===============================================================================

    def _load_structure(self):
        if not hasattr(self, 'struc_man'):
            if not self.args.non_interactive and self.args.input_structure_path is None:
                self.args.input_structure_path = input("Enter input structure path (PDB, mmcif): ")
            self.struc_man = StructureManager()
            self.struc_man.loadStructure(self.args.input_structure_path, 'force', False, self.args.debug)
            print ("Structure", self.args.input_structure_path, "loaded")

    def _get_structure(self):
        return self.struc_man.get_structure()

    def _save_structure(self):
        if not self.args.non_interactive and self.args.output_structure_path is None:
            self.args.output_structure_path = input("Enter output structure path: ")
        self.struc_man.saveStructure(self.args.output_structure_path)

    def json(self):
        json_writer = JSONWriter()
        
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
