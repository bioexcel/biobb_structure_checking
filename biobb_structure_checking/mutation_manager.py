"""
  Module to manage mutations
"""

import re
import sys

import biobb_structure_checking.model_utils as mu

ADD = 'Add'
DEL = 'Del'
MOV = 'Mov'

class MutationManager():
    """ Class to manage list of mutations. """
    def __init__(self, id_list, chain_ids):
        self.mutation_list = []
        self.chain_ids = chain_ids
        if 'file:' in id_list:
            #Load from file
            id_list = id_list.replace('file:', '')
            print(f'Reading mutation list from file {id_list}')
            self.id_list = [
                line.replace('\n', '').replace('\r', '') for line in open(id_list, 'r')
            ]
        else:
            self.id_list = id_list.replace(' ', '').split(',')
        #convert to list ids to MutationSet instances and make the list
        self.mutation_list = [MutationSet(mut_id, self.chain_ids) for mut_id in self.id_list]

    def prepare_mutations(self, struc, stop_on_error=True):
        """  Check mutation_list and unroll chains/models."""
        for mut in self.mutation_list:
            mut.prepare(struc, stop_on_error) 
            
    def apply_mutations(self, mutation_map, residue_lib, remove_h='mut'):
        """ Perform the requested mutations."""
        mutated_res = []
        for mut in self.mutation_list:
            mutated_res += mut.apply(mutation_map, residue_lib, remove_h)
        return mutated_res

    def __str__(self):
        return ','.join(self.id_list)

#===============================================================================
class MutationSet():
    """ Class to manage instances of a mutation. """
    def __init__(self, mut_id, chain_ids):
        if ':' not in mut_id:
            mut_id = '*:' + mut_id

        self.id = mut_id.upper()
        self.chain_ids = chain_ids
        self.chain, mut = mut_id.split(':')

        mut_comps = re.match('([A-z]*)([0-9]*)([A-z]*)', mut)
        self.old_id = mut_comps.group(1).upper()
        self.new_id = mut_comps.group(3).upper()
        self.res_num = mut_comps.group(2)

        self.id = ''.join([self.chain, ":", self.old_id, self.res_num, self.new_id])
        self.mutations = []

    def prepare(self, struc, stop_on_error=True):
        """ Check which mutations are possible. """
        mut_ok = 0
        for model in struc.get_models():
            if self.chain == '*':
                chains = self.chain_ids.keys()
            else:
                chains = [self.chain]

            for chn in chains:
                old_id = mu.valid_residue_check(self.old_id, self.chain_ids[chn])
                new_id = mu.valid_residue_check(self.new_id, self.chain_ids[chn])
                if int(self.res_num) in model[chn]:
                    res = model[chn][int(self.res_num)]
                    if res.get_resname() == old_id:
                        self.mutations.append({
                            'model':model.get_id(),
                            'chain':chn,
                            'type': self.chain_ids[chn],                    
                            'residue':old_id,
                            'new_id':new_id,
                            'resobj': res
                        })
                        mut_ok += 1
                    else:
                        print(f'#WARNING: Unknown residue {chn}:{old_id}{self.res_num}')
                else:
                    print(f'#WARNING: Unknown residue {chn}:{old_id}{self.res_num}')

        if not mut_ok and stop_on_error:
            sys.exit(f'#ERROR: no mutations available for {self.id}')

    def apply(self, mut_map, res_lib, remove_h):
        """ Perform the individual mutations on the set. """
        mutated_res = []
        for mut in self.mutations:
            res = mut['resobj']
            #struc[mut['model']][mut['chain']][mut['residue']]
            rname = res.get_resname().replace(' ', '')
            # Deleting H
            if remove_h == 'mut':
                mu.remove_H_from_r(res, verbose=True)
            # checking side chain
            side_atoms = []
            bck_atoms = []
            for atm in res.get_atoms():
                if atm.id in mut_map[rname]['side_atoms']:
                    side_atoms.append(atm.id)
                else:
                    bck_atoms.append(atm.id)

            if mut['type'] == mu.PROTEIN:
                for at_id in ['N', 'CA', 'C']:
                    if at_id not in bck_atoms:
                        sys.exit(f'#ERROR: Backbone atoms missing for {mu.residue_id(res)}, aborting')
            else:
                for at_id in ["C1'", "O4'", "C4'"]:
                    if at_id not in bck_atoms:
                        sys.exit(f'#ERROR: Backbone atoms missing for {mu.residue_id(res)}, aborting')

            missing_ats = [
                at_id
                for at_id in mut_map[rname]['side_atoms']
                if at_id not in side_atoms
            ]
            print(f"Replacing {mu.residue_id(res)} into {mut['new_id']}")
            in_rules = []
            extra_adds = []
             # Deleting atoms
            if DEL in mut_map[rname][mut['new_id']]:
                for at_id in mut_map[rname][mut['new_id']][DEL]:
                    print(f'  Deleting {at_id}')
                    if at_id in side_atoms:
                        mu.delete_atom(res, at_id)
                    else:
                        print(f'#WARNING: atom {at_id} already missing in {mu.residue_id(res)}')
                    in_rules.append(at_id)
            # Renaming ats
            if MOV in mut_map[rname][mut['new_id']]:
                for rule in mut_map[rname][mut['new_id']][MOV]:
                    old_at, new_at = rule.split("-")
                    print(f'  Renaming {old_at} to {new_at}')
                    if old_at in side_atoms:
                        mu.rename_atom(res, old_at, new_at)
                    else:
                        print(f'#WARNING: atom {old_at} missing in {mu.residue_id(res)}')
                        extra_adds.append(new_at)
                    in_rules.append(old_at)

            # Adding atoms (new_id required as r.resname is still the original)
            # Adding missing atoms that keep name
            for at_id in mut_map[rname]['side_atoms']:
                if at_id not in in_rules and at_id in missing_ats:
                    print(f'  Adding missing atom {at_id}')
                    mu.build_atom(res, at_id, res_lib, mut['new_id'])
            for at_id in extra_adds:
                print(f'  Adding new atom {at_id}')
                mu.build_atom(res, at_id, res_lib, mut['new_id'])
            if ADD in mut_map[rname][mut['new_id']]:
                for at_id in mut_map[rname][mut['new_id']][ADD]:
                    print(f'  Adding new atom {at_id}')
                    mu.build_atom(res, at_id, res_lib, mut['new_id'])
            #Renaming residue
            res.resname = mut['new_id']
            mutated_res.append(res)
        print("")
        return mutated_res

    def __str__(self):
        return self.id
