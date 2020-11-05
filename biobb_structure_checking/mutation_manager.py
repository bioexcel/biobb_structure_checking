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
    def __init__(self, id_list):
        self.mutation_list = []
        print(id_list)

        if 'file:' in id_list:
            #Load from file
            id_list = id_list.replace('file:', '')
            print('Reading mutation list from file {}'.format(id_list))
            self.id_list = [
                line.replace('\n', '').replace('\r', '') for line in open(id_list, 'r')
            ]
        else:
            self.id_list = id_list.replace(' ', '').split(',')
        #convert to list ids to MutationSet instances and make the list
        self.mutation_list = list(map(MutationSet, self.id_list))

    def prepare_mutations(self, struc, stop_on_error=True):
        """  Check mutation_list and unroll chains/models."""
        for mut in self.mutation_list:
            mut.prepare(struc, stop_on_error)

    def apply_mutations(self, mutation_map, residue_lib, remove_h='mut'):
        """ Perform the requestet mutations."""
        mutated_res = []
        for mut in self.mutation_list:
            mutated_res += mut.apply(mutation_map, residue_lib, remove_h)
        return mutated_res

    def __str__(self):
        return ','.join(self.id_list)

#===============================================================================
class MutationSet():
    """ Class to manage instances of a mutation. """
    def __init__(self, mut_id):
        if ':' not in mut_id:
            mut_id = '*:' + mut_id

        self.id = mut_id.upper()

        self.chain, mut = mut_id.split(':')

        mut_comps = re.match('([A-z]*)([0-9]*)([A-z]*)', mut)

        self.old_id = mu.protein_residue_check(mut_comps.group(1))
        self.new_id = mu.protein_residue_check(mut_comps.group(3))
        self.res_num = mut_comps.group(2)

        self.id = ''.join([self.chain, ":", self.old_id, self.res_num, self.new_id])
        self.mutations = []

    def prepare(self, struc, stop_on_error=True):
        """ Check which mutations are possible. """
        mut_ok = 0
        for model in struc.get_models():
            if self.chain == '*':
                chains = []
                for chn in model.get_list():
                    chains.append(chn.get_id())
            else:
                chains = [self.chain]

            for chn in chains:
                if int(self.res_num) in model[chn]:
                    res = model[chn][int(self.res_num)]
                    if res.get_resname() == self.old_id:
                        self.mutations.append({
                            'model':model.get_id(),
                            'chain':chn,
                            'residue':res.get_id(),
                            'new_id':self.new_id,
                            'resobj': res
                        })
                        mut_ok += 1
                    else:
                        print(
                            '#WARNING: Unknown residue {}:{}{}'.format(
                                chn, self.old_id, self.res_num
                            )
                        )
                else:
                    print(
                        '#WARNING: Unknown residue {}:{}{}'.format(
                            chn, self.old_id, self.res_num
                        )
                    )

        if not mut_ok and stop_on_error:
            sys.exit(
                '#ERROR: no mutations available for {}'.format(self.id),
            )

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
            for at_id in ['N', 'CA', 'C']:
                if at_id not in bck_atoms:
                    sys.exit(
                        '#ERROR: Backbone atoms missing for {}, aborting'.format(
                            mu.residue_id(res)
                        )
                    )
            missing_ats = [
                at_id
                for at_id in mut_map[rname]['side_atoms']
                if at_id not in side_atoms
            ]
            print("Replacing " + mu.residue_id(res) + " into " + self.new_id)
            in_rules = []
            extra_adds = []
            # Renaming ats
            if MOV in mut_map[rname][self.new_id]:
                for rule in mut_map[rname][self.new_id][MOV]:
                    old_at, new_at = rule.split("-")
                    print('  Renaming {} to {}'.format(old_at, new_at))
                    if old_at in side_atoms:
                        mu.rename_atom(res, old_at, new_at)
                    else:
                        print(
                            '#WARNING: atom {} missing in {}'.format(
                                old_at, mu.residue_id(res)
                            )
                        )
                        extra_adds.append(new_at)
                    in_rules.append(old_at)
            # Deleting atoms
            if DEL in mut_map[rname][self.new_id]:
                for at_id in mut_map[rname][self.new_id][DEL]:
                    print('  Deleting {}'.format(at_id))
                    if at_id in side_atoms:
                        mu.delete_atom(res, at_id)
                    else:
                        print(
                            '#WARNING: atom {} already missing in {}'.format(
                                at_id, mu.residue_id(res)
                            )
                        )
                    in_rules.append(at_id)
            # Adding atoms (new_id required as r.resname is still the original)
            # Adding missing atoms that keep name
            for at_id in mut_map[rname]['side_atoms']:
                if at_id not in in_rules and at_id in missing_ats:
                    print('  Adding missing atom {}'.format(at_id))
                    mu.build_atom(res, at_id, res_lib, self.new_id)
            for at_id in extra_adds:
                print('  Adding new atom {}'.format(at_id))
                mu.build_atom(res, at_id, res_lib, self.new_id)
            if ADD in mut_map[rname][self.new_id]:
                for at_id in mut_map[rname][self.new_id][ADD]:
                    print('  Adding new atom {}'.format(at_id))
                    mu.build_atom(res, at_id, res_lib, self.new_id)
            #Renaming residue
            res.resname = self.new_id
            mutated_res.append(res)
        print("")
        return mutated_res

    def __str__(self):
        return self.id
