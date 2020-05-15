"""
 Module to read AMBER residue lib files
"""

import re
import sys

class ResidueLib():
    """ Class to manage a residue library taken from AMBER prep files.
    """
    def __init__(self, library_path):
        try:
            lib_file_h = open(library_path, "r")

        except IOError:
            print("ERROR: unable to open residue library" + library_path, file=sys.stderr)
            sys.exit(1)

        line = lib_file_h.readline()
        line = lib_file_h.readline()

        self.residues = {}
        at_group = False
        ch_group = False
        im_group = False
        res = ResidueDef()
        for line in lib_file_h:
            line = line.replace("\n", "").replace("\r", "")
            if line == '':
                continue

            elif line == "DONE":
                self.residues[res.id] = res
                at_group = False
                ch_group = False
                im_group = False
                res = ResidueDef()

            elif line == 'CHARGE':
                at_group = False
                ch_group = True

            elif line == 'IMPROPER':
                ch_group = False
                im_group = True

            elif re.match(' (...)  INT', line):
                resid_str = re.match(' (...)  INT', line)
                res.id = resid_str.group(1)
                at_group = True

            elif at_group:
                data = line.split()
                if len(data) < 11:
                    continue
                atm = AtomDef(data)
                for i in range(0, 3):
                    if atm.link[i] > 0:
                        atm.link_ats[i] = res.ats[atm.link[i]].id
                res.ats.append(atm)

            elif ch_group:
                for chrg in line.split():
                    res.charges.append(chrg)

            elif im_group:
                res.improper.append(line)

    def get_atom_def(self, res_id, at_id):
        """ Gets an atom definition given residue and atom names."""
        resid_def = self.residues[res_id]
        i = 1
        while resid_def.ats[i].id != at_id and i < len(resid_def.ats):
            i = i + 1
        if resid_def.ats[i].id == at_id:
            return resid_def.ats[i]
        return None

class ResidueDef():
    """ Class to pack a residue definition."""
    def __init__(self):
        self.id = ''
        self.name = ''
        self.ats = ['']
        self.improper = []
        self.charges = []

class AtomDef():
    """ Class to pack an atom definition."""
    def __init__(self, data):
        self.id = data[1]
        self.type = data[2]
        self.branch = data[3]
        self.link = [int(data[4]), int(data[5]), int(data[6])]
        self.link_ats = ['', '', '']
        self.geom = [float(data[7]), float(data[8]), float(data[9])]
        self.coord = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
        self.chrg = float(data[10])

    def __str__(self):
        return self.id
