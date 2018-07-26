"""
 Module to read AMBER residue lib
"""

import re
import sys

class ResidueLib():

    def __init__(self, library_path):

        try:
            fh = open(library_path, "r")

        except IOError:
            print ("ERROR: unable to open residue library" + library_path, file=sys.stderr)
            sys.exit(1)

        line = fh.readline()
        line = fh.readline()

        self.residues = {}
        at_group = False
        ch_group = False
        im_group = False
        res = ResidueDef()
        for line in fh:
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
                at = AtomDef(data)
                res.ats.append(at)

            elif ch_group:
                for c in line.split():
                    res.charges.append(c)

            elif im_group:
                res.improper.append(line)

class ResidueDef():

    def __init__(self):
        self.id = ''
        self.name = ''
        self.ats = ['']
        self.improper = []
        self.charges = []

class AtomDef():

    def __init__(self, data):
        self.id = data[1]
        self.type = data[2]
        self.branch = data[3]
        self.link = [int(data[4]), int(data[5]), int(data[6])]
        self.geom = [float(data[7]), float(data[8]), float(data[9])]
        self.coord = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
        self.ch = float(data[10])

    def __str__(self):
        return self.id
