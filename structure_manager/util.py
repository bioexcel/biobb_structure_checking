"""
  Utility functions to manipulate structures
"""

import math
import numpy as np
from numpy import cos
from numpy import pi
from numpy import sin
from numpy.linalg import norm

# TODO: replace by Bio.PDB equivalent
one_letter_residue_code = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
    'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

three_letter_residue_code = {
    'A':'ALA', 'C': 'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
    'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'
}

# Residue manipulation =======================================================
def removeHFromRes (r, verbose=False):
    H_list = []
    for at in r.get_atoms():
        if at.element == 'H':
            H_list.append(at.id)
    for at_id in H_list:
        if verbose:
            print ("  Deleting atom " + at_id)
        r.detach_child(at_id)

def residueid(r, models=False):
    rid = r.get_resname() + " " \
        + str(r.get_parent().id) \
        + str(r.id[1])
    if models:
        rid += "/" + str(r.get_parent().get_parent().id)
    
    return rid

def atomid(at, models=False):
    return residueid(at.get_parent(), models) + "." + at.id

def residueCheck(r):
    r = r.upper()
    rid = ''
    if r in three_letter_residue_code.keys():
        rid = three_letter_residue_code[r]
    elif r in one_letter_residue_code.keys():
        rid = r
    else:
        print ('#ERROR: unknown residue id ' + r)
        sys.exit(1)

    return rid

# Atom building ===============================================================
def buildCoordsOther(r, res_lib, new_res, at_id):

    resid_def = res_lib.residues[new_res]
    i = 1
    while resid_def.ats[i].id != at_id and i < len(resid_def.ats):
        i = i + 1
    if resid_def.ats[i].id == at_id:
        return buildCoords(
                           r[resid_def.ats[resid_def.ats[i].link[0]].id].get_coord(),
                           r[resid_def.ats[resid_def.ats[i].link[1]].id].get_coord(),
                           r[resid_def.ats[resid_def.ats[i].link[2]].id].get_coord(),
                           resid_def.ats[i].geom
                           )
    else:
        print ("#ERROR: Unknown target atom")
        sys.exit(1)

def buildCoordsCB(r): # Get CB from Backbone

    return buildCoords(
                       r['CA'].get_coord(),
                       r['N'].get_coord(),
                       r['C'].get_coord(),
                       [1.5, 115.5, -123.]
                       )

def buildCoords(avec, bvec, cvec, geom):

    dst = geom[0]
    ang = geom[1] * pi / 180.
    tor = geom[2] * pi / 180.0

    v1 = avec-bvec
    v2 = avec-cvec

    n = np.cross(v1, v2)
    nn = np.cross(v1, n)

    n /= norm(n)
    nn /= norm(nn)

    n *= -sin(tor)
    nn *= cos(tor)

    v3 = n + nn
    v3 /= norm(v3)
    v3 *= dst * sin(ang)

    v1 /= norm(v1)
    v1 *= dst * cos(ang)

    return avec + v3 - v1

# Structure utils =============================================================
def calcRMSdAll (st1, st2):
    ats1 = []
    ats2 = []

    for at in st1.get_atoms():
        ats1.append(at)
    for at in st2.get_atoms():
        ats2.append(at)

    rmsd = 0

    i = 0
    while i < len(ats1)and i < len(ats2):
        d = ats1[i]-ats2[i]
        rmsd = rmsd + d * d / len(ats1)
        i = i + 1

    return (math.sqrt(rmsd))

def get_all_rr_distances(r1, r2, with_h=False):
    dist_mat = []
    for at1 in r1.get_atoms():
        if at1.element == 'H' and not with_h:
            continue
        for at2 in r2.get_atoms():
            if at2.element == 'H' and not with_h:
                continue
            if at1.serial_number < at2.serial_number:
                dist_mat.append ([at1, at2, at1-at2])
    return dist_mat

def same_residue (at1, at2):
    return at1.get_parent() == at2.get_parent()

def same_model(r1, r2):
    return r1.get_parent().get_parent() == r2.get_parent().get_parent()

def same_chain(r1, r2):
    return r1.get_parent() == r2.get_parent() and same_model(r1, r2)

def seq_consecutive(r1, r2):
    resnum1 = r1.id[1]
    resnum2 = r2.id[1]
    return same_chain(r1, r2) and abs(resnum1-resnum2) == 1


