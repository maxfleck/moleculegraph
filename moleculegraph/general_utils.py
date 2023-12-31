import numpy as np
import glob
import csv
from rapidfuzz.distance.Levenshtein import normalized_distance
import json
from .molecule_utils import sort_force_fields, make_graph

"""
loose collection of functions in connection with the molecule class.
Mainly heading towards generate stuff for mapping for our inhouse projects.
Might be useless for external users.

subjects:
- read
- write
- angles, distances...
- forces
"""

def unite_atoms( elements, atoms, masses, bond_list, UA=True ):
    """
    generates united atom coordinates from a list of atoms.

    TO DO:
    - the atoms and masses objects are not obvious...

    Parameters:
    - elements:
        - list of atoms elements.
    - atoms:
        - list of atoms coordinates.
    - masses:
        - masses of the atoms.
    - bond_matrix:
        - bond matrix of the corresponding molecule
        - ...distance matrix works too
        
    Returns:
    - new_matrix:
        - list of united atoms with coordinates.
    - new_masses:
        - list of corresponding united atom masses.
    - new_bonds:
        - list of corresponding united atom bonds.
        - you can generate a new bond or distance matrix with this list.
    """
    keep = []
    keep_bonds = []
    for i, p in enumerate(bond_list):
        p0 = int(np.min(p))
        p1 = int(np.max(p))
        name0 = elements[p0][0]
        name1 = elements[p1][0]
        if (name0 == "H" and name1 == "C" and UA) or (
            name1 == "H" and name0 == "C" and UA
        ):
            if name0 == "H" and name1 == "C":
                p0 = int(np.max(p))
                p1 = int(np.min(p))
                name0 = elements[p0][0]
                name1 = elements[p1][0]
            coos_a = np.array(atoms[p0]).astype(float)
            coos_b = np.array(atoms[p1]).astype(float)
            mass_a = masses[p0]
            mass_b = masses[p1]
            dummy = (coos_a * mass_a + coos_b * mass_b) / (mass_a + mass_b)
            atoms[p0][0] = dummy[0]
            atoms[p0][1] = dummy[1]
            atoms[p0][2] = dummy[2]
            masses[p0] = mass_a + mass_b
            masses[p1] = -1
            keep.append(p0)
        else:
            keep.append(p0)
            keep.append(p1)
            keep_bonds.append(i)

    keep = np.array(np.sort(np.unique(keep))).astype(int)

    new_atoms = atoms[keep]
    new_bonds = bond_list[keep_bonds]
    new_elements = elements[keep]

    p1 = np.where(masses != -1)
    new_masses = masses[p1]
    return new_elements, new_atoms, new_masses, new_bonds


def assign_coos_via_distance_mat(coos_list, distance_matrix, reference):
    """
    NOT MATURE !!! (but sometimes useful)
    
    Assigns coos to suit a reference based on a distance matrix relying to the coos.
    Reference and distance matrix/ coos belong to the same molecule type but are sorted in
    different ways. Sorted rows (or cols) of a distance matrix belong to the same atom when
    they are equal.

    TO DO:
    - add check based on levensthein to omit errors through branches!!!
    - gernalizable for more than only coos???
    - bad names :(

    Parameters:
    - coos_list:
        - list of coordinates (or anything else???!!!).
    - distance_matrix:
        - distance matrix which belongs to the coos_list.
    - reference:
        - distance matrix which belongs to the reference you want to apply the coos to.

    Returns:
    - new_coos_list:
        - list of coordinates fitting the reference.
    - idx:
        - indexes to translate sth. to reference.
    """
    distance_matrix_sort = np.sort(distance_matrix, axis=1)
    reference_sort = np.sort(reference, axis=1)
    idx = []
    for ref in reference_sort:
        for i, row in enumerate(distance_matrix_sort):
            if np.array_equal(ref, row) and i not in idx:
                idx.append(i)
                break
    idx = np.array(idx).astype(int)

    print(
        """\n\nWARNING:
    Assign_coos_via_distance_mat is not mature yet.
    In branched molecules errors are conceivable because atom elements are not checked.\n
    In structurally symmetric molecules with non-symmetric atome types errors might occur \n
    Double-check your results!!! \n \n"""
    )

    return coos_list[idx], idx

def pair_from_row(row):
    """
    Gets pair potential from Joachims MC Code pair_potentials file.
    Charges are added later ;)

    Parameters:
    - row:
        - splitted line i.e. list of line in pair_potentials file.
    
    Returns:
    - pair potential dictionary
    """
    pair = {}
    pair["name"] = str(row[0])
    pair["mass"] = float(row[1])
    pair["epsilon"] = float(row[2])
    pair["sigma"] = float(row[3])
    pair["m"] = float(row[4])
    pair["cut"] = float(row[5])
    pair["charge"] = 0.0
    return pair


def pair_of_h():
    """
    Builds pair potential dummy for hydrogen atoms.
    Charges are added later ;)

    Parameters:
        none
    
    Returns:
    - hydrogen pair potential dictionary
    """
    pair = {}
    pair["name"] = "H"
    pair["mass"] = 1.0
    pair["epsilon"] = 0.0
    pair["sigma"] = 1.0
    pair["m"] = 12.0
    pair["cut"] = 0.0
    pair["charge"] = 0.0
    return pair


def get_charge(row, pair_dict):
    """
    Assigns charge to pair potential dict.
    Uses the normalized_levenshtein algorithm :)

    Parameters:
    - row:
        - splitted charge line i.e. list of line in pair_potentials file.
    - pair_dict:
        - dictionary containing pair potentials
    
    Returns:
    - key:
        - dict key i.e. atom name the charge belongs to.
    - charge:
        - float value of the charge
    """
    keys = list(pair_dict.keys())
    dummy = np.zeros(len(keys))
    for i, key in enumerate(keys):
        dummy[i] = normalized_distance(row[0][1:], key)
    pointer = np.squeeze(np.where(dummy == np.amin(dummy)))
    charge = float(row[2])
    if pointer.shape:
        print("ERROR: several atoms found for one charge:")
        print("charge:", charge)
        print("atoms:",keys[pointer])
        return None, charge
    else:
        return keys[pointer], charge


def read_pair_potentials(path):
    """
    Reads pair potentials from Joachims MC Code pair_potentials file.

    Parameters:
    - path:
        - path to pair_potentials file.
            
    Returns:
    - pair_dict:
        - dictionary containing pair potentials
    """
    keylist = ["!", "#", "models:", "types", "References"]
    pair_dict = {}
    flag = -1
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")
        for row in spamreader:
            row = [x for x in row if x]
            if row:
                if "VdW-site" in row:
                    flag = 1
                elif "coulomb-site" in row:
                    flag = 2
                elif row[0] in keylist:
                    break
                elif flag == 1 and row:
                    pair_dict[row[0]] = pair_from_row(row)
                elif flag == 2 and row:
                    # print(row)
                    if row[0].split("_")[0] == "cH":
                        pair_dict[row[0]] = pair_of_h()
                        pair_dict[row[0]]["charge"] = float(row[2])
                    else:
                        p, ch = get_charge(row, pair_dict)
                        # print(p,ch)
                        pair_dict[p]["charge"] = ch
    return pair_dict


def read_bond_potentials(path):
    """
    Reads bond potentials from Joachims MC Code bond_potentials file.

    Parameters:
    - path:
        - path to bond_potentials file.
        
    Returns:
    - bond_dict:
        - dictionary containing bond potentials
    """
    keylist = ["!", "#", "models:", "types"]
    bond_dict = {}
    flag = 0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")
        for row in spamreader:
            row = [x for x in row if x]
            if row:
                if row[0] == "model":
                    flag = 1
                elif flag == 1 and row[0] in keylist:
                    break
                elif flag == 1:
                    dummy = sort_force_fields(row[1:3])
                    name = make_graph(dummy)
                    bond_dict[name] = {}
                    bond_dict[name]["list"] = row[1:3]
                    bond_dict[name]["type"] = int(row[0])
                    # bond_dict[name]["len"] = float(row[3])
                    # bond_dict[name]["spring"] = float(row[4])
                    bond_dict[name]["p"] = [float(r) for r in row[3:]]

    return bond_dict


def read_angle_potentials(path):
    """
    Reads angle potentials from Joachims MC Code angle_potentials file.

    Parameters:
    - path:
        - path to angle_potentials file.
        
    Returns:
    - angle_dict:
        - dictionary containing angle potentials
    """
    keylist = ["!", "#", "models:", "types", "Note", "Note:"]
    angle_dict = {}
    flag = 0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")
        for row in spamreader:
            row = [x for x in row if x]
            if row:
                if row[0] == "model":
                    flag = 1
                elif flag == 1 and row[0] in keylist:
                    break
                elif flag == 1:
                    dummy = sort_force_fields(row[1:4])
                    name = make_graph(dummy)
                    angle_dict[name] = {}
                    angle_dict[name]["list"] = row[1:4]
                    angle_dict[name]["type"] = int(row[0])
                    # angle_dict[name]["angle"] = float(row[4])
                    # angle_dict[name]["p"] = float(row[5])
                    angle_dict[name]["p"] = [float(r) for r in row[4:]]

    return angle_dict


def read_torsion_potentials(path):
    """
    Reads torsion potentials from Joachims MC Code torsion_potentials file.

    Parameters:
    - path:
        - path to torsion_potentials file.
        
    Returns:
    - torsion_dict:
        - dictionary containing torsion potentials
    """
    keylist = ["!", "#", "models:", "types", "Note:"]
    torsion_dict = {}
    flag = 0
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")
        for row in spamreader:
            row = [x for x in row if x]
            if row:
                if row[0] == "model":
                    flag = 1
                elif flag == 1 and row[0] in keylist:
                    break
                elif flag == 1:
                    dummy = sort_force_fields(row[1:5])
                    name = make_graph(dummy)
                    torsion_dict[name] = {}
                    torsion_dict[name]["list"] = row[1:5]
                    torsion_dict[name]["type"] = int(row[0])
                    torsion_dict[name]["p"] = [float(r) for r in row[5:]]

    return torsion_dict


def read_xyz(path, energy=False):
    """
    Reads xyz file. Supports tubomole xyz with energy in header

    Parameters:
    - path:
        - path to xyz file.
    - energy:
        - returns energy if true.
        - use this with turbomol QM-results.
    
    Returns:
    - xyz:
        - list of dicts, keys: "atom": atom name and "xyz": coordinates.
    - energy:
        - (optional) if energy = True (bad style ?)
    """
    data = []
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ")
        for row in spamreader:
            data.append([r for r in row if r])
    n = int(data[0][0])
    xyz = {}
    for i, d in enumerate(data[2 : n + 2]):
        xyz[i] = {}
        xyz[i]["atom"] = d[0]
        xyz[i]["xyz"] = np.array([float(x) for x in d[1:4]])

    if energy:
        energy = float(data[1][2])
        return xyz, energy

    return xyz


def assign_CHx(xyz):
    """
    Builds united atom groups fromm all atom xyz.
    Works with turbomol results and read_xyz def.

    Parameters:
    - xyz:
        - list of atom dicts, keys: "atom": atom name and "xyz": coordinates.
        
    Returns:
    - xyz_CHx:
        - list of united atom dicts, keys: "atom": atom name and "xyz": coordinates.
    """
    xyz_CHx = {}
    xkeys = sorted(xyz.keys())
    ii = 0
    x = np.min(xkeys)
    flag = 0
    while x <= np.max(xkeys):
        if flag == 0 and xyz[x]["atom"] == "C":
            coos = np.array(xyz[x]["xyz"])
            flag = 1
            x += 1
        elif flag > 0:
            if xyz[x]["atom"] == "H":
                coos = np.column_stack((coos, xyz[x]["xyz"]))
                x += 1
                flag += 1
            else:
                no = flag - 1
                xyz_CHx[ii] = {}
                xyz_CHx[ii]["atom"] = "CH" + str(no)
                xyz_CHx[ii]["xyz"] = np.sum(
                    coos * np.array([15] + [1] * no), axis=1
                ) / (15 + 1 * no)
                coos = []
                flag = 0
                ii += 1
        else:
            xyz_CHx[ii] = xyz[x]
            x += 1
            ii += 1
    print(
        """\n\nWARNING:
    Very basic tool. Only works with well-sorted coordinates. Your better off doing it by hand. At least check your results!!! \n \n"""
    )
    return xyz_CHx


def kill_CHx(xyz):
    """
    Kills C-bonded hydrogens fromm all atom xyz.
    Works with turbomol results and read_xyz def.

    Parameters:
    - xyz:
        - list of atom dicts, keys: "atom": atom name and "xyz": coordinates.
    
    Returns:
    - xyz_CHx:
        - list of atom dicts, keys: "atom": atom name and "xyz": coordinates.
    """
    xyz_CHx = {}
    xkeys = sorted(xyz.keys())
    ii = 0
    x = np.min(xkeys)
    flag = 0
    while x <= np.max(xkeys):
        if flag == 0 and xyz[x]["atom"] == "C":
            coos = np.array(xyz[x]["xyz"])
            flag = 1
            x += 1
        elif flag > 0:
            if xyz[x]["atom"] == "H":
                x += 1
                flag += 1
            else:
                no = flag - 1
                xyz_CHx[ii] = {}
                xyz_CHx[ii]["atom"] = "C"
                xyz_CHx[ii]["xyz"] = coos
                coos = []
                flag = 0
                ii += 1
        else:
            xyz_CHx[ii] = xyz[x]
            x += 1
            ii += 1
    print(
        """\n\nWARNING:
    Very basic tool. Only works with well-sorted coordinates. Your better off doing it by hand. At least check your results!!! \n \n"""
    )
    return xyz_CHx


def to_xyz(xyz, path):
    """
    Writes xyz file.
    Works with turbomol results and read_xyz def.

    Parameters:
    - xyz:
        - list of atom dicts, keys: "atom": atom name and "xyz": coordinates.
    - path:
        - path to xyz file.
        
    Returns:
        nothing
    """
    f = open(path, "w")
    f.write(str(len(xyz)) + "\n")
    f.write("\n")
    for x in xyz:
        line = "    ".join([xyz[x]["atom"][0]] + [str(y) for y in xyz[x]["xyz"]])
        f.write(line + "\n")
    f.close()
    return


def distance(x1, x2):
    """
    Returns spatial distance.

    Parameters:
    - x1, x2:
        - vectors.
        
    Returns:
    - spatial distance between x1 and x2.
    """
    return np.linalg.norm(x1 - x2)


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    Parameters:
    - vector:
        - vector.
        
    Returns:
    - unit vector
    """
    return vector / np.linalg.norm(vector)


def angle_between(x1, x2, x3):
    """
    Returns the angle in radians between vectors x1, x2, x3.

    Parameters:
    - x1, x2, x3:
        - vectors.
            
    Returns:
    - rad angle
    """
    v1 = unit_vector(x1 - x2)
    v2 = unit_vector(x3 - x2)

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def dihedral(x0, x1, x2, x3):
    """
    Returns the dihedral angle in radians between vectors x0, x1, x2, x3.
    Praxeolitic formula -> 1 sqrt, 1 cross product

    Parameters:
    - x0, x1, x2, x3:
        - vectors.
    
    Returns:
    - deg angle
    """
    b0 = -1.0 * (x1 - x0)
    b1 = x2 - x1
    b2 = x3 - x2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    # print(b1)
    # print(np.linalg.norm(b1))
    if np.linalg.norm(b1) > 0:
        b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)
