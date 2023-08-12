import moleculegraph
import numpy as np
import toml
import unittest


"""
gly test:
"""
gly = "[OH][C][b1][=O][C][NH2]"

molecule = moleculegraph.molecule(gly)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize("gly.png")
print(molecule.angle_keys)
print(molecule.torsion_keys)

gly_expected = ["[C][OH]", "[C][C]", "[C][NH2]", "[=O][C]"]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(gly_expected)), "gly lists not equal"
)

"""
alanine test:
"""

ala = "[OH][C][b1][=O][C][b1][C][NH2]"

molecule = moleculegraph.molecule(ala)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)

ala_expected = ["[C][OH]", "[=O][C]", "[C][C]", "[C][C]", "[C][NH2]"]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(ala_expected)), "ala lists not equal"
)

"""
arginine test:
"""

arg = "[OH][C][b1][=O][C][b1][NH2][C][C][C][NH][b1][NH][H2N]"

molecule = moleculegraph.molecule(arg)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)

arg_expected = [
    "[C][OH]",
    "[=O][C]",
    "[C][C]",
    "[C][NH2]",
    "[C][C]",
    "[C][C]",
    "[C][C]",
    "[C][NH]",
    "[NH][NH]",
    "[H2N][NH]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(arg_expected)),
    "arginine lists not equal",
)

"""
phenylalanine test:
"""

phe = "[OH][C][b1][=O][C][b1][NH2][C][C1][C2][C3][C4][C5][C6][r6]"

molecule = moleculegraph.molecule(phe)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)

phe_expected = [
    "[C][OH]",
    "[=O][C]",
    "[C][C]",
    "[C][NH2]",
    "[C][C]",
    "[C1][C]",
    "[C1][C2]",
    "[C1][C6]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][C5]",
    "[C5][C6]",
]
tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(phe_expected)),
    "phenylalanine lists not equal",
)

"""
tyrosine test:
"""

tyr = "[OH][C][b1][=O][C][b1][NH2][C][C1][C2][C3][C4][b1][OH][C5][C6][r6]"

molecule = moleculegraph.molecule(tyr)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)

tyr_expected = [
    "[C][OH]",
    "[=O][C]",
    "[C][C]",
    "[C][NH2]",
    "[C][C]",
    "[C1][C]",
    "[C1][C2]",
    "[C1][C6]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][OH]",
    "[C4][C5]",
    "[C5][C6]",
]
tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(tyr_expected)),
    "tyrosine lists not equal",
)

"""
tryptophane test:
"""

trp = (
    "[OH][C][b1][=O][C][b1][NH2][C][C1][C2][C3][b4][CC1][CC2][CC3][CC4][r6][NH][C5][r5]"
)

molecule = moleculegraph.molecule(trp)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

trp_expected = [
    "[C][OH]",
    "[=O][C]",
    "[C][C]",
    "[C][NH2]",
    "[C][C]",
    "[C1][C]",
    "[C1][C2]",
    "[C1][C5]",
    "[C2][C3]",
    "[C2][CC4]",
    "[C3][CC1]",
    "[C3][NH]",
    "[CC1][CC2]",
    "[CC2][CC3]",
    "[CC3][CC4]",
    "[C5][NH]",
]
tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(trp_expected)), "lists not equal"
)

"""
CH3 tests:
"""

CH3 = "[C][b1][H][b1][H][b1][H][C][C][C]"

molecule = moleculegraph.molecule(CH3)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

CH3_expected = ["[C][C]", "[C][C]", "[C][C]", "[C][H]", "[C][H]", "[C][H]"]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(CH3_expected)),
    "CH3 case 0 lists not equal",
)

CH3 = "[C][b1][H][b1][H][b1][H]"

molecule = moleculegraph.molecule(CH3)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

CH3_expected = ["[C][H]", "[C][H]", "[C][H]"]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(CH3_expected)),
    "CH3 case 1 lists not equal",
)

CH2 = "[C][C][b1][H][H][b1][H][C]"

molecule = moleculegraph.molecule(CH2)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

CH3_expected = ["[C][C]", "[C][H]", "[C][H]", "[C][H]", "[H][H]"]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)),
    sorted(list(CH3_expected)),
    "CH2 case xxx lists not equal",
)
