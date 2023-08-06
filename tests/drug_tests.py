import moleculegraph
import numpy as np
import toml
import unittest

"""
MDMA
"""

mdma = "[CH3][C][b2][HN][CH3][C][C1][C2][C3][C4][b3][O1][C][O2][r5][C5][C6][r6]"

molecule = moleculegraph.molecule(mdma)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize(saveto="mdma.pdf")
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

mdma_expected = [
    "[C][CH3]",
    "[C][HN]",
    "[C][C]",
    "[CH3][HN]",
    "[C1][C]",
    "[C1][C2]",
    "[C1][C6]",
    "[C2][C3]",
    "[C3][C4]",
    "[C3][O2]",
    "[C4][O1]",
    "[C4][C5]",
    "[C][O1]",
    "[C][O2]",
    "[C5][C6]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(mdma_expected)), "lists not equal"
)

"""
Ketamin
"""
keta = "[CH3][NH][C1][b5][C2][b1][=O][C3][C4][C5][C6][r6][CC1][CC2][b1][Cl][CC3][CC4][CC5][CC6][r6]"

molecule = moleculegraph.molecule(keta)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

keta_expected = [
    "[CH3][NH]",
    "[C1][NH]",
    "[C1][C2]",
    "[C1][C6]",
    "[C1][CC1]",
    "[=O][C2]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][C5]",
    "[C5][C6]",
    "[CC1][CC2]",
    "[CC1][CC6]",
    "[CC2][Cl]",
    "[CC2][CC3]",
    "[CC3][CC4]",
    "[CC4][CC5]",
    "[CC5][CC6]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(keta_expected)), "lists not equal"
)

"""
Ketamin Mutation for testing reasons
"""
mut_keta = "[CH3][NH][C1][b6][C2][b1][=O][C3][C4][C5][C6][r6][=O][CC1][CC2]"  # "[b1][Cl][CC3][CC4][CC5][CC6][r6]"

molecule = moleculegraph.molecule(mut_keta)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize(saveto="mut.pdf")
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

mut_keta_expected = [
    "[CH3][NH]",
    "[C1][NH]",
    "[C1][C2]",
    "[C1][C6]",
    "[C1][CC1]",
    "[=O][C2]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][C5]",
    "[C5][C6]",
    "[=O][C6]",
    "[CC1][CC2]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(mut_keta_expected)), "lists not equal"
)


"""
Ketamin Fold for testing reasons
"""
print("roots outside molecule expected")
keta = "[CH3][NH][C][C1][b5][C2][r16][C3][C4][C5][C6][r6][C][CC1][b5][CC2]"  # "[b1][Cl][CC3][CC4][CC5][CC6][r6]"

molecule = moleculegraph.molecule(keta)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)

lllong_expected = [
    "[C1][C2]",
    "[C1][C6]",
    "[C1][C]",
    "[C1][C]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][C5]",
    "[C5][C6]",
    "[CC1][CC2]",
    "[CC1][C]",
    "[CH3][NH]",
    "[C][NH]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(lllong_expected)), "lists not equal"
)

"""
Ketamin Fold for testing reasons no 2
"""
print("roots outside molecule expected")
keta = "[CH3][NH][C][C1][C2][b4][C3][C4][C5][C6][r6][C][CC1][b5][CC2]"  # "[b1][Cl][CC3][CC4][CC5][CC6][r6]"

molecule = moleculegraph.molecule(keta)
print(np.array([molecule.atom_numbers.astype(float)]))
print(molecule.distance_matrix)
# print( molecule.bond_matrix )
print(molecule.bond_list)
molecule.visualize()
print(molecule.angle_keys)
print(molecule.torsion_keys)
print(molecule.bond_keys)
lllong_expected2 = [
    "[CH3][NH]",
    "[C][NH]",
    "[C1][C]",
    "[C1][C2]",
    "[C2][C]",
    "[CC1][C]",
    "[CC1][CC2]",
    "[C2][C3]",
    "[C3][C4]",
    "[C4][C5]",
    "[C5][C6]",
    "[C1][C6]",
]

tc = unittest.TestCase()
tc.assertListEqual(
    sorted(list(molecule.bond_keys)), sorted(list(lllong_expected2)), "lists not equal"
)
