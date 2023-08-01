"""
moleculegraph

A graph representation for molecules and everything else for data mapping
"""

__version__ = "0.0.0"
__author__ = "Maximilian Fleck"
__credits__ = "ITT Uni Stuttgart"


import numpy as np
import csv
from .molecule_utils import *
from .funcs import *

# from rapidfuzz.string_metric import normalized_levenshtein
import json


class molecule:
    """
    write some random stuff here... lalalala
    """

    def __init__(self, mol):
        """
        initializes molecules based on graph-list representation.
        Future: allow different input types i.e. string, list...

        Main results of the initialization are:
            - distance matrix of the molecule.
            - atom, bond, angle, torsion lists.
        ...you can use those to:
            - map force-field information onto your molecule.
            - generate inputs for templating trajectories and coordinates
            - use group contribution approaches
            - process coordinates from the PubChem-API
        ...furthermore:
            this graph-representation is quite general and not limited to molecules.
            You might use it to represent your favorite baking recipes and optimize them.
        ...GOALS:
            this reprsentation is easily translatable to SMILEYS and SMILES-Codes.
            SILEYS is powerful when it comes to chemical space exploration *zwinkersmiley*

        TO DO (??):
        - function pointer i.e. r -> _r and b -> _b or something like this.
        - OR: allow function dict ;)

        Args:
            mol:
                -list representation of a molecule.
        Returns:
            nothing
        """
        if isinstance(mol, str):
            mol = mol[1:-1].split("][")

        self.molecule = np.array(mol)
        self.f = np.zeros(len(mol))  # shows function
        self.n = -1 * np.ones(len(mol))  # shows atomnumber
        self.i = np.arange(len(mol))  # shows index
        self.len = len(self.i)
        n = 0
        for i, m in enumerate(mol):
            if m[0] == "b":
                self.f[i] = int(m[1:])
            elif m[0] == "r":
                # dum = -int(m[1:])
                self.f[i] = -int(m[1:])  # self.get_atom_index()
            else:
                self.n[i] = n
                n += 1

        # contains indexes of atoms
        # index of atom no n: self.atom_indexes[n]
        self.atom_indexes = self.i[self.f == 0]
        self.atom_names = self.molecule[self.f == 0]
        self.atom_number = len(self.atom_indexes)
        self.atom_numbers = np.arange(self.atom_number)

        # get bond list
        self.get_neighbour_list()
        # get branch ends
        n, count = np.unique(self.bond_list, return_counts=True)
        pp = np.squeeze(np.where(count == 1))
        if np.sum(pp) > 0:
            self.branch_end_numbers = n[pp]
            self.branch_end_indexes = self.atom_indexes[self.branch_end_numbers]
        else:
            self.branch_end_numbers = np.empty(0)
            self.branch_end_indexes = np.empty(0)
        # get ring closures
        self.ring_close_indexes = np.squeeze(np.where(self.f < 0)) - 1
        if np.sum(self.ring_close_indexes.shape) > 0:
            self.ring_close_numbers = self.n[self.ring_close_indexes]
            self.ring_root_numbers = self.n[self.ring_root_indexes]
        else:
            self.ring_close_numbers = np.empty(0)
            self.ring_root_numbers = np.empty(0)

        # if self.neighbour_list.size > 0:
        self.get_distance_matrix()
        # else:
        #    print("either small molecule or error")

        """ heres an inconsistency: 
            indexes refer to indexes in the graph.
            numbers refer to atom nubers etc. i.e. indexes excluding functionals like b=branch and r=ring
            ... here numbers are called indexes.
            ... is this ok? Is the whole concept intuitive?
        """
        dummy = unique_sort(self.atom_names, return_inverse=True)
        (
            self.unique_atom_keys,
            self.unique_atom_indexes,
            self.unique_atom_inverse,
        ) = dummy
        self.unique_atom_names = self.atom_names[self.unique_atom_indexes]
        self.unique_atom_numbers = self.atom_numbers[self.unique_atom_indexes]
        assert np.array_equal(
            self.unique_atom_keys[self.unique_atom_inverse], self.atom_names
        ), "unique atom in inverse doesnt work"

        dummy = unique_sort(self.bond_keys, return_inverse=True)
        (
            self.unique_bond_keys,
            self.unique_bond_indexes,
            self.unique_bond_inverse,
        ) = dummy
        self.unique_bond_names = self.bond_names[self.unique_bond_indexes]
        self.unique_bond_numbers = self.bond_list[self.unique_bond_indexes]
        if len(self.bond_keys) > 0:
            assert np.array_equal(
                self.unique_bond_keys[self.unique_bond_inverse], self.bond_keys
            ), "unique bond in inverse doesnt work"

        dummy = unique_sort(self.angle_keys, return_inverse=True)
        (
            self.unique_angle_keys,
            self.unique_angle_indexes,
            self.unique_angle_inverse,
        ) = dummy
        self.unique_angle_names = self.angle_names[self.unique_angle_indexes]
        self.unique_angle_numbers = self.angle_list[self.unique_angle_indexes]
        if len(self.angle_keys) > 0:
            assert np.array_equal(
                self.unique_angle_keys[self.unique_angle_inverse], self.angle_keys
            ), "unique angle in inverse doesnt work"

        dummy = unique_sort(self.torsion_keys, return_inverse=True)
        (
            self.unique_torsion_keys,
            self.unique_torsion_indexes,
            self.unique_torsion_inverse,
        ) = dummy
        self.unique_torsion_names = self.torsion_names[self.unique_torsion_indexes]
        self.unique_torsion_numbers = self.torsion_list[self.unique_torsion_indexes]
        if len(self.torsion_keys) > 0:
            assert np.array_equal(
                self.unique_torsion_keys[self.unique_torsion_inverse], self.torsion_keys
            ), "unique angle in inverse doesnt work"

        self.unique_atom_pair_names = []
        self.unique_atom_pair_keys = []
        for i, a0 in enumerate(self.unique_atom_names):
            for j, a1 in enumerate(self.unique_atom_names):
                if i > j:
                    self.unique_atom_pair_names.append(sort_force_fields([a0, a1]))
                    self.unique_atom_pair_keys.append(
                        make_graph(self.unique_atom_pair_names[-1])
                    )
        self.unique_atom_pair_names = np.array(self.unique_atom_pair_names)
        self.unique_atom_pair_keys = np.array(self.unique_atom_pair_keys)
        return

    def get_neighbour_list(self):
        """
        generates a neighbour list from an initialized molecule

        Args:
            self:
                -initialized molecule.
        Returns:
            neighbour_list:
                -list of neighboured atoms containing their atom numbers.
            idx_neighbour_list:
                -list of neighboured atoms containing their atom indexes.
        """

        idx_neighbour_list = []
        idx_branch_ends = np.empty(0)
        idx_branch_starts = np.empty(0)
        idx_ring_close = np.empty(0)
        idx_branch_funs = np.empty(0)

        branches = np.zeros(len(self.i))
        branch_origins = np.zeros(len(self.i))
        branch_count = 1

        if self.atom_number < 2:
            self.idx_neighbour_list = np.array([])
            self.neighbour_list = np.array([])
            self.idx_bond_list = np.array([])
            self.bond_list = np.array([])
            return np.array([]), np.array([])

        """
        find branches and assign numbers.
        get main path through molecule.
        build molecule from this...
        """
        inv_f = self.f[::-1]
        branch_no = 1
        inv_f_p = np.atleast_1d(np.squeeze(np.where(inv_f != 0)))
        inv_branches = inv_f.copy()
        inv_branches[inv_branches != 0] = -1
        for i, ifi in zip(inv_f_p, inv_f[inv_f_p]):
            if ifi < 0:  # ring
                inv_branches[i] = -1
            if ifi > 0:  # branch :)
                pp = np.atleast_1d(np.squeeze(np.where(inv_branches[:i] == 0)))
                if len(pp) > ifi:
                    p_branch = pp[-int(ifi) :]
                    inv_branches[p_branch] = branch_no
                    inv_branches[i] = branch_no
                    branch_no += 1
        branches = reenumerate_branches(inv_branches[::-1])
        print(branches)

        main_path_p = np.where(branches == 0)
        main_path_i = self.i[main_path_p]

        # get bond list of main path
        bond_list = bond_list_from_simple_path(main_path_i)
        # print( "branches",branches )

        """
        build bond list
        """
        f_p = np.atleast_1d(np.squeeze(np.where(self.f != 0)))
        self.ring_root_indexes = []
        print("f_p , self.f[f_p]", f_p, self.f[f_p])
        for i, fi in zip(f_p, self.f[f_p]):
            print(branches)
            if fi < 0:
                idx = np.max(
                    get_diff_in_bond_lists(f_p[f_p < i], np.arange(i))
                )  # index of last atom
                branches[i] = get_branch_root(
                    branches[idx], branches
                )  # connect ring to acompanying branch
                subgraph = graph_from_bonds(bond_list)
                path_to_main = get_shortest_path(subgraph, 0, idx)

                bond_list_to_main = bond_list_from_simple_path(path_to_main)
                graph_to_main = graph_from_bonds(bond_list_to_main)
                root = get_root(graph_to_main, idx, int(np.abs(fi)))
                if root < 0:
                    print("root of ring outside molecule. ignore.")
                else:
                    bond_list = np.concatenate([bond_list, [[root, idx]]])
                    self.ring_root_indexes.append(root)
            elif fi > 0:
                link = get_branch_root(branches[i], branches)
                branch = branches[i]
                p_subchain_i = np.atleast_1d(np.squeeze(np.where(branches == branch)))[
                    1:
                ]
                if len(p_subchain_i) < fi:
                    print("end of branch outside molecule. ignore.")
                else:
                    p_subchain_i = np.concatenate([[link], p_subchain_i])
                    subchain_bond_list = bond_list_from_simple_path(p_subchain_i)
                    bond_list = np.concatenate([bond_list, subchain_bond_list])
            else:
                print("fi = 0... this should not happen :P")
        print("fin branches", branches)
        self.idx_neighbour_list = bond_list
        self.idx_bond_list = bond_list
        self.ring_root_indexes = np.array(self.ring_root_indexes)

        self.neighbour_list = np.vectorize(self.get_atom_no)(self.idx_neighbour_list)
        self.bond_list = self.neighbour_list
        return self.neighbour_list, self.idx_neighbour_list

    def get_bond_matrix(self):
        """
        generates bond matrix from an initialized molecule

        Args:
            self:
                -initialized molecule.
        Returns:
            bond_matrix:
                - bond matrix of you molecule
        """
        self.bond_matrix = get_bond_matrix(self.neighbour_list, self.atom_number)
        return self.bond_matrix.copy()

    def get_distance_matrix(self):
        """
        generates distance matrix from an initialized molecule.
        Furthermore different angle and torsion lists are generated.

        Args:
            self:
                -initialized molecule.
        Returns:
            distance_matrix:
                - distance matrix of you molecule
        """
        self.distance_matrix, self.bond_lists = get_distance_matrix(
            self.neighbour_list, self.atom_number
        )
        self.bond_list = self.bond_lists[0]
        self.bond_names = np.array(
            [sort_force_fields(self.atom_names[a]) for a in self.bond_list]
        )
        self.bond_keys = np.array([make_graph(t) for t in self.bond_names])
        if len(self.bond_lists) > 1:
            self.angle_list = np.array(self.bond_lists[1])
            self.angle_names = np.array(
                [sort_force_fields(self.atom_names[a]) for a in self.angle_list]
            )
            self.angle_keys = np.array([make_graph(t) for t in self.angle_names])
        else:
            self.angle_list = np.array([])
            self.angle_names = np.array([])
            self.angle_keys = np.array([])
        if len(self.bond_lists) > 2:
            self.torsion_list = np.array(self.bond_lists[2])
            self.torsion_names = np.array(
                [sort_force_fields(self.atom_names[a]) for a in self.torsion_list]
            )
            self.torsion_keys = np.array([make_graph(t) for t in self.torsion_names])
        else:
            self.torsion_list = np.array([])
            self.torsion_names = np.array([])
            self.torsion_keys = np.array([])
        return self.distance_matrix.copy()

    def get_atom_index(self, no):
        """
        returns index of atom.
        Number: atom number
        Index:  index in list representation of molecule.
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            no:
                -number of the atom
        Returns:
            index of atom nuber "no"
        """
        return int(self.atom_indexes[int(no)])  # int(np.squeeze(np.where(self.n==no)))

    def get_atom_no(self, idx):
        """
        returns index of atom.
        Number: atom number
        Index:  index in list representation of molecule.
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            idx:
                -idx of the atom
        Returns:
            number of atom with index
        """
        idx = int(idx)
        if self.f[idx] >= 0:
            no = int(np.squeeze(np.where(self.atom_indexes == idx)))
            assert no == self.n[idx], "atom numbers dont work"
            return no
        else:
            assert self.n[idx] == -1, "atom numbers dont work"
            print("Error in get_atom_no: not an atom")
            return False

    def get_next_atom_idx(self, idx):
        """
        returns index of the next atom after index idx.
        Index:  index in list representation of molecule.
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            idx:
                -index in list representation of molecule
        Returns:
            index of the next atom after index idx
        """
        return int(np.min(self.atom_indexes[(self.atom_indexes - idx) > 0]))

    def get_last_atom_idx(self, idx):
        """
        returns index of the atom before index idx.
        Index:  index in list representation of molecule.
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            idx:
                -index in list representation of molecule
        Returns:
            index of the atom before index idx
        """
        return int(np.max(self.atom_indexes[(self.atom_indexes - idx) < 0]))

    def get_next_atom_no(self, idx):
        """
        returns number of the atom after index idx.
        Index:  index in list representation of molecule.
        Number: atom number
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            idx:
                -index in list representation of molecule
        Returns:
            number of the atom after index idx
        """
        return int(self.n[self.get_next_atom_idx(idx)])

    def get_last_atom_no(self, idx):
        """
        returns number of the atom before index idx.
        Index:  index in list representation of molecule.
        Number: atom number
        Number and index are differing due to functionals in the list.

        Args:
            self:
                -initialized molecule.
            idx:
                -index in list representation of molecule
        Returns:
            number of the atom before index idx
        """
        return int(self.n[self.get_last_atom_idx(idx)])

    def map_molecule_via_atom_names(self, data):
        """
        maps information from a dictionary onto the molecule.
        Use this for molecule information i.e. coordinates, group contributions,...

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            data:
                -dictionary with: keys == atom names !!!
        Returns:
            mapped data
        """
        return [*map(data.get, self.atom_names)]

    def map_unique_via_atom_names(self, data):
        """
        Maps information from a dictionary onto unique elements of the molecule.
        Preserves the order of the original list.
        Use this for force-field information i.e. potentials,...

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            data:
                -dictionary with: keys == atom names !!!
        Returns:
            mapped unique data
        """
        print(
            "map_unique_via_atom_names is outtdated... dont use\n shouldnt cause errors btw"
        )
        return [*map(data.get, self.unique_atom_names)]

    def map_molecule(self, goal, data):
        """
        maps information from a dictionary onto a goal.
        Use this for molecule information i.e. coordinates, group contributions,...

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            goal:
                -list of keys
            data:
                -dictionary with keys in goal
        Returns:
            mapped data
        """
        return [*map(data.get, goal)]

    def map_unique(self, goal, data):
        """
        Maps information from a dictionary onto unique elements a goal.
        Preserves the order of the original list.
        Use this for force-field information i.e. potentials,...

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            goal:
                -list of keys
            data:
                -dictionary with keys in goal
        Returns:
            mapped unique data
        """
        print("map_unique is outtdated... dont use\n shouldnt cause errors btw")
        return [*map(data.get, unique_sort(goal))]

    def map_nested(self, goal, data):
        """
        Maps information from a dictionary onto nested elements of a goal.
        Preserves the order of the original list.
        Use this for advanced force-field information i.e. mixing rules,...

        same as:
        [ self.map_unique( m, data) for m in goal ]

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            goal:
                -with list list of keys
            data:
                -dictionary with keys in goal
        Returns:
            mapped unique data
        """
        return [[*map(data.get, m)] for m in goal]

    def get_distance(self, x, y):
        """
        Returns distance between atoms no x and y
        0: same atom i.e. x==y
        1: bond
        2: angle
        3: torsion

        TO DO:
        - checks for easy debugging

        Args:
            self:
                -initialized molecule.
            x,y:
                -atom numbers
        Returns:
            distance
        """
        return self.distance_matrix[int(x)][int(y)]

    def get_sorted_key(self, x):
        """
            Sorts a list alphabetically.
            Very important for dictionary keys in the molecule class!!!

        Args:
            x:
                -string list to sort.
        Returns:
            key (string) and sorted list
        """
        x_sort = sort_force_fields(x)
        return make_graph(x_sort), x_sort

    def visualize(self, saveto=""):
        """
        plots and visualizes molecule
        Args:
            self:
                -initialized molecule.
        Kwargs:
            saveto:
                - optional, location to save pic of graph to,
                - add file format at the end ;)
        Returns:
            nothing
        """
        if self.atom_numbers.size > 1:
            plot_graph(self.atom_names, self.bond_list, saveto=saveto)
        else:
            print("one bead only")
        return
