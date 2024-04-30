import numpy as np
import csv
import json
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem

"""
functions to treat moleculestring:
"""
def split_molstring(molstring):
    return molstring[1:-1].split("][")


def get_syntax_from_numbers(mol_array, 
                            branch_operators=["b"], 
                            ring_operators=["r"]):
    """
    - gets syntactic elemnts from a splitted molstring
    - rings and brances are marked with letters followed
      by a number giving the size
    """
    n=0
    ff = np.zeros(len(mol_array))
    nn = -1 * np.ones(len(mol_array))
    for i,m in enumerate( mol_array ):
        if re.sub(r"\d+","", m) in branch_operators:
            ff[i] = int(re.sub("[^0-9]","", m ))
        elif re.sub(r"\d+","", m) in ring_operators:
            ff[i] = -int(re.sub("[^0-9]","", m ))
        else:
            nn[i] = n
            n+=1
    return ff ,nn

"""
collection of functions needed by the molecule class.

- sorting functionalities
- generate bond lists and matrices
- plotting
- smiles interfacing
- tools to generate equivalent graphstrings

"""      

def sort_force_fields(x):
    """
    sorts force field i.e. atom lists alphabetically.
    Very important for dictionary keys in the molecule class!!!

    Parameters:
    - x:
        - string list to sort.
    
    Returns:
    - sorted list
    """

    # old... might cause errors
    # p = np.argsort( [ x[0], x[-1] ] )

    # new... should omit errors
    p = np.argsort(["".join(x), "".join(np.flip(x))])
    if p[0] <= p[1]:
        return np.array(x)
    else:
        return np.flip(x)


def sort_graph_key(x):
    """
    sorts string i.e. the encapsulated atom lists alphabetically.
    Very important for dictionary keys in the molecule class!!!

    Parameters:
    - x:
        -string to sort.
    
    Returns:
    - sorted list
    """

    x = x[1:-1].split("][")
    x = sort_force_fields(x)
    return "[" + "][".join(x) + "]"

def read_json(path):
    """
    well...

    Parameters:
    - path:
        -path to json file.
    
    Returns:
    - dictionary, list of dictionaries (or whatever is in the file)
    """
    with open(path) as json_file:
        return json.load(json_file)


def unique_sort(x, return_inverse=False):
    """
    Returns unique elements of a list array.
    
    Important:
    
    not the list is sorted but the indexes
    i.e. this function preserves the order
    of the original list.

    Parameters:
    - X:
        - list array.
    - return_inverse:
        - boolean
        
    Returns:
    - unique list array
    ( and index and inverse if return_inverse = True )
    """
    _, index = np.unique(x, return_index=True)
    p = np.sort(index)
    p = np.array(p).astype(int)
    unique_names = np.array(x)[p]
    inverse = np.array([np.squeeze(np.where(unique_names == xx)) for xx in x])
    if return_inverse:
        return unique_names, p, inverse
    else:
        return unique_names


def get_bond_list(distance_matrix, distance=1):
    """
    generates a list of atom pairs fullfilling a certain
    distance criterium.

    Parameters:
    - distance_matrix:
        - distance matrix.
    - distance:
        - distance to look for
            - 1: bond
            - 2: angle
            - 3: torsion/dihedral
            
    Returns:
    - bond list
    """
    blist = np.array(np.where(distance_matrix == distance)).T
    return np.unique(blist, axis=0)


def get_distances(y, xlist, size):
    """
    used in get_distance_list.

    Parameters:
    - y:
        - hmmm.
    - xlist:
        - hmmm
    - size:
        - hmmm
        
    Returns:
    - hmmm
    """
    md = []
    dummy = np.array([])
    for x in xlist:
        if x[-1] == y[0]:
            dummy = np.concatenate((x, y[1:]))
        elif x[0] == y[-1]:
            dummy = np.concatenate((y, x[1:]))
        elif x[-1] == y[-1]:
            x = np.flip(x)
            dummy = np.concatenate((y, x[1:]))
        elif x[0] == y[0]:
            x = np.flip(x)
            dummy = np.concatenate((x, y[1:]))
        if np.unique(dummy).size == size:
            md.append(dummy)
    return np.unique(np.array(md), axis=0)


def get_distance_list(xlist, ylist, size):
    """
    used in get_distance_matrix.

    Parameters:
    - y:
        - hmmm.
    - xlist:
        - hmmm
    - size:
        - hmmm
        
    Returns:
    - hmmm
    """
    dlist = []
    for y in ylist:
        dummy = get_distances(y, xlist, size)
        if dummy.size > 0:
            for d in dummy:
                dlist.append(d)
    dlist = np.unique(np.array(dlist), axis=0)
    if dlist.size > 1:
        _, p = np.unique(np.sort(dlist, axis=1), axis=0, return_index=True)
        return dlist[p]
    else:
        return np.array([])


def get_distance_matrix(neighbour_list, atom_number):
    """
    generates distance matrix from a list of neighboured atoms.

    Parameters:
    - neighbour_list:
        - list of neighboured atoms.
    - atom_number:
        - number of atoms in the analyzed molecule.
    
    Returns:
    - distance_matrix:
        - distance matrix of a molecule
        - contains all distances
    - bond_lists:
        - list of bond lists
        - every bond list contains bonds of one distance
            - bond_lists[0]: neighbour list
            - bond_lists[1]: angle list
            - bond_lists[2]: torsion/dihedral list
            - ...
    """
    distance_matrix = get_bond_matrix(neighbour_list, atom_number)
    # print(distance_matrix)
    bond_lists = [neighbour_list]
    for distance in range(2, atom_number):
        dlist = get_distance_list(bond_lists[-1], neighbour_list, distance + 1)
        if dlist.size > 0:
            bond_lists.append(dlist.astype(int))
            for d in dlist:
                i = int(d[0])
                j = int(d[-1])
                if distance_matrix[i][j] == 0 or distance_matrix[i][j] > distance:
                    distance_matrix[i][j] = distance
                if distance_matrix[j][i] == 0 or distance_matrix[j][i] > distance:
                    distance_matrix[j][i] = distance
        else:
            break
    # print(distance_matrix)
    return distance_matrix, bond_lists


def get_bond_matrix(neighbour_list, atom_number):
    """
    generates bond matrix from a list of neighboured atoms.
    redundant if distance matrix is available.

    Parameters:
    - neighbour_list:
        - list of neighboured atoms.
    - atom_number:
        - number of atoms in the analyzed molecule.
    Returns:
        bond_matrix:
            - bond matrix of a molecule
                - bond: 1
                - everything else: 0
    """
    bond_matrix = np.zeros((atom_number, atom_number))
    for b in neighbour_list:
        i = int(b[0])
        j = int(b[1])
        bond_matrix[i][j] = 1
        bond_matrix[j][i] = 1
    return bond_matrix

def plot_graph(
    atom_names,
    bond_list,
    saveto="",
    options={
        "node_size": 1200,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 5,
        "width": 3,
        "with_labels": True,
        "alpha": 0.3,
        "font_size": 12,
    } ):
    """
    plots graph...

    Parameters:
    - atom_names:
        - list or np.array of atom names
    - bond_list:
        - list of bonds
    - saveto:
        - optional, location to save pic of graph to,
        - add file format at the end ;)
    - options:
        - dict
        - networkx visualization options
            
    Returns:
        nothing
    """
    graph = nx.Graph()
    for b0, b1 in bond_list:
        graph.add_edge(b0, b1)

    labels = {}
    for i, name in enumerate(atom_names):
        labels[i] = name

    pos = nx.kamada_kawai_layout(graph)

    nx.draw_networkx(graph, pos, labels=labels, **options)
    nx.draw_networkx_labels(graph, pos, labels, font_size=12, font_color="black")

    # Set margins for the axes so that nodes aren't clipped
    ax = plt.gca()
    ax.margins(0.20)
    plt.axis("off")
    plt.tight_layout()
    if saveto:
        plt.savefig(saveto)
    plt.show()
    return


def get_bond_list_from_smiles(smiles):
    """
    gets a bond list from a smiles code using rd kit

    Parameters:
    - smiles:
        - string, smiles code

    Returns:
    - bond_list:
        - np.array of list of bonds
    """
    m = Chem.MolFromSmiles(smiles)
    bond_list = []
    for b in m.GetBonds():
        dummy = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond_list.append(dummy)
    bond_list = np.sort(np.array(bond_list))
    return bond_list, m


def graph_from_bonds(bond_list):
    """
    gets networkx graph object from a bond list

    Parameters:
    - bond_list:
        - np.array of list of bonds

    Returns:
    - graph:
        - networkx graph object
    """
    graph = nx.Graph()
    for b0, b1 in bond_list:
        graph.add_edge(b0, b1)

    return graph


def visualize_graph(
    graph,
    names,
    options={
        "node_size": 1200,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 5,
        "width": 3,
        "with_labels": True,
        "alpha": 0.3,
        "font_size": 12,
    },
):
    """
    visualizes graph...

    Parameters:
    - graph:
        - networkx graph object
    - names:
        - list or np.array of atom names
    - options:
        - dict
        - networkx visualization options
    
    Returns:
        nothing
    """

    labels = {}

    for i, name in zip(graph.nodes, names):
        labels[i] = name

    pos = nx.kamada_kawai_layout(graph)

    nx.draw_networkx(graph, pos, labels=labels, **options)
    nx.draw_networkx_labels(graph, pos, labels, font_size=12, font_color="black")
    return


def longest_simple_paths(graph, source, target):
    """
    uses networkx to get longest simple paths
    i.e. paths without double paths
    between source and target

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
    - target:
        - int, atom number to end at
            
    Returns:
    - list of longest paths

    """
    longest_paths = []
    longest_path_length = 0
    for path in nx.all_simple_paths(graph, source=source, target=target):
        if len(path) > longest_path_length:
            longest_path_length = len(path)
            longest_paths.clear()
            longest_paths.append(path)
        elif len(path) == longest_path_length:
            longest_paths.append(path)
    return longest_paths


def get_longest_path(graph, source=0):
    """
    uses networkx to get longest simple path
    starting from source

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
            
    Returns:
     - np.array, path (first that was found)

    """
    longest_path = []
    longest_path_end = 0
    for iii in graph.nodes:
        dummy = np.atleast_2d(longest_simple_paths(graph, source=source, target=iii))[0]
        if len(dummy) > len(longest_path):
            longest_path = dummy
            longest_path_end = iii
    return np.array(longest_path)


def get_shortest_paths(graph, source, target):
    """
    uses networkx to get shortest paths
    between source and target

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
    - target:
        - int, atom number to end at
    
    Returns:
    - np.array, list of shortest paths

    """
    shortest_paths = []
    shortest_path_length = 1e6
    for path in nx.shortest_simple_paths(graph, source=source, target=target):
        if len(path) < shortest_path_length:
            shortest_path_length = len(path)
            shortest_paths.clear()
            shortest_paths.append(path)
        elif len(path) == shortest_path_length:
            shortest_paths.append(path)
    return shortest_paths


def get_shortest_path(graph, source, target):
    """
    uses networkx to get shortest paths
    between source and target

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
    - target:
        - int, atom number to end at
            
    Returns:
    - np.array, shortest path (first that was found)

    """
    shortest_path = []
    shortest_path_length = 1e6
    for path in nx.shortest_simple_paths(graph, source=source, target=target):
        if len(path) < shortest_path_length:
            shortest_path_length = len(path)
            shortest_path = path.copy()
    return np.array(shortest_path)


def get_shortest_nontrivial_path(graph, source, target):
    """
    uses networkx to get shortest nontrivial paths
    i.e. longer than one bond between source and target
    useful to detect ring structures

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
    - target:
        - int, atom number to end at
            
    Returns:
    - np.array, shortest nontrivial path (first that is found)

    """
    shortest_path = []
    shortest_path_length = 1e6
    for path in nx.shortest_simple_paths(graph, source=source, target=target):
        if len(path) < shortest_path_length and len(path) > 2:
            shortest_path_length = len(path)
            shortest_path = path.copy()
    return np.array(shortest_path)


def get_root(graph, source, range_back):
    """
    gets root of a substructure i.e. branch

    Parameters:
    - graph:
        - networkx graph object
    - source:
        - int, atom number to start from
    - range_back:
        - int, atom number to go back to find root
            
    Returns:
    - int, atom number to of root

    """
    nodes = np.array(graph.nodes)
    for target in np.sort(nodes[nodes < source]):
        shortest_path = get_shortest_path(graph, source, target)
        if len(shortest_path) == range_back:
            return target
    return -1


def bond_list_from_simple_path(path):
    """
    builds a bond list from a simple path
    i.e. path without branches or rings

    Parameters:
    - path:
        - np.array, simple path
            
    Returns:
    - np.array, bond list

    """
    return np.sort(np.array([path[:-1], path[1:]]).T)


def get_diff_in_bond_lists(bond_list_1, bond_list_2):
    """
    gets difference between two bond lists

    Parameters:
    - bond_list_1:
        - np.array, bond list
    - bond_list_2:
        - np.array, bond list
            
    Returns:
    - np.array, bond list containing diffs
    
    """
    dd = np.concatenate([bond_list_1, bond_list_2])
    unique_bonds, count = np.unique(dd, axis=0, return_counts=True)
    unique_bonds = unique_bonds[np.where(count == 1)]
    return np.sort(unique_bonds)


def get_next_index(main_path, remaining_bonds):
    """
    gets index of next substructure

    Parameters:
    - main_path:
        - networkx graph object
    - remaining_bonds:
        - np.array, bond list
            
    Returns:
    - int, atom number of next index
    """
    connections = np.intersect1d(main_path, np.unique(remaining_bonds))
    ii = 1e10
    for c in connections:
        p = np.squeeze(np.where(main_path == c))
        if p < ii:
            ii = main_path[p]
    return ii


def get_graphstring(bond_list, names, source=0):
    """
    generates a graphstring from a bond list and atom names
    uses the longest path from source to an end as main path

    Parameters:
    - bond_list:
        - np.array, bond list
    - names:
        - np.array, atom names
    - source:
        - int, start of main path, default=0
            
    Returns:
    - str, graphstring to use with moleculegraph
    """
    graph = graph_from_bonds(bond_list)
    main_path = get_longest_path(graph, source=source)
    return get_graphstring_set_main(graph, bond_list, names, main_path)


def get_graphstring_set_main(graph, bond_list, names, main_path):
    """
    generates a graphstring from a bond list and atom names
    uses the longest path from source to an end as main path

    Parameters:
    - graph:
        - networkx graph object
    - bond_list:
        - np.array, bond list
    - names:
        - np.array, atom names
    - main_path:
        - np.array, main path to build graph from
            
    Returns:
    - str, graphstring to use with moleculegraph
    """
    funs = np.zeros(main_path.shape)
    fun_ranges = np.zeros(main_path.shape)

    main_path_bond_list = bond_list_from_simple_path(main_path)
    remaining_bonds = get_diff_in_bond_lists(bond_list, main_path_bond_list)

    while True:
        subgraph = graph_from_bonds(remaining_bonds)

        idx = get_next_index(main_path, remaining_bonds)
        # print("iii",p, idx)

        subpath = get_longest_path(subgraph, source=idx)
        subpath_bond_list = bond_list_from_simple_path(subpath)
        match = np.intersect1d(main_path, subpath)

        if len(match) == 2 and len(subpath) == 2:
            print("ring")
            subpath = get_shortest_nontrivial_path(graph, match[0], match[1])
            # print("sub",subpath)
            i = np.squeeze(np.where(main_path == match[0]))
            j = np.squeeze(np.where(main_path == match[1]))
            # print("match ij",match, i, j)
            iinsert = np.max((i, j)) + 1
            # print(iinsert)

            main_path = np.insert(main_path, iinsert, [-1])
            fun_ranges = np.insert(fun_ranges, iinsert, [len(subpath)])
            funs = np.insert(funs, iinsert, [-1])

        elif len(match) == 1:
            print("branch")
            i = np.squeeze(match)
            # print(i,subpath)
            if subpath[0] != i:
                subpath = subpath[::-1]
            subpath = subpath[1:]
            # print(subpath)
            subfuns = np.concatenate([[1], np.zeros(subpath.shape)])
            subfun_ranges = np.concatenate([[len(subpath)], np.zeros(subpath.shape)])
            subpath = np.concatenate([[-1], subpath])

            iinsert = np.squeeze(np.where(main_path == i)) + 1

            main_path = np.insert(main_path, iinsert, subpath)
            fun_ranges = np.insert(fun_ranges, iinsert, subfun_ranges)
            funs = np.insert(funs, iinsert, subfuns)

        else:
            print("ERROR")
            return None

        main_path_bond_list = np.concatenate([main_path_bond_list, subpath_bond_list])
        remaining_bonds = get_diff_in_bond_lists(bond_list, main_path_bond_list)
        if len(remaining_bonds) == 0:
            break

    molecule_list = []
    for i, (fun, rang) in zip(main_path, zip(funs, fun_ranges)):
        if fun == 0:
            molecule_list.append(names[i])
        elif fun == 1:
            molecule_list.append("b" + str(int(rang)))
        elif fun == -1:
            molecule_list.append("r" + str(int(rang)))

    molstring = "[" + "][".join(molecule_list) + "]"

    return molstring


def reenumerate_branches(branches):
    """
    reenumerates branches due to occurence

    Parameters:
    - branches:
        - list with branch numbers
    
    Returns:
    - branches:
        - rearranged list with branch numbers
    """
    dummy = branches.copy()
    _, uniq_idx = np.unique(dummy, return_index=True)
    for bno, xx in enumerate(np.sort(uniq_idx)):
        branches[dummy == dummy[xx]] = bno
    return branches


def get_branch_root(i, branches):
    """
    gets root of a branch

    Parameters:
    - i: int
        - number of branch to get root for
    - branches: np.array
        - array with branch numbers
            
    Returns: 
    - int number of root branch
    """
    if i == 0:
        return 0
    p = np.atleast_1d(np.squeeze(np.where(branches == i)))
    for ii in np.unique(branches)[::-1]:
        pp = np.atleast_1d(np.squeeze(np.where(branches == ii)))
        if pp[0] < p[0] and pp[-1] > p[-1]:
            return pp[np.max(np.atleast_1d(np.squeeze(np.where(pp < p[0]))))]

    print(
        "WARNING... root not found. Set to zero (check visualization of ur molecule pls)"
    )
    return None


def make_graph(stringlist):
    """
    builds graph string from string list
    
    Parameters:
    - stringlist:
        - list of strings containing names and
        - add file format at the end ;)
    
    Returns:
        nothing
    """
    return "[" + "][".join(stringlist) + "]"
