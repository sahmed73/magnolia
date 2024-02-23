# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 17 19:29:26 2023
"""

#### Module for analyze atom_connectivity using graph theory
## List of functions ##
# 1. longest_chain
# 2. find_corresponding_node
# 3. draw_rdkit2D

import networkx as nx
from rdkit import Chem
import matplotlib.pyplot as plt
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io
import warnings
from collections import deque

def longest_chain(molecule_adjList):
    """
    Calculate the longest chain in a molecule represented as an adjacency list.
    The molecule can be represented either as a NetworkX graph or
    as a dictionary. Assumes that the graph is a tree (i.e., no cycles).

    Parameters:
    molecule_adjList (nx.Graph or dict): The adjacency list representation of
                                         the molecule.

    Returns:
    tuple: A tuple containing a list of nodes that are the starting points of 
           the longest chain(s), and an integer representing the length of the
           longest chain.
    """
    
    # Check the type of molecule_adjList and convert to a uniform format (dict)
    if isinstance(molecule_adjList, nx.Graph):
        temp_dict      = dict(molecule_adjList.adjacency()) # dict of dict
        adjacency_list = {x: list(y.keys()) for x, y in temp_dict.items()}
    elif isinstance(molecule_adjList, dict):
        adjacency_list = molecule_adjList.copy()
    else:
        raise TypeError("molecule_adjList must be either a nx.Graph or a dict")
        
    # Helper function to perform BFS and find the max depth
    def bfs_max_depth(source):
        visited   = []
        queue     = [(source,1)] # node and depth (depth started from 1)
        max_depth = 1
        while queue:
            parent, depth = queue.pop()
            children      = adjacency_list[parent]
            max_depth     = max(max_depth,depth)
            visited.append(parent)
            for child in children:
                if child not in visited:
                    queue.append((child,depth+1))
        return max_depth
    
    # identify the leaf (a node that does not have any children)
    leaves  = []
    for parent, children in adjacency_list.items():
        if len(children)==1:
            leaves.append(parent)
    
    # Compute the length of the longest chain starting from each leaf
    longest = 0
    source_leaves = []
    for leaf in leaves:
        current_length = bfs_max_depth(leaf)
        if current_length > longest:
            longest = current_length
            source_leaves = [leaf]
        elif current_length == longest:
            source_leaves.append(leaf)
     
    return source_leaves,longest

def number_nodes_bfs(graph, start_node):
    """
    Number nodes in a graph using breadth-first search.

    :param graph: NetworkX graph.
    :param start_node: The starting node for BFS.
    :return: Dictionary of node numbers.
    """
    node_number = {}
    current_number = 1
    queue = deque([start_node])
    visited = set()

    while queue:
        node = queue.popleft()
        if node in visited:
            continue
        visited.add(node)

        node_number[node] = current_number
        current_number += 1

        for neighbor in graph.neighbors(node):
            if neighbor not in visited:
                queue.append(neighbor)

    return node_number

def number_nodes_dfs(graph, start, visited=None, counter=1):
    if visited is None:
        visited = {}  # To keep track of visited nodes and their numbers
    visited[start] = counter
    # print(f'Carbon {start} is numbered: {counter}')
    for neighbor in graph[start]:
        if neighbor not in visited:
            counter += 1
            number_nodes_dfs(graph, neighbor, visited, counter)
    return visited



def find_corresponding_node(graph1, graph2, node_in_graph1):
    """
    Finds the corresponding node in graph2 for a given node in graph1 based on
    an isomorphism. Accepts input as either NetworkX graph objects or
    dictionary representations of graphs.
    Raises a ValueError if the graphs are not isomorphic.

    Parameters:
    graph1 (nx.Graph or dict): The first graph or its dict representation.
    graph2 (nx.Graph or dict): The second graph or its dict representation.
    selected_node: The node in graph1 for which to find the corresponding
                   node in graph2.

    Returns:
    The corresponding node in graph2 if an isomorphism exists, otherwise raises ValueError.
    """
    # Convert dictionaries to NetworkX graphs if necessary
    if isinstance(graph1, dict):
        graph1 = nx.Graph(graph1)
    if isinstance(graph2, dict):
        graph2 = nx.Graph(graph2)

    # Check for isomorphism and find the mapping
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph1, graph2)

    if matcher.is_isomorphic():
        iso_mapping = next(matcher.isomorphisms_iter())
        return iso_mapping.get(node_in_graph1)
    else:
        raise ValueError("The graphs are not isomorphic.")

def draw_rdkit2D_from_smiles(smiles, highlight_atoms=None, ax=None):
    # Create RDKit mol object
    rdkit_mol = Chem.MolFromSmiles(smiles)

    # Set up drawing
    drawer = rdMolDraw2D.MolDraw2DCairo(1200, 1200)

    # Validate highlight_atoms indices and filter out invalid indices
    valid_highlight_atoms = []
    if highlight_atoms is not None:
        num_atoms = rdkit_mol.GetNumAtoms()
        for atom in highlight_atoms:
            if atom < num_atoms:
                valid_highlight_atoms.append(atom)
            else:
                warnings.warn(f"Atom index {atom} in highlight_atoms exceeds the number of atoms in the molecule ({num_atoms}). Skipping this atom.")

    # Highlight specified atoms (only valid ones)
    if valid_highlight_atoms:
        highlight_colors = {atom: (0, 1, 0) for atom in valid_highlight_atoms}  # RGB for green
        drawer.DrawMolecule(rdkit_mol, highlightAtoms=valid_highlight_atoms, highlightAtomColors=highlight_colors)
    else:
        drawer.DrawMolecule(rdkit_mol)

    drawer.FinishDrawing()

    # Convert the drawing to an image
    image_data = drawer.GetDrawingText()
    image = Image.open(io.BytesIO(image_data))

    # Display the image in the provided ax or create a new one
    if ax is None:
        fig, ax = plt.subplots()
    ax.imshow(image)
    ax.axis('off')  # Hide the axis

def create_mol_from_graph(graph, atypes, atomsymbols):
    mol = Chem.RWMol()

    # Add atoms
    atom_index_map = {}
    for node_id in graph.nodes():
        atype = atypes.get(node_id, 1)  # Default to the first atom type if not specified
        atom_symbol = atomsymbols[atype - 1]  # Adjust index for zero-based Python lists
        atom = Chem.Atom(atom_symbol)
        atom_idx = mol.AddAtom(atom)
        atom_index_map[node_id] = atom_idx

    # Add bonds
    for start, end in graph.edges():
        if start in atom_index_map and end in atom_index_map:
            mol.AddBond(atom_index_map[start], atom_index_map[end], Chem.rdchem.BondType.SINGLE)

    Chem.SanitizeMol(mol)
    return mol

# Draw the molecule
def draw_rdkit2D_from_graph(graph, atypes, atomsymbols, highlight_atoms=None, ax=None):
    rdkit_mol = create_mol_from_graph(graph, atypes, atomsymbols)

    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    if highlight_atoms is not None:
        highlight_colors = {atom: (0, 1, 0) for atom in highlight_atoms}
        drawer.DrawMolecule(rdkit_mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
    else:
        drawer.DrawMolecule(rdkit_mol)
    drawer.FinishDrawing()

    image_data = drawer.GetDrawingText()
    image = Image.open(io.BytesIO(image_data))
    if ax is None:
        fig, ax = plt.subplots()
    ax.imshow(image)
    ax.axis('off')

def draw_molecule_asGraph(G,atomsymbols,elabel=True):
    ## Input: Graph with node attr 'atom_type', edge attr 'bond_order'
    ##              get the atomConnectivity as graph using 
    ##              bfp.parsebondfile_asGraph function
    ## Draw graph with 'CHO' node and bond_order as edge label
    
    # Color mapping
    atom_type_color = {
        1: 'white',  # atom_type 1, white with black edge
        2: 'gray',   # atom_type 2, gray with black edge
        3: 'red'     # atom_type 3, red with black edge
    }
    
    # Apply the color mapping to each node
    node_colors = [atom_type_color[attr['atom_type']] for node, attr in G.nodes(data=True)]
    
    # Node labels: {node_id: atom_type}
    node_labels = {node: atomsymbols[attr['atom_type']-1] for node, attr in G.nodes(data=True)}
    
    
    # Position nodes using a layout
    pos = nx.spring_layout(G)
    
    # Draw the graph
    nx.draw(G, pos, with_labels=True, labels=node_labels,
            node_color=node_colors, edgecolors='black',node_size=500)
    
    if elabel:
        # Edge labels: {(u, v): bond_order}
        edge_labels = {(u, v): attr['bond_order'] for u, v, attr in G.edges(data=True)}
        # Draw edge labels (bond orders)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
            
def get_moleculeGraph(molAdjList,atomtypes,bondorders):
    pass