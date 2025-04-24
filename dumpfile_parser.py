# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:32:10 2023

@author: Shihab
Dump File Parser
"""
import math
from functools import wraps
import time
import sys
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN


def function_runtime(f):
    @wraps(f)
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = f(*args,**kwargs)
        stop_time = time.time()
        runtime = stop_time-start_time
        minute = int(runtime/60)
        second = runtime%60
        if minute==0:
            msg = "Execution time of {}: {:0.1f} sec".format(f.__name__,second)
        else:
            msg = "Execution time of {}: {} min {:0.1f} sec".format(f.__name__,minute,second)
        print(msg)
        return result
    return wrapper

@function_runtime
def parsedumpfile_OLD(dumpfile,**kwargs):
    with open(dumpfile) as df:
        # We are parsing text file line  by line.
        # So it is not possible to get the next line.
        # That's why we are declaring some flags. Flags will enable us to 
        # collect a data that is flagged to be collected at earlier itteration.
        # Output: 
            # a dictionary of dictionary:
                # first key  = step
                # second key = atom id
                # value      = [atomtype, position = (x,y,z)]
                
        ###-get kwargs-####
        # here.............
        
        state    = None # timestep, atoms
    
        dumpdata   = {}
        atomtypes  = {}
        position = {}
        for line in df:
            splitted = line.split()  
            
            # flagging
            if splitted[-1]=='TIMESTEP':
                state = 'timestep'
                continue
            if line.startswith('ITEM: ATOMS'):
                column_heads = splitted[2:]
                state = 'atoms'
                continue
                
            # get values
            if state=='timestep':
                step = int(splitted[0])
                state = None
                dumpdata[step] = {}
                
            if state=='atoms':
                atomid = int(splitted[0])
                dumpdata[step][atomid] = splitted[1:]
            
    
    return dumpdata

@function_runtime
def parsedumpfile(dumpfile, Nfreq=1, **kwargs):
    """
    Parses a LAMMPS dump file and extracts atomic data for selected timesteps.

    This function reads the file line by line and constructs a dictionary of
    DataFrames, each representing atomic data for a specific timestep. Only
    every `Nfreq`-th timestep is stored, based on the given sampling frequency.

    Parameters
    ----------
    dumpfile : str
        Path to the LAMMPS dump file to parse.
    Nfreq : int, optional
        Sampling frequency. Only every Nfreq-th timestep is parsed. Default is 1.
    **kwargs : dict
        Placeholder for future keyword arguments (currently unused).

    Returns
    -------
    dumpdata : dict of pandas.DataFrame
        Dictionary where keys are timestep integers and values are DataFrames
        containing atom-level data for that timestep.

    Notes
    -----
    - Columns like 'id', 'mol', and 'type' are cast to `int`.
    - Positions (`x`, `y`, `z`) and other numerical values are cast to `float`.
    - If a value cannot be cast to `int` or `float`, it remains as a string.
    - Column headers are taken from the `ITEM: ATOMS` section of the dump file.
    """

    with open(dumpfile) as df:
        state    = None # timestep, atoms, skip
        step     = None
        step_count = 0 # for Nfreq
        
        dumpdata   = {}
        ATOM = []
        for line in df:
            splitted = line.split()  
            
            # flagging
            if splitted[-1]=='TIMESTEP':
                step_count+=1
                state = 'timestep'
                continue
            
            if line.startswith('ITEM: ATOMS'):
                column_heads = splitted[2:]
                for i in range(len(column_heads)):
                    if column_heads[i]=='xu': column_heads[i]='x'
                    if column_heads[i]=='yu': column_heads[i]='y'
                    if column_heads[i]=='zu': column_heads[i]='z'
                state = 'atoms'
                continue
                
            # get values
            if state=='timestep':
                if step is not None and (step_count-1)%Nfreq==0:
                    dumpdata[step] = pd.DataFrame(ATOM, columns=column_heads)
                ATOM = []
                step = int(splitted[0])
                state = None
                
            if state=='atoms':
                for i in range(len(splitted)):
                    try:
                        if i < 3:
                            try:
                                splitted[i] = int(splitted[i])
                            except ValueError:
                                splitted[i] = float(splitted[i])
                        else:
                            splitted[i] = float(splitted[i])
                    except ValueError:
                        print('--')
                        pass  # leave as original string


                ATOM.append(splitted.copy())
    if (step_count-1)%Nfreq==0:
        dumpdata[step] = pd.DataFrame(ATOM, columns=column_heads)
    return dumpdata

def add_element_symbols(dumpdata, datafile):
    import magnolia.lammps_datafile_parser as ldfp
    masses, atoms, bonds = ldfp.parse_lammps_data(datafile)
    
    for step, snapshot in dumpdata.items():
        snapshot = snapshot.merge(masses[['type', 'symbol']], on='type', how='left')
        dumpdata[step] = snapshot
    
    return dumpdata

@function_runtime
def creat_graph_network(dumpdata_with_symbols, bond_tolerance=0.1):
    network = {}
    covalent_radii = {'H': 0.31, 'C': 0.76, 'O': 0.66}

    for step, snapshot in dumpdata_with_symbols.items():
        coords = snapshot[['x', 'y', 'z']].values
        symbols = snapshot['symbol'].values
        ids = snapshot['id'].values

        radii = np.array([covalent_radii[s] for s in symbols])
        tree = cKDTree(coords)

        # Step 1: broad query
        candidate_pairs = tree.query_pairs(r=(radii.max() * 2 + bond_tolerance))

        # Step 2: precise filtering by actual distance and pairwise cutoff
        bonds = []
        for i, j in candidate_pairs:
            dist = np.linalg.norm(coords[i] - coords[j])
            cutoff = radii[i] + radii[j] + bond_tolerance
            if dist <= cutoff:
                # Use atom IDs, not just indices
                bonds.append((ids[i], ids[j]))

        G = nx.Graph()
        G.add_nodes_from(ids)     # Make sure all atoms are in the graph
        G.add_edges_from(bonds)   # Add valid bonds
        network[step] = G

    return network

def cluster_molecules(coords, eps=2.0, min_samples=1):
    """
    Cluster atoms into molecules based on spatial proximity using DBSCAN.

    Parameters
    ----------
    coords : ndarray of shape (n_atoms, 3)
        Cartesian coordinates of atoms.
    eps : float, optional
        Maximum distance between two atoms to be considered in the same cluster.
        Should be similar to a bond distance threshold (e.g., ~1.5–2.5 Å).
    min_samples : int, optional
        Minimum number of atoms to form a cluster. Default is 1 for molecules.

    Returns
    -------
    labels : ndarray of shape (n_atoms,)
        Cluster label for each atom. Same order as input `coords`.
        Each unique label corresponds to a detected molecule.
    n_clusters : int
        Total number of clusters (molecules) detected.
    """
    
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    labels = clustering.labels_

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    return labels, n_clusters



def compute_msd(positions):
    """
    positions: numpy array of shape (N, 3) where N is the number of time steps
    returns: 1D array of MSD values for each time lag τ
    """
    N = len(positions)
    msd = np.zeros(N)

    for tau in range(1, N):
        diffs = positions[tau:] - positions[:-tau]
        squared_displacements = np.sum(diffs**2, axis=1)
        msd[tau] = np.mean(squared_displacements)

    return msd

def distance_tracker(dumpdata,atom1,atom2):  
    # Calculate distance between two given atoms atom1 and atom2 in each step
    # 
    distances = {}     
    for step in dumpdata.keys():
        coord_atom1 = dumpdata[step][atom1][1] # index 1 contain position
        coord_atom2 = dumpdata[step][atom2][1] # index 1 contain position
         
        # Calculate distance between two points using math module
        # New in version 3.8
        distances[step]=math.dist(coord_atom1, coord_atom2)
        
    return distances
