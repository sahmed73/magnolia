# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Oct  6 20:23:08 2023
"""

import numpy as np
import pandas as pd
import networkx as nx
from functools import wraps
import time

#%% runtime wrapper
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
#%% get_steps
@function_runtime
def get_steps(bondfilepath):
    steps = []
    with open(bondfilepath) as bf:
        for line in bf:
            if 'Timestep' in line:
                step = int(line.strip().split()[-1]) 
                steps.append(step)
    return np.array(steps)
#%% get_atomConnectivity
@function_runtime
def get_atomConnectivity(bondfilepath,timestep,**kwargs):
    connectivity = {}
    
    with open(bondfilepath) as bf:
        breakFlag = False
        
        for line in bf:
            if 'Timestep' in line:
                step = int(line.strip().split()[-1])
                # break if required data collected for the given timestep
                if step!=timestep and breakFlag:
                    break
            
            if step==timestep:
                breakFlag=True
                splitted = line.strip().split()
                if splitted[0].isnumeric():
                    parent  = int(splitted[0])
                    nb      = int(splitted[2]) #number of bond
                    # all connected atom
                    children= np.array((splitted[3:3+nb])).astype(int)
                    connectivity[parent]=children  
    return connectivity
#%% get_atomtypes
def get_atomtypes(bondfilepath,timestep,**kwargs):
    atypes = {}
    with open(bondfilepath) as bf:
        for line in bf:
            if 'Timestep' in line:
                step = int(line.strip().split()[-1])
            
            if step==timestep:
                splitted = line.strip().split()
                if splitted[0].isnumeric():
                    parent  = int(splitted[0])
                    typ     = int(splitted[1]) # atom type
                    atypes[parent]=typ  
    return atypes
#%% get_moleculeID
def get_moleculeID(bondfilepath,timestep,**kwargs):
    molids = {}
    with open(bondfilepath) as bf:
        for line in bf:
            if 'Timestep' in line:
                step = int(line.strip().split()[-1])
            
            if step==timestep:
                splitted = line.strip().split()
                if splitted[0].isnumeric():
                    parent  = int(splitted[0])
                    nb      = int(splitted[2]) # number of bonds
                    molid   = int(splitted[nb+3]) # molecule id
                    molids[parent]=molid
    return molids
#%% get_atomCharge
def get_atomCharge(bondfilepath,timestep,**kwargs):
    charges = {}
    with open(bondfilepath) as bf:
        for line in bf:
            if 'Timestep' in line:
                step = int(line.strip().split()[-1])
            
            if step==timestep:
                splitted = line.strip().split()
                if splitted[0].isnumeric():
                    parent  = int(splitted[0])
                    charge  = float(splitted[-1])
                    charges[parent]=charge 
    return charges
#%% get_molecules
def get_molecules(atomConnectivity):
    graph     = nx.Graph(atomConnectivity)
    molecules = nx.connected_components(graph)
    return molecules
#%% get_cycles        
def get_cycles(atomConnectivity):
    graph    = nx.Graph(atomConnectivity)
    cycles   = nx.cycle_basis(graph)
    return cycles
#%% parse_bondfile
@function_runtime
def parse_bondfile(bondfilepath,Nevery=1):
    neighbours = {}
    steps    = get_steps(bondfilepath)[::Nevery]
    for frame,step in enumerate(steps):
        print("Frame",frame,":",step)
        neighbours[step]=get_atomConnectivity(bondfilepath, step)
    return neighbours
#%% --main--
bondfilepath=r"C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Antioxidants\ABCDE\A\5A_30_O2\Equilibrate\Sim-1"+"\\bonds.reaxc"

steps = get_steps(bondfilepath)
n_frame   = len(steps)
frame     = -1
timestep  = steps[frame]
print('Number of frame:',n_frame)
print(timestep)
# connectivity = get_atomConnectivity(bondfilepath, timestep)
# atypes       = get_atomtypes(bondfilepath, timestep)
# molecules    = get_molecules(connectivity)
# cycles       = get_cycles(connectivity)
# molids       = get_moleculeID(bondfilepath, timestep)
# atomCharges  = get_atomCharge(bondfilepath, timestep)
# print(connectivity)
# print(atypes)
# print(sum(1 for _ in molecules))
# print(cycles)
# print(molids)
# print(atomCharges)
# neighbours = parse_bondfile(bondfilepath)
# print(neighbours)
#%%
import magnolia.bondfile_parser as bfp

bonddata = bfp.parsebondfile(bondfilepath)
neighbours = bonddata['neighbours']

