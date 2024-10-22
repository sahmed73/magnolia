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
def parsedumpfile(dumpfile,**kwargs):
    with open(dumpfile) as df:
        # We are parsing text file line  by line.
        # So it is not possible to get the next line.
        # That's why we are declaring some flags. Flags will enable us to 
        # collect a data that is flagged to be collected at earlier itteration.
        # Output data:
            # dict of pd.DataFrames()
                
        ###-get kwargs-####
        # here.............
        ##
        '''
            id: int
            mol: int
            type: int
            x,y,z : float
            now all are 'object', need to change it to respective type
        '''
        
        state    = None # timestep, atoms
        step     = None
        
        dumpdata   = {}
        atomtypes  = {}
        position = {}
        ATOM = []
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
                if step is not None:
                    dumpdata[step] = pd.DataFrame(ATOM, columns=column_heads)
                ATOM = []
                step = int(splitted[0])
                state = None
                
            if state=='atoms':
                 ATOM.append(splitted.copy())
    
    dumpdata[step] = pd.DataFrame(ATOM, columns=column_heads)
    return dumpdata


                
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

