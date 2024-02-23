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
def parsedumpfile(dumpfile,**kwargs):
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
        fso      = kwargs.get('fso',False)
        fs_flag  = 0
        
        get_timestep    = False
        get_coordinates = False
        
        dumpdata   = {}
        atomtypes  = {} # need only one timestep
        position = {}   # dynamic position (Step-->atomid-->position)
        for line in df:
            splitted = line.split()            
            if splitted[-1]=='TIMESTEP':
                get_timestep    = True
                get_coordinates = False
                fs_flag +=1
                continue
            
            if fso and fs_flag>1:
                break
                
            if get_timestep:
                step = int(splitted[0])
                position[step] = {}
                get_timestep = False
                continue
                
            if line.find('ITEM: ATOMS')!=-1:
                get_coordinates = True
                continue
            
            if get_coordinates:
                atomid,atype,*coordinates = splitted
                atomid   = int(atomid)
                atype = int(atype)
                coordinates = tuple(map(float,coordinates))
                position[step][atomid] = coordinates
                atomtypes[atomid] = atype
        
        dumpdata['position'] = position  #dictionary of dictionary 
        dumpdata['atypes']   = atomtypes #dictionary (atomid-->atomtypes)
    
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