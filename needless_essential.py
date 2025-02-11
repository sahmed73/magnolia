# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 06:59:23 2023

@author: arup2
"""
import numpy as np  # Import NumPy

def print_runtime(sec): #convert second into M minutes S seconds
    minutes   = int(sec/60)
    sec = round(sec%60,1)
    if minutes==0:
        print('Runtime: {} sec'.format(sec))
    else:
        print('Runtime: {} min {} sec'.format(minutes, sec))
   

def randstr(low=100000,high=999999):
    import random
    rand = random.randint(low, high)
    return str(rand)

def getlab(key):
    key = key.lower()
    if key in ['pe','poteng','potential_energy']:
        result = 'Potential Energy (Kcal/mol)'
    elif key in ['ke','kineng']:
        result = 'Kinetic Energy (Kcal/mol)'
    elif key in ['eng','toteng','energy']:
        result = 'Total Energy (Kcal/mol)'
    elif key in ['density']:
        result = 'Density (g/cm$^3$)'
    elif key in ['time']:
        result = 'Time (ps)'
    elif key in ['temp','temperature']:
        result = 'Temperature (K)'
    elif key in ['volume','vol']:
        result = 'Volume ($\AA^3$)'
    elif key in ['press','pressure']:
        result = 'Pressure (atm)'
    else:
        print('No Key Found')
        result = key
    
    return result

def atomids2expression(atoms, key='ParticleIdentifier'): # atom ids to ovito expression selection
    
    def inner(atomss):
        expression = ''
        for atom in atomss:
            expression+=f'|| {key}=='+str(atom)
        
        return expression

    if isinstance(atoms[0], int): # [374, 375, ....]
        result = inner(atoms)
    elif isinstance(atoms[0], str): # ['374', '375', ....]
        result = inner(atoms)
    elif isinstance(atoms[0], tuple): # [(374, 375), (.., ..), ....]
        aa = []
        for a in atoms:
            aa.extend(a)
        result = inner(aa)
    
    return result


def atomids2expression_v2(atoms):  # atom ids to ovito expression selection
    result = ''  # Initialize result to handle cases where atoms might not match any condition
    
    def inner(atomss):
        expression = ''
        for index, atom in enumerate(atomss):  # Use enumerate for efficiency and correctness
            if index == len(atomss) - 1:
                expression += f'ParticleIdentifier=={atom}'
            else:
                expression += f'ParticleIdentifier=={atom} || '
        return expression

    if isinstance(atoms, np.ndarray):  # Check if input is a NumPy array
        if atoms.ndim == 1:
            result = inner(atoms)
        elif atoms.ndim == 2:
            # Flatten the array if it's two-dimensional
            flattened_atoms = atoms.flatten()
            result = inner(flattened_atoms)
    elif isinstance(atoms, list) and atoms:  # Check if it's a non-empty list
        if isinstance(atoms[0], (int, str)):  # Check for int or str directly
            result = inner(atoms)
        elif isinstance(atoms[0], tuple):
            aa = [atom for a in atoms for atom in a]  # Flatten list of tuples
            result = inner(aa)
        elif isinstance(atoms[0], list):  # Handle list of lists
            flattened_atoms = [atom for sublist in atoms for atom in sublist]
            result = inner(flattened_atoms)
    # Optionally, handle unexpected or empty input more explicitly
    # else:
    #     result = ''  # or raise an exception
    
    return result



        
        