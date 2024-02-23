# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 06:59:23 2023

@author: arup2
"""


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
    elif key in ['press','Press','pressure','Pressure']:
        result = 'Pressure (atm)'
    else:
        print('No Key Found')
    
    return result

def atomids2expression(atoms): # atom ids to ovito expression selection
    
    def inner(atomss):
        expression = ''
        for atom in atomss:
            if atomss.index(atom)==len(atomss)-1:
                expression+='ParticleIdentifier=='+str(atom)
            else:
                expression+='ParticleIdentifier=='+str(atom)+'|| '
        
        return expression

    if isinstance(atoms[0], int):
        result = inner(atoms)
    elif isinstance(atoms[0], str):
        result = inner(atoms)
    elif isinstance(atoms[0], tuple):
        aa = []
        for a in atoms:
            aa.extend(a)
        result = inner(aa)
    
    return result

        
        