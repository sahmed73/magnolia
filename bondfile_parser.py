# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 23:56:22 2023 (Version 1.0)
Version 2.0 released on June 1, 2023
last modified: 12/10/2023

@author: shihab

This is my own LAMMPS post-processing library
"""

import networkx as nx
from networkx.algorithms import isomorphism
from collections import Counter
import re
import sys
import math
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import random
import os
import pickle
import json
from functools import wraps
import numpy as np
from scipy.optimize import curve_fit, newton
from rdkit import Chem
from sklearn.cluster import KMeans
from collections.abc import Iterable

# =============================================================================
## Dacorator functions:
 #   1. function_runtime 
 #   2. loadpickle_or_execute
 #   3. loadjson_or_execute
## Module functions:
 #   4. get_neighbours
 #   5. parsebondfile
 #   6. parsebondfile_asGraph
 #   7. merge_bondfiles
 #   8. merge_bondfiles_OLD
 #   9. text_isequal
 #   10. get_molecules
 #   11. get_molecular_formula
 #   12. get_SpeciesCountAtEveryTimestep
 #   13. stepwise_species_count
 #   14. stepwise_species_count_v2
 #   15. expression_selector
 #   16. pathway_tracker
 #   17. step2picosecond
 #   18. get_nearestindex
 #   19. sort_molecular_formula
 #   20. make_molecular_formula_latex
 #   21. compute_molecular_weight
 #   22. plot_species_heatmap
 #   23. plot_species_heatmap_v2
 #   24. atomic_weight
 #   25. get_speciesVStemp
 #   26. get_onset
 #   27. bondorder_evolution
 #   28. get_nbondsVStime
 #   29. species_to_molecule
 #   30. cumulative_nspecies
 #   31. onset_plot
 #   32. get_species_count
 #   33. count_functinalGroup
 #   34. atomConnectivity2smiles_OLD
 #   35. moleculeGraph2smiles
 #   36. get_bondtypes
 #   37. draw_molecule_asGraph
# =============================================================================
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

def loadpickle_or_execute(pickle_path,function, *args, **kwargs):
    if os.path.exists(pickle_path):
        # Load data from pickle
        start = time.time()
        print('Loading data from pickle instead of {}...'.format(function.__name__))
        with open(pickle_path, 'rb') as pf:
            data = pickle.load(pf)
            stop = time.time()
            _min = int((stop-start)/60)
            _sec = (stop-start)%60
            if _min==0:
                print('Loaded Succesfully!! Loading time: {:0.1f} sec'.format(_sec))
            else:
                print('Loaded Succesfully!! Loading time: {} min {:0.1f} sec'.format(_min,_sec))
            print('-'*60)
    else:
        # Execute the provided function to get the data
        print('Loading data from {}...'.format(function.__name__))
        data = function(*args,**kwargs)
        print('Loaded Succesfully!!!')
        print('-'*60)
        # Save data to pickle for future use
        with open(pickle_path, 'wb') as pf:
            pickle.dump(data, pf)
    return data

def loadjson_or_execute(pickle_path,function, *args, **kwargs):
    if os.path.exists(pickle_path):
        # Load data from pickle
        print('Loading data from JSON instead of {}...'.format(function.__name__))
        print('Loaded Succesfully!!!')
        print('-'*60)
        with open(pickle_path, 'r') as pf:
            data = json.load(pf)
    else:
        # Execute the provided function to get the data
        print('Loading data from {}...'.format(function.__name__))
        data = function(*args,**kwargs)
        print('Loaded Succesfully!!!')
        print('-'*60)
        # Save data to pickle for future use
        with open(pickle_path, 'w') as pf:
            json.dump(data, pf)
    return data

# ---------------------List of Neighbours-----------------------------------
@function_runtime
def get_neighbours(bondfile, **kwargs):
    #10 times faster than version 1
    
    #-------keyward arguments-------------
    bo     = kwargs.get('bo',False) # to get bondorders
    mtypes = kwargs.get('mtypes',False) # to get molecule ids
    charge = kwargs.get('charge',False) # get charge of atoms
    nlp    = kwargs.get('nlp',False) # get the number of lone pair
    #-------------------------------------
    
    neighbours = {}
    atomtypes  = {}
    if bo: bondorders    = {}
    if mtypes: molecule_types    = {}
    if charge: charge = {}
    if nlp: nlp = {}
    
    with open(bondfile) as bf:
        prev_natoms = 0
        natoms_flag = False
        warning_flag = False
        first_warning_ignore = True
        atom_counter = 0
        
        for line in bf:
            splitted = line.split()
            if line.find('Timestep')!=-1:
                step = int(splitted[-1])
                neighbours[step]={}
                atomtypes = {}
                if bo: bondorders[step]  = {}
                if charge : charge[step] = {}
                if nlp: nlp[step] = {}
                
            if line.find('Number of particles')!=-1:
                current_natoms = int(splitted[-1])
                
                if atom_counter!=current_natoms:
                    if first_warning_ignore:
                        first_warning_ignore=False
                    else:
                        warning_flag=True
                    
                if natoms_flag and current_natoms!=prev_natoms:
                    print('User warning from get_neighbours function: Lost atom warning!!')
                atom_counter = 0   
                prev_natoms  = current_natoms #new change
            if splitted != [] and splitted[0].isnumeric():
                atom_counter            +=1
                parent                   = int(splitted[0])
                parent_type              = int(splitted[1])
                number_of_children       = int(splitted[2])
                children                 = list(map(int,splitted[3:3+number_of_children]))
                neighbours[step][parent] = children
                atomtypes[parent]        = parent_type
                
                if mtypes:
                    molecule_types[parent] = int(splitted[3+number_of_children])
                if bo:
                    bo_id = 4+number_of_children
                    bo_values = list(map(float,splitted[bo_id:bo_id+number_of_children]))
                    bondorders[step][parent]={}
                    for i,child in enumerate(children):
                        bondorders[step][parent][child] = bo_values[i]
                if charge:
                    charge[step][parent] = float(splitted[-1])
                if nlp:
                    nlp[step][parent] = float(splitted[-2])
        
        if warning_flag:
            print('User warning from get_neighbours function: Repeated atom information detected in some timesteps!')
            
        result = (neighbours,atomtypes,)
        if bo: result+=(bondorders,)
        if mtypes: result+=(molecule_types,)
        
        return result

# ---------------------List of Neighbours-----------------------------------
@function_runtime
def parsebondfile(bondfilepath, cutoff=0.3,**kwargs):
    print(f"bond order cutoff {cutoff} is used")
    #10 times faster than version 1
    
    #-------keyward arguments-------------
    bo              = kwargs.get('bo',False) # to get bondorders
    mtypes          = kwargs.get('mtypes',False) # to get molecule ids
    charge          = kwargs.get('charge',False) # get charge of atoms
    abo             = kwargs.get('abo',False)
    nlp             = kwargs.get('nlp',False) # get the number of lone pair
    mols            = kwargs.get('mols',False) # get mtypes wise molecues
    ALL             = kwargs.get('ALL',False) # get everything
    pkl             = kwargs.get('pkl',None) # Fatser the process (directory)
    firststep = kwargs.get('firststep',False) # first step only
    #-------------------------------------
    if firststep:
        pkl = None
        print('If firststep is True, pkl automatically set to None')
        
    if pkl is not None:
        ALL = True
        print('If pkl is not None, ALL automatically set to True')
    #-------------------------------------
    
    
    def parse():
        bonddata = {}
        neighbours = {}
        atomtypes  = {}
        # bondorders[step][atom_id1][atom_id2]=bo between atom1 and atom2
        if bo or ALL: bondorders    = {}
        if mtypes or ALL: molecule_types    = {}
        if charge or ALL: charges = {}
        if nlp or ALL: NLP = {}
        if abo or ALL: ABO = {}
        if mols or ALL: molecules = {}
        if firststep:
            fs_flag = False
        
        with open(bondfilepath) as bf:
            prev_natoms = 0
            natoms_flag = False
            warning_flag = False
            first_warning_ignore = True
            atom_counter = 0
            
            for line in bf:
                splitted = line.split()
                if line.find('Timestep')!=-1:
                    step = int(splitted[-1])
                    
                    if firststep:
                        if fs_flag:
                            break
                        fs_flag = True
                        
                    if bo or ALL: bondorders[step]  = {}
                    if charge or ALL: charges[step] = {}
                    if nlp or ALL: NLP[step] = {}
                    if abo or ALL: ABO[step] = {}
                    
                    neighbours[step]={}
                    atomtypes = {}
                    
                if line.find('Number of particles')!=-1:
                    current_natoms = int(splitted[-1])
                    
                    if atom_counter!=current_natoms:
                        if first_warning_ignore:
                            first_warning_ignore=False
                        else:
                            warning_flag=True
                        
                    if natoms_flag and current_natoms!=prev_natoms:
                        print('User warning from get_neighbours function: Lost atom warning!!')
                    atom_counter = 0   
                    prev_natoms  = current_natoms #new change
                    
                if splitted != [] and splitted[0].isnumeric():
                    atom_counter       +=1
                    parent             = int(splitted[0])
                    parent_type        = int(splitted[1])
                    number_of_children = int(splitted[2])
                    children           = list(map(int,splitted[3:3+number_of_children]))
                    
                    # bond orders
                    bo_id1 = 4+number_of_children
                    bo_idn = bo_id1+number_of_children
                    bo_values = list(map(float,splitted[bo_id1:bo_idn]))
                    
                    ## cutoffs
                    updated_children = []
                    for child,bond_order in zip(children,bo_values):
                        if bond_order>=cutoff:
                            updated_children.append(child)
                    
                    ## assigning values
                    neighbours[step][parent] = updated_children
                    atomtypes[parent]        = parent_type
                    
                    if mtypes or ALL:
                        molecule_types[parent] = int(splitted[3+number_of_children])
                    if bo or ALL:
                        bondorders[step][parent]={}
                        for i,child in enumerate(children):
                            bondorders[step][parent][child] = bo_values[i]
                        
                    if charge or ALL:
                        charges[step][parent] = float(splitted[-1])
                    if nlp or ALL:
                        NLP[step][parent] = float(splitted[-2])
                    if abo or ALL:
                        ABO[step][parent] = float(splitted[-3])
                    if mols or ALL:
                        mol_type = int(splitted[3+number_of_children])
                        if mol_type in molecules.keys():
                            molecules[mol_type].add(parent)
                        else:
                            molecules[mol_type] = set()
            
            if warning_flag:
                print('User warning from get_neighbours function: Repeated atom information detected in some timesteps!')
                
            bonddata['neighbours'] = neighbours
            bonddata['atypes'] = atomtypes
            if bo or ALL: bonddata['bondorders']  = bondorders
            if mtypes or ALL: bonddata['mtypes']  = molecule_types
            if charge or ALL: bonddata['charge']  = charges
            if nlp or ALL: bonddata['nlp'] = NLP
            if abo or ALL: bonddata['abo'] = ABO
            if mols or ALL: bonddata['molecules'] = molecules
            
            if firststep:
                for key, data in bonddata.items():
                    if key in ['neighbours','bondorders']:
                        for step, info in data.items():
                            bonddata[key]=info
                bonddata['fstep'] = step
            
            return bonddata
    
    if pkl is not None:
        if pkl == 'yes':
            path   = path = bondfilepath[:bondfilepath.rfind('.')]+'.pickle'
        else:
            path = pkl
            
        if os.path.exists(path):
            # Load data from pickle
            print('Loading data from pickle...')
            with open(path, 'rb') as pf:
                result = pickle.load(pf)
        else:
            # Execute the provided function to get the data
            print('No pickle file exist. Data is dumping for further use.')
            result = parse()
            with open(path, 'wb') as pf:
                pickle.dump(result, pf)
    
    else:
        result = parse()
        
    return result

#%%
@function_runtime
def parsebondfile_asGraph(bondfilepath):
    atomConnectivity = {} ## dict of nx graphs
            
    with open(bondfilepath) as bf:
        prev_natoms = 0
        natoms_flag = False
        repeated_atom_warning_flag = False
        first_warning_ignore = True
        atom_counter = 0
        hashline_flag = False
        
        ## predefining is necessary
        step = -1
        atomConnectGraph = nx.Graph()
        
        for line in bf:
            splitted = line.split()
            
            if hashline_flag and 'Timestep' in splitted:
                if not nx.is_empty(atomConnectGraph):
                    atomConnectivity[step]=atomConnectGraph
            
            if 'Timestep' in line:
                step             = int(splitted[-1])
                atomConnectGraph = nx.Graph()
                
                
            if line.find('Number of particles')!=-1:
                current_natoms = int(splitted[-1])
                
                if atom_counter!=current_natoms:
                    if first_warning_ignore:
                        first_warning_ignore=False
                    else:
                        repeated_atom_warning_flag=True
                    
                if natoms_flag and current_natoms!=prev_natoms:
                    print('User warning from get_neighbours function: Lost atom warning!!')
                atom_counter = 0   
                prev_natoms  = current_natoms #new change
                
            if splitted != [] and splitted[0].isnumeric():
                atom_counter    +=1
                parent           = int(splitted[0])
                parent_type      = int(splitted[1])
                n_children       = int(splitted[2])
                children         = list(map(int,splitted[3:3+n_children]))
                
                bo_id = 4+n_children
                bo_values = list(map(float,splitted[bo_id:bo_id+n_children]))
                
                atomConnectGraph.add_node(parent, atom_type=parent_type)
                for child, bo in zip(children,bo_values):
                    atomConnectGraph.add_edge(parent, child, bond_order=bo)

            if len(splitted)==1 and splitted[0]=='#':
                hashline_flag = True
            else:
                hashline_flag = False
                    
        if repeated_atom_warning_flag:
            print('User warning from get_neighbours function: Repeated atom information detected in some timesteps!')
    
    return atomConnectivity    
    
#%%---------get molecule from list neighbour-----------------------
def get_molecules(neigh):
    '''
    Parameters
    ----------
    neigh : Dictionary
        DESCRIPTION.
        neighbours of a single timestep or neighbours[step]

    Returns
    -------
    molecules : Set
        DESCRIPTION.
    '''
    graph = nx.Graph(neigh)
    molecules = nx.connected_components(graph)    
    return molecules

#%%-----------get chemical formula-----------------------
def get_molecular_formula(molecule,atomtypes,atomsymbols):
    '''
    Parameters
    ----------
    molecule : any itterable
        DESCRIPTION.
    atomtypes : TYPE
        DESCRIPTION.
    atomsymbols : TYPE, optional
        DESCRIPTION. The default is 'HCO'.

    Returns
    -------
    formula : TYPE
        DESCRIPTION.

    '''
    species = map(lambda x: atomtypes[x],molecule)
    counter = [0]*len(atomsymbols)
    for s in species:
        counter[s-1]+=1
    
    formula = ''
    for i in range(len(atomsymbols)):
        if counter[i]!=0:
            if counter[i]!=1:
                formula += atomsymbols[i]+str(counter[i])
            else:
                formula += atomsymbols[i]
    return formula

#%%---Get All Speceis Count At every Timestep--------
@function_runtime
def get_SpeciesCountAtEveryTimestep(neighbours,atomtypes,atomsymbols,exception=[],step2ps=False,step2psargs=[],minps=-math.inf,maxps=math.inf):
    '''
    Parameters
    ----------
    neighbours : Dictionary
        DESCRIPTION.
    atomtypes : Dictionary
        DESCRIPTION.
    atomsymbols : any iterable
        DESCRIPTION.
    exception : list (any iterable), optional
        DESCRIPTION. The default is [].
    step2ps : bool, optional
        DESCRIPTION. The default is False.
    step2psargs : list (any iterable), optional
        DESCRIPTION. The default is [].
    minps : TYPE, optional
        DESCRIPTION. The default is -math.inf.
    maxps : TYPE, optional
        DESCRIPTION. The default is math.inf.

    Returns
    -------
    stepwise_species_count : Dictionary
        DESCRIPTION.

    '''
    #local function
    def get_species(molecules,atomtypes,atomsymbols): #take the molecule file from get_molecute()
        #atomlist=['H','C','O'], same order as TYPE
        allspecies = []
        for molecule in molecules:
            species = get_molecular_formula(molecule, atomtypes,atomsymbols)
            allspecies.append(species)
        allspecies_count = Counter(allspecies)
        
        return allspecies_count

    ## if step2ps=True --> Convert timesteps to piccoseconds ##
    _steps  = list(neighbours.keys())
    s2ps = _steps.copy()
    s2ps_dict = {}
    
    if step2ps:
        if not step2psargs:
            sys.exit('Value error: \'step2psargs\' is empty')
        del s2ps
        s2ps = step2picosecond(_steps, *step2psargs)
    
    s2ps_dict = dict(zip(_steps,s2ps))
    ##-------------------------------------------------------##
    
    stepwise_species_count = {}
    for step,neigh in neighbours.items():
        molecules = get_molecules(neigh)
        species_count = get_species(molecules, atomtypes,atomsymbols)
        
        ## Delete species from the Exception list ##
        for ex in exception: del species_count[ex]
        ##------------------------------------------##
        if minps<=s2ps_dict[step]<maxps:
            stepwise_species_count[s2ps_dict[step]]=species_count
                
    return stepwise_species_count     
#%%
@function_runtime
def stepwise_species_count(neighbours,atomtypes,atomsymbols,step2ps=None):
    # if requires ps instead of steps set step2ps = timestep
    _swsc = {} # stepwise specis count _swsc[step][species]=count
    if step2ps:
        steps = list(neighbours.keys())
        ps    = step2picosecond(steps, step2ps)
    for i,items in enumerate(neighbours.items()):
        step,neigh = items
        molecules = get_molecules(neigh)
        _sc = {} # species count
        for molecule in molecules:
            species = get_molecular_formula(molecule, atomtypes,atomsymbols)
            if species in _sc.keys():
                _sc[species]+=1
            else:
                _sc[species] = 1
                
        if step2ps:
            _swsc[ps[i]]=_sc
        else:
            _swsc[step]=_sc
    return _swsc    
#%%
@function_runtime
def stepwise_species_count_v2(atomtypes,atomsymbols,*args,**kwargs):
    # if requires ps instead of steps set step2ps = timestep
    step2ps = kwargs.get('step2ps',None)
    kind    = kwargs.get('kind','sum')
    
    def inner(neighbours):
        # if requires ps instead of steps set step2ps = timestep
        _swsc = {} # stepwise specis count _swsc[step][species]=count
        if step2ps:
            steps = list(neighbours.keys())
            ps    = step2picosecond(steps, step2ps)
        for i,items in enumerate(neighbours.items()):
            step,neigh = items
            molecules = get_molecules(neigh)
            _sc = {} # species count
            for molecule in molecules:
                species = get_molecular_formula(molecule, atomtypes,atomsymbols)
                if species in _sc.keys():
                    _sc[species]+=1
                else:
                    _sc[species] = 1
                    
            if step2ps:
                _swsc[ps[i]]=_sc
            else:
                _swsc[step]=_sc
        return _swsc
    
    result = []
    for arg in args:
        result.append(inner(arg))
    
    output = {}
    for out in result:
        for step, species_count in out.items():
            if step not in output.keys():
                output[step]={}
            for species,count in species_count.items():
                if species not in output[step].keys():
                    output[step][species] = count
                else:
                    output[step][species] +=count
    if kind=='sum':
        pass
    elif kind=='mean':
        n = len(args)
        print(n)
        for step, species_count in output.items():
            for species,count in species_count.items():
                    cc = int(count/n)
                    output[step][species] = cc+1 if cc!=count/n else cc
    else:
        print('User Error from stepwise_species_count: kind key error')
    return output
        
#%%-----------Expression Selector------------------------------
def expression_selector(neighbours,atomtypes,file='expression_selector.txt',atomsybols='HCO',heading=None):
    steps = neighbours.keys()
    with open(file, 'w') as f:
        if heading:
            f.write(heading+'\n'+'-------------------------------'+'\n')
        for step in steps:
            molecules = get_molecules(neighbours[step])
            f.write('Timestep='+str(step)+'\n')
            
            for molecule in molecules:
                formula = get_molecular_formula(molecule, atomtypes,atomsybols)
                f.write('Timestep='+str(step)+'\t')
                f.write('Molecule: '+formula+'\n')
                f.write('Molecule: '+str(molecule)+'\n')
                for atom in molecule:
                    f.write('ParticleIdentifier=='+str(atom)+'|| ')
                f.write('\n\n')
                
#%%----------Pathway Tracker Of a singler molecule----------------------
def pathway_tracker(seek_molecule,neighbours,atomtypes,file='pathway_tracker.txt',atomsybols='HCO'):
    seek_molecule_formula = get_molecular_formula(seek_molecule, atomtypes)
    steps = neighbours.keys()
    file = '{}'.format(seek_molecule)+seek_molecule_formula+'_'+file
    
    with open(file, 'w') as f:
        f.write('This is a pathway tracker for {}'.format(seek_molecule_formula)+'\n')
        f.write('{}'.format(seek_molecule))
        f.write('-----------------------------------------------------------------\n\n')
        pathway = []
        for step in steps:
            molecules = get_molecules(neighbours[step])
            for molecule in molecules:
                formula = get_molecular_formula(molecule,atomtypes,'HCO')
                intersection = set(seek_molecule) & set(molecule)
                if intersection and formula not in pathway:
                    pathway.append(formula)
                    f.write('Timestep='+str(step)+'\t')
                    f.write('Molecule: '+formula+'\n')
                    for atom in molecule:
                        f.write('ParticleIdentifier=='+str(atom)+'|| ')
                        
                    f.write('\n{}\n\n'.format(molecule))
    return pathway

#%%--------------------------------------------------------
def step2picosecond(steps,*args):
    timestep =[]
    if len(args)==1:
        tstep1, = args
        current_time = steps[0]*tstep1/1000
        previous_step = steps[0]
        
        for current_step in steps:
            current_time += (current_step-previous_step)*tstep1/1000  
            previous_step = current_step
            timestep.append(current_time)
        return timestep
    
    elif len(args)==3:
        tstep1,tstep2,limit = args
        current_time = steps[0]*tstep1/1000
        previous_step = steps[0]
        
        for current_step in steps:
            if current_step<=limit: 
                current_time += (current_step-previous_step)*tstep1/1000
            else:        
                current_time += (current_step-previous_step)*tstep2/1000
            previous_step = current_step
            timestep.append(current_time)
        return timestep
    else:
        sys.exit('Error: The length of the argument in the step2picosecond() function is not specified correctly.\n')

def step2ps(steps,timestep):
    # timestep in fs
    # steps can be np array
    ps = steps*timestep/1000
    return ps
#%%--sort chemical formula in this order: C,H,O---------------
def sort_molecular_formula(molecular_formula,order=['C','H','O']):
    """
    Sorts the elements in a molecular formula (or a list of formulas) according
    to a specified order.
    
    Parameters:
    molecular_formula (str or iterable of str): The molecular formula(s)
                                                to be sorted.
    order (list of str): The order in which elements should be sorted within
                         the formula. Default is ['C', 'H', 'O'].

    Returns:
    str or list of str: Sorted molecular formula(s).
    """
    
    def do_sort(species):
        """
        Helper function to sort a single molecular formula.
        """
        item = re.findall('[A-Z][a-z]?|\d+|.', species)+['']
        elements = []
        for i in range(len(item)-1):
            if item[i].isalpha() and item[i+1].isnumeric():
                elements.append((item[i],item[i+1]))
            elif item[i].isalpha() and not item[i+1].isnumeric():
                elements.append((item[i],''))
                
        # Check if all elements are in the specified order list
        symbols = [x for x, _ in elements]
        if not set(symbols).issubset(set(order)):
            raise ValueError('Some elements are not in the order list! '
                             'Please specify the order list accordingly.')
        
        for i in range(len(elements)):
            elements[i]+=(order.index(elements[i][0]),)
        elements.sort(key=lambda x:x[2])
        sorted_chem_formula = ''
        for element in elements:
            sorted_chem_formula+=element[0]+element[1]
        return sorted_chem_formula
          
    
    # Handle both single string and iterable of strings cases
    if isinstance(molecular_formula, str):
        # molecular_formula is a single string
        result = do_sort(molecular_formula)
    else:
        try:
            iterator = iter(molecular_formula)
            if all(isinstance(item, str) for item in iterator):
                # Input is an iterable of strings
                result = []
                for species in molecular_formula:
                    result.append(do_sort(species))
            else:
                raise TypeError("Iterable must contain only strings")
                
        except TypeError:
            # Input is neither a single string nor an iterable of strings
            raise TypeError("Input must be a string or an iterable of strings")
        
    return result
   
#%%----------------make_latex-----------------------------------------
def make_molecular_formula_latex(molecular_formula,sort=False):
    """
    Converts a molecular formula or a list of formulas to LaTeX format.

    Parameters:
    molecular_formula (str or iterable of str): The molecular formula(s)
                                                to be converted.
    sort (bool): Whether to sort the elements in the formula(s) according to
                 a predefined order.

    Returns:
    str or list of str: Molecular formula(s) in LaTeX format.
    """
    
    if sort:
        molecular_formula = sort_molecular_formula(molecular_formula)
    
    if isinstance(molecular_formula, str):
        formula = [molecular_formula]
    else:
        formula = list(molecular_formula).copy()
    
    
    
    latex_formula = ['']*len(formula)
    for i in range(len(formula)):
        formula[i]+='$'
        for j in range(len(formula[i])-1):
            if formula[i][j].isalpha() and formula[i][j+1].isnumeric():
                latex_formula[i]+=formula[i][j] + '_{'
            elif formula[i][j].isnumeric() and not formula[i][j+1].isnumeric():
                latex_formula[i]+=formula[i][j] + '}'
            else:
                latex_formula[i]+=formula[i][j]
        latex_formula[i] = '$'+latex_formula[i]+'$'
    
    if len(latex_formula)==1:
        return latex_formula[0]
    else:
        return latex_formula

#%%-------------------------------------------------------------
def compute_molecular_weight(species,exclude=[]):
    atomic_info    = atomic_weight(key='all')
    item = re.findall('[A-Z][a-z]?|\d+|.', species)
    item.append('q') #fake
    molecular_weight = 0.0
    for i in range(len(item)-1):
        now,nxt = item[i],item[i+1]
        if now in exclude:
            continue
        if now.isalpha() and nxt.isalpha():
            if now in atomic_info:
                molecular_weight+=atomic_info[now]
            else:
                sys.exit('Error: No such element \'{}\' in the compute_molecular_weight() function'.format(now))
        elif now.isalpha() and nxt.isnumeric():
            if now in atomic_info:
                molecular_weight+=atomic_info[now]*float(nxt)
            else:
                sys.exit('Error: No such element \'{}\' in the compute_molecular_weight() function'.format(now))
    return molecular_weight

#%%-------------------Heatmap of Species--------------------------------
@function_runtime
def plot_species_heatmap(atomtypes,atomsymbols,*neighbours_list,**kwargs):
    # keyward agruments:
        # nspecies (int)   = number of species to be plot (default value = 20)
        # ts (int)         = timestep (default value = 0.25)
        # savedir (string) = directory to save (default value = '')
        # titile (string)  = title to be shown on the top (default value = '')
        # order (list)     = order of the species. Takes list of species (default = mean abundance)
        # pikle (string)   = load from pickle. It takes pickle_path (default = None)
        # skipts (float)   = skip some data from the first (default = None)
        # topspec (list)   = species that appears on the top (default = None)
        # 
        
    nspecies    = kwargs.get('nspecies', 20)
    ts          = kwargs.get('ts', 0.25)
    savedir     = kwargs.get('savedir', '')
    title       = kwargs.get('title', '')
    order       = kwargs.get('order',None)
    pickle      = kwargs.get('pickle',None)
    skipts      = kwargs.get('skipts',None) # skip timestep
    topspec     = kwargs.get('topspec',None)
    ignor       = kwargs.get('ignor',None)
    kind        = kwargs.get('kind','sum') # kind = sum,mean
    log         = kwargs.get('log',True)
    exclude     = kwargs.get('exclude',None)
    fontsize    = kwargs.get('fontsize', 12)
    figsize     = kwargs.get('figsize',(15,8))
    
    if pickle:
        data = loadpickle_or_execute(pickle, stepwise_species_count_v2, atomtypes, atomsymbols, *neighbours_list, step2ps=ts)
    else:
        data = stepwise_species_count_v2(atomtypes, atomsymbols,*neighbours_list,step2ps=ts,kind=kind)
    
    if log:
        df = pd.DataFrame(data).fillna(1)
    else:
        df = pd.DataFrame(data).fillna(0)
    
    if order:
        df = df.loc[order,:]
    else:
        df.loc[:,'sum'] = df.sum(axis=1)
        df = df.sort_values(by='sum',ascending=False)
        df = df.drop(['sum'],axis=1)
        
    #############-exclude-##################
    if exclude is not None:
        df = df.drop(index=exclude)
    ######################################
    
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################
    
    #############-ignor-##################
    if ignor is not None:  # ignor = (main_species,'close')
        if isinstance(ignor, tuple):
            main_species, kind = ignor
            main = compute_molecular_weight(main_species)
            
            H = compute_molecular_weight('H')
            O = compute_molecular_weight('O')
            criteria = []
            
            for a in [0,1,2,3]:
                for b in [0,1,2,3]:
                    if a==0 and b==0:
                        continue
                    cri = a*H+b*O
                    criteria.append(cri)
            
            for chem in df.index:
                current = compute_molecular_weight(chem)            
                for cri in criteria:
                    if abs(abs(main-current)-cri)<=0.01:
                        print(chem,'ignored')
                        df = df.drop(chem)
                        continue
        elif isinstance(ignor, str):
            df = df.drop([ignor],axis=0)
    #####################################
    
    
    df = df.head(nspecies)
    ##############-topspec-###############
    if topspec is not None:
        desired_order = []
        for spec in topspec:
            if spec in df.index:
                desired_order.append(spec)
        for spec in df.index:
            if spec not in desired_order:
                desired_order.append(spec)
        df = df.reindex(index=desired_order)
    ######################################
        
    print('Species ORDER:',list(df.index))
    df.index = make_molecular_formula_latex(df.index,sort=True)
    _, ax1 = plt.subplots(figsize=figsize)
    if log:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                     cbar_kws={'label': 'Number of species'},
                     norm=LogNorm(),xticklabels=1000)
    else:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                     cbar_kws={'label': 'Number of species'},
                     xticklabels=1000)
    ax.set_xlabel('Time (ps)',fontsize=fontsize+1)
    # ax.set_ylabel('Species')
    plt.tick_params(left=False,bottom=False)
    plt.yticks(rotation=0,fontsize=fontsize)
    fig = ax.get_figure()
    title += '-{}\nTop {} species\n'.format(kind,nspecies)
    ax.set_title(title)
    ax.set_xticklabels([0,250,500,750,1000],fontsize=fontsize)
    
    
    # cbar
    ax.figure.axes[-1].yaxis.label.set_size(fontsize+1) # cbar label size
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=fontsize)
    
    saveplot = 'species_heatmap_'+str(random.randint(100000, 999999))
    if savedir: savedir += '\\'+saveplot
    else: savedir += saveplot
    fig.savefig(savedir, dpi=400, bbox_inches='tight')
    
#%%-------Heatmap of Species from multiple neighbours-------
@function_runtime
def plot_species_heatmap_v2(neighbours,atomtypes,atomsymbols,**kwargs):
    # keyward agruments:
        # nspecies (int)   = number of species to be plot (default value = 20)
        # ts (int)         = timestep (default value = 0.25)
        # savedir (string) = directory to save (default value = '')
        # titile (string)  = title to be shown on the top (default value = '')
        # order (list)     = order of the species. Takes list of species (default = mean abundance)
        # pikle (string)   = load from pickle. It takes pickle_path (default = None)
        # skipts (float)   = skip some data from the first (default = None)
        # topspec (list)   = species that appears on the top (default = None)
        
    nspecies    = kwargs.get('nspecies', 20)
    ts          = kwargs.get('ts', 0.25)
    savedir     = kwargs.get('savedir', '')
    title       = kwargs.get('title', '')
    order       = kwargs.get('order',None)
    pickle      = kwargs.get('pickle',None)
    skipts        = kwargs.get('skipts',None) # skip timestep
    topspec     = kwargs.get('topspec',None)
    ignor       = kwargs.get('ignor',None)
    kind        = kwargs.get('kind','sum') # kind = 'sum' or 'mean'
    log         = kwargs.get('log',True)
    exclude     = kwargs.get('exclude',None)
    
    if pickle:
        data = loadpickle_or_execute(pickle, stepwise_species_count, neighbours,atomtypes, atomsymbols, step2ps=ts)
    else:
        data = stepwise_species_count(neighbours, atomtypes, atomsymbols,
                                      step2ps=ts)
    ###########-fill with###############
    if log:
        df = pd.DataFrame(data).fillna(1)
    else:
        df = pd.DataFrame(data).fillna(0)
    ####################################
    
    if order:
        df = df.loc[order,:]
    else:
        df.loc[:,'sum'] = df.sum(axis=1)
        df = df.sort_values(by='sum',ascending=False)
        df = df.drop(['sum'],axis=1)
    
    #############-exclude-##################
    if exclude is not None:
        df = df.drop(index=exclude)
    ######################################
    
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################
    
    #############-ignor-##################
    if ignor is not None:  # ignor = (main_species,'close')
        if isinstance(ignor, tuple):
            main_species, kindd = ignor
            main = compute_molecular_weight(main_species)
            
            H = compute_molecular_weight('H')
            O = compute_molecular_weight('O')
            criteria = []
            
            for a in [0,1,2,3]:
                for b in [0,1,2,3]:
                    if a==0 and b==0:
                        continue
                    cri = a*H+b*O
                    criteria.append(cri)
            
            for chem in df.index:
                current = compute_molecular_weight(chem)            
                for cri in criteria:
                    if abs(abs(main-current)-cri)<=0.01:
                        print(chem,'ignored')
                        df = df.drop(chem)
                        continue
        elif isinstance(ignor, str):
            df = df.drop([ignor],axis=0)
    #####################################
    
    
    df = df.head(nspecies)    
    ##############-topspec-###############
    if topspec is not None:
        desired_order = []
        for spec in topspec:
            if spec in df.index:
                desired_order.append(spec)
        for spec in df.index:
            if spec not in desired_order:
                desired_order.append(spec)
        df = df.reindex(index=desired_order)
    ######################################
        
    print('Species ORDER:',list(df.index))
    df.index = make_molecular_formula_latex(df.index,sort=True)
    #df.loc['$O_{2}$',:] = df.loc['$O_{2}$',:] + 200
    _, ax1 = plt.subplots(figsize=(15,8))
    if log:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                         cbar_kws={'label': 'Number of species'},
                         norm=LogNorm(),xticklabels=1000)
    else:
        ax = sns.heatmap(df,cmap='jet',ax=ax1 ,
                         cbar_kws={'label': 'Number of species'},
                         xticklabels=1000)
        
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Species')
    plt.tick_params(left=False,bottom=False)
    plt.yticks(rotation=0)
    fig = ax.get_figure()
    title += '---{}\nTop {} species'.format(kind,nspecies)
    ax.set_title(title)
    plt.show()
    
    saveplot = 'species_heatmap_'+str(random.randint(100000, 999999))
    if savedir: savedir += '\\'+saveplot
    else: savedir += saveplot
    fig.savefig(savedir, dpi=400, bbox_inches='tight')
#%%-----Get Atomic Mass --------------------------------
def atomic_weight(*args,**kwargs):
    # args   = element string
    # kwargs = {key:'all'}
    key = kwargs.get('key',None)
    MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182,
                      'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,'F': 18.9984032,
                      'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
                      'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
                      'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415,
                      'Cr': 51.9961, 'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934,
                      'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216,
                      'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62,
                      'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063,
                      'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411,
                      'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6, 'I': 126.90447,
                      'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547,
                      'Ce': 140.116, 'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36,
                      'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5,
                      'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421,
                      'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49,
                      'Ta': 180.9479, 'W': 183.84, 'Re': 186.207,
                      'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
                      'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833,
                      'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824,
                      'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197,
                      'Ra': 226.0254, 'Ac': 227.0278, 'Th': 232.03806,
                      'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482,
                      'Pu': 244.0642, 'Am': 243.0614, 'Cm': 247.0703,
                      'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829,
                      'Fm': 257.0951, 'Md': 258.0951, 'No': 259.1009,
                      'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271,
                      'Bh': 270, 'Hs': 269, 'Mt': 278, 'Ds': 281,
                      'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289,
                      'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,'': 0}

    if len(args)==1:
        element = args[0]
        return MM_of_Elements[element]
    elif key=='all':
        return MM_of_Elements
    else:
        sys.exit('User Error from atomic_weight: Unknown kwargs or number of args')
        
#%%-----ONSET Finding----------------------
@function_runtime
def get_speciesVStemp(species,neighbour,atomtypes,atomsymbols,**kwargs):
    ########## Soft Warnings ################
    if 'ramp' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Ramping rate not specified by the user. The default value of 1 K/ps has been set')
    if 'ts' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Timestep not specified by the user. The default value of 0.25 fs has been set')
    if 'it' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Initial temperature not specified by the user. The default value of 0 K has been set')
        
    ######### Getting kwargs #################    
    ramp         = kwargs.get('ramp',1) # temp ramp rate (Default: 1 K/ps)
    it           = kwargs.get('it',0) # initial temp (Default: 0 K)
    ts           = kwargs.get('ts',0.25)
    method       = kwargs.get('method','mf')
                   # calculate onset based on method
                   # mf = molecular formula
                   # mw = molecular weight
    lim        = kwargs.get('lim',0)
                   # if method='mw', it will count species if the mw
                   # within the 'lim'
    
    if method=='mf':
        count = {}
        for step,neigh in neighbour.items():
            count[step]=0
            molecules = get_molecules(neigh)
            for molecule in molecules:
                mf = get_molecular_formula(molecule, atomtypes, atomsymbols)
                if mf==species:
                    #print(mf)
                    count[step]+=1
            #print(step,count[step])

    elif method=='mw':
        count = {}
        mw_species = compute_molecular_weight(species)
        for step,neigh in neighbour.items():
            count[step]=0
            molecules = get_molecules(neigh)
            for molecule in molecules:
                mf = get_molecular_formula(molecule, atomtypes, atomsymbols)
                mw = compute_molecular_weight(mf)
                if abs(mw-mw_species)<=lim:
                    count[step]+=1
                    #print(mw,mw_species,end='\t')
    else:
        sys.exit('Method not found!!!')
    
    steps          = list(count.keys())
    sc             = np.array(list(count.values())) # species count
    ps             = np.array(step2picosecond(steps, ts))
    temp           = it + ramp*ps
    return temp,sc

#%%
def get_onset(temp,sc,**kwargs):
    imc          = kwargs.get('imc',50) # initial molecular count
    ig           = kwargs.get('ig',[2200,107,0,50]) # initial guess for curve fit
    show_fit     = kwargs.get('show_fit',False)
    
    #-----Fit Function----------------------
    def fit(x, A, B, C, D):   
        y = C+((D-C)/(1+np.exp((x-A)/B)))
        return y

    def inv_fit(y,A,B,C,D):
        if D-y<0:
            print('Warning from get_onset: Negative D value adjusted!!')
            D = 49.98
        x = A + B*np.log((D-y)/(y-C))
        return x
    #-----------------------------------------
    
    p, cov = curve_fit(fit, temp, sc,ig,maxfev=2000)
    if show_fit :
        print('Parameters: ',p)
        print('Covariance',cov)
    sc_fit         = fit(temp,*p)
    onset          = inv_fit(imc-1, *p) # onset temp at which 1 molecule
    
    # R^2 Value
    sc_mean = np.mean(sc)
    tss = np.sum((sc - sc_mean)**2)
    rss = np.sum((sc - sc_fit)**2)
    r2 = 1 - (rss / tss)
    print('R square value: ', r2)
    
    return onset, sc_fit
#%%
def bondorder_evolution(bondorders,bondlist,**kwargs):
    if 'ts' not in kwargs.keys():
        print('User Warning from get_speciesVStemp: Timestep not specified by the user. The default value of 0.25 fs has been set')
    ## Getting kwargs ##
    ts           = kwargs.get('ts',0.25)
    title        = kwargs.get('title',"")
    savedir      = kwargs.get('savedir', '')
    skipts       = kwargs.get('skipts',None) # skip timestep
    sort         = kwargs.get('sort',None) #ascending,descending or index(ts)
    figsize      = kwargs.get('figsize',(6,10))
    fontsize     = kwargs.get('fontsize',10)
    plot         = kwargs.get('plot','yes')
    ps2temp      = kwargs.get('ps2temp',None) # (it,ramp)
                        # it   = initial temp
                        # ramp = ramping rate
    
    steps = list(bondorders.keys())
    bo_evolution = {}
    ps = np.array(step2picosecond(steps, ts))
    xlabel = ''#'Time (ps)'
    
    if ps2temp:
        xlabel = 'Temrature (K)'
        it,ramp = ps2temp
        ps = it + ps*ramp

    for index, atompair in enumerate(bondlist):
        atom1,atom2 = atompair
        bo_evolution[index+1]={}
        for step,p in zip(steps,ps):
            bo = bondorders[step][atom1].get(atom2,0)
            bo_evolution[index+1][p]=bo
            
    #############-sort-##################       
    df = pd.DataFrame(bo_evolution)
    df = df.transpose()    
    if sort is not None:
        if sort == 'ascending':
            df.loc[:,'sum'] = df.sum(axis=1)
            df = df.sort_values(by='sum')
            df = df.drop(['sum'],axis=1)      
        elif sort == 'descending':
            df.loc[:,'sum'] = df.sum(axis=1)
            df = df.sort_values(by='sum',ascending=False)
            df = df.drop(['sum'],axis=1)
        elif type(sort) in [float,int]:
            df = df.sort_values(by=sort)
        else:
            print('User Error from bondorder_evolution: sort key not found')
            sys.exit()
    #####################################
    df.index = range(1,len(bondlist)+1)
    #############-skipts-##################
    if skipts is not None:
        df = df.loc[:,skipts:]
        df.columns = df.columns-skipts
    ######################################
    
    if plot=='yes':
        _, ax1 = plt.subplots(figsize=figsize)
        hm = sns.heatmap(df,cmap='jet',cbar_kws={ 'label': 'Bond order',},xticklabels=1000,ax=ax1,vmin=0.0,vmax=3.0)
        hm.set_xlabel(xlabel,fontsize=fontsize+2)
        # hm.set_ylabel('Bonds ',fontsize=fontsize)
        plt.yticks(rotation=0)
        
        #############-Y Tick Labels-###########
        hm.set_yticklabels([])
        
        ## set xticks location and labels ####
        # hm.set_xticks([795*x for x in range(5)])
        print([int(x) for x in hm.get_xticks()])
        hm.set_xticklabels([0,250,500,750,1000])
        plt.xticks(rotation=90,fontsize=fontsize+2)
        ###########
        
        hm.set_title(title,fontsize=fontsize-5)
        # Set the font size of the y-label
        hm.set_yticklabels(hm.get_yticklabels(), fontsize=fontsize)
    
    
        ######## colorbar settings #########
        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(labelsize=fontsize+2)
        hm.figure.axes[-1].yaxis.label.set_size(15+4)
        ####################################
        
        if savedir:
            fig = hm.get_figure()
            savedir = savedir+'\\bo_evolution_'+str(random.randint(100000, 999999))
            fig.savefig(savedir, dpi=400, bbox_inches='tight')
    
    return df

#%%
@function_runtime
def get_nbondsVStime(bondorders, atomtypes, bondlist,**kwargs):
    # plot the number of bonds vs time
    ## Getting kwargs and setting up default values ##
    cutoff = kwargs.get('cutoff',1)
    tol    = kwargs.get('tol',0.20)
    ts     = kwargs.get('ts',0.25) # timestep
    plot   = kwargs.get('plot',None) # values: 'yes'
    skipts = kwargs.get('skipts',None) # skip timestep
    cumu   = kwargs.get('cumu',False) 
    
    ## Getting args
    
    nbonds = []
    cumu_nbonds = []
    cumu_count = 0
    steps  = []
    for step,bo in bondorders.items():
        steps.append(step)
        bondcount = 0
        for bond in bondlist:
            atom1, atom2 = bond
            try:
                bovalue = bo[atom1][atom2]
            except:
                bovalue = 0
                
            # if bovalue!=0: print(step,bovalue,end=' | ')
            if abs(bovalue-cutoff)<=tol:
                bondcount+=1
                cumu_count+=1
        nbonds.append(bondcount)
        cumu_nbonds.append(cumu_count)
        
    ps = step2picosecond(steps, ts)
    
    if skipts is not None:
        startid = ps.index(skipts)
        ps      = np.array(ps[startid+1:])
        nbonds  = np.array(nbonds[startid+1:])
        if cumu:
            cumu_nbonds  = np.array(cumu_nbonds[startid+1:])
        ps = ps - skipts
        
    # ploting
    if plot == 'yes':
        plt.plot(ps,nbonds,color='black')
        plt.xlabel('Time (ps)')
        plt.ylabel('Number of bonds')
    
    if cumu:
        return ps, cumu_nbonds
    else:
        return ps, nbonds
        
    
#%%
@function_runtime
def species_to_molecule(bonddata,atomsymbols, species,**kwargs): # species-->molecule
    # speceis is chemcical formula of a molecule
    # molecule here is the list of atom ids that refer the molecule
    
    ###########-get kwargs-###############
    source    = kwargs.get('source',False)
    dump      = kwargs.get('dump',None)
    shortinfo = kwargs.get('shortinfo',False)
    ######################################
    def atomids2expression(atoms): # atom ids to ovito expression selection
        atoms = list(atoms)
        expression = ''
        for atom in atoms:
            if atoms.index(atom)==len(atoms)-1:
                expression+='ParticleIdentifier=='+str(atom)
            else:
                expression+='ParticleIdentifier=='+str(atom)+'|| '
        
        return expression
    
    
    
    neighbours = bonddata['neighbours']
    # bondorders = bonddata['bondorders']
    atypes  = bonddata['atypes']
    mtypes     = bonddata['mtypes']
    molecules  = bonddata['molecules']
    
    seeklist    = []
    seekstep    = []
    if dump is not None:
        f = open(dump,'w')
    if shortinfo:
        sf = open(dump[:dump.rfind('.')]+'_shortinfo.txt','w')
    for step,neigh in neighbours.items():
        PRINT = []
        graph = nx.Graph(neigh)
        connected = nx.connected_components(graph)
        for component in connected:
            componentformula = get_molecular_formula(component, atypes, atomsymbols)
            if componentformula == species and component not in seeklist:                
                seeklist.append(component)
                seekstep.append(step)
                # print(step)
                # print(componentformula,component)
                # print(atomids2expression(component))
                if dump is not None:
                    print(step,file=f)
                    print(componentformula,component,file=f)
                    print(atomids2expression(component),file=f)
                if shortinfo:
                    print(step,file=sf)
                    print(componentformula+' =',end=' ',file=sf)
                if source:
                    src = {mtypes[x] for x in component}
                    for j,s in enumerate(src):
                        # print('Source {}:'.format(s))
                        # print(atomids2expression(molecules[s]))
                        # print('-------------------')
                        if dump is not None:
                            print('Source {}:'.format(s),file=f)
                            print(atomids2expression(molecules[s]),file=f)
                            print('-------------------',file=f)
                        if shortinfo:
                            PRINT.append('M-{}'.format(s))
                # print('-'*200)
                if dump is not None: print('-'*200,file=f)
                if shortinfo:
                    print(*PRINT,sep=' + ',file=sf)
                    print(*PRINT,sep=' + ')
    f.close()
    return seeklist,seekstep

#%%
@function_runtime
def cumulative_nspecies(neighbours,atypes,atomsymbols,specieslist,**kwargs):
    ###########-get kwargs-###############
    skipts    = kwargs.get('skipts',None)
    ts        = kwargs.get('ts',0.25)
    ######################################
    
    count = [0]*len(specieslist)
    steps = list(neighbours.keys())
    ccount = [[] for i in range(len(specieslist))]
    checklist = [[] for i in range(len(specieslist))]
    for step,neigh in neighbours.items():
        molecules = get_molecules(neigh)
        for molecule in molecules:
            formula = get_molecular_formula(molecule, atypes, atomsymbols)
            for i,species in enumerate(specieslist):
                if species==formula and molecule not in checklist[i]:
                    checklist[i].append(molecule)
                    count[i]+=1
        for i in range(len(ccount)):
            ccount[i].append(count[i])
        
    ps     = step2picosecond(steps, ts)
    
    if skipts is not None:
        index  = ps.index(skipts)
        ps     = ps[index:]
        for i in range(len(ccount)):
            ccount[i] = ccount[i][index:]
    
    ps = np.array(ps)-min(ps)
    print(len(ccount))
    ccount = [np.array(x) for x in ccount]  
    
    return ps,ccount
#%%
def onset_plot(path,whole,atomsymbols,timestep,temp_ramp,initial_temp,**kwargs):
    
    ## getting kwargs
    color    = kwargs.get('color',None)  # assign random color
    sim_path = kwargs.get('sim_path',['Sim-1','Sim-2','Sim-3'])
    ax       = kwargs.get('ax',None)
    
    
    df = pd.DataFrame()
    ## Iterate over number of simulations
    for sim in sim_path:
        bondfilepath = path+'\\'+sim+'\\bonds.reaxc'
    
        ## Geting neighbours for each simulations
        bonddata   = parsebondfile(bondfilepath)
        neighbours = bonddata['neighbours']
        atypes     = bonddata['atypes']
        
        ## Getting temperatures using steps, timesteps, temp_ramp, initial_temp
        steps      = np.array(list(neighbours.keys()))
        time       = steps*timestep/1000  # in piccosecond
        temp       = initial_temp + time*temp_ramp
        df['temp'] = temp
        
        ##  Getting number of whole
        nwhole     = []
        for step, neigh in neighbours.items():
            count = 0
            molecules = get_molecules(neigh)
            for molecule in molecules:
                species = get_molecular_formula(molecule, atypes, atomsymbols)
                if species == whole:
                    count+=1
            nwhole.append(count)
        df[sim]=nwhole
    
    ## Getting the Upper and Lower bound
    nsim = len(sim_path)
    upper_bound = df.iloc[:, -nsim:].max(axis=1)
    lower_bound = df.iloc[:, -nsim:].min(axis=1)
    
    ## plotting the fill-between plot
    x       = df['temp']
    if ax is None: fig, ax = plt.subplots()
    ax.fill_between(x, lower_bound, upper_bound, color=color,alpha=0.2)
    
    ## default fit function
    def function(x, A, B, C, D):   
        y = C+((D-C)/(1+np.exp((x-A)/B)))
        return y
    ## creating average sim data out of df
    y         = df.iloc[:,-nsim:].mean(axis=1)
    
    ## curve fit
    popt, cov = curve_fit(function, x, y,p0=[2200,107,0,25])
    y_fit     = function(x, *popt)
    
    ## plot fit values
    ax.plot(x,y_fit,color=color)
    
    ## getting onset using newton-rahpson or secant
    y_target  = df.iloc[:100,-1].mean()
    print(y_target)
    onset     = newton(lambda x: function(x,*popt)-y_target, x0=300)
    print('Onset: ', onset)
    
    ## plotting a onset indicating verticle dashed line
    ax.plot([onset]*2,[0,y_target],'--',color=color)
    
    return fig,ax

#%%
@function_runtime
def get_species_count(bondfilepath,atomsymbols, cutoff: float = 0.3,
                      timestep = None, restart_time: bool = False,
                      temp: tuple = None, frame: bool = False):
    
    '''
    input: ---, temp = (initial_temp, ramp_rate)
    output: a panda dataframe of specie timeseries
    '''
    
    # Show errors
    if temp is not None:
        if not isinstance(temp, Iterable) or len(temp) != 2:
            raise ValueError("The 'temp' parameter must be an iterable with exactly two elements: (initial_temp, ramp_rate)")

    bonddata   = parsebondfile(bondfilepath,cutoff=cutoff)
    neighbours = bonddata['neighbours']
    atypes     = bonddata['atypes']
    
    count = {}
    for step, neigh in neighbours.items():
        molecules = get_molecules(neigh)
        count[step] = {}
        for molecule in molecules:
            species = get_molecular_formula(molecule, atypes, atomsymbols)
            if species in count[step]:
                count[step][species]+=1
            else:
                count[step][species]=1
        
    df = pd.DataFrame(count).fillna(0).T
    df.index.name = 'Timestep'
    
    if timestep is not None:
        df['Time'] = df.index*timestep/1000
    
    if restart_time:
        df['Time'] = df['Time'] - df['Time'].min()
    
    if temp is not None:
        time = df.index*timestep/1000
        initial_temp, ramp_rate = temp
        df['Temperature'] = initial_temp+ramp_rate*time
    
    if frame:
        df = df.reset_index(drop=True)  # Reset to default 0-based index first
        df.index = df.index + 1         # Make it 1-based
        df.index.name = 'Frame'
    
    return df

#%%
## doubt
@function_runtime
def count_functinalGroup(neighbours,atypes,seek):
    group = {'OH': [3,[1,2]],
             'COOH': [2,[1,3,3]],
             'Keto': [2,[2,2,3]],
             'Aldy': [2,[1,2,3]]}
    
    count = [0]*len(seek)
    for i, ss in enumerate(seek):
        if ss not in group:
            target_parent_type    = group[ss][0]
            target_children_types = group[ss][1]
        else:
            print('hi')
            raise('{} functional group not found!'.format(ss))
        
        checklist = []
        for step, neigh in neighbours.items():
            for parent, children in neigh.items():
                children_types = sorted([atypes[x] for x in children])
                match = atypes[parent]==target_parent_type and \
                            children_types == sorted(target_children_types)    
                if match and parent not in checklist:
                    count[i]+=1
                    checklist.append(parent)
    return count

#%%
def atomConnectivity2smiles_OLD(atomConnectivity,atypes,atomic_num):
    # atomConnectivity is a networkx Graph (to be speciec subgraph, same
    # shit though).  
    # This converts connected atom graphs, which are obtained from the atom
    # neighbors in a single frame to SMILES. Essentially, the connectivity 
    # represents a subgraph of each molecule, extracted from its neighbors.
    
    mol = Chem.RWMol()
    
    # Tips:
    # If the atom IDs in your atomConnectivity graph are not sequential
    # starting from 0 or if they are random, you'll need to create a mapping
    # between these IDs and the indices of atoms in the RDKit molecule.
    # This is because RDKit expects atom indices to start from 0 and
    # increment sequentially.
    # No worries! I have done this already. atomid2index is the mapping
    
    
    ## adding atoms
    atomid2index = {}
    for node in atomConnectivity.nodes():
        atomicNumber = atomic_num[atypes[node]-1]
        if atomicNumber==1:
            continue
        atom_index = mol.AddAtom(Chem.Atom(atomicNumber))
        atomid2index[node] = atom_index
        
    ## adding bonds
    for u, v, bond_type in atomConnectivity.edges(data='bond_type'):
        uan = atomic_num[atypes[u]-1]
        van = atomic_num[atypes[v]-1]
        if uan==1 or van==1:
            continue
        
        if bond_type==1:
            rdkitBondType = Chem.BondType.SINGLE
        elif bond_type==2:
            rdkitBondType = Chem.BondType.DOUBLE
        elif bond_type==3:
            rdkitBondType = Chem.BondType.TRIPLE
        else:
            print('User ERROR: bond type is not recognized')
            print('Setting the bond type as single')
            rdkitBondType = Chem.BondType.SINGLE
            
        mol.AddBond(atomid2index[u], atomid2index[v], rdkitBondType)
    
    mol = mol.GetMol()
    smiles = Chem.MolToSmiles(mol)
    return smiles

#%%
def moleculeGraph2smiles(moleculeGraph, atomic_num,
                         n_clusters, plot_cluster=False,
                         bo_analysis=True, atom_types=None):
    # moleculeGraph is a networkx Graph (to be speciec subgraph, same
    # shit though) of a single molecule having node attr as atom_type and
    # edge attr as bond_order. This converts molecule graphs to SMILES.
    import warnings
    warnings.filterwarnings("ignore")
    
    mol = Chem.RWMol()
    
    # Tips:
    # If the atom IDs in your moleculeGraph graph are not sequential
    # starting from 0 or if they are random, you'll need to create a mapping
    # between these IDs and the indices of atoms in the RDKit molecule.
    # This is because RDKit expects atom indices to start from 0 and
    # increment sequentially.
    # No worries! I have done this already. atomid2index is the mapping
    
    #######################################################################
    # auto-defining the bond type using unsupervised clustering algorithm
    #######################################################################
    if bo_analysis:
        bo_list = []
        for u,v,bond_order in moleculeGraph.edges(data='bond_order'):
            bo_list.append(bond_order)
        bo_array = np.array(bo_list).reshape(-1, 1)
        # Apply KMeans clustering
    
        kmeans = KMeans(n_clusters=n_clusters,n_init=10,random_state=0)
        kmeans.fit(bo_array)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        
        ref_bond_orders = [1.0,1.5,2.0]
        
        
        # Assign labels to clusters
        cluster_labels = {}
        for i, center in enumerate(centroids):
            # Find which reference bond order is closest to this cluster center
            closest_ref_label = np.argmin(np.abs(ref_bond_orders - center)) + 1
            cluster_labels[i] = closest_ref_label
        
        
        # Labeling each bond order
        new_labels = [cluster_labels[label] for label in kmeans.labels_]
        bond_order2type = dict(zip(bo_list,new_labels))
        
        # Visualize the clusters
        # Create scatter plots for each type
        colors = ['r','b','g','grey','orange','k','cyan']
        if plot_cluster:
            for i in range(n_clusters):
                plt.scatter(np.arange(bo_array.size)[labels == i],
                            bo_array[labels==i],s=60,edgecolor='k',
                            color=colors[i])
            plt.show()
        ######################################################################
        ######################################################################
        
        ## adding atoms
        atomid2index = {}
        for node, atom_type in moleculeGraph.nodes(data='atom_type'):
            atomicNumber = atomic_num[atom_type-1]
            if atomicNumber==1: # we are neglecting Hydrogen in SMILE
                continue
            atom_index = mol.AddAtom(Chem.Atom(atomicNumber))
            atomid2index[node] = atom_index
            
        ## adding bonds
        for u, v, bond_order in moleculeGraph.edges(data='bond_order'):
            bond_type = bond_order2type[bond_order]
            u_atom_type = moleculeGraph.nodes[u]['atom_type']
            v_atom_type = moleculeGraph.nodes[v]['atom_type']
            uan = atomic_num[u_atom_type - 1]
            van = atomic_num[v_atom_type - 1]
            if uan == 1 or van == 1:
                continue
            
            if bond_type == 1:
                rdkitBondType = Chem.BondType.SINGLE
            elif bond_type == 2:
                rdkitBondType = Chem.BondType.DOUBLE
            elif bond_type == 3:
                rdkitBondType = Chem.BondType.TRIPLE
            else:
                print('User ERROR: bond type is not recognized')
                print('Setting the bond type as single')
                rdkitBondType = Chem.BondType.SINGLE
                
            mol.AddBond(atomid2index[u], atomid2index[v], rdkitBondType)
    
    else: # if bo_analysis == False
        atomid2index = {}
        for node in moleculeGraph.nodes:
            atomicNumber = atomic_num[atom_types[node]-1]
            if atomicNumber==1: # we are neglecting Hydrogen in SMILE
                continue
            atom_index = mol.AddAtom(Chem.Atom(atomicNumber))
            atomid2index[node] = atom_index
            
        ## adding bonds
        for u, v in moleculeGraph.edges:
            u_atom_type = atom_types[u]
            v_atom_type = atom_types[v]
            uan = atomic_num[u_atom_type - 1]
            van = atomic_num[v_atom_type - 1]
            if uan == 1 or van == 1:
                continue
            mol.AddBond(atomid2index[u], atomid2index[v], Chem.BondType.SINGLE)
    
    mol = mol.GetMol()
    smiles = Chem.MolToSmiles(mol)
    return smiles

#%%
def get_bondtypes(bo,n_clusters,plot=False):   
    ## getting all bond orders in a np array
    ## bo is the bond order of of atomconnectivity of each frame
    
    bo_list = []
    for atom1, inner_dict in bo.items():
        for atom2, bond_order in inner_dict.items():
            if bond_order not in bo_list:
                bo_list.append(bond_order)
    
    bo_array = np.array(bo_list).reshape(-1, 1)
    # Apply KMeans clustering

    kmeans = KMeans(n_clusters=n_clusters,n_init=10)
    kmeans.fit(bo_array)
    
    # Re-labeling the cluster. min_bo=Label-1, mid_bo=Label-2, max_bo=Label-3
    # One problem. Only work if single-double-triple, this order maintained
    # For example: single-triple = won't work; double-triple = won't work
    # Only single = will work, single-double = will work
    labels = kmeans.labels_
    centroids = kmeans.cluster_centers_
    
    # Pair up cluster labels with their centroids and sort
    sorted_clusters = sorted(enumerate(centroids), key=lambda x: x[1])
    
    # Create a mapping from old to new labels
    label_mapping = {old_label: new_label+1 for new_label, (old_label, _) in enumerate(sorted_clusters)}
    
    # Apply the new labels to your data
    new_labels = [label_mapping[label] for label in labels]
    
    # Visualize the clusters
    # Create scatter plots for each type
    colors = ['r','b','g','grey','orange','k','cyan']
    if plot:
        for i in range(n_clusters):
            plt.scatter(np.arange(bo_array.size)[labels == i],
                        bo_array[labels==i],s=60,edgecolor='k',
                        color=colors[i])
        plt.show()
    
    return bo_array, np.array(new_labels)

#%%
def assign_speciesID(neighbors,life=False):
    speciesID = {}
    life_frame = {}

    ## asign a species id to each and every species
    uniqueID = 1
    for frame, (step, neigh) in enumerate(neighbors.items()):
        molecules = get_molecules(neigh)
        for molecule in molecules:
            frozen_molecule = frozenset(molecule)
            if frozen_molecule not in speciesID:
                speciesID[frozen_molecule] = uniqueID
                uniqueID+=1
                life_frame[frozen_molecule] = (frame,0)
            else:
                life_frame[frozen_molecule][1] = frame
    
    if life:
        return speciesID,life_frame
    else:
        return speciesID
#%%
# 3/14/2024
# Take the bondfile and return a dict: key=bond_index, value=(atom1,atom2)#bond
# molecule_length to scan only the particular molecule
@function_runtime
def map_isomer_bonds(bondfilepath, ref=None, molecule_length=None):
    # ref is the reference molecule. If None ref=first_molecule
    bonddata = parsebondfile(bondfilepath,firststep=True) #only the first step
    neigh    = bonddata['neighbours']
    atypes   = bonddata['atypes']
    molecules= get_molecules(neigh)
    
    # get list of molecule_graph from the first steps
    neigh_graph = nx.Graph(neigh) 
    mol_graphs  = []
    for molecule in molecules:
        if molecule_length is None or len(molecule)==molecule_length:
            mol_graph = neigh_graph.subgraph(molecule)
            mol_graphs.append(mol_graph)
    # ref is a nx.Graph
    if ref is None:
        ref=mol_graphs[0]
    
    #####
    isomer_bonds_list = []
    # isomer_bonds_list = [e1,e2,.....,en], n=#of molecule
    # e1 = [(e1_u1,e1_v1), (e1_u2,e1_v2), ......, #of bonds]
    # e2 = [(e2_u1,e2_v1), (e2_u2,e2_v2), ......]
    # where (e1_u1,e1_v1) and (e2_u1,e2_v1) are isomer bond
    
    ref_bonds = [(u, v) for u, v in ref.edges()]
    for mol_graph in mol_graphs:
        matcher = isomorphism.GraphMatcher(ref, mol_graph)
        if matcher.is_isomorphic():
            mapping = matcher.mapping
            isomer_bonds = [(mapping[u], mapping[v]) for u, v in ref_bonds]
            isomer_bonds_list.append(isomer_bonds)
        else:
            raise ValueError('Some molecule_graphs are not isomers')
    
    # returning the ref so that you can use this again to call this function    
    return isomer_bonds_list, ref

























