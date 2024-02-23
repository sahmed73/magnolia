# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:24:40 2023
last checked on 12/10/2023

@author: arup2
"""
import pandas as pd

def get_species_count(speciesfile):
    with open(speciesfile,'r') as sf:
        species = {}
        for line in sf:
            if "Timestep" in line:
                headers       = line.strip().split()[1:] # not taking '#'
                species_name  = headers[3:] # skip Timestep, No_moles, No_Specs
            else:
                values        = [int(x) for x in line.strip().split()]
                timestep      = values[0]
                species_count = values[3:]
                species[timestep]={}
                for key,value in zip(species_name,species_count):
                    if key in species:
                        species[timestep][key]+=value
                    else:
                        species[timestep][key] = value
    
    df = pd.DataFrame(species).fillna(0)
    df.index.name = 'Timestep'
    return df