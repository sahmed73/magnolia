# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:24:40 2023
last checked on 12/10/2023

@author: arup2
"""
import pandas as pd

def get_species_count(speciesfile, timestep=None, restart_time=False,
                      Nfreq=1, first_timestep=False):
    with open(speciesfile,'r') as sf:
        species = {}
        for line in sf:
            if "Timestep" in line:
                headers       = line.strip().split()[1:] # not taking '#'
                species_name  = headers[3:] # skip Timestep, No_moles, No_Specs
            else:
                values        = [int(x) for x in line.strip().split()]
                step      = values[0]
                species_count = values[3:]
                species[step]={}
                for key,value in zip(species_name,species_count):
                    if key in species:
                        species[step][key]+=value
                    else:
                        species[step][key] = value
                        
                if first_timestep: break
    
    df = pd.DataFrame(species).fillna(0).T
    df.index.name = 'Timestep'
    
    if timestep is not None:
        df['Time'] = df.index*timestep/1000
    
    if restart_time:
        df['Time'] = df['Time'] - df['Time'].min()
    
    return df[::Nfreq]