# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  5 20:01:53 2025
"""

import pandas as pd
import re
import periodictable

def parse_lammps_data(filepath):
    atom_lines = []
    bond_tuples = []
    mass_lines = []
    
    with open(filepath, 'r') as file:
        lines = file.readlines()

    section = None
    for line in lines:
        line = line.strip()

        if not line or line.startswith("#"):
            continue  # skip empty/comment lines

        # Detect section headers
        if "Masses" in line:
            section = "masses"
            continue
        elif "Atoms" in line:
            section = "atoms"
            continue
        elif "Bonds" in line:
            section = "bonds"
            continue
        elif line.endswith("Coeffs") or "Velocities" in line:
            section = None  # skip non-relevant sections
            
        # Parse relevant sections
        if section == "masses":
            if line[0].isdigit():
                parts = re.split(r'\s+', line, maxsplit=2)
                atom_type = int(parts[0])
                mass = float(parts[1])
                mass_lines.append([atom_type, mass])
            else: # blank lines are skipped at the top so no worries
                section = None
                
        elif section == "atoms":
            if line[0].isdigit():
                parts = line.split()
                atom_id = int(parts[0])
                mol_id   = int(parts[1])
                atom_type = int(parts[2])
                x, y, z = map(float, parts[3:6])
                atom_lines.append([atom_id, mol_id, atom_type, x, y, z])
            else: # blank lines are skipped at the top so no worries
                section = None

        elif section == "bonds":
            if line[0].isdigit():
                parts = line.split()
                atom1 = int(parts[2])
                atom2 = int(parts[3])
                bond_tuples.append((atom1, atom2))
            else: # blank lines are skipped at the top so no worries
                section = None

    # Convert atoms to DataFrame
    atom_df = pd.DataFrame(atom_lines, columns=['id', 'mol-id', 'type', 'x', 'y', 'z'])
    mass_df = pd.DataFrame(mass_lines, columns=['type', 'mass'])
    mass_df['symbol'] = mass_df['mass'].apply(get_symbol_from_mass)

    return mass_df, atom_df, bond_tuples



def get_symbol_from_mass(mass, tol=0.01):
    closest = None
    min_diff = float("inf")

    for element in periodictable.elements:
        if element.number == 0:
            continue  # skip "None" element
        diff = abs(element.mass - mass)
        if diff < tol and diff < min_diff:
            min_diff = diff
            closest = element

    if closest:
        return closest.symbol
    else:
        return None
