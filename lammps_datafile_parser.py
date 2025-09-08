# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Apr  5 20:01:53 2025
"""

import pandas as pd
import re
import periodictable
from pathlib import Path

def read_lammps_full_datafile(filepath, box_size = False):
    atom_lines = []
    bond_tuples = []
    mass_lines = []
    
    if box_size: box = []
    
    with open(filepath, 'r') as file:
        lines = file.readlines()

    section = None
    from periodictable import elements
    for line in lines:
        line = line.strip()

        if not line or line.startswith("#"):
            continue  # skip empty/comment lines
        
        # Box size
        if box_size:
            if "xlo" in line and "xhi" in line:
                stripped = line.split()
                box.append((float(stripped[0]), float(stripped[1])))
            
            if "ylo" in line and "yhi" in line:
                stripped = line.split()
                box.append((float(stripped[0]), float(stripped[1])))
                
            if "zlo" in line and "zhi" in line:
                stripped = line.split()
                box.append((float(stripped[0]), float(stripped[1])))
                
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
                # Split into numbers and optional comment
                if "#" in line:
                    main, comment = line.split("#", 1)
                    comment = comment.strip()
                else:
                    main, comment = line, ""
                parts = re.split(r'\s+', main.strip())
                atom_type = int(parts[0])
                mass = float(parts[1])
                # Guess symbol from mass
                symbol = ""
                min_diff = float('inf')
                for el in elements:
                    try:
                        diff = abs(el.mass - mass)
                        if diff < min_diff and diff < 1.0:  # 1 amu tolerance
                            symbol = el.symbol
                            min_diff = diff
                    except Exception:
                        continue
                mass_lines.append([atom_type, mass, symbol, comment])
            else:
                section = None

        elif section == "atoms":
            if line[0].isdigit():
                parts = line.split()
                atom_id = int(parts[0])
                mol_id   = int(parts[1])
                atom_type = int(parts[2])
                charge    = float(parts[3])
                x, y, z = map(float, parts[4:7])
                atom_lines.append([atom_id, mol_id, atom_type, charge, x, y, z])
            else:
                section = None

        elif section == "bonds":
            if line[0].isdigit():
                parts = line.split()
                atom1 = int(parts[2])
                atom2 = int(parts[3])
                bond_tuples.append((atom1, atom2))
            else:
                section = None

    # Convert to DataFrames
    atom_df = pd.DataFrame(atom_lines, columns=['id', 'mol-id', 'type', 'charge', 'x', 'y', 'z'])
    mass_df = pd.DataFrame(mass_lines, columns=['type', 'mass', 'symbol', 'comment'])
    
    if box_size:
        return box, mass_df, atom_df, bond_tuples
    
    return mass_df, atom_df, bond_tuples


def write_modified_datafile(original_path, atom_df, fill_missing_fields=True):
    original_path = Path(original_path)
    modified_path = original_path.with_name(original_path.stem + "_modified" + original_path.suffix)

    with open(original_path, 'r') as f:
        lines = f.readlines()

    # Recover ix, iy, iz, comment from original Atoms section if requested
    if fill_missing_fields:
        original_data = {}
        in_atoms = False
        for line in lines:
            if line.strip().startswith("Atoms"):
                in_atoms = True
                continue
            if in_atoms:
                if not line.strip():
                    continue
                if line.strip()[0].isalpha():
                    break  # end of Atoms section
                raw, *comment_part = line.split("#", 1)
                parts = raw.strip().split()
                if len(parts) < 10:
                    continue  # skip malformed line
                atom_id = int(parts[0])
                original_data[atom_id] = {
                    "ix": int(parts[7]),
                    "iy": int(parts[8]),
                    "iz": int(parts[9]),
                    "comment": comment_part[0].strip() if comment_part else ""
                }

        # Fill missing fields if not in atom_df
        for col in ["ix", "iy", "iz", "comment"]:
            if col not in atom_df.columns:
                atom_df[col] = atom_df["id"].map(lambda i: original_data.get(i, {}).get(col, "" if col == "comment" else 0))

    # Write the modified data file
    output_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith("Atoms"):
            output_lines.append(line)  # "Atoms # full"
            i += 1
            while i < len(lines) and (lines[i].strip() == "" or lines[i].strip().startswith("#")):
                output_lines.append(lines[i])
                i += 1
            # Skip old Atoms section data
            while i < len(lines) and lines[i].strip() and lines[i].strip()[0].isdigit():
                i += 1

            # Write new Atoms section
            for _, row in atom_df.iterrows():
                formatted = (
                    f"{int(row['id']):>5d} {int(row['mol-id']):>5d} {int(row['type']):>4d} "
                    f"{row['charge']:>10.6f} {row['x']:>15.9f} {row['y']:>15.9f} {row['z']:>15.9f} "
                    f"{int(row['ix']):>2d} {int(row['iy']):>2d} {int(row['iz']):>2d}"
                )
                if pd.notna(row["comment"]) and row["comment"]:
                    formatted += f"  # {row['comment']}"
                output_lines.append(formatted + "\n")
            output_lines.append("\n")
            continue

        output_lines.append(line)
        i += 1

    with open(modified_path, "w") as f:
        f.writelines(output_lines)

    print(f"Modified Atoms section written to: {modified_path.name}")
    return modified_path



    

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
