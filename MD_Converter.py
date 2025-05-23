# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 19:42:53 2023

@author: Shihab
"""
# Install 'periodictable' module using:
#    conda install periodictable (Anaconda User)
#    pip install periodictable (Non-Anaconda User)

from periodictable import elements
import sys
import warnings
import magnolia.bondfile_parser as bfp # for bond2speciecfile
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# This Function Convert .xsd file to LAMMPS Data (.data) formate.
# The output lammps data file I used for reaxff simulation
# lmpdata file contain:
    # Number of atoms
    # Simulation box dimension
    # Masses
    # atom_id mol_id atom_type charge x y z
# Charges generally set to zero. Because reaxff calculate charges
    
def xsd2lmpdata(xsdfile, atomtype_order=None):
    # This function takes xsd file directory+xsd_filename as first argument.
    # It takes the element list as the optional second argument
    
    # I tried it with xsd file, generated by Modules > Amorphous Cell
    ##  If you draw a molecule in BIoVIA Material Studio
        # and export it as xsd, it will not work
    
    if xsdfile[xsdfile.rfind('.'):] != '.xsd':
        raise TypeError('Input must be a .xsd file!!')

    lammps_datafile = xsdfile[:xsdfile.rfind('.')]+'.data'

    xsd = open(xsdfile, 'r')
    datafile = open(lammps_datafile, 'w')
    #################################
    atom_types = {}
    atom_symbol = []
    atomic_mass = []
    pos = []
    atom_id = []
    molecule_id = []
    mol = 1

    for line in xsd.readlines():
        if line.find('Molecule ID=') != -1:
            mol += 1

        if line.find('AVector') != -1:
            # X-Lenght of the simulation box
            start_id = line.find('AVector')+9
            end_id = line.find('"', start_id)
            X_lenght = float(line[start_id:end_id].split(',')[0])

        if line.find('BVector') != -1:
            # Y-Lenght of the simulation box
            start_id = line.find('BVector')+9
            end_id = line.find('"', start_id)
            Y_lenght = float(line[start_id:end_id].split(',')[1])

        if line.find('CVector') != -1:
            # Z-Lenght of the simulation box
            start_id = line.find('CVector')+9
            end_id = line.find('"', start_id)
            Z_lenght = float(line[start_id:end_id].split(',')[2])

        if line.find('Atom3d ID') != -1:
            # molecule id
            molecule_id.append(mol)

            # position
            start_id = line.find('XYZ')
            if start_id != -1:
                end_id = line.index('"', start_id+5)
                XYZ = list(map(float, line[start_id+5:end_id].split(',')))
                pos.append(XYZ)

            # atom type
            start_id = line.find('Components')
            if start_id != -1:
                end_id = line.index('"', start_id+12)
                element = line[start_id+12:end_id]
                atom_symbol.append(element)

    # De-Normalize the Coordinates
    box_length = [X_lenght, Y_lenght, Z_lenght]
    for i in range(len(pos)):
        for j in range(len(pos[i])):
            pos[i][j] = pos[i][j]*box_length[j]

    # Assign Atom Type
    unique_atom_symbol = list(set(atom_symbol))
    if atomtype_order:
        if set(atom_symbol) != set(atomtype_order):
            warnings.warn(
                'atomtype_order does not match with the elements in the .xsd file\nAssigning default atom type...')
        else:
            unique_atom_symbol = atomtype_order
    atom_types  = dict(zip(unique_atom_symbol,range(1,1+len(unique_atom_symbol))))

    # arranege data for printing purpose
    Masses = []
    for i in range(len(unique_atom_symbol)):
        element = unique_atom_symbol[i]
        atomic_mass = elements.symbol(element).mass
        Masses.append((atom_types[element],atomic_mass,element))
        
    Masses_print = []
    
    for item in Masses:
        f = '{}   {} # {}'.format(*item)
        Masses_print.append(f)
    
    header = '''# LAMMPS data file converted from Material Studio .xsd file # Created by SA
{}  atoms
{}  atom types
0.000  {:0.3f}   xlo xhi
0.000  {:0.3f}   ylo yhi
0.000  {:0.3f}   zlo zhi

Masses
'''.format(len(pos), len(unique_atom_symbol), *box_length)
    
    print('Masses\n\n')
    datafile.writelines(header+'\n')
    #printing Masses
    for line in Masses_print:
        print(line)
        datafile.writelines(line+'\n')
    
    datafile.writelines('\nAtoms   #full\n\n')
    
    atom_id = list(range(1,len(pos)+1))
    
    ####Error###
    if not len(pos)==len(atom_id)==len(molecule_id)==len(atom_symbol):
        sys.exit('ERROR: number of coordiates, atom id, molecule id and atom symbols are not same!!\nDETAILS:\n  Number of coordinates: {}\n  Number of atom id:{}\n  Number of molecule id:{}\n  Number of atom symbol:{}'.format(len(pos),len(atom_id),len(molecule_id),len(atom_symbol)))
    ###########

    for ID, molid, sym, p in zip(atom_id, molecule_id, atom_symbol, pos):
        sep = '\t'
        charge = 0
        text = list(map(str, [ID, molid,atom_types[sym], charge, *p]))
        wl = sep.join(text)+'\n'
        datafile.writelines(wl)
    xsd.close()
    datafile.close()

# it convert the fix reaxff/bonds output to fix reaxff/species output
# reaxff/bonds example: fix 1 all reaxff/bonds 100 bonds.reaxff
# reaxff/species example: fix 1 all reaxff/species 10 10 100 species.out    
def bond2speciecfile(bondfile, atomsymbols, cutoff):
    df = bfp.get_species_count(bondfile,atomsymbols,cutoff=cutoff).T
    speciesfile = bondfile[:bondfile.rfind('\\')]+f'\\species_{cutoff}.out'
    with open(speciesfile,'w') as sf:
        for header in df.columns:
            column = df.loc[:,header]
            # Count the number of species with non-zero values
            column_as_int = column.round(0).astype(int)
            non_zero_species = column_as_int[column_as_int > 0]
            no_moles = non_zero_species.sum()
            no_species = non_zero_species.count()        
            
            # Writing the header line
            string = "# Timestep\tNo_Moles\tNo_Specs\t"
            # line = string+"\t".join(non_zero_species.index) + "\n"
            line = string + "\t".join(str(x) for x in non_zero_species.index) + "\n"
            # print(line,end='')
            sf.write(line)
            # Writing the species counts line
            
            species_counts = "\t".join(map(str, non_zero_species.values))
            line = f"{header}\t{no_moles}\t{no_species}\t{species_counts}\n"
            # print(line,end='')
            sf.write(line)


def pdb2lmpdata(pdbfile, atomtype_order=None, box_size=None, box_incr=None,
                adjust_coord=True):
    # atom_types is a dictionary key=symbol, value=atom_id
    # box_size=[0., 48., 0., 48., 0., 48.] # xlo, xhi, ylo, yhi, zlo, zhi
    # if box_size is None then it will assign the min, max coord as dimension
    # box_incr=[0., 2. , 0., 2.,  0. , 2.] # example
    # adjust_coord True: atom coord trans so that min atom pos change to 0
    
    '''
    COLUMNS        DATA TYPE       FIELD         DEFINITION (short)
    ----------------------------------------------------------------------------
      1 - 6        Record name     "ATOM  "
      7 -11        Integer         Atom serial number (atom_id)
    13 - 16        Atom            Atom name (atom_name)
    17             Character       Alternate location indicator
    18 - 20        Residue name    Residue name (res_name)
    22             Character       Chain identifier
    23 - 26        Integer         Residue sequence number (res_id)
    27             AChar           Code for insertion of residues
    31 - 38        Real(8.3)       Orthogonal coordinates for X (X)
    39 - 46        Real(8.3)       Orthogonal coordinates for Y (Y)
    47 - 54        Real(8.3)       Orthogonal coordinates for Z (Z)
    55 - 60        Real(6.2)       Occupancy
    61 - 66        Real(6.2)       Temperature factor
    77 - 78        LString(2)      Element symbol (sym)
    79 - 80        LString(2)      Charge on the atom (charge)
    '''
    
    start = None
    with open(pdbfile, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                end = i
                if start is None:
                    start = i
    
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22),
                (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60),
                (60, 66), (76, 78), (78, 80)]
    
    names = ['Record name', 'atom_id', 'atom_name',
             'Alternate location indicator',
             'res_name', 'Chain identifier', 'res_id',
             'Insertion code', 'X', 'Y', 'Z', 'Occupancy',
             'Temperature factor', 'sym', 'charge']

    # Read the PDB file
    df_pdb = pd.read_fwf(pdbfile, colspecs=colspecs, header=None, names=names,
                     skiprows=start, nrows=end-start+1)
    
    lmpdata_header = ['atom_id', 'res_id', 'sym', 'charge', 'X', 'Y', 'Z']
    df_lmpdata = df_pdb[lmpdata_header].fillna(0)
    
    # Cahrge handling: convert 1- charges to -1
    def convert_charge(charge):
        if isinstance(charge, str) and charge[-1] in ['-', '+']:
            return charge[::-1]
        return charge
    df_lmpdata['charge'] = df_lmpdata['charge'].apply(convert_charge)
    
    # atom type mapping
    if atomtype_order is None: # alphabatcal order
        atomtype_order = sorted(df_lmpdata['sym'].unique())
    
    atom_type_map = dict(zip(atomtype_order,range(1,1+len(atomtype_order))))
    
    
    df_lmpdata['sym'] = df_lmpdata['sym'].map(atom_type_map)
    df_lmpdata.rename(columns={'sym': 'atom_type'}, inplace=True)
    
    if box_size is None:
        xlo, ylo, zlo = np.round(df_lmpdata[['X','Y','Z']].min())
        xhi, yhi, zhi = np.round(df_lmpdata[['X','Y','Z']].max())
        if box_incr is not None:
            temp_box = np.array([xlo, xhi, ylo, yhi, zlo, zhi])
            xlo, xhi, ylo, yhi, zlo, zhi = temp_box+np.array(box_incr)
    else:
        if box_incr is None:  
            xlo, xhi, ylo, yhi, zlo, zhi = box_size
        else:
            xlo, xhi, ylo, yhi, zlo, zhi = [x+y for x,y in zip(box_size,box_incr)]
            
    # translate the atoms to the middle
    if box_incr is not None:
        xmin, xmax, ymin, ymax, zmin, zmax = box_incr
        trans = np.array([(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2])
        df_lmpdata[['X', 'Y', 'Z']] = df_lmpdata[['X', 'Y', 'Z']] + trans
        
    # adjust coordinate: make min coordinate 0 using translation
    if adjust_coord:
        df_lmpdata['X']-=xlo
        df_lmpdata['Y']-=ylo
        df_lmpdata['Z']-=zlo
        
        temp_box = np.array([xlo, xhi, ylo, yhi, zlo, zhi])
        trans    = np.array([xlo,xlo,ylo,ylo,zlo,zlo])
        # update xlo, xhi, ylo, yhi, zlo, zhi
        xlo, xhi, ylo, yhi, zlo, zhi = temp_box-trans
    
    n_atom_types = len(atom_type_map)
    natoms, _    = df_lmpdata.shape
    
    header = f'''# LAMMPS data file converted PDB file # Created by SA
    
{natoms}  atoms
{n_atom_types}  atom types

{xlo}  {xhi}   xlo xhi
{ylo}  {yhi}   ylo yhi
{zlo}  {zhi}   zlo zhi

Masses

'''
    masses = ''
    for symbol, atom_type in atom_type_map.items():
        element = elements.symbol(symbol)
        mass = f'{atom_type} {element.mass} # {element}\n'
        masses+=mass
    
    print(header+masses)
    output = pdbfile[:pdbfile.rfind('.')]+'.data' 
    string = df_lmpdata.to_string(header=None, index=None)
    with open(output, 'w') as out:
        out.write(header)
        out.write(masses+'\n')
        out.write('Atoms\n\n')
        out.write(string)
    return df_pdb

def smiles2pdb(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_with_h = Chem.AddHs(mol)

    # Generate 3D coordinates for the molecule with hydrogens
    AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_with_h)

    # Convert the molecule to a PDB format
    pdb_block = Chem.MolToPDBBlock(mol_with_h)
    
    return pdb_block

def write_xyzfile(filename, atom_symbols, positions):
    if os.path.exists(filename):
        print('A xyz file with the same name already exists. Skipping!')
        return
    
    lines=[f"{len(atom_symbols)}\n"]
    positions = np.array(positions)
    for atom, position in zip(atom_symbols, positions):
        lines.append(f"{atom} {position[0]:.6f} {position[1]:.6f} {position[2]:.6f}")
    
    # Step 4: Write to a file
    with open(filename, "w") as f:
        f.write("\n".join(lines))
    
    print("XYZ file generated.")
