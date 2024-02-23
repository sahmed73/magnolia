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
import magnolia.bondfile_parser as bfp # only for bond_to_speciecfile

# This Function Convert .xsd file to LAMMPS Data (.data) formate.



def xsd2data(xsdfile, atomtype_order=None):
    # This function takes xsd file directory+xsd_filename as first argument.
    # It takes the element list as the optional second argument

    if xsdfile[xsdfile.rfind('.'):] != '.xsd':
        sys.exit('Input must be a .xsd file!!')

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
    
def bond_to_speciecfile(bondfile, atomsymbols, cutoff):
    df = bfp.get_species_count(bondfile,atomsymbols,cutoff=cutoff)
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
            line = string+"\t".join(non_zero_species.index) + "\n"
            # print(line,end='')
            sf.write(line)
            # Writing the species counts line
            
            species_counts = "\t".join(map(str, non_zero_species.values))
            line = f"{header}\t{no_moles}\t{no_species}\t{species_counts}\n"
            # print(line,end='')
            sf.write(line)
    