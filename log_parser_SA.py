'''Created by Shihab Ahmed'''

##### Edit history #######
# Shihab Ahmed: 04:06:2024      | Define thermo_panda function
# Shihab Ahmed: 02:21:2025      | Updated serial input handling to support 'start:end' notation

import pandas as pd

def thermo_panda(logfile, serial,
                 start_string = 'Per MPI',
                 end_string   = 'Loop time',
                 zero_ref : str = None,
                 verbose: bool = True,
                 warn: bool = True):
    '''
    Parses LAMMPS log file to extract thermodynamic data.
    
    Parameters:
    logfile : str
        Path to the LAMMPS log file.
    serial : int, list, or str
        Specifies which blocks of thermo to extract (1-based):
        - int: Extract a single serial.
        - list: Extract multiple serials.
        - 'start:end': Extract a range of serials.
        - ':': Extract all serials.
        - 'start:': Extract from start to last.
        - ':end': Extract from first to end.
    start_string : str, default='Per MPI'
        String that marks the start of a thermo data block.
    end_string : str, default='Loop time'
        String that marks the end of a thermo data block.
    zero_ref : str, optional
        Specifies which variables to apply zero-referencing to.
        Options include:
            - "energy": Apply to all energy (TotEng, PotEng, KinEng).
            - any thermo parameter e.g. 'Time'
            - Use '+' to combine multiple options, e.g. 'Time+energy'.
    
    Returns:
    pandas.DataFrame
        Thermodynamic data extracted from the log file.
    '''
    
    timestep=None # starting None value
    start_lines  = []
    end_lines    = []
    with open(logfile, 'r') as file:
        for i, line in enumerate(file):
            if start_string in line:
                start_lines.append(i)
            if end_string in line:
                end_lines.append(i)
            
            # getting timestep
            if line.strip().startswith('timestep'):
                try:
                    timestep=float(line.split('#')[0].split()[1])
                except:
                    pass
                else:
                    if verbose:
                        print(f'timestep found from the log file: {timestep}')
    
    if verbose:
        print(f'Total number of serials: {len(start_lines)}')
    
    if timestep is None:
        timestep=1.0
        if warn: print("Warning: No 'timestep' found; set to 1.0 fs")
                
    if len(start_lines)!=len(end_lines):
        if warn: print('Warning: Log file is incomplete')
        end_lines.append(i)
        
    feed = list(zip(start_lines,end_lines))  
    
    if isinstance(serial,str):
        if ':' in serial:
            start, end = serial.split(':')
            start = 1 if start == '' else int(start)
            end = len(feed) if end == '' else int(end)
            
            # negative indexing
            if start<0: start = len(feed)+start+1
            if end<0: end=len(feed)+end+1
            
            # range error
            if start<1 or start>len(feed):
                raise ValueError(f"Serial {serial.split(':')[0]} out of range!")
            if end<1 or end>len(feed):
                raise ValueError(f"Serial {serial.split(':')[1]} out of range!")
            
            serial = list(range(start,end+1))
            if verbose: print(f'Serial asked for {serial}')
        else:
            serial = int(serial)
    
    if isinstance(serial,int):      
        if abs(serial)>len(feed):
            raise ValueError(f"Serial {serial} out of range!")
            
        if serial<0: serial = len(feed)+serial+1
        if verbose: print(f'Serial asked for {serial}')
            
        start, end = feed[serial-1]
        thermo = pd.read_csv(logfile,sep=r'\s+',skiprows=start+1,
                           nrows=end-start-2)
        
    elif isinstance(serial,list):
        if any(s > len(feed) or s < 0 for s in serial):
            raise ValueError(f"Serial {serial} out of range!")
        
        thermo_list = []  # Create an empty list to store dataframes

        for s in serial:
            start, end = feed[s-1]
            thermo_part = pd.read_csv(logfile, sep=r'\s+', skiprows=start+1,
                                      nrows=end-start-2)
            thermo_list.append(thermo_part)  
        
        thermo = pd.concat(thermo_list, ignore_index=True)
            
    else:
        raise TypeError("Serial must be an int, a list of ints, or a string in the format 'start:end'")
    
    # Keep only rows where all values are numeric, then drop Nan
    thermo = thermo.apply(pd.to_numeric, errors='coerce')
    thermo = thermo.dropna()
    
    ps = thermo['Step']*timestep/1000
    thermo['Time'] = ps
        
    if zero_ref:
        zrefs = zero_ref.split('+')
        for z in zrefs:
            if z.lower()=='energy':
                energy_columns = ['PotEng', 'TotEng', 'KinEng']
                for col in energy_columns:
                    if col in thermo.columns:
                        thermo[col] = thermo[col] - thermo[col].min()
                    else:
                        if warn: print(f"Warning: {z} in not in the thermo list")
            else:
                if z in thermo:
                    thermo[z] = thermo[z] - thermo[z].min() 
                else:
                    if warn: print(f"Warning: {z} in not in the thermo list")

    return thermo
