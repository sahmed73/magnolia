'''Created by Shihab Ahmed'''

##### Edit history #######
# Shihab Ahmed: 04:06:2024      | Define thermo_panda function

import pandas as pd

def thermo_panda(logfile, serial,
                 start_string = 'Per MPI',
                 end_string   = 'Loop time',
                 zero_ref : str = None):
    '''
    zero_ref: str
    Specifies which variables to apply zero-referencing to.
    Options include:
        - "energy": Apply zero-referencing to all energy (tot, pe, ke) columns.
        - "time": Apply zero-referencing to time values.
        - "energy+time" or "time+energy": Apply zero-referencing to both energy and time.
        - "user_defined_thermo_parameter": Apply zero-referencing.
        use '+' to use multiple options
    The values can be combined using an underscore (+) in any order.
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
                    print(f'timestep found from the log file: {timestep}')
    
    if timestep is None:
        timestep=1.0
        print("Warning: No 'timestep' found; set to 1.0 fs")
                
    if len(start_lines)!=len(end_lines):
        print('Warning: Log file is incomplete')
        end_lines.append(i)
        
    feed = list(zip(start_lines,end_lines))  
    
    if serial=='all':
        thermo_list = []  # Create an empty list to store dataframes
        for f in feed:
            start, end = f
            # Read each serial and append it to the list of thermos
            thermo_part = pd.read_csv(logfile, sep=r'\s+', skiprows=start+1,
                                      nrows=end-start-2)
            thermo_list.append(thermo_part)  # Append the read data to the list
        
        # Combine all the individual thermos into one dataframe
        thermo = pd.concat(thermo_list, ignore_index=True)
        
    elif isinstance(serial,int):
        if serial>len(feed):
            raise ValueError(f'Only {len(feed)} serial exists but found {serial}')
            
        start, end = feed[serial-1]
        thermo = pd.read_csv(logfile,sep=r'\s+',skiprows=start+1,
                           nrows=end-start-2)
        
    elif isinstance(serial,list):
        if any(s > len(feed) for s in serial):
            raise ValueError('Found a serial number greater than '
                             f'the total number of serials ({len(feed)}).')
        
        # check what if serials are not consecutive
        # under construction
        
        thermo_list = []  # Create an empty list to store dataframes

        for s in serial:
            start, end = feed[s-1]
            # Read each serial and append it to the list of thermos
            thermo_part = pd.read_csv(logfile, sep=r'\s+', skiprows=start+1,
                                      nrows=end-start-2)
            thermo_list.append(thermo_part)  # Append the read data to the list
        
        # Combine all the individual thermos into one dataframe
        thermo = pd.concat(thermo_list, ignore_index=True)
            
    else:
        raise TypeError("Serial must be a single int or an array of int or 'all'")
    
    # Time column
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
                        print(f"{z} in not in the thermo list")
            else:
                if z in thermo:
                    thermo[z] = thermo[z] - thermo[z].min() 
                else:
                    print(f"{z} in not in the thermo list")

    return thermo