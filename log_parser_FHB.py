'''Created by Fakhrul Hasan Bhuiyan'''

##### Edit history #######

# Fakhrul H Bhuiyan: 10.13.2020 | Finalized the first version of the module
# Fakhrul H Bhuiyan: 10.28.2020 | Fixed issues: Same thermodynamic keyword as header
#                               | Incomplete line in the data block
#                               | Added a new flag: 'Per MPI rank memory allocation'
# Shihab Ahmed: 04:06:2024      | Define thermo_panda function

from warnings import showwarning, warn
import numpy as np
from scipy.optimize import curve_fit
import sys
import pandas as pd

def thermo_dict(filename, serial):
    with open(filename, 'r') as f:
        start = 0
        thermo_data = {}

        for step, line in enumerate(f):
            #print (line)
            if line.startswith('Memory usage per processor') or line.startswith('Per MPI rank memory allocation'):   # looks for this line in the file
                #print('start')
                start += 1      # seed for start taking in data
                data_row = 0    # this is equal to the line number of the data
                thermo_headers = []     # thermo information in log files
                continue
            if start == serial:     # if the code reaches the desired portion of the log file
                data_row += 1
                data = line.split()
                if data_row == 1:   # in the first line you have the headers
                    for i in (data):
                        if i not in thermo_headers:
                            thermo_headers.append(i)    # creating a list of str of headers
                            thermo_data[i] = []         # Creating dictionary keys named after each thermodynamic property
                        else:
                            # There are identical thermodynamic keywords as header
                            name = i + '_' + str(data.index(i) + 1)
                            thermo_headers.append(name)
                            thermo_data[name] = []
                            showwarning('Log file has same Thermodynamic keyword ' + str(i),
                                        UserWarning, 'log_parser_FHB.py', 31)

                    #print ('len', len(thermo_data), thermo_data)

                elif data_row > 1:      # this is where the data information starts in the log file
                    #print (start)
                    for i in range(0, len(data)):
                        try:
                            data[i] = float(data[i])  # making the str to float

                            if len(data) != len(thermo_headers):    # checking if the data line is complete
                                showwarning('Line number: ' + str(step+1) + ' in the log file is incomplete.',
                                            UserWarning, 'log_parser.py', 46)
                                #print (data)
                                break
                            else:
                                thermo_data[thermo_headers[i]].append(
                                    data[i])  # putting in the data info in respective header's dict
                        except:
                            start += 1
                            #print ('breaking')
                            break  # if data lines end, you get str again, this terminates taking any more data




    #print (thermo_data['Step'])
    ##print (thermo_data['KinEng'])
    #for x in thermo_data: print (x)
    f.close()
    
    from colorama import Fore, Style
    print(Fore.CYAN + '#### From log_parser: Serial = {}, Total Steps = {} ####'.format(serial,len(thermo_data['Step'])),end='')
    print(Style.RESET_ALL)


    return thermo_data      # Returns the whole dictionary


##############
#   Created by Shihab from here
#   Date: 4/18/2023
##############



def thermo_dict_v2(filename, serial):               
    ## Local Function That can be called ##
    def thermo_serial(filename,serial):
        with open(filename, 'r') as f:
            start = 0
            thermo_data = {}
    
            for step, line in enumerate(f):
                #print (line)
                if line.startswith('Memory usage per processor') or line.startswith('Per MPI rank memory allocation'):   # looks for this line in the file
                    #print('start')
                    start += 1      # seed for start taking in data
                    data_row = 0    # this is equal to the line number of the data
                    thermo_headers = []     # thermo information in log files
                    continue
                if start == serial:     # if the code reaches the desired portion of the log file
                    data_row += 1
                    data = line.split()
                    if data_row == 1:   # in the first line you have the headers
                        for i in (data):
                            if i not in thermo_headers:
                                thermo_headers.append(i)    # creating a list of str of headers
                                thermo_data[i] = []         # Creating dictionary keys named after each thermodynamic property
                            else:
                                # There are identical thermodynamic keywords as header
                                name = i + '_' + str(data.index(i) + 1)
                                thermo_headers.append(name)
                                thermo_data[name] = []
                                showwarning('Log file has same Thermodynamic keyword ' + str(i),
                                            UserWarning, 'log_parser_FHB.py', 31)
    
                        #print ('len', len(thermo_data), thermo_data)
    
                    elif data_row > 1:      # this is where the data information starts in the log file
                        #print (start)
                        for i in range(0, len(data)):
                            try:
                                data[i] = float(data[i])  # making the str to float
    
                                if len(data) != len(thermo_headers):    # checking if the data line is complete
                                    showwarning('Line number: ' + str(step+1) + ' in the log file is incomplete.',
                                                UserWarning, 'log_parser.py', 46)
                                    #print (data)
                                    break
                                else:
                                    thermo_data[thermo_headers[i]].append(
                                        data[i])  # putting in the data info in respective header's dict
                            except:
                                start += 1
                                #print ('breaking')
                                break  # if data lines end, you get str again, this terminates taking any more data
        return thermo_data

    
    if type(serial) in [list,tuple]:
        result = {}
        for s in serial:
            thermo = thermo_serial(filename,s)
            
            for key,value in thermo.items():
                if key in result.keys():
                    result[key] +=thermo[key]
                else:
                    result[key] = thermo[key]
         
    else:
        result = thermo_serial(filename, serial)
            
            
    from colorama import Fore, Style
    print(Fore.CYAN + '#### From log_parser: Serial = {}, Total Steps = {} ####'.format(serial,len(result['Step'])),end='')
    print(Style.RESET_ALL)


    return result      # Returns the whole dictionary


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
        raise ValueError("No 'timestep' value is found")
                
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
            
        print('Code is under construction..')
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
                        
            elif z.lower()=='time':
                if 'Time' in thermo.columns:
                    thermo['Time'] = thermo['Time'] - thermo['Time'].min()
            
            else:

                warn(f"The 'zero_ref' parameter '{z}' is not recognized. "
                    "Please use 'energy', 'time', or a combination of them.",
                    UserWarning
                )
                print()

    return thermo