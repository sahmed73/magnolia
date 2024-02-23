'''Created by Fakhrul Hasan Bhuiyan'''

##### Edit history #######

# Fakhrul H Bhuiyan: 10.13.2020 | Finalized the first version of the module
# Fakhrul H Bhuiyan: 10.28.2020 | Fixed issues: Same thermodynamic keyword as header
#                               | Incomplete line in the data block
#                               | Added a new flag: 'Per MPI rank memory allocation'

from warnings import showwarning
import numpy as np
from scipy.optimize import curve_fit
import sys

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

def tempramp(thermo,timestep,**kwargs):
    #getting kwargs
    Print = kwargs.get('Print',False)
    
    
    from magnolia.bondfile_parser import step2picosecond
    def ramp_fit(x,m,c):
        y = m*x+c
        return y
    
    steps = np.array(thermo['Step'])
    ps = np.array(step2picosecond(steps, timestep))
    temp = np.array(thermo['Temp'])
    p, cov = curve_fit(ramp_fit, ps, temp)
    if Print:
        print('Ramping Rate (m):',round(p[0],2))
        print('Initial Temperature (c):',round(p[1],2))
    return p


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