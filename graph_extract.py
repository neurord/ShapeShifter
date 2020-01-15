'''<Extracts data values from .swc files>
    Copyright (C) <2019>  <Jonathan Reed>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.'''

#Usage:  python -i graph_extract.py --path /name/of/path 
#        where path is name of path to directory containing archive directories
#Output: morphology_extract.py

#George Mason University
#Jonathan Reed
#Jan 15, 2020

import numpy as np
import btmorph
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter

import argparse
import os
import glob
import re

'''Function Definitions for Fit Equations'''
def func0(x,m,b):                                     
    return m*x+b

def func1(x,z,a,b):
    return b + a*(x**z)

def func_merge(x1,m,b,x2,z,a,max_val):
    return m*x1*(x1<max_val)+a*(x2**z)+b#(x2> max_val)

'''
def func2(x,a,b):
    return b + a*(x**3)

def func3(x,a,b): 
    return b + a*np.log(x)

def func4(x,a,b):
    return b + a*np.log10(x)
'''

#want to send in list of possible functions to test at once:

#may want to combine plot functions together and call separately with optional variable

#simply have to send in arguments of label of data intended for plot

'''Coding Functions'''

'''Will Flatten Nested List to Single (Non-Nested) List'''
def flatten(container):    
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

'''Will Plot 2-D (Simple), 3-D, or Data Across Archives (Merge)'''             
def plot(data,to,extras = None):                                                   #extras will be dictionary to call by key i.e. params will be list of params.keys()
    if to == 'Simple': #'2-D'
        fig = plt.figure()
        plt.plot(data['values'][0],data['values'][1], 'o')
        if extras:
            comp_type = extras['labels']['comp_type']; param_x = extras['labels']['param_x']; param_y = extras['labels']['param_y']
            if len(data) > 1:
                line_label = extras['line_label']     
                plt.plot(data['line'][0],data['line'][1], color = 'orange', label = line_label)
            plt.legend()
            plt.xlabel(param_x)
            plt.ylabel(param_y)
            #print(extras.keys())
            #plt.title(labels[2])
        plt.show()
        #plt.savefig(title + '.png')                                                             
        #print('File Created : ' + str(title) + '.png')
        
    elif to == 'Merge':  
        for param in params:
            fig = plt.figure()#figsize = (14,8))
            for num,comp_type in enumerate(data):
                ax1 = fig.add_subplot(len(data.keys()),1,num+1)
                for archive in data[comp_type]:                   #may want to select Residual instead as option
                    ax1.plot(data[comp_type][archive][param], data[comp_type][archive]['RADIUS'],'o',label = archive)
                    ax1.set_title(comp_type)
                    plt.legend()
            plt.xlabel(str(param))
            plt.ylabel('Radius')
            #title = param + ' to RADIUS'
            #plt.title(title)
            #plt.savefig(title + '.png')                      
            #print('File Created : ' + str(title) + '.png')
            #plt.close()
        plt.show()
        
    elif to == '3d':
        combinations = []
        for x in params:
            for y in params:
                if x != y and x != 'RADIUS' and y != 'RADIUS':
                    temp = [x,y]
                    if not combinations:
                        combinations.append(temp)
                    if combinations:
                        if temp not in combinations and temp[::-1] not in combinations:
                            combinations.append(temp)
        for comb in combinations:
            fig = plt.figure(figsize = (14,8))
            for num,comp_type in enumerate(data):
                #fig = plt.figure(figsize = (14,8))
                ax1 = fig.add_subplot(len(data.keys()),1,num+1, projection = '3d')
                for archive in data[comp_type]:
                    dict_loc = data[comp_type][archive]             #may want to select Residual instead as option
                    ax1.scatter(dict_loc[comb[0]],dict_loc[comb[1]],dict_loc['RADIUS'],'o', label = archive)
                ax1.set_xlabel(comb[0])
                ax1.set_ylabel(comb[1])
                ax1.set_zlabel('RADIUS')
                ax1.set_title(comp_type)
                plt.legend()
            title = comb[0] + ' ' + comb[1] + ' ' + 'RADIUS'
            #title = comp_type + ': ' + comb[0] + ' ' + comb[1] + ' ' + 'RADIUS'
            #plt.title(title)
            plt.savefig(title + '.png')                      
            print('File Created : ' + str(title) + '.png')
            #plt.show()
            plt.close()

'''Pearson's r Correlation for Feature Relationships'''

def corr(data,keys): 
    corr_dict = {}
    for comp_type in data.keys():
        corr_dict[comp_type] = pd.DataFrame(index = keys, columns = keys)
        for rparam in keys:
            for cparam in keys:
                combx = []; comby = []
                for archive in data[comp_type].keys():                   
                    combx.append(list(data[comp_type][archive][rparam])) 
                    comby.append(list(data[comp_type][archive][cparam])) 
                combx = list(flatten(combx))
                comby = list(flatten(comby))
                corr_dict[comp_type].loc[rparam,cparam] = round(pearsonr(combx,comby)[0],4)
    return corr_dict

'''Descrete Functions (top) to Selected Data; Outputs Variables and Indicates Possible Fit Fine-Tuning'''

def fit(data,labels,fit_func):
    
    func = fit_func['func']; x_max = fit_func['x_max']; x_min = fit_func['x_min']; x_range = fit_func['x_range']; line_label = fit_func['line']
    comp_type = labels['comp_type']; param_x = labels['param_x']; param_y = labels['param_y']; keys = labels['params']
    counts = {}; select_data = {}; select_data['in'] = []; select_data['out'] = []; outbound = {}

    '''Splits Data if Boundaries from Function Call'''
    for archive in data[comp_type]:                
        count = []
        for num,val in enumerate(data[comp_type][archive][param_x]):
            x = data[comp_type][archive][param_x][num]
            y = data[comp_type][archive][param_y][num]
            point = [x,y]
            if x_max != None and x_min != None:              
                if val > x_max and val < x_min:                
                    count.append(num)                          
                    select_data['out'].append(point)
                else:
                    select_data['in'].append(point)
            elif x_max != None:
                if val > x_max:
                    count.append(num)
                    select_data['out'].append(point)
                else:
                    select_data['in'].append(point)
            elif x_min != None:
                if val < x_min:
                    count.append(num)
                    select_data['out'].append(point)
                else:
                    select_data['in'].append(point)
                    
            elif x_max == None and x_min == None:
                select_data['in'].append(point)
                    
        counts[archive] = count

    for data_set in select_data.keys():                    #reorganizes data into xlist and ylist of values
        select_data[data_set] = zip(*select_data[data_set])
        select_data[data_set] = [list(i) for i in select_data[data_set]]
        
    if counts[archive]:                                          
        plot({'values':[select_data['out'][0],select_data['out'][1]]}, 'Simple', {'labels':labels})
        temp_dict = {}#;outbound = {}
        temp_dict[comp_type] = {}; outbound[comp_type] = {}         #if values found outside of bounds
        for archive in data[comp_type]:                             #save to new dictionary with same structure
            temp_dict[comp_type][archive] = {}
            outbound[comp_type][archive] = {}
            for num,param in enumerate(zip(keys)):
                outbound[comp_type][archive][param[0]] = []
                temp_dict[comp_type][archive][param[0]] = []
                for num,val in enumerate(data[comp_type][archive][param[0]]):
                    if num in counts[archive]:
                        outbound[comp_type][archive][param[0]].append(val)
                    else:
                        temp_dict[comp_type][archive][param[0]].append(val) #update sent data to include only 'in' bounds
        data = temp_dict

    '''Calculate Variables for Curve Fitting'''
    popt, pcov = optimize.curve_fit(func,select_data['in'][0],select_data['in'][1])  #find values of function variables
    print('popt',popt)
    print('pcov',pcov)

    plot({'values':[select_data['in'][0],select_data['in'][1]],'line':[x_range,func(x_range, *popt)]}, 'Simple', {'labels':labels, 'line_label':line_label})

    '''Calculate Residuals within Dictionary Organization'''
    res_label = 'Residual_' + str(param_x)
    for archive in data[comp_type]:
        prediction = []
        data[comp_type][archive][res_label] = []
        for x in data[comp_type][archive][param_x]:
            predict_y = func(x,*popt)
            prediction.append(predict_y)
        residuals = []
        for y_val, y_predict in zip(data[comp_type][archive][param_y],prediction): 
            residuals.append(y_val - y_predict)
        data[comp_type][archive][res_label] = residuals

    xlist = []; ylist = []
    labels['param_x'] = 'Radius'
    labels['param_y'] = 'Residuals'
    for archive in data[comp_type]:
        for val in zip(data[comp_type][archive][param_y], data[comp_type][archive][res_label]):
            xlist.append(val[0])
            ylist.append(val[1])

    plot({'values': [xlist,ylist]}, 'Simple', {'labels':labels})

    return data, outbound, select_data, [line_label,popt,pcov]

'''Start of Working Code'''
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   

args = parser.parse_args()                      
path = args.path                                

'''Locates and Organizes Archive Data'''
root, dirs, files = list(os.walk(path))[0]                    #all _extract files require identical parameters
archive_dict = {}; header = {}                                
str_list = ['PARENT','X','Y','Z','CHILD','TYPE','NUM_ENDS']   #parameters not useful or unique for Radius comparision

for d1 in dirs:
    fullpath = path + d1 + '/*CNG_extract.txt'
    data_list = []; apical_list = []; basal_list = []
    print('Working on Directory : ', str(d1))
    for fname in glob.glob(fullpath):                               #locates and loops through _extract files within fullpath
        with open(fname) as f:
            for line in f:
                if line.strip():                                    #removes any empty lines
                    if '*C' in line and not header:                 #finds header line and separately saves from data
                        if 'XYZ' in line:                           
                            line = re.sub('XYZ','X; Y; Z', line)    #find and fix header for each file correctly
                        line = line.split('; ')                     
                        for num, val in enumerate(line):            
                            if '\n' in val:
                                val = re.sub('\n','',val)
                            if val == '*CHILD':
                                header['CHILD'] = num
                            else:
                                header[val] = num
                    elif line[0] != '*' and line[0] != ' ' and line[0] != '/n':  #data values    
                        temp_line = line.split()
                        if len(temp_line) == len(header):
                            for point, val in enumerate(temp_line):        #Child, Type, Parent use DESCRETE numerical labelling from .swc
                                if point != header['CHILD'] and point != header['TYPE'] and point != header['PARENT']:
                                    temp_line[point] = float(temp_line[point])
                            data_list.append(temp_line)                   

    for line in data_list:
        if line[header['TYPE']] == '4':                     #splits by node type (Apical or Basal if present)
            apical_list.append(line)                        #currently recognizes only Apical and Basal types 
        elif line[header['TYPE']] == '3':
            basal_list.append(line)

    if not 'Basal' in archive_dict.keys():                  #archive_dict will hold all possible .swc values as list within dictionary
        archive_dict['Basal'] = {}                          #'Basal' and 'Apical' (if present) as dictionaries within archive_dict
    basal_list = zip(*basal_list)                           #zip shifts column/row to call all parameter data by parameter number as tuple
    archive_dict['Basal'][d1] = basal_list
    if apical_list:
        if not 'Apical' in archive_dict.keys():
            archive_dict['Apical'] = {}                         
        apical_list = zip(*apical_list)                        
        archive_dict['Apical'][d1] = apical_list

'''Plot Parameters to Radius For Initial Comparison'''
params = {key:header[key] for key in header if key not in str_list}   #new dictionary to only contain useful parameters

param_data = {}                                                       #param_data will hold useful parameter values as nested dictionary
for comp_type in archive_dict:
    param_data[comp_type] = {}
    for archive in archive_dict[comp_type]:
        param_data[comp_type][archive] = {}
        for param in params:
            param_data[comp_type][archive][param] = archive_dict[comp_type][archive][header[param]]

#plot(param_data,'Merge', {'params':params.keys()}) 

#plot(param_data,'3d', {'params':params.keys()})

'''Parameter Correlation (Pearson's r) to Radius For Linear Relationships'''
initial_corr = corr(param_data, params.keys())

#could try and transform the data to 'linearize' relationships


'''Initiate Fit and Residuals'''

Apical_PR, Apical_PRout, Apical_PRval, Apical_PRline = fit(param_data, {'comp_type':'Apical', 'param_x':'PARENT_RAD', 'param_y':'RADIUS','params':params.keys()},{'func':func0, 'x_range':np.arange(0,7,1), 'x_min':None, 'x_max':7, 'line':'y= mx + b'})

#for parent_rad maybe use btmorph to see if points are directly connected to soma
#would have to keep filenames and unique compartment order maintained if so...


'''
>>> temp_dict = {}
>>> temp_dict['labels'] = {}
>>> temp_dict['labels']['comp_type'] = 'Apical'
>>> temp_dict['labels']['param_x'] = 'PATH_TO_END'
>>> temp_dict['labels']['param_y'] = 'RADIUS'
>>> temp_dict = {'labels': {'comp_type': 'Apical', 'param_y': 'RADIUS', 'param_x': 'PATH_TO_END'}}
{'labels': {'comp_type': 'Apical', 'param_y': 'RADIUS', 'param_x': 'PATH_TO_END'}}
>>> plot({'values':[Apical_PRout['Apical']['Groen']['PATH_TO_END'], Apical_PRout['Apical']['Groen']['RADIUS']]}, 'Simple', temp_dict)
'''

#maybe check remaining data points within function
#huh how to combine residuals to remaining parameters. may have to maintain archives when pushing to function.
#plot_fit(param_data['Apical']['PARENT_RAD'], param_data['Apical']['RADIUS'], func0, np.arange(0,7,1), x_max = 7) #y-x labels send in key for x,y
    #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func1, np.arange(0,np.max(param_data['Apical']['BRANCH_LEN']),1))
    #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func2, np.arange(0,np.max(param_data['Apical']['BRANCH_LEN']),1))
    #find_min = param_data['Apical']['BRANCH_LEN']  
    #find_min = list(set(find_min))  #for some reason we have to remove duplicates before sorting...
    #find_min.sort()
    #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func4, np.arange(find_min[1],np.max(param_data['Apical']['BRANCH_LEN']),1), x_max = None, x_min = np.min(param_data['Apical']['BRANCH_LEN']))
    #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func4, np.arange(find_min[1],np.max(param_data['Apical']['BRANCH_LEN']),1), x_max = None, x_min = np.min(param_data['Apical']['BRANCH_LEN']))
    #I wonder if I can make a function to fit optimize for certain parameters with good correlation (linear and nonlinear) values


#Additional Types of Plotting Data
#plot mulitple comp type (if present) with single archive


#TODO
#Implement Random Forest with parameters ideal for Radius estimation (want to lower the amount of features needed BUT still have an accurate equation for estimate)

