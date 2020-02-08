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

'''Only works in Python 3 because certain packages statsmodels and collections'''
#George Mason University
#Jonathan Reed
#Feb. 7, 2020

import numpy as np
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
#import collections
from collections import Counter
from collections.abc import Mapping
#import math
#import sklearn
#from sklearn import linear_model

import argparse
import os
import glob
import re
import statsmodels.formula.api as smf

def forward_selected(data, response):
    """Linear model designed by forward selection.

    Parameters:
    -----------
    data : pandas DataFrame with all possible predictors and response

    response: string, name of response column in data

    Returns:
    --------
    model: an "optimal" fitted statsmodels linear model
           with an intercept
           selected by forward selection
           evaluated by adjusted R-squared
    """
    remaining = set(data.columns)
    remaining.remove(response)
    selected = []
    current_score, best_new_score = 0.0, 0.0
    while remaining and current_score == best_new_score:
        scores_with_candidates = []
        for candidate in remaining:
            formula = "{} ~ {} + 1".format(response,
                                           ' + '.join(selected + [candidate]))
            score = smf.ols(formula, data).fit().rsquared_adj
            scores_with_candidates.append((score, candidate))
        scores_with_candidates.sort()
        best_new_score, best_candidate = scores_with_candidates.pop()
        if current_score < best_new_score:
            remaining.remove(best_candidate)
            selected.append(best_candidate)
            current_score = best_new_score
    formula = "{} ~ {} + 1".format(response,
                                   ' + '.join(selected))
    model = smf.ols(formula, data).fit()
    return model

'''Function Definitions for Fit Equations'''
def func0(x,m,b):                                     
    return m*x+b

def func1(x,z,a,b):
    return b + a*(x**z)

def func_merge(x,m,b,z,a,c,x_max):
    return (b+m*x)*(x<x_max)+(c+a*(x**z))*(x>x_max)#+b (x2> max_val)

def func2(x,a,b):
    return b + a*(x**3)

def func3(x,a,b): 
    return b + a*np.log(x)

def func4(x,a,b):
    return b + a*np.log10(x)

def func5(a,x,c):
    return a*(1-np.exp(x/c))
    #return (a-b**(c/x))

def piecewise_linear(x, x0, k1, k2, b):#(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + b, lambda x:k2*x + b])
#want to send in list of possible functions to test at once:

'''Coding Functions'''

'''Will Flatten Nested List to Single (Non-Nested) List'''
def flatten(container):    
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

def save_png(png, title):
    png.savefig(title)
    print('File Created : ' + str(title))
    png.close()

'''Will Plot 2-D (Simple), 3-D, or Data Across Archives (Merge)'''
def plot(data,to,extras = None):                                                   
    plt.ion()  
    if to == 'Simple': #'2-D'                           #so i am thinking that if i send in requests, i would have to find correct key
        param_x = extras[0]; param_y = extras[1]        #in dictionary to plot...
        fig = plt.figure(figsize = (14,8))              #instead of just sending in the separated data...and then unwrapping it??
        for archive in data:
            plt.plot(data[archive][param_x], data[archive][param_y], 'o', label = archive)
        if len(extras) > 2:
            line_data = extras[2][0]; line_label = extras[2][1]  #could add correct x-y labels in extras instead 0 1
            plt.plot(line_data[0], line_data[1], color = 'orange', label = line_label)
        plt.xlabel(param_x)
        plt.ylabel(param_y)
        plt.legend()
        #save_png(plt, str(param_x) + ' vs ' + str(param_y) + '.png') 
        
    elif to == 'Merge':
        if extras[0] == 'Separate':                   #will plot separate archives for direct/indirect soma connection
            for connect in param_data:
                for param in extras[1]:
                    fig = plt.figure(figsize = (14,8))
                    for num,comp_type in enumerate(data[connect]):
                        ax1 = fig.add_subplot(len(data[connect].keys()),1,num+1)
                        for archive in data[connect][comp_type]:                  
                            ax1.plot(data[connect][comp_type][archive][param], data[connect][comp_type][archive]['RADIUS'],'o',label = archive)
                            ax1.set_title(comp_type)
                            plt.legend()
                    plt.xlabel(str(param))
                    plt.ylabel('Radius')
                    save_png(plt, str(connect) + ' ' + str(param) + ' vs ' + 'RADIUS' + '.png')

        if extras[0] == 'Combine':                    #will plot combined archives for direct/indirect soma connection
            direct = {}; indirect = {}
            for comp_type in extras[2]:
                direct[comp_type] = {}; indirect[comp_type] = {}
                for param in extras[1]:
                    direct[comp_type][param] = []; indirect[comp_type][param] = []
                    for archive in extras[3]:
                        direct[comp_type][param].append(param_data['Direct'][comp_type][archive][param])
                        indirect[comp_type][param].append(param_data['Indirect'][comp_type][archive][param])
                    direct[comp_type][param] = list(flatten(direct[comp_type][param]))
                    indirect[comp_type][param] = list(flatten(indirect[comp_type][param]))
            
            for param in extras[1]:
                fig = plt.figure(figsize = (14,8))
                for num,comp_type in enumerate(extras[2]):
                    ax1 = fig.add_subplot(len(extras[2]),1,num+1)
                    ax1.plot(direct[comp_type]['PARENT_RAD'], direct[comp_type][param], 'o', label = 'Direct')
                    ax1.plot(indirect[comp_type]['PARENT_RAD'], indirect[comp_type][param], 'o', label = 'Indirect')
                    ax1.set_title(comp_type)
                plt.legend()
                plt.xlabel('Parent Radius (um)')
                plt.ylabel(str(param))
                save_png(plt, 'Combined'  + ' ' + str(param) + ' vs ' + 'Parent Radius' + '.png')
                        
    elif to == '3d':                                  #will plot 3d view of parameters to radius
        combinations = []
        for x in extras:
            for y in extras:
                if x != y and x != 'RADIUS' and y != 'RADIUS':
                    temp = [x,y]
                    if not combinations:
                        combinations.append(temp)
                    if combinations:
                        if temp not in combinations and temp[::-1] not in combinations:
                            combinations.append(temp)
                            
        for connect in data:
            for comb in combinations:
                fig = plt.figure(figsize = (14,8))
                for num,comp_type in enumerate(data[connect]):
                    #fig = plt.figure(figsize = (14,8))
                    ax1 = fig.add_subplot(len(data[connect].keys()),1,num+1, projection = '3d')
                    for archive in data[connect][comp_type]:
                        dict_loc = data[connect][comp_type][archive]             #may want to select Residual instead as option
                        ax1.scatter(dict_loc[comb[0]],dict_loc[comb[1]],dict_loc['RADIUS'],'o', label = archive)
                    ax1.set_xlabel(comb[0])
                    ax1.set_ylabel(comb[1])
                    ax1.set_zlabel('RADIUS')
                    ax1.set_title(comp_type)
                    plt.legend()
                save_png(plt, connect + ' ' + comb[0] + ' ' + comb[1] + ' ' + 'RADIUS' + '.png')

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
    #may want to fit multiple x values for possible multivariante equation...
    
    func = fit_func['func']; x_max = fit_func['x_max']; x_min = fit_func['x_min']; x_range = fit_func['x_range']; line_label = fit_func['line']
    comp_type = labels['comp_type'];
    #if len(labels['param_x']) > 1:
        #for num, param in eumerate(labels['param_x']):
            #how to organize the multiple x's recieved here....
    param_x = labels['param_x']; param_y = labels['param_y']; keys = labels['params']
    counts = {}; select_data = {}; select_data['in'] = []; select_data['out'] = []; outbound = {}

    '''Splits Data if Boundaries from Function Call'''
    for archive in data[comp_type]:                
        count = []; test_count = 0
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
            elif x_max != None and x_min == None:
                if val > x_max:
                    count.append(num)
                    select_data['out'].append(point)
                else:
                    select_data['in'].append(point)
            elif x_min != None and x_max == None:
                if val < x_min:
                    count.append(num)
                    select_data['out'].append(point)
                else:
                    select_data['in'].append(point)
                    
            elif x_max == None and x_min == None:
                test_count = test_count + 1
                select_data['in'].append(point)
                  
        counts[archive] = count

    for data_set in select_data.keys():                    #reorganizes data into xlist and ylist of values
        select_data[data_set] = list(zip(*select_data[data_set])) #can probably relabel for param_x and param_y sent into function
        select_data[data_set] = [list(i) for i in select_data[data_set] if i]

    if counts[archive]:                                          
        #plot({'values':[select_data['out'][0],select_data['out'][1]]}, 'Simple', {'labels':labels})
        plot(select_data, 'Simple', [0, 1]) #can leave this here: plots data that is outside of bounds
        temp_dict = {}
        temp_dict[comp_type] = {}; outbound[comp_type] = {}         #if values found outside of bounds
        for archive in data[comp_type]:                             #save to new dictionary with same structure
            temp_dict[comp_type][archive] = {}
            outbound[comp_type][archive] = {}
            for num,param in enumerate(list(zip(keys))):
                outbound[comp_type][archive][param[0]] = []
                temp_dict[comp_type][archive][param[0]] = []
                for num,val in enumerate(data[comp_type][archive][param[0]]):
                    if num in counts[archive]:
                        outbound[comp_type][archive][param[0]].append(val)
                    else:
                        temp_dict[comp_type][archive][param[0]].append(val) #update sent data to include only 'in' bounds
        data = temp_dict
        
    else:
        del select_data['out']
        
    '''Calculate Variables for Curve Fitting'''
    popt, pcov = optimize.curve_fit(func,select_data['in'][0],select_data['in'][1])  #find values of function variables
    print('popt',popt)
    print('pcov',pcov)
    
    #plot({'values':[select_data['in'][0],select_data['in'][1]],'line':[x_range,func(x_range, *popt)]}, 'Simple', {'labels':labels, 'line_label':line_label})
    plot(select_data,'Simple', [0,1,[[x_range,func(x_range, *popt)],line_label]])
    
    
    '''Calculate Residuals within Dictionary Organization'''
    res_label = 'Residual_' + str(param_x)
    for archive in data[comp_type]:
        prediction = []
        data[comp_type][archive][res_label] = []
        for x in data[comp_type][archive][param_x]:
            predict_y = func(x,*popt)
            prediction.append(predict_y)
        residuals = []
        for y_val, y_predict in list(zip(data[comp_type][archive][param_y],prediction)):
            residuals.append(y_val - y_predict)
        data[comp_type][archive][res_label] = residuals
    '''
    xlist = []; ylist = []
    labels['param_x'] = 'Radius'
    labels['param_y'] = 'Residuals'
    for archive in data[comp_type]:
        for val in zip(data[comp_type][archive][param_y], data[comp_type][archive][res_label]):
            xlist.append(val[0])
            ylist.append(val[1])
    '''
    plot(data[comp_type], 'Simple', ['RADIUS', res_label])

    for select_param in data[comp_type][archive].keys():
        plot(data[comp_type], 'Simple', [select_param, res_label])

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
    basal_list = list(zip(*basal_list))                           #zip shifts column/row to call all parameter data by parameter number as tuple
    archive_dict['Basal'][d1] = basal_list
    if apical_list:
        if not 'Apical' in archive_dict.keys():
            archive_dict['Apical'] = {}                         
        apical_list = list(zip(*apical_list))
        archive_dict['Apical'][d1] = apical_list

'''Plot Parameters to Radius For Initial Comparison'''
params = {key:header[key] for key in header if key not in str_list}   #new dictionary to only contain useful parameters

'''Separate Data by Connection to Soma (either 0 -> directly connected, 1 -> all others)'''
param_data = {}                          #contains useful parameter values (as defined by params) as nested dictionary

param_data['Direct'] = {}                #Directly connected to soma
param_data['Indirect'] = {}              #All other nodes (not soma points)

#Setup Dictionaries before sorting data
for connect in param_data:
    for comp_type in archive_dict:
        param_data[connect][comp_type] = {}
        for archive in archive_dict[comp_type]:
            param_data[connect][comp_type][archive] = {}
            for param in params.keys():
                param_data[connect][comp_type][archive][param] = []

#Fill Dictionaries by Connection to Soma (as 0, others as 1) for all parameters
for comp_type in archive_dict:
    for archive in archive_dict[comp_type]:
        for num,val in enumerate(archive_dict[comp_type][archive][header['DIRECT_SOMA']]):
            for param in params.keys():
                if val == 0:
                    param_data['Direct'][comp_type][archive][param].append(archive_dict[comp_type][archive][header[param]][num])
                else:
                    param_data['Indirect'][comp_type][archive][param].append(archive_dict[comp_type][archive][header[param]][num])

#plot(param_data, 'Merge', ['Separate', params.keys()]) #will plot based on archives and connection to soma
#plot(param_data, 'Merge', ['Combine', params.keys(), param_data['Direct'], dirs]) #will plot based on comp_type and connection to soma

#plot(param_data,'3d', params.keys())

'''Parameter Correlation (Pearson's r) to Radius For Linear Relationships'''
initial_corr = {}
excluded = ['DIRECT_SOMA','NODE_ORDER','NODE_COUNT']
for i in param_data:
    if i == 'Direct':
        initial_corr[i] = corr(param_data[i], [i for i in params.keys() if i not in excluded])
    else:
        initial_corr[i] = corr(param_data[i], [i for i in params.keys() if i != 'DIRECT_SOMA'])

with open('initial_corr_dict.txt', 'a') as outfile:
    for i in initial_corr:
        for j in initial_corr[i]:
            outfile.write(str(i) + ' - ' + str(j) + '\n')
            outfile.write('\n')
            outfile.write(str(initial_corr[i][j]) + '\n')
            outfile.write('\n')
    outfile.close()              
    
#plot(#transformedx, transformedy)

'''Testing Linear Regression Model from sklearn'''
'''OldNebish cannot access module sklearn --> works on Nebish but cannot seem to install package; maybe requires super-user access'''
#clf = linear_model.LinearRegression()
#clf.fit([[getattr(t, 'x%d' % i) for i in range(len(params.keys()))] for archive in param_data[comp_type]], [archive.y for archive in param_data[comp_type]]) #is this the list of data values?
'''Initiate Fit and Residuals'''

#Apical_PR, Apical_PRout, Apical_PRval, Apical_PRline = fit(param_data, {'comp_type':'Apical', 'param_x':'PARENT_RAD', 'param_y':'RADIUS','params':params.keys()},{'func':func0, 'x_range':np.arange(0,7,1), 'x_min':None, 'x_max':7, 'line':'y= mx + b'})

#Apical_PR, Apical_PRout, Apical_PRval, Apical_PRline = fit(param_data, {'comp_type':'Apical', 'param_x':'PARENT_RAD', 'param_y':'RADIUS','params':params.keys()},{'func':func_merge, 'x_range':np.arange(0,7,1), 'x_min':None, 'x_max':7, 'line':'y = mx if(x1<x_max) + ax**z if(x2>x_max) + b'})
'''
xmax = 0
for archive in dirs:
    temp = max(param_data['Apical'][archive]['BRANCH_LEN'],key = lambda item:item)
    #temp = np.max(param_data['Apical'][archive]['BRANCH_LEN'])
    if temp > xmax:
        xmax = temp
'''

'''Testing Pandas Dataframe for Stepwise Regression'''
'''
for param in params.keys():
    temp = []
    for archive in param_data['Apical']:
        temp.append(param_data['Apical'][archive][param])
    temp = list(flatten(temp))
    temp_dict[param] = temp
'''
'''
temp_dict = {}
for param in params.keys():   #converts to DataFrame
    temp = []
    for comp_type in param_data.keys():   #lose archive organization or could save number at which append happens...
        if not comp_type in temp_dict.keys():
            temp_dict[comp_type] = {}
        for archive in param_data[comp_type]:
            temp.append(param_data[comp_type][archive][param])
        temp = list(flatten(temp))
        temp_dict[comp_type][param] = temp
'''
'''
Framed = {}
#memory error...
for comp_type in temp_dict.keys():
    Framed[comp_type] = forward_selected(pd.DataFrame(temp_dict[comp_type]), 'RADIUS')
#try single archives at a time
'''
#try stepwise regression --> add features at a time and keep if R is increased
#plot residuals to remaining parameters

#temp_frame = pd.DataFrame(temp_dict)

#model = forward_selected(temp_frame, 'RADIUS')
#print(model.model.formula)
#print(model.rsquared_adj)

'''It looks like all parameters would be considered significant to estimate diameter...'''

#Try and plot only the (0) and (1) for Direct Soma to see if there is a significant difference between their supposed equations
'''Data Transformations'''
#transform = ['param_log', 'param_exp', 'param_sqrt']

#Apical_PR, Apical_PRout, Apical_PRval, Apical_PRline = fit(param_data, {'comp_type':'Apical', 'param_x':'PARENT_RAD', 'param_y':'RADIUS','params':params.keys()},{'func':func0, 'x_range':np.arange(0,7,1), 'x_min':None, 'x_max':7, 'line':'y = mx + b'})

#possibly loop over the parameters to the new residual and possible pick best function --> careful for bounds
#Apical_PR_P, Apical_PRout_P, Apical_PRval_P, Apical_PRline_P = fit(Apical_PR, {'comp_type':'Apical', 'param_x':'PATH_TO_END', 'param_y':'Residual_PARENT_RAD','params':Apical_PR['Apical']['Groen'].keys()},{'func':func0, 'x_range':np.arange(0,250,1), 'x_min':None, 'x_max':250, 'line':'y = mx + b'})
'''
#could probably make function to determine max/min number within data list
xmax = 0
for archive in dirs:
    temp = max(param_data['Apical'][archive]['BRANCH_LEN'],key = lambda item:item)
    #temp = np.max(param_data['Apical'][archive]['BRANCH_LEN'])
    if temp > xmax:
        xmax = temp
'''
#Apical_BL, Apical_BLout, Apical_BLval, Apical_BLline = fit(param_data, {'comp_type':'Apical', 'param_x':'BRANCH_LEN', 'param_y':'RADIUS','params':params.keys()},{'func':func1, 'x_range':np.arange(0,xmax,1), 'x_min':0, 'x_max':None, 'line':'b + a*(x**z)'})

#Apical_BL, Apical_BLout, Apical_BLval, Apical_BLline = fit(param_data, {'comp_type':'Apical', 'param_x':'BRANCH_LEN', 'param_y':'RADIUS','params':params.keys()},{'func':func2, 'x_range':np.arange(0,xmax,1), 'x_min':0, 'x_max':None, 'line':'b + a*(x**3)'})

#Apical_BL, Apical_BLout, Apical_BLval, Apical_BLline = fit(param_data, {'comp_type':'Apical', 'param_x':'BRANCH_LEN', 'param_y':'RADIUS','params':params.keys()},{'func':func3, 'x_range':np.arange(0,xmax,1), 'x_min':0, 'x_max':None, 'line':'b + a*np.log(x)'})

#Apical_BL, Apical_BLout, Apical_BLval, Apical_BLline = fit(param_data, {'comp_type':'Apical', 'param_x':'BRANCH_LEN', 'param_y':'RADIUS','params':params.keys()},{'func':func4, 'x_range':np.arange(0,xmax,1), 'x_min':0, 'x_max':None, 'line':'b + a*np.log10(x)'})

#Apical_BL, Apical_BLout, Apical_BLval, Apical_BLline = fit(param_data, {'comp_type':'Apical', 'param_x':'BRANCH_LEN', 'param_y':'RADIUS','params':params.keys()},{'func':func5, 'x_range':np.arange(0,xmax,1), 'x_min':0, 'x_max':None, 'line':'a*(1-exp(x/c))'})

#Apical_PR, Apical_PRout, Apical_PRval, Apical_PRline = fit(param_data, {'comp_type':'Apical', 'param_x':'PARENT_RAD', 'param_y':'RADIUS','params':params.keys()},{'func':func_merge, 'x_range':np.arange(0,7,1), 'x_min':None, 'x_max':7, 'line':'y = (b+m*x)*(x<x_max)+(c+a*(x**z))*(x>x_max)'})

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

#THINGS TO TRY
#Separating the Direct Soma points (0) and (1) to see if they have separate equations to estimate radius
#transform the data to see any nonlinear relationships (if linear after the transformation)
#Check to see why direct_soma nodes have differing values of Soma Radius...(or their parent radius)

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

