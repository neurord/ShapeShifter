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

#Usage:  python graph_extract.py --path /name/of/path --parameters PATH_TO_END BRANCH_LEN
#        where path is name of path to directory containing archive directories
#        where parameters requried as capitalized in code
#Output: morphology_extract.py

#George Mason University
#Jonathan Reed
#Oct 31, 2019

import numpy as np
#import scipy
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
import pickle
from collections import Counter

import argparse
import os
import glob
import re

#import sympy as sy
def func(x,m,b): 
    return m*x+b

def func1(x,z,a,b):
    return b + a*(x**z)

def func2(x,a,b):
    return b + a*(x**3)

def func3(x,a,b): #odd error where we are dividing by 0 (even when I think I am not sending in 0 values here)
    return b + a*np.log(x)

def func4(x,a,b):
    return b + a*np.log10(x)

'''
def makefunc(num,form):
    func = 'func' + num
    #extract the equation variables as list
    
    def func(#send in separated list values, but varied in len...)
    
for num,form in enumerate(list):
makefunc(num,form)

#this would initiate the 
#want to send in list of possible functions to test at once:
'''
#may want to combine plot functions together and call separately with optional variable
def plot_label(xlabel, ylabel, title):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(title + '.png')                       #works better but still missing the basal graphs if both basal and apical data exist in particular dictionary
    #plt.clf()                                         #more explicit way to refer to figures
    plt.close()                                       #close window or remove figure from window...
    #plt.show()
    print('File Created : ' + str(title) + '.png')

def plot_3d(archive_dict, comp_type, header, x, y, z):
    fig = plt.figure(figsize = (14,8))
    ax = fig.add_subplot(111, projection='3d')
    for archive in archive_dict[comp_type]:
        ax.scatter(archive_dict[comp_type][archive][header[x]],archive_dict[comp_type][archive][header[y]],archive_dict[comp_type][archive][header[z]],'o', label = archive)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_zlabel(z)
    plt.title(comp_type)                        #i wonder if I can save multiple views of the 3d plots since we cannot manipulate outside of python --> trying pickle
    plt.legend()                                                                                                     #seems to print everything, not individual files
    title = comp_type + ' : ' + str(x) + ' ' + str(y) + ' ' + str(z)                                          
    plt.savefig(title + '.png')                       #works better but still missing the basal graphs if both basal and apical data exist in particular dictionary
    #plt.savefig(title + '.pickle')
    #pickle.dump(fig, open(title + '.pickle', 'wb'))
    print('File Created : ' + str(title) + '.png')
    #print('File Created : ' + str(title) + '.pickle')
    #plt.show()
    plt.close() 

def plot_fit(xdata, ydata, fit_func, x_line, x_max = None, x_min = None):
    #do fit work here for each of the plots we are focusing on
    data_list = []
    if x_max != None and x_min != None:
        data_list = [val for val in zip(xdata,ydata) if val[0] < x_max and val[0] > x_min]
    elif x_max != None:
        data_list = [val for val in zip(xdata,ydata) if val[0] < x_max]
    elif x_min != None:
        data_list = [val for val in zip(xdata,ydata) if val[0] > x_min]
    elif x_max == None and x_min == None:
        data_list = [val for val in zip(xdata,ydata)]

    data_list = zip(*data_list) #reorganizes to list of x and y values within data_list

    popt, pcov = optimize.curve_fit(fit_func,data_list[0],data_list[1])         #may have new issue as no longer optimizing for parameter estimation...
    print('popt',popt)
    print('pcov',pcov)
    
    plt.figure(figsize = (14,8))
    plt.plot(data_list[0],data_list[1], 'o', label = 'Data')   
    plt.plot(x_line, fit_func(x_line, *popt), color = 'orange', label = 'Fitted Function')   #plots fitted function with data it represents
    plt.legend()                                                                             #x_line is element sent in to be simple x-values for function (0,x_max,1)
    plt.show() 

    prediction = []          #predicts values for all xdata for residual calculation
    '''
    for x in xdata:          #may need specific range within x_min to x_max
        if x_max != None and x_min != None:
            predict_y = fit_func(x, *popt) if x < x_max and x > x_min else 0 
        elif x_max != None:
            predict_y = fit_func(x, *popt) if x < x_max else 0    #if x out of range of accepted values prediction = 0
        elif x_min != None:                                               #may need to change if particular data with BOTH x_max and x_min
            predict_y = fit_func(x, *popt) if x > x_min else 0
        elif x_max == None and x_min == None:
            predict_y = fit_func(x, *popt)
        prediction.append(predict_y)
    '''
    for x in data_list[0]:                        #only potential issue is we are not plotting all possible points (since log 0 giving error)
        predict_y = fit_func(x, *popt)            #plotting shortened data_list[0], not all original xdata sent to function
        prediction.append(predict_y)
        
    residuals = []
    for y_val, y_predict in zip(data_list[1],prediction):
        residuals.append(y_val - y_predict)
    plt.figure(figsize = (14,8))
    #plt.plot([i[1] for i in data_list], residuals, 'o') #possible outlier, why would single value of parent radius > 40??
    plt.plot(data_list[0],residuals,'o')
    #plt.title('Residuals')
    #plt.plot(x_line, func(x_line, popt[0], popt[1]), color = 'orange', label = 'Fitted Function')
    plt.show()
    
    
['m*x+b','b + a*(x**3)','b + a*np.log(x)','b + a*np.log10(x)']
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   
parser.add_argument("--graphing", type = str, nargs = '+', choices = {'archives','comp_type','3d','fit'})
parser.add_argument("--view") #, type = str, choices = {'3d'})

args = parser.parse_args()                      
path = args.path                                
graphing = args.graphing
view = args.view
''' 
if view:
    fullpath = path + '/*.pickle'
    for fname in glob.glob(fullpath):
        with open(fname, 'rb') as file:
            figx = pickle.load(file)
            figx.show()
            #figx = pickle.load(open(fname, 'rb'))
            #figx.close()
'''
'''
    for param in params:
        for num,comp_type in enumerate(archive_dict):
    #Axes3D.plot_surface()look for seqeuntial color map
'''

'''Extracts data from .CNG_extract.swc files moved to archive directories under main brain region directory (i.e. 'Groen' and 'Jaffe' directories under Hippocampus)'''
if graphing:
    root, dirs, files = list(os.walk(path))[0]                           #make sure all extracted files have same parameters before graphing                                            
    archive_dict = {}; header = {}; params = {}                          #params will be filled once str_list params are removed from header.keys()
    str_list = ['PARENT','X','Y','Z','CHILD','TYPE','NUM_ENDS']                     #identifiers of particular compartments not useful for Radius comparision


    for d1 in dirs:
        fullpath = path + d1 + '/*CNG_extract.txt'
        data_list = []; apical_list = []; basal_list = []
        print('Working on Directory : ', str(d1))
        for fname in glob.glob(fullpath):
            with open(fname) as f:
                for line in f:
                    if line.strip():                                    #removes any empty lines
                        if '*C' in line and not header:
                            if 'XYZ' in line:                           
                                line = re.sub('XYZ','X; Y; Z', line)    #find and fix header for each file correctly
                            line = line.split('; ')                     
                            for num, val in enumerate(line):            
                                if '\n' in val:
                                    val = re.sub('\n','',val)
                                if val == '*CHILD':
                                    header['CHILD'] = num
                                else:
                                    header[val] = num                          #once this is complete dictionary of param locations is created
                        elif line[0] != '*' and line[0] != ' ' and line[0] != '/n':                #for rest of data line by line
                            temp_line = line.split()
                            if len(temp_line) == len(header):
                                for point, val in enumerate(temp_line):        #Child, Type, Parent use DESCRETE numerical labelling
                                    if point != header['CHILD'] and point != header['TYPE'] and point != header['PARENT']:
                                        temp_line[point] = float(temp_line[point])
                                data_list.append(temp_line)                   

        for line in data_list:
            if line[header['TYPE']] == '4':                     #splits by node type (Apical or Basal)
                apical_list.append(line)
            elif line[header['TYPE']] == '3':
                basal_list.append(line)

        if not 'Basal' in archive_dict.keys():
            archive_dict['Basal'] = {}                                   #sets compartment(cell) type as nested dictionary within complete archive dictionary
        basal_list = zip(*basal_list)                           #zip will shift column/row organization to follow order of params (as tuple)
        archive_dict['Basal'][d1] = basal_list                  #each archive (of certain region i.e. Hippocampal archives) will be under Basal or Apical in archive_dict
        if apical_list:
            if not 'Apical' in archive_dict.keys():
                archive_dict['Apical'] = {}                         
            apical_list = zip(*apical_list)                        
            archive_dict['Apical'][d1] = apical_list

    #header is dictionary of all data parameters found in file
    params = {key:header[key] for key in header if key not in str_list}  

    #plot mulitple archives with single comp type 

    if 'archives' in graphing:
        for param in params:                  #make y value for plot specified (input) possibly function to plot against radius and separate residuals
            for num,comp_type in enumerate(archive_dict):
                figure(num = num, figsize = (14,8))
                for archive in archive_dict[comp_type]:           
                    label = archive + ' ' + comp_type
                    plt.plot(archive_dict[comp_type][archive][header[param]], archive_dict[comp_type][archive][header['RADIUS']], 'o', label = label)
                title = str(comp_type) + ' ' + str(param)
                plot_label(str(param),'Radius',title)

    #plot mulitple comp type (if present) with single archive

    if 'comp_type' in graphing:    
        for param in params:           
            #maybe initialize here axis = fig.axes
            #ax = fig.axis
            for num,archive in enumerate(dirs):
                #would have to initialize fig,ax = plt.subplot()
                figure(num = num, figsize = (14,8))
                #fig,ax = plt.subplot(nrows = len(archive_dict.keys()))
                for comp_type in archive_dict.keys()[::-1]:  #only two comp_types of Apical and Basal, will plot Basal after Apical for easier viewing
                    #plt.figure(num)
                    label = comp_type
                    #ax[num].plot(archive_dict[comp_type][archive][header[param]], archive_dict[comp_type][archive][header['RADIUS']], 'o', label = label)        #if i wanted to plot basal and apical in subplot, same figure...
                    plt.plot(archive_dict[comp_type][archive][header[param]], archive_dict[comp_type][archive][header['RADIUS']], 'o', label = label)
                title = str(archive) + ' ' + str(param) + ' vs. ' + 'RADIUS' #str(root) + ' ' + str(param) #str(cell_type)
                plot_label(str(param),'Radius',title)

    if '3d' in graphing:
        z = 'RADIUS'
        for comp_type in archive_dict:
            param_comb = []                                    #list to hold parameter combinations already plotted (no overlap desired)
            for x,y in product(params.keys(),params.keys()):   #iterates over different combinations if unique and if contains Radius
                if x != y and x != z and y != z:
                    param_list = [x,y,z]
                    if param_comb:
                        count = 0
                        for comb in param_comb:
                            if Counter(param_list) == Counter(comb):
                                count = 1
                        if count != 1:                        #if combination of parameters not found in list 3d plot
                            param_comb.append(param_list)
                            plot_3d(archive_dict, comp_type, header, x, y, z)
                    if not param_comb:                        #initiates plotting with first parameter combination sent to list
                        param_comb.append(param_list)
                        plot_3d(archive_dict, comp_type, header, x, y, z)

    if 'fit' in graphing:         #pearson's r correlation for parameters (only linear relationship)

        param_data = {}

        for comp_type in archive_dict.keys():        #saves param data into intergrated list (from archives) param_dict --> comp_type --> param_list --> param
            param_data[comp_type] = {}          
            for param in params.keys():
                if not param in param_data[comp_type]:
                    param_data[comp_type][param] = []
                for archive in archive_dict[comp_type].keys():     
                    param_data[comp_type][param] = list(archive_dict[comp_type][archive][header[param]]) + list(param_data[comp_type][param])  #option to create concatenation dictionary

        corr_dict = {}
        for comp_type in param_data.keys():
            corr_dict[comp_type] = pd.DataFrame(index = param_data[comp_type].keys(), columns = param_data[comp_type].keys())
            for rparam in param_data[comp_type].keys():                            #how to concatenate each archive for single corr table
               for cparam in param_data[comp_type].keys():                         #later this is will refer to already made xlist and ylist (instead of rewriting ever iteration)
                   corr_dict[comp_type].loc[rparam,cparam] = round(pearsonr(param_data[comp_type][rparam], param_data[comp_type][cparam])[0], 4)  #in or out of loop for temp or permanent x and y list
        for comp_type in param_data.keys():                                        #since this is Pearson's R it is showing strength of LINEAR correlation between params
            print(comp_type)
            print(corr_dict[comp_type])

        #plot_fit(param_data['Apical']['PARENT_RAD'], param_data['Apical']['RADIUS'], func, np.arange(0,7,1), x_max = 7)
        #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func1, np.arange(0,np.max(param_data['Apical']['BRANCH_LEN']),1))
        #plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func2, np.arange(0,np.max(param_data['Apical']['BRANCH_LEN']),1))
        find_min = param_data['Apical']['BRANCH_LEN']  
        find_min = list(set(find_min))  #for some reason we have to remove duplicates before sorting...
        find_min.sort()
        plot_fit(param_data['Apical']['BRANCH_LEN'], param_data['Apical']['RADIUS'], func4, np.arange(find_min[1],np.max(param_data['Apical']['BRANCH_LEN']),1), x_max = None, x_min = np.min(param_data['Apical']['BRANCH_LEN']))

        #I wonder if I can make a function to fit optimize for certain parameters with good correlation (linear and nonlinear) values


#TODO
#Implement Random Forest with parameters ideal for Radius estimation (want to lower the amount of features needed BUT still have an accurate equation for estimate)
