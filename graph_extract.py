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
#Sep 3, 2019

import numpy as np
#import scipy
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D

import argparse
import os
import glob
import re

def func(x,m,b): 
    return m*x+b

def plot_func(xlabel, ylabel, title, var = None):
    if var == 'label':
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        #plt.savefig(title + '.png')                       #works better but still missing the basal graphs if both basal and apical data exist in particular dictionary
        #plt.clf()                                         #more explicit way to refer to figures
        plt.show()
       
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   
parser.add_argument("--graphing", type = str, nargs = '+', choices = {'archives','comp_type','fit'})

args = parser.parse_args()                      
path = args.path                                
graphing = args.graphing

'''Extracts data from .CNG_extract.swc files moved to archive directories under main brain region directory (i.e. 'Groen' and 'Jaffe' directories under Hippocampus)''' 
root, dirs, files = list(os.walk(path))[0]                           #make sure all extracted files have same parameters before graphing                                            
archive_dict = {}; header = {}; params = {}                          #params will be filled once str_list params are removed from header.keys()
str_list = ['PARENT','X','Y','Z','CHILD','TYPE']                     #identifiers of particular compartments not useful for Radius comparision


for d1 in dirs:
    fullpath = path + d1 + '/*CNG_extract.txt'
    data_list = []; apical_list = []; basal_list = []
    print('Working on Directory : ', str(d1))
    for fcount,fname in enumerate(glob.glob(fullpath)):
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
            plot_func(str(param),'Radius',title,'label')

#plot mulitple comp type (if present) with single archive

if 'comp_type' in graphing:    
    for param in params:           
        #maybe initialize here axis = fig.axes
        #ax = fig.axis
        for num,archive in enumerate(dirs):
            #would have to initialize fig,ax = plt.subplot()
            figure(num = num, figsize = (14,8))
            #fig,ax = plt.subplot(nrows = len(archive_dict.keys()))
            for comp_type in archive_dict.keys()[::-1]:
                #plt.figure(num)
                label = comp_type
                #ax[num].plot(archive_dict[comp_type][archive][header[param]], archive_dict[comp_type][archive][header['RADIUS']], 'o', label = label)        #if i wanted to plot basal and apical on same plot...
                plt.plot(archive_dict[comp_type][archive][header[param]], archive_dict[comp_type][archive][header['RADIUS']], 'o', label = label)
            title = str(archive) + ' ' + str(param) + ' vs. ' + 'RADIUS' #str(root) + ' ' + str(param) #str(cell_type)
            plot_func(str(param),'Radius',title,'label')

#could also make function to store values of residuals as well...just into other dictionaries I would send into function



if 'fit' in graphing:
   
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

    #if I were to limit the data, I would first have to desginate where to split, and then pair x-y data together so that order and length of x and y list is mantained
    #for parent radius, this would only allow the fit to be for the data that is linear and not the extraneous points

    xdata = param_data['Apical']['PARENT_RAD']                    #so, this may be improved if I choose another user argument of 'fit' and allow them to select the parameter of choice
    #xdata = param_data['Basal']['PARENT_RAD']
    ydata = param_data['Apical']['RADIUS']
    #ydata = param_data['Basal']['RADIUS']
    data_list = []

    x_max = 7 #x_max is the maximum accepted value for linear relationship, not = max(xdata)
    for x,y in zip(xdata, ydata):
        if x < x_max:                   #this is point on parent_rad where there is no longer linear relationship to radius
            data_list.append((x,y))

    popt, pcov = optimize.curve_fit(func,xdata,ydata)
    print('popt',popt)
    print('pcov',pcov)
    plt.figure(figsize = (14,8))
    plt.plot([i[0] for i in data_list],[i[1] for i in data_list], 'o', label = 'Data')   #i would have to order the xdata being sent into func
    print(max(xdata))
    x_line = np.arange(0, x_max, 1)
    #xpredict = [i[0] for i in data_list]  #trying to do residual list, but probably in separate plot than here
    #xpredict.sort()
    plt.plot(x_line, func(x_line, popt[0], popt[1]), color = 'orange', label = 'Fitted Function')
    plt.legend()
    plt.show()

    prediction = []
    for x in xdata:
        predict_y = func(x, popt[0], popt[1]) if x < x_max else 0
        prediction.append(predict_y)
    #print(len(prediction))
    #print(len(xdata))
    
    #x_predict = [i[0] for i in data_list]
    #y_predict = []                             #actually residuals are y vals - pred vals
    residuals = []
    #for x in x_predict:
        #y_predict.append(func(x, popt[0], popt[1]))
    #for y_val, y_pred in zip([i[1] for i in data_list],y_predict):
        #residuals.append(y_val - y_pred)
    for y_val, y_predict in zip(ydata,prediction):
        residuals.append(y_val - y_predict)
    #print(len(xdata))
    #print(len(ydata))
    #print(len(prediction))
    #print(len(residuals))
    plt.figure(figsize = (14,8))
    #plt.plot([i[1] for i in data_list], residuals, 'o')
    plt.plot(xdata,residuals,'o')
    #plt.plot(x_line, func(x_line, popt[0], popt[1]), color = 'orange', label = 'Fitted Function')
    plt.show()
    
    #so this is working; how would I limit the remaining data so that it would have equal lengths and val order for other parameter comparisons to residuals....

    #I wonder if I can make a function to fit optimize for certain parameters with good correlation (linear and nonlinear) values
    #I also wonder if I can limit the data entering optimize.fit (since parent rad is good up to a certain parent_rad value)

    #maybe if correlation (linear) is good for particular combination do something else with that in code
    #i.e. if corr value is above certain acceptance level of good relationship





'''

if len(Radius_a):
        #plot both apical and basal data
        plt.plot(HoSt_a, Radius_a, 'o', label = 'apical dendrite')
        plt.plot(HoSt_b, Radius_b, 'o', label = 'basal dendrite')
    else:
        #plot only basal data as no apical data exists in particular archive
        plt.plot(HoSt_b, Radius_b, 'o', label = 'basal dendrite')
    #Axes3D.plot_surface()look for seqeuntial color map

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if len(Radius_a):
        ax.scatter(End_Distance_a, Radius_a, HoSt_a, c='r', marker = 'o')
        ax.scatter(End_Distance_b, Radius_b, HoSt_b, c='b', marker = 'o')
    #else:
    ax.scatter(Parent_Rad_b, Radius_b, End_Distance_b, c='g', marker = 'o')
    ax.set_xlabel('Parent_Rad')
    ax.set_ylabel('Radius')
    ax.set_zlabel('End_Distance')
    plt.show()

    #either send data to separate for loop(plot entire data set) or function 
    #plt.ylabel('Radius')
    #plt.xlabel('HS')
    
    #plt.legend()
    #plt.savefig(morph_file + 'HS' + ".png", format="PNG")
    #print(morph_file + 'HS' + ".png")
    #plt.show()
  
    #3D plotting

    #X,Y,Z = [][][]

    X = 
    Y = 
    Z = 
    ax.scatter(X, Y, Z, marker = 'o')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    plt.plot(Degree, Radius, 'o')
    plt.ylabel('Radius')
    plt.xlabel('Degree')

    print(len(Radius))
    print(len(Distance))
    print(len(Degree))
    print(len(Order))
    print(len(Parent_Rad))
    print(len(End_Distance))
    print(len(Ends))
    print(len(HoSt))

    plt.savefig(morph_file + 'degree' + ".png", format="PNG")
    print(morph_file + 'degree' + ".png")
    plt.show()
    
    
    if Ends == Degree:
        print('y')
    else:
        print('n')
    
'''
