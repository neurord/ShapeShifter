'''<Extracts data values from .swc files>
    Copyright (C) <2020>  <Jonathan Reed>

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
#        where path is name of path to directory containing archive directories with archive files


'''Must use Python 3 for seemless usage in statsmodels package'''
#George Mason University
#Jonathan Reed
#Mar. 1, 2020

import numpy as np
import math
from scipy import optimize
from scipy.stats import pearsonr
from scipy.stats import linregress
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm

import argparse
import os
import glob
import re

'''Function Definitions for Fit Equations'''

def inverse_func(X,m1,m2,b):  #x2 currently sent in as 1/(value within Path to End)
    x1,x2 = X
    return m1*x1 + m2*x2 + b  #(where 1 --> Parent Radius and 2 --> Path to End)

def func0(x,m,b):                                     
    return m*x+b


'''Flatten Nested List to Single Non-Nested List'''
def flatten(container):    
    for i in container:                 #useful for appending lists to list and returning values
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

'''Save plot as png'''            
def save_png(png, title):               #default for all plots instead of plotting to window
    png.savefig(title)
    print('File Created : ' + str(title))
    png.close()

    
'''Main Plot Function for Basic, Fit or 3d Color Plots'''
def simple_plot(data, var_labels, plot_type, fit_line = None, _alpha = None): 
    #plt.ion()                                     #default commented out to avoid extensive numbers of plots to screen
    fig = plt.figure(figsize = (14,8))
    if plot_type == 'Color3d':                 #plot 3 variables, with variable 2 using colorbar
        fig, ax = plt.subplots(figsize = (14,8))
        scat = ax.scatter(data[var_labels[0]], data[var_labels[1]], c=data[var_labels[2]], s=100, marker='o', label=var_labels[2])
        fig.colorbar(scat) 
        title = plot_type + ' Plot: ' + var_labels[0] + ' vs ' + var_labels[1] + ' vs ' + var_labels[2]
    else:                                      #plot 2 variables either within or outside fit function
        if _alpha:
            plt.plot(data[var_labels[0]], data[var_labels[1]], 'o', alpha = _alpha, label = plot_type) #can set alpha = 0.03 to see overlapping values on plot
        else:
            plt.plot(data[var_labels[0]], data[var_labels[1]], 'o', label = plot_type) 
        if plot_type == 'Fit':                 #plot fitted equation with data values
            plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) #residual currently only line
        title = plot_type + ' Plot: ' + var_labels[0] + ' vs ' + var_labels[1]        #could be useful if points for certain plots
    plt.legend()
    plt.xlabel(var_labels[0])
    plt.ylabel(var_labels[1])
    plt.title(title)
    save_png(plt, title + '.png')

'''Initial Plot Function comparing Archive Data'''
def merge_plot(data, to, extras):
    #plt.ion()                             #default commented out to avoid extensive numbers of plots to screen
    if to == 'Separate':                   #plot distinct archives based on connection to soma and compartment-type
        for connect in data:               #separate plots for connection to soma, subplots for compartment-type
            for param in extras:           
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

    if to == 'Combine':                    #plot archives together based on connection to soma and compartment-type
        for param in extras[0]:            #separate subplots for Apical/Basal compartment type
            fig = plt.figure(figsize = (14,8))
            for num,comp_type in enumerate(extras[1]):
                ax1 = fig.add_subplot(len(extras[1]),1,num+1)
                for connect in data:
                    ax1.plot(data[connect][comp_type]['PARENT_RAD'], data[connect][comp_type][param], 'o', label = connect)
                ax1.set_title(comp_type)  
            plt.legend()                   #default is 'Parent Radius' to visualize difference of nodes directly connected to soma 
            plt.xlabel('Parent Radius (um)') 
            plt.ylabel(str(param))        
            save_png(plt, 'Combined'  + ' ' + str(param) + ' vs ' + 'Parent Radius' + '.png')

'''Alternative Plot Function for Visualizing Radius with 2 Features'''
def _3d_plot(data,extras):                 #plot features to Radius for 3-d visualization                     
    combinations = []                      #create combinations of features as variables for plotting
    for x in extras:
        for y in extras:               
            if x != y and x != 'RADIUS' and y != 'RADIUS':
                temp = [x,y]
                if not combinations:                        #removes any duplicate combinations of features to plot
                    combinations.append(temp)
                if combinations:
                    if temp not in combinations and temp[::-1] not in combinations:
                        combinations.append(temp)

    for connect in data:                   #create 3d plots from unique feature combinations 
        for comb in combinations:
            fig = plt.figure(figsize = (14,8))
            for num,comp_type in enumerate(data[connect]):
                ax1 = fig.add_subplot(len(data[connect].keys()),1,num+1, projection = '3d')
                for archive in data[connect][comp_type]:
                    dict_loc = data[connect][comp_type][archive]             
                    ax1.scatter(dict_loc[comb[0]],dict_loc[comb[1]],dict_loc['RADIUS'],'o', label = archive)
                ax1.set_xlabel(comb[0])                      #3rd variable default to Radius
                ax1.set_ylabel(comb[1])
                ax1.set_zlabel('RADIUS')
                ax1.set_title(comp_type)
                plt.legend()
            save_png(plt, connect + ' ' + comb[0] + ' ' + comb[1] + ' ' + 'RADIUS' + '.png')

                
'''Pearson's r Correlation for Feature Relationships'''
def corr(data,keys):
    corr_df = pd.DataFrame(index = keys, columns = keys) #correlation of all features to each other
    for rparam in keys:                                  #helpful for initial relation to Radius and multicollinearity
        for cparam in keys:
            corr_df.loc[rparam,cparam] = round(pearsonr(data[rparam],data[cparam])[0],4)
    return corr_df                                       #compare with results from scipy.stats.linregress


'''Fit Data to Selected Equation and Obtain Residuals'''
def initial_fit(data,xlabel,ylabel,func):
    
    '''Fit Data to Choosen Function'''
    if len(xlabel) > 1:       #if xlabel has more than single variable for fit
        temp = []             #fit multiple x-variables to y
        for i in xlabel: 
            temp.append(data[i])  
        x = tuple(temp)           
    else:
        x = data[xlabel[0]]
    y = data[ylabel[0]]                             #assumes only fit to single y-variable
    popt, pcov = optimize.curve_fit(func[0],x,y)    #fits equation to choosen function if possible
    print('This is fit estimate for : ', func[1])
    print('popt',popt)
    print('pcov',pcov)
    print('Where x = ', xlabel, ' and y = ', ylabel) 

    '''Plot Data to Function and Find Residuals'''
    for i in xlabel:          #plot x-variables individually to y and fitted function
        temp_max = round(max(data[i])) 
        temp_min = round(min(data[i]))
        x_range = np.arange(temp_min, temp_max + 1, 1)
        simple_plot(data, [i, ylabel[0]], 'Fit', [x_range, func[0](x_range, *popt), func[1]])
            #fails here if sending in >1 x variable with too many values to unpack error
        res_label = str(i) + '_Res'  
        predictions = [func[0](x,*popt) for x in data[i]]
        data[res_label] = [yval - est for yval,est in zip(data[ylabel[0]],predictions)]
        print('Residual Plot : ', i, ' to ', ylabel[0])
        simple_plot(data, [ylabel[0], res_label], 'Residual')
        
    #for param in data:                           #plot features to residuals to reveal possible trends not explained by estimate
    #    simple_plot(data, [param, res_label], 'Residual')
        
    return data, [popt,pcov]
    
def ols_fit(data,xlabel,ylabel,constant = None):
    temp_df = pd.DataFrame(data)
    X = temp_df[xlabel]
    Y = temp_df[ylabel]
    if constant:
        X = sm.add_constant(X)       #adds intercept if desired #<-- could make this if statement
    model = sm.OLS(Y,X).fit()     
    print(model.summary())
    return model,model.predict(X)


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
    archive_dict['Basal'][d1] = list(zip(*basal_list))  
    if apical_list:
        if not 'Apical' in archive_dict.keys():
            archive_dict['Apical'] = {}                         
        archive_dict['Apical'][d1] = list(zip(*apical_list))

        
'''Separate Data by Connection to Soma (either 0 -> directly connected, 1 -> all others)'''
params = {key:header[key] for key in header if key not in str_list}   #new dictionary to hold useful features to Radius
param_data = {}                                                       #maintains archive file orgin
param_data['Direct'] = {}                #directly connected to soma
param_data['Indirect'] = {}              #all other nodes (excluding soma points)

#Fill Dictionaries by Connection to Soma (as 0, others as 1) for all parameters
for comp_type in archive_dict:                         
    param_data['Direct'][comp_type] = {}; param_data['Indirect'][comp_type] = {}
    for archive in archive_dict[comp_type]:
        param_data['Direct'][comp_type][archive] = {}; param_data['Indirect'][comp_type][archive] = {}
        for param in params.keys():
            param_data['Direct'][comp_type][archive][param] = []; param_data['Indirect'][comp_type][archive][param] = []
        for num,val in enumerate(archive_dict[comp_type][archive][header['DIRECT_SOMA']]):
            for param in params.keys():                  #if directly connected (0) to soma --> Direct, else (1) --> Indirect
                if val == 0:
                    param_data['Direct'][comp_type][archive][param].append(archive_dict[comp_type][archive][header[param]][num])
                else:
                    param_data['Indirect'][comp_type][archive][param].append(archive_dict[comp_type][archive][header[param]][num])

#Combine archive data together for trends towards Radius
comb_data = {}; comb_df = {}
for connect in param_data:
    comb_data[connect] = {}; comb_df[connect] = {}         #dictionary organization with data as list (data) or dataframe (df)
    for comp_type in param_data[connect]:   
        if not comp_type in comb_data[connect].keys():
            comb_data[connect][comp_type] = {}; comb_df[connect][comp_type] = {}
        for param in params.keys():   
            temp = []
            for archive in param_data[connect][comp_type]:                    #appends archive file data as single list
                temp.append(param_data[connect][comp_type][archive][param])   
            temp = list(flatten(temp))
            comb_data[connect][comp_type][param] = temp
        comb_df[connect][comp_type] = pd.DataFrame(comb_data[connect][comp_type]) 

        
'''Parameter Correlation (Pearson's r) to Radius For Linear Relationships'''
initial_corr = {}
excluded = ['DIRECT_SOMA','NODE_ORDER','NODE_COUNT']       #not useful for nodes directly connected to soma
for connect in comb_data:
    initial_corr[connect] = {}
    for comp_type in comb_data[connect]:
        if connect == 'Direct':
            initial_corr[connect][comp_type] = corr(comb_data[connect][comp_type], [i for i in params.keys() if i not in excluded])
        else:
            initial_corr[connect][comp_type] = corr(comb_data[connect][comp_type], [i for i in params.keys() if i != 'DIRECT_SOMA'])
'''
with open('initial_corr_dict.txt', 'a') as outfile:
    for connect in initial_corr:
        for comp_type in initial_corr[connect]:
            outfile.write(str(connect) + ' - ' + str(comp_type) + '\n')
            outfile.write('\n')
            outfile.write(str(initial_corr[connect][comp_type]) + '\n')
            outfile.write('\n')
    outfile.close()              
'''

'''Ideal Usage of PYTHON Statistical Tools'''

'''Check Initial Correlations between Features and Radius (Linear)'''
#print(initial_corr)

'''Plot Features to Visualize Relationship to Radius'''
#merge_plot(comb_data, 'Combine', [params.keys(),archive_dict.keys()])
#merge_plot(param_data, 'Separate', params.keys())                      #separate keeps archive designation
#for param in params.keys():
    #simple_plot(comb_data['Indirect']['Apical'], [param,'RADIUS'],'Basic')
#_3d_plot(param_data, params.keys())

'''Initiate Fit and Residuals for Best Feature to Estimate Radius'''
#DA_PR, DA_PR_fit = initial_fit(comb_data['Direct']['Apical'], ['PARENT_RAD'],['RADIUS'], [func0,'y = mx + b'])
#IA_PR, IA_PR_fit = initial_fit(comb_data['Indirect']['Apical'], ['PARENT_RAD'], ['RADIUS'], [func0,'y = mx + b'])
#for param in IA_PR.keys():
    #simple_plot(IA_PR,[param,'PARENT_RAD_Res'],'Residual') #prints residuals if commented out in fit function
    
'''Utilize statsmodels.ols for Multiple Regression and send in possible Transformed Feature Values'''
#how to represent infinity appropriately...
#trans_PTE = [1/i if i != 0 else 0 for i in comb_data['Indirect']['Apical']['PATH_TO_END']]  #how do we effectively approach infinite values or address properly
#trans_PTE = [np.log(i) if i != 0 else 0 for i in comb_data['Indirect']['Apical']['PATH_TO_END']] #(-np.finfo('d').max) for i in comb_data['Indirect']['Apical']['PATH_TO_END']]

#trans_PTE = [1/np.exp(i) if i < np.log(np.finfo('d').max) else 0 for i in comb_data['Indirect']['Apical']['PATH_TO_END']] #high values 
#trans_PTE = [1/(1+np.exp(i)) if i < np.log(np.finfo('d').max) else 1 for i in comb_data['Indirect']['Apical']['PATH_TO_END']] #high values

#trans_ND = [1/i if i != 0 else 0 for i in comb_data['Indirect']['Apical']['NODE_DEGREE']]
#trans_ND = [np.log(i) if i != 0 else 0 for i in comb_data['Indirect']['Apical']['NODE_DEGREE']]
#trans_ND = [1/np.exp(i) if i < np.log(np.finfo('d').max) else 0 for i in comb_data['Indirect']['Apical']['NODE_DEGREE']]
#trans_ND = [1/(1+np.exp(i)) if i < np.log(np.finfo('d').max) else 1 for i in comb_data['Indirect']['Apical']['NODE_DEGREE']]

#fit_data = {}                      #will hold selected features to fit to Radius

#fit_data['trans_PTE'] = trans_PTE
#fit_data['trans_ND'] = trans_ND
#for param in params:
    #fit_data[param] = comb_data['Indirect']['Apical'][param]
#model,predictions = ols_fit(fit_data,['PARENT_RAD','trans_ND'],'RADIUS') #can send in multiple X values to estimate Y
#model,predictions = ols_fit(fit_data,['PARENT_RAD','trans_ND'],'RADIUS','constant')

'''Plot Transformed Values with Radius or Residuals to better Select Equation Features'''
#fit_data['residuals'] = [yval - est for yval,est in zip(fit_data['RADIUS'],predictions)]
#simple_plot(fit_data, ['RADIUS',residuals], 'Residual', _alpha = 0.3)

#for param in fit_data:
    #simple_plot(fit_data, ['PATH_TO_END', residuals, param], 'Color3d')



