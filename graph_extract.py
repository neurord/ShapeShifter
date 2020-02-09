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
#Feb. 9, 2020

import numpy as np
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from collections.abc import Mapping
#from sklearn.linear_model import LinearRegression
#import seaborn as sns

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

#def piecewise_linear(x, x0, k1, k2, b):#(x, x0, y0, k1, k2):
    #return np.piecewise(x, [x < x0], [lambda x:k1*x + b, lambda x:k2*x + b])

def combo_func(X):#, a, b, c):
    x1,x2 = X          #where x1 and x2 are unpacked in func
    #return np.log(a) + b*np.log(x) + c*np.log(y)
    return m1*x1 + m2*x2 + b


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
def plot(data,to,extras,fit_line = None):                                                   
    plt.ion()  
    if to == 'Simple': #'2-D'                           
        fig = plt.figure(figsize = (14,8))
        if not fit_line:
            for archive in data:
                plt.plot(data[archive][extras[0]], data[archive][extras[1]], 'o', label = archive)
        elif 'Res' in fit_line:       #Gives FUTURE_WARNING for each elif statement...
            plt.plot(data[extras[0]], data[extras[1]], 'o', label = 'Residuals')
        elif 'Fit' in fit_line:
            plt.plot(data[extras[0]], data[extras[1]], 'o', label = 'Data')
            plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2])
        elif 'Post'in fit_line:
            plt.plot(data[extras[0]], data[extras[1]], 'o', label = 'Data')
        plt.xlabel(extras[0])
        plt.ylabel(extras[1])
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
    corr_df = pd.DataFrame(index = keys, columns = keys)
    for rparam in keys:
        for cparam in keys:
            corr_df.loc[rparam,cparam] = round(pearsonr(data[rparam],data[cparam])[0],4)
    return corr_df


'''Descrete Functions (top) to Selected Data; Outputs Variables and Indicates Possible Fit Fine-Tuning'''
def fit(data,labels,func):
    '''Calculate Variables for Curve Fitting'''
    '''Calculate Residuals within Dictionary Organization'''
    if len(labels[0]) > 1: 
        temp = []
        for i in labels[0]: #this will be the x params selected for function
            temp.append(data[i])  #might only print out each letter, not the extent of labels[0]
        x = tuple(temp)           #<--- this will be useful if sending in multiple X to function
    else:
        x = data[labels[0][0]]
    y = data[labels[1][0]]
    popt, pcov = optimize.curve_fit(func[0],x,y)
    print('This is fit estimate for : ', func[1])
    print('popt',popt)
    print('pcov',pcov)
    xlist = [str(i) for i in labels[0]]
    print('Where x = ', xlist[0], ' and y = ', labels[1][0]) 
    
    for i in labels[0]: #will keep for now as it will print individual select X to Radius
        temp_max = round(max(data[i]))  #Plots Fitted Equation to Data
        temp_min = round(min(data[i]))
        x_range = np.arange(temp_min, temp_max + 1, 1)
        plot(data, 'Simple', [i, labels[1][0]], [x_range, func[0](x_range, *popt), func[1], 'Fit'])
        
        res_label = str(i) + '_Res'  #Finds Residuals and saves as new parameter within data
        prediction = []; residuals = []
        data[res_label] = []
        for x in data[i]:
            predict_y = func[0](x,*popt)
            prediction.append(predict_y)
        for y_val, y_predict in list(zip(data[labels[1][0]],prediction)):
            residuals.append(y_val - y_predict)
        data[res_label] = residuals
        print('Residual Plot : ', i, ' to ', labels[1][0])
        plot(data, 'Simple', [labels[1][0], res_label], ['Res'])

    #for select_param in data[comp_type][archive].keys():
        #plot(data[comp_type], 'Simple', [select_param, res_label])

    return data, [func1,popt,pcov]
    

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

        
'''Separate Data by Connection to Soma (either 0 -> directly connected, 1 -> all others)'''
params = {key:header[key] for key in header if key not in str_list}   #new dictionary to only contain useful parameters
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

#Combine archive data together for initial forward selection towards Radius
comb_dict = {}; dataf = {}
for connect in param_data:
    if not connect in comb_dict.keys():
        comb_dict[connect] = {}; dataf[connect] = {}
    for comp_type in param_data[connect]:   
        if not comp_type in comb_dict[connect].keys():
            comb_dict[connect][comp_type] = {}; dataf[connect][comp_type] = {}
        for param in params.keys():   
            temp = []
            for archive in param_data[connect][comp_type]:
                temp.append(param_data[connect][comp_type][archive][param])
            temp = list(flatten(temp))
            comb_dict[connect][comp_type][param] = temp
        dataf[connect][comp_type] = pd.DataFrame(comb_dict[connect][comp_type]) #will hold dataframe organized values in same hierarchy...

        
'''Parameter Correlation (Pearson's r) to Radius For Linear Relationships'''
initial_corr = {}
excluded = ['DIRECT_SOMA','NODE_ORDER','NODE_COUNT']
for connect in comb_dict:
    initial_corr[connect] = {}
    for comp_type in comb_dict[connect]:
        if connect == 'Direct':
            initial_corr[connect][comp_type] = corr(comb_dict[connect][comp_type], [i for i in params.keys() if i not in excluded])
        else:
            initial_corr[connect][comp_type] = corr(comb_dict[connect][comp_type], [i for i in params.keys() if i != 'DIRECT_SOMA'])

with open('initial_corr_dict.txt', 'a') as outfile:
    for connect in initial_corr:
        for comp_type in initial_corr[connect]:
            outfile.write(str(connect) + ' - ' + str(comp_type) + '\n')
            outfile.write('\n')
            outfile.write(str(initial_corr[connect][comp_type]) + '\n')
            outfile.write('\n')
    outfile.close()              

'''Plot Parameters to Radius For Initial Comparison'''
'''Plot Types'''
   #Merge - Separate to plot Direct/Indirect separately each archive to parameters
   #Merge - Combine to plot Direct/Indirect together w/o archive specific to parameters
   #Simple to plot selected key in selected dictionary
   #3d to plot all parameters to Radius to see 2 parameter influence on Radius

#plot(param_data, 'Merge', ['Separate', params.keys()]) 
#plot(param_data, 'Merge', ['Combine', params.keys(), archive_dict.keys(), dirs]) 

'''
for i in params.keys():
    if 'Apical' in archive_dict.keys():
        plot(param_data['Direct']['Apical'], 'Simple', [i, 'RADIUS'])
        plot(param_data['Indirect']['Apical'], 'Simple', [i, 'RADIUS'])
    if 'Basal' in archive_dict.keys():
        plot(param_data['Direct']['Basal'], 'Simple', [i, 'RADIUS'])
        plot(param_data['Indirect']['Basal'], 'Simple', [i, 'RADIUS'])
'''
#plot(param_data,'3d', params.keys())

#comb_dict quite useful if we do not need to pay attention to archive orgin...
#seaborn package may be useful for residuals
'''
x = dataf['Direct']['Apical']['PARENT_RAD']
y = dataf['Direct']['Apical']['RADIUS']
fig = plt.figure()
sns.residplot(x,y)
plt.xlabel('Parent Radius')
plt.ylabel('Radius')
'''
#Residual plot looks same if i were to show parent_rad function...

'''
Framed = {}
for connect in param_data:
    Framed[connect] = {}
    for comp_type in param_data[connect]:
        Framed[connect][comp_type] = forward_selected(pd.DataFrame(comb_dict[connect][comp_type]), 'RADIUS')

for connect in Framed:
    for comp_type in Framed[connect]:
        print(Framed[connect][comp_type].model.formula)
        print(Framed[connect][comp_type].rsquared_adj)
'''


#try stepwise regression --> add features at a time and keep if R is increased
#plot residuals to remaining parameters

#plot(#transformedx, transformedy)

'''Testing Linear Regression Model from sklearn'''
'''
mlr = LinearRegression()
mlr.fit(dataf['Indirect']['Apical'][['PARENT_RAD','NODE_ORDER']],dataf['Indirect']['Apical']['RADIUS'])
print(mlr.intercept_)
print(mlr.coef_)
'''
#check to see if mlr(LinearRegression()) has some significance value for model fit

'''Initiate Fit and Residuals'''
DA_PR, DA_PR_fit = fit(comb_dict['Direct']['Apical'], [['PARENT_RAD'],['RADIUS']], [func0,'y = mx + b'])
IA_PR, IA_PR_fit = fit(comb_dict['Indirect']['Apical'], [['PARENT_RAD'],['RADIUS']], [func0,'y = mx + b'])

IA_PR_Res_corr = corr(IA_PR, [i for i in IA_PR.keys() if i != 'DIRECT_SOMA'])

'''
for param in IA_PR:
    plot(IA_PR, 'Simple', [param, 'PARENT_RAD_Res'], ['Post'])
'''
#could try and transform data now...

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

'''Data Transformations'''

#THINGS TO TRY
#Separating the Direct Soma points (0) and (1) to see if they have separate equations to estimate radius
#transform the data to see any nonlinear relationships (if linear after the transformation)
#Check to see why direct_soma nodes have differing values of Soma Radius...(or their parent radius)

#TODO
#Implement Random Forest with parameters ideal for Radius estimation (want to lower the amount of features needed BUT still have an accurate equation for estimate)

