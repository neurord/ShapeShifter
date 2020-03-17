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
#Mar. 12, 2020

import numpy as np
import math
from scipy import optimize
from scipy.stats import pearsonr
from scipy.stats import linregress
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.offsetbox import AnchoredText

import statsmodels.api as sm
import random
import seaborn as sns
import argparse
import os
import glob
import re

'''Function Definitions for Fit Equations'''
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

'''Will Split N sequences into 'equal' sizes'''
def split_seq(seq, size):
        newseq = []
        splitsize = 1.0/size*len(seq)
        for i in range(size):
                newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
        return newseq
    
'''Main Plot Function for Basic, Fit or 3d Color Plots'''
def simple_plot(data, var, _alpha, title, legend = None, plot_type = None, labels = None, fit_line = None, add = None, where = None): 
    #plt.ion()                                     #default commented out to avoid extensive numbers of plots to screen
    fig, ax = plt.subplots(figsize = (14,8))
    plt.rc('font', size = 18)
    
    if add:
        at = AnchoredText(add,
                  prop=dict(size=14), frameon=True,
                  loc='upper center',
                  )
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") #maybe look for ways to remove border around annotation
        ax.add_artist(at)
        
    if plot_type == 'Color3d':                     #plot 3 variables(x,y,z), with variable 2 as 'z' using colorbar
        scat = ax.scatter(data[var[0]], data[var[1]], c=data[var[2]], s=100, marker='o', alpha = _alpha, label=legend)
        fig.colorbar(scat)
    elif plot_type == 'Fit':                       #plot initial fitted equation to data values
        plt.plot(data[var[0]], data[var[1]], 'o', alpha = _alpha)
        plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) 
        legend = True
    else:
        if len(data) > 1:                          #plot multiple data sets with same axes
            for num in range(len(data)): 
                if legend:
                    plt.plot(data[num][var[0]], data[num][var[1]], 'o', alpha = _alpha[num], label = legend[num])
                else:
                    plt.plot(data[num][var[0]], data[num][var[1]], 'o', alpha = _alpha)
        else:
            if legend:                             #plot single data set
                plt.plot(data[var[0]], data[var[1]], 'o', alpha = _alpha, label = legend)
            else:
                plt.plot(data[var[0]], data[var[1]], 'o', alpha = _alpha)
    if legend:
        if where:
            plt.legend(loc = where)
        else:
            plt.legend()
    if labels:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    else:
        plt.xlabel(var[0])           #default will label axes with data label
        plt.ylabel(var[1])
    plt.title(title)
    save_png(plt,title + '.png')


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
    return corr_df                                      

'''Fit Data to Selected Equation and Obtain Residuals'''
def initial_fit(data,xlabel,ylabel,func):
    x = data[xlabel[0]]
    y = data[ylabel[0]]                             #assumes only fit to single y-variable
    popt, pcov = optimize.curve_fit(func[0],x,y)    #fits equation to choosen function if possible
    print('This is fit estimate for : ', func[1])
    print('popt',popt)
    print('pcov',pcov)
    print('Where x = ', xlabel, ' and y = ', ylabel) 
 
    '''Plot Data to Function and Find Residuals'''
    for i in xlabel:
        temp_max = round(max(data[i])) 
        temp_min = round(min(data[i]))
        x_range = np.arange(temp_min, temp_max + 1, 1)
        simple_plot(data, [i, ylabel[0]], 1.0, 'Fitted Equation ' + i + ' to ' + ylabel[0], fit_line = [x_range, func[0](x_range, *popt), func[1]])
        res_label = str(i) + '_Res'  
        predictions = [func[0](x,*popt) for x in data[i]]
        simple_plot(data, [ylabel[0], predictions], 1.0, 'Predicted vs. Actual Plot : ' + i)
        
        data[res_label] = [yval - est for yval,est in zip(data[ylabel[0]],predictions)]
        print('Residual Plot : ', i, ' to ', ylabel[0])
        simple_plot(data, [ylabel[0], res_label], 1.0, 'Residual Plot : ' + i)
        
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

def save_file(coeff,filename):
    new_file = open(filename,'w')
    for i in coeff:
        new_file.write(str(i) + '\n')
    
'''Start of Working Code'''


parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   
args = parser.parse_args()                      
path = args.path                                

'''Locates and Organizes Archive Data'''
root, dirs, files = list(os.walk(path))[0]                    #all _extract files require identical parameters
archive_dict = {}; header = {}                                
file_list = []; complete_list = []

for d1 in dirs:
    fullpath = path + d1 + '/*CNG_extract.txt'
    data_list = []; apical_list = []; basal_list = []
    print('Working on Directory : ', str(d1))
    for fname in glob.glob(fullpath):                               #locates and loops through _extract files within fullpath
        temp_name = re.search(d1 + '/(.*).txt', fname)
        file_list.append(temp_name.group(1))
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
                            temp_line.append(temp_name.group(1))
                            temp_line.append(d1)
                            data_list.append(temp_line)
                            complete_list.append(temp_line)
                            
    '''Separation of Data By Apical/Basal Compartment Type & Archive Orgin'''
    for line in data_list:                                  #data_list has all compartment values
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

'''Separation of Data By Apical/Basal Compartment Type'''
complete_dict = {}  
temp_list = list(zip(*complete_list))
header['FNAME'] = len(header)        #useful for splitting data into training and test sets
header['ARCHIVE'] = len(header)
for num,param in enumerate(header):
    complete_dict[param] = temp_list[num]

if '4' in complete_dict['TYPE']:       #checks if there is at least one apical compartment
    apical_dict = {}; basal_dict = {}  #creates new dictionaries of only values within apical or basal compartment types
    for param in header:
        apical_dict[param] = []
        basal_dict[param] = []
    for num,i in enumerate(complete_dict['TYPE']):
        if i == '4':                   #will only keep apical and basal compartments separately, not axonal if present
            for param in header:
                apical_dict[param].append(complete_dict[param][num])
        elif i == '3':
            for param in header:
                basal_dict[param].append(complete_dict[param][num])
                
if not '4' in complete_dict['TYPE']:   #checks if there is no apical compartments
    temp_dict = {}
    for param in header:
        temp_dict[param] = []
    for num,i in enumerate(complete_dict['TYPE']):
        if i == '3':                   #will only keep basal compartments, not axonal if present
            for param in header:
                temp_dict[param].append(complete_dict[param][num])
    basal_dict = temp_dict


'''Separate Data by Connection to Soma (either 0 -> directly connected, 1 -> all others)'''
str_list = ['PARENT','X','Y','Z','CHILD','NUM_ENDS']   #parameters not useful or unique for Radius comparision
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
            if param != 'FNAME' and param != 'ARCHIVE' and param != 'TYPE':
                temp = []
                for archive in param_data[connect][comp_type]:                    #appends archive file data as single list
                    temp.append(param_data[connect][comp_type][archive][param])   
                temp = list(flatten(temp))
                comb_data[connect][comp_type][param] = temp
        comb_df[connect][comp_type] = pd.DataFrame(comb_data[connect][comp_type]) 

'''Separation by Connection to Soma'''
'''
new_basal = {}; new_apical = {}
for i in param_data:
    new_basal[i] = {}

for i in new_basal:
    for param in params:
        new_basal[i][param] = []
        for archive in dirs:
            new_basal[i][param].append(param_data[i]['Apical'][archive][param])
        new_basal[i][param] = list(flatten(new_basal[i][param]))

simple_plot([new_basal['Indirect'],new_basal['Direct']],['PARENT_RAD','RADIUS'],[1,1],'Apical Hippocampus Compartments_SC', legend = ['Continuing','Initial'], labels = ['Parent Radius','Radius'])
'''

'''Separation by Archive'''
'''
new_basal = {}; new_apical = {}
for i in dirs:
    new_basal[i] = {}

for i in new_basal:
    for param in params:
        new_basal[i][param] = []
        for connect in param_data:
            new_basal[i][param].append(param_data[connect]['Basal'][i][param])
        new_basal[i][param] = list(flatten(new_basal[i][param]))

#simple_plot([new_basal[dirs[0]],new_basal[dirs[1]]],['PARENT_RAD','RADIUS'],[1,1],'Basal Hippocampus Compartments', legend = list(dirs), labels = ['Parent Radius','Radius'])
'''
        
'''Parameter Correlation (Pearson's r) to Radius For Linear Relationships'''
'''
initial_corr = {}
excluded = ['DIRECT_SOMA','NODE_ORDER','NODE_COUNT','FNAME','ARCHIVE','TYPE']       #not useful for nodes directly connected to soma

for connect in comb_data:
    initial_corr[connect] = {}
    for comp_type in comb_data[connect]:
        if connect == 'Direct':
            initial_corr[connect][comp_type] = corr(comb_data[connect][comp_type], [i for i in params.keys() if i not in excluded])
        else:
            initial_corr[connect][comp_type] = corr(comb_data[connect][comp_type], [i for i in params.keys() if i != 'DIRECT_SOMA' and i != 'FNAME' and i != 'ARCHIVE' and i != 'TYPE'])
'''
#do we need to limit which parameters go into correlation

'''
test_dict = {key:val for key,val in comb_data['Direct']['Basal'].items() if key != 'DIRECT_SOMA' and key != 'NODE_ORDER' and key != 'NODE_COUNT'}
columns = ['Branch Length','Horton Strahler', 'Node Degree', 'Parent Radius', 'Path Distance', 'Path to End', 'Radius'] #for initial comps
#columns = ['Branch Length','Horton Strahler', 'Node Count', 'Node Degree', 'Node Order', 'Parent Radius', 'Path Distance', 'Path to End', 'Radius'] #for continuing comps
test_df = pd.DataFrame(test_dict)
try_corr = test_df.corr()
fig,ax = plt.subplots()
plt.rc('font', size = 15)
#fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
sns.set(font_scale = 1.4)
ax = sns.heatmap(try_corr,vmin = -1,vmax = 1, center = 0,cmap = sns.diverging_palette(20, 220, n=256),square = True,annot = True,cbar_kws={'label':"Pearson's R"})
#ax.set_xticklabels(test_df.keys(), rotation=45, horizontalalignment='right')
ax.set_yticklabels(columns)
ax.set_xticklabels(columns, rotation=45, horizontalalignment='right')
plt.title('Initial Basal Compartments')
plt.show()
'''

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
    #simple_plot(comb_data['Indirect']['Apical'], [param,'RADIUS'], 1.0, 'Basic Plot : ' + param + ' to Radius')
    #complete_dict does not remove axon compartments
#_3d_plot(param_data, params.keys())

'''Initiate Fit and Residuals for Best Feature to Estimate Radius'''
#DA_PR, DA_PR_fit = initial_fit(comb_data['Direct']['Apical'], ['PARENT_RAD'],['RADIUS'], [func0,'y = mx + b'])
#IA_PR, IA_PR_fit = initial_fit(comb_data['Indirect']['Apical'], ['PARENT_RAD'], ['RADIUS'], [func0,'y = mx + b'])
#for param in IA_PR.keys():
    #simple_plot(IA_PR,[param,'PARENT_RAD_Res'], 1.0, 'Residual Plot : ' + param + ' to Parent Radius Residual')
    
'''Utilize statsmodels.ols for Multiple Regression and send in possible Transformed Feature Values'''
#fit_data = {}                      #will hold selected features to fit to Radius
#for param in params:
    #fit_data[param] = comb_data['Indirect']['Apical'][param]


#model,predictions = ols_fit(fit_data,['PARENT_RAD','trans_ND'],'RADIUS') #can send in multiple X values to estimate Y
#model,predictions = ols_fit(fit_data,['PARENT_RAD','trans_ND'],'RADIUS','constant')


#fit_data['comb_feat2'] = [1/(1 + np.exp(ND*PTE)) if ND*PTE < np.log(np.finfo('d').max) else 1 for ND,PTE in zip(fit_data['NODE_DEGREE'],fit_data['PATH_TO_END'])]

#fit_data[new_try] = [1/j if j >= 1 else 1 for j in fit_data[i]]
#fit_data[new_try] = [np.log(j) if j >= 1 else 0 for j in fit_data[i]]
#fit_data[new_try] = [1/np.exp(j) if j < np.log(np.finfo('d').max) else 0 for j in fit_data[i]] #high values gives overflow error: larger values will evaluate towards 0
#fit_data[new_try] = [1/(1+np.exp(j)) for j in fit_data[i]]
#fit_data[new_try] = [np.sqrt(j) for j in fit_data[i]]

#for i in fit_data1:
    #maxed = max(fit_data[i])
    #new_try = str(i) + '_try'
    #fit_data[new_try] = [(dp * maxed)/((dp/2) * val) for dp,val in zip(fit_data['PARENT_RAD'],fit_data[i])]
    #fit_data[new_try] = [(dp * maxed)/((dp/2) + val) for dp,val in zip(fit_data['PARENT_RAD'],fit_data[i])]
    #fit_data[new_try] = [(dp * val)/((dp/2) * val) for dp,val in zip(fit_data['PARENT_RAD'],fit_data[i])]
    #fit_data[new_try] = [(dp * val)/((dp/2) + val) for dp,val in zip(fit_data['PARENT_RAD'],fit_data[i])]
    #fit_data[new_try] = [np.log(val) if val >= 1 else 0 for val in fit_data[i]]
    #model,predictions = ols_fit(fit_data,['PARENT_RAD',new_try],'RADIUS')
    #model,predictions = ols_fit(fit_data,[new_try],'RADIUS')
    #model,predictions = ols_fit(fit_data,['PARENT_RAD',new_try],'RADIUS','constant')
    #model,predictions = ols_fit(fit_data,[new_try],'RADIUS','constant')

#for i in params:
    #new_try = str(i) + '_try'
    #fit_data[new_try] = [1/j if j >= 1 else 1 for j in fit_data[i]]
    #model,predictions = ols_fit(fit_data,['PARENT_RAD',new_try],'RADIUS')


'''Plot Transformed Values with Radius or Residuals to better Select Equation Features'''
#fit_data['residuals'] = [yval - est for yval,est in zip(fit_data['RADIUS'],predictions)]
#simple_plot(fit_data, ['RADIUS','residuals'], 'Residual', _alpha = 0.3)

#for param in fit_data:
    #simple_plot(fit_data, ['PATH_TO_END', residuals, param], 'Color3d')

'''Separate Training and Testing Data and Determine Single Equation Coefficients'''

#random.seed(14)

file_copies = list(file_list) 
random.Random(14).shuffle(file_copies) #how would you split up data from single file...
#random.seed(14)
#random.shuffle(file_copies)
split_files = split_seq(file_copies,2)
training_files = split_files[0]
testing_files = split_files[1]


#num_comps = len(basal_dict['CHILD'])
'''
comp_list = [i for i in range(len(basal_dict['CHILD']))]
copy_list = list(comp_list)
random.Random(14).shuffle(copy_list)
split_files = split_seq(copy_list,2)
training_files = split_files[0]
testing_files = split_files[1]
'''
training_basal = {}
testing_basal = {}

for param in header.keys():
    training_basal[param] = []
    testing_basal[param] = []
for num,i in enumerate(basal_dict['FNAME']):
    #for num in copy_list:
    if i in training_files:
        for param in header.keys():
            training_basal[param].append(basal_dict[param][num])
    else:
        for param in header.keys():
            testing_basal[param].append(basal_dict[param][num])
        
train_basal = training_basal
test_basal = testing_basal
#train_apical['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_apical['DIRECT_SOMA'],train_apical['PARENT_RAD'])]
#train_apical['Cont_PR'] = [SC*PR for SC,PR in zip(train_apical['DIRECT_SOMA'],train_apical['PARENT_RAD'])]
#train_apical['Init_logPD'] = [(1-SC)*np.log(PD) for SC,PD in zip(train_apical['DIRECT_SOMA'],train_apical['PATH_DIS'])]
train_basal['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PARENT_RAD'])]
train_basal['Cont_PR'] = [SC*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PARENT_RAD'])]
model,predictions = ols_fit(train_basal,['Init_PR','Cont_PR'],'RADIUS')
train_coef = model.params
save_file(train_coef,'Basal Ganglia Coeff')
train_basal['predictions'] = predictions
#test_apical['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*SC*PR + train_coef[2]*SC*np.log(PD) for SC,PR,PD in zip(test_apical['DIRECT_SOMA'],test_apical['PARENT_RAD'],test_apical['PATH_DIS'])]
test_basal['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*SC*PR for SC,PR in zip(test_basal['DIRECT_SOMA'],test_basal['PARENT_RAD'])]
f_stat = '{:.2e}'.format(model.fvalue)
r_val = '{:.4}'.format(model.rsquared_adj)
additions = 'F-Statistic : ' + f_stat + ' Adjusted R-Squared : ' + r_val
simple_plot([train_basal,test_basal],['RADIUS','predictions'],[1,0.3],'Parent Radius Predicted vs. Actual: Basal Basal Ganglia', legend = ['Training','Testing'], add = additions,labels = ['Radius','predictions'])


