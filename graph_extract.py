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

#####statsmodels requires python3#####

#Usage:  python3 -i morph_feature_analysis.py --path /path/to/folder/with/archives 
#        where path is name of path to folder containing archive folder(s) with .CNG_extract.txt files
#        .CNG_extract.txt files created from morph_feature_extract.py

#George Mason University
#Jonathan Reed
#September 18, 2020

import numpy as np
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from itertools import chain
import statsmodels.api as sm
import random
import seaborn as sns
import argparse
import os
import glob
import re

'''Save plot as png'''            
def save_png(png, title):               #default to save plots instead of plot to window
    png.savefig(title)
    print('File Created : ' + str(title))
    png.close()
    
'''Main Plot Function for Basic, Fit or 3d Color Plots'''
def main_plot(data, var, title, ax_titles = None, labels = None, save_as = None, plot_type = None, fit_line = None, add = None, where = None):
    plt.rc('font', size = 28)                          #default plot and font sizes
    fig, ax = plt.subplots(figsize = (20,10))          
    if add:                                            #additional information to plot
        at = AnchoredText(add,
                          prop=dict(size=20), frameon=True,
                          loc='upper center')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") 
        ax.add_artist(at)
       
    if plot_type == 'Color3d':                     #plot 3 variables(x,y,z) with variable 'z' using colorbar
        scat = ax.scatter(data[var[0]], data[var[1]], c=data[var[2]], s=100, marker='o', label = labels)
        fig.colorbar(scat)
        
    elif plot_type == 'Fit':                       #plot initial fitted equation to data values
        plt.plot(data[var[0]], data[var[1]], 'o')
        plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) 
    
    else:
        for num in range(len(data)):
            if labels:
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', label = labels[num], markersize = 15)
            else:
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', markersize = 15)
            
    if labels or fit_line:
        plt.legend(loc = where) if where else plt.legend()
    plt.xlabel(ax_titles[0]) if ax_titles else plt.xlabel(var[0])
    plt.ylabel(ax_titles[1]) if ax_titles else plt.ylabel(var[1])
    plt.title(title)                 #save_as can be separate from plot title
    save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')

'''Plot Parameter Correlation to Diameter with Pearson's R'''
def corr(data,selection,title,save_as = None):
    dframe = pd.DataFrame(data,columns = selection.keys())  #dataframe reorders alphabetically columns by default 
    dcorr = dframe.corr()                                   #set to order columns as in select_params - passed as 'selection'
    fig,ax = plt.subplots(figsize = (20,15))
    sns.set(font_scale = 2.8)
    plt.rc('font', size = 26)
    ax = sns.heatmap(dcorr,vmin = -1,vmax = 1, center = 0,cmap = sns.diverging_palette(20, 220, n=256),square = True,annot = True,cbar_kws={'shrink':1,'label':"Pearson's R"})
    ax.set_yticklabels([sub['short'] for sub in selection.values()],fontsize = 30)
    ax.set_xticklabels([sub['short'] for sub in selection.values()],rotation = 0,horizontalalignment = 'right',fontsize = 30,ha = 'center')
    ax.figure.axes[-1].yaxis.label.set_size(26)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 26)
    plt.title(title)   
    #save_png(plt, title)
    save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')
    
'''Function Definitions for Fit Equations'''
#fit selected variables to equations
#utilized with initial_fit function
def func0(x,m,b):                                     
    return m*x+b

'''Fit Data to Selected Equation and Obtain Residuals'''
#First method to plot predicted values and find residuals
def initial_fit(data, var, func):
    x = var[0]
    y = var[1]                             #assumes only fit to single y-variable
    popt, pcov = optimize.curve_fit(func[0],data[x],data[y])    #fits equation to choosen function if possible
    print('This is fit estimate for : ', func[1])
    print('popt',popt)
    print('pcov',pcov)
    print('Where x = ', x, ' and y = ', y) 
 
    '''Plot Data to Function and Find Residuals'''
    for i in data[x]:
        temp_max = round(max(data[i])) 
        temp_min = round(min(data[i]))
        x_range = np.arange(temp_min, temp_max + 1, 1)        #plot original feature values to line equation
        main_plot(data, [i, y], 'Fitted Equation ' + i + ' to ' + y, fit_line = [x_range, func[0](x_range, *popt), func[1]])
        res_label = str(i) + '_Res'  
        predictions = [func[0](xval,*popt) for xval in data[i]]     #plot original feature values to predicted feature values
        main_plot(data, [y, predictions], 'Predicted vs. Actual Plot : ' + i)
                                                              #plot original feature values to residuals
        data[res_label] = [yval - est for yval,est in zip(data[y],predictions)]
        print('Residual Plot : ', i, ' to ', y)
        main_plot(data, [y, res_label], 'Residual Plot : ' + i)
        
    return data, [popt,pcov]

'''OLS Regression to fit Variables'''
#Second PREFERRED Method to predict variable values
def ols_fit(data,xlabel,ylabel,constant = None,extended = None):
    temp_df = pd.DataFrame(data)   #xlabel can be multiple names of parameters to fit towards Diameter
    X = temp_df[xlabel]
    Y = temp_df[ylabel]
    if constant:                   #by default no constant or y-intercept fitted to equation
        X = sm.add_constant(X)      
    model = sm.OLS(Y,X).fit()
    if extended:
        print(model.summary())     #by default no extended output printed within function
    return model,model.predict(X)

def regress_val(model):
    f_stat = '{:.2e}'.format(model.fvalue)
    r_val = '{:.4}'.format(model.rsquared_adj)
    aic = '{:.4}'.format(model.aic)
    bic = '{:.4}'.format(model.bic)
    cond_num = '{:.4}'.format(model.condition_number)
    vals = [f_stat,r_val,cond_num,aic,bic]
    return vals

'''Save equation coefficients and update file'''
def add_coeff(temp_list, temp_dict, comp_type):
    for i in temp_dict:
        to_add = [comp_type + '_' + i, temp_dict[i]]   #adds feature name and relation to soma designation
        temp_list.append(to_add)
    return temp_list

'''Flatten nested lists of values into single list of values'''
def flatten(container):    
    for i in container:                 
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

'''Split N sequences into 'equal' sizes'''
def split_seq(seq, size):               #splits file_list into equal amounts for testing-training datasets
        newseq = []
        splitsize = 1.0/size*len(seq)
        for i in range(size):
                newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
        return newseq

'''Save data to file'''    
def save_file(coeff,filename):
    with open(filename + '.txt','a') as f:
        for i in coeff:
            f.write(str(i) + ' ')
        f.write('\n')

        
'''Start of Working Code'''
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   
args = parser.parse_args()                      
path = args.path                                

'''Locates Archive Data'''
root, dirs, files = list(os.walk(path))[0]           #all CNG_extract.txt files require identical parameters and length data
if len(dirs) == 0:
    if os.path.basename(path) == '':
        path = root.rstrip('/')                      #fix for trailing '/' in path
    dirs = [os.path.basename(path)]
    root = os.path.dirname(path) + '/'

archive_dict = {}; header = {}                                
file_list = []; complete_list = []

'''Read in Morphology Data in Folder with Single/Multiple Archive(s)'''
for d1 in dirs:                                      #can accept individual archives or multiple archives within single directory           
    fullpath = root + d1 + '/*CNG_extract.txt'                      
    data_list = []; apical_list = []; basal_list = []               
    print('Working on Directory : ', str(d1))
    for fname in glob.glob(fullpath):                               #locates and loops through _extract files in path
        temp_name = re.search(d1 + '/(.*).txt', fname)
        file_list.append(temp_name.group(1))
        with open(fname) as f:
            for line in f:
                if line.strip():                                    #removes any empty lines
                    if '*C' in line and not header:                 #finds feature names as header line and separately saves from data
                        if 'XYZ' in line:                           
                            line = re.sub('XYZ','X; Y; Z', line)    #fix for header length by separating XYZ into X,Y,Z 
                        line = line.split('; ')                     
                        for num, val in enumerate(line):            
                            if '\n' in val:                         #fix for trailing '\n' attached to final parameter in header
                                val = re.sub('\n','',val)
                            if val == '*CHILD':                     #fix for leading '*' to indicate text in .swc file
                                header['CHILD'] = num
                            else:
                                header[val] = num
                            
                    elif line[0] != '*' and line[0] != ' ' and line[0] != '/n':  #organizes remaining parameter data values    
                        temp_line = line.split()
                        for point, val in enumerate(temp_line):                  #Child, Type, Parent are DESCRETE values defined by .swc
                            if point != header['CHILD'] and point != header['TYPE'] and point != header['PARENT']:
                                temp_line[point] = float(temp_line[point])
                        temp_line.extend([temp_line[header['RADIUS']]*2,temp_line[header['PARENT_RAD']]*2,temp_name.group(1),d1]) 
                        data_list.append(temp_line)              #adding new values of Diameter, Parent Diameter, File name, Archive name
                        complete_list.append(temp_line)          #complete_list will be used to test equations changing radius
                                                                                                                                    
    '''Initial Data Separation by .swc Compartment Types (Basal/Apical if present) and Archive'''
    if len(temp_line) > len(header):                             
        for i in ['DIAMETER','PARENT_DIA','FNAME','ARCHIVE']:    #header with new added parameters
            header[i] = len(header)
    
    for line in data_list:                                  #'3' - basal and '4' - apical .swc compartment type
        if line[header['TYPE']] == '4':                     #currently does not keep axon data
            apical_list.append(line)                        
        elif line[header['TYPE']] == '3':
            basal_list.append(line)                        

    if not 'Basal' in archive_dict.keys():                  #assumes basal compartments present in all .swc iles
        archive_dict['Basal'] = {}                          #archive_dict will organize data by archive and by compartment type
    archive_dict['Basal'][d1] = {}                          #individual data values will be stored in parameter lists 
    for param,vals in zip(header.keys(),list(zip(*basal_list))):
        archive_dict['Basal'][d1][param] = vals
        
    if apical_list:                                         #add apical values if present in .swc file
        if not 'Apical' in archive_dict.keys():
            archive_dict['Apical'] = {}
        archive_dict['Apical'][d1] = {} 
        for param,vals in zip(header.keys(),list(zip(*apical_list))):
            archive_dict['Apical'][d1][param] = vals
    
'''Organize Data by Parent Connection - Initial, Branch Children, Continuing Compartments'''
comp_types = [i for i in archive_dict.keys()]
comp_list = ['Initial','BP_Child','Continuing','All']
complete_dict = {}; sep_data = {}; comb_data = {}  #create dictionaries of separate or combined archive morphologies
for i in comp_types:
    complete_dict[i] = {}
    for connect in comp_list:                                
        if len(sep_data.keys()) != len(comp_list) and len(comb_data.keys()) != len(comp_list):
            sep_data[connect] = {}; comb_data[connect] = {}
        sep_data[connect][i] = {}; comb_data[connect][i] = {}
        for param in header.keys():
            comb_data[connect][i][param] = []
            if len(complete_dict[i].keys()) != len(header.keys()):
                complete_dict[i][param] = []
            for d1 in dirs:
                if connect == comp_list[0]:
                    complete_dict[i][param].extend(archive_dict[i][d1][param])    #fill in values for complete_list only once
                if len(sep_data[connect][i].keys()) != len(dirs):                 #create new sep_data dictionaries if not present
                    sep_data[connect][i][d1] = {}
                sep_data[connect][i][d1][param] = []         
                if connect == comp_list[-1]:                                      #fill in values after dictionaries created
                    for num,val in enumerate(archive_dict[i][d1]['PAR_CONNECT']):
                        sep_data['All'][i][d1][param].append(archive_dict[i][d1][param][num])
                        comb_data['All'][i][param].append(archive_dict[i][d1][param][num])
                        for cnum,ctype in enumerate(comp_list):                                         #'Initial'    - 0
                            if val == cnum:                                                             #'BP_Child'   - 1
                                sep_data[ctype][i][d1][param].append(archive_dict[i][d1][param][num])   #'Continuing' - 2
                                comb_data[ctype][i][param].append(archive_dict[i][d1][param][num])

'''Select parameters to analyze with Compartment Diameter'''
select_params = {'NODE_ORDER': {'short':'IB','long':'Initial Branch Order'},
                 'PATH_TO_END': {'short':'LP','long':'Longest Path to Terminal'},
                 'PARENT_DIA': {'short':'PD','long':'Parent Diameter'},
                 'PATH_DIS': {'short':'PS','long':'Path To Soma'},
                 'NODE_DEGREE': {'short':'TB','long':'Terminal Branch Order'},
                 'BRANCH_LEN': {'short':'TD','long':'Total Dendritic Length'},
                 'DIAMETER': {'short':'D','long':'Diameter'}}

initial_params = {key:value for (key,value) in select_params.items() if key != 'NODE_ORDER'}
#may have to make initial_params to not include Node_Order for plots as all same value

'''Initial Plots and Parameter Analysis''' #save back in folder with fullpath os.cwd

#could possibly save in separate location
if not dirname in os.listdr(rootdir):
    os.mndir(rootdir + dirname) #check for '/'
    #rootdir + 'images'
    #rootdir as os.getcwd()
    #os.basename(fullpath) 
    
for i in comp_types:
    for param in select_params:
        if param != 'DIAMETER':
            d1_data = []; label_list = []
            plabel = select_params[param]['long']
            title = i + ' ' + plabel + ' to Diameter'
            saving = 'All ' + i + ' ' + plabel + ' to Diameter'

            '''Plot Parameters to Diameter'''
            main_plot([comb_data['All'][i]], [param,'DIAMETER'], title, ax_titles = [plabel,'Diameter'], save_as = saving)
            
            '''Plot Parameters by Parent Connection'''
            main_plot([comb_data['Continuing'][i],comb_data['BP_Child'][i],comb_data['Initial'][i]], [param,'DIAMETER'], title, ax_titles = [plabel,'Diameter'], labels = ['Continuing','Branch Children', 'Initial'], where = 'upper right', save_as = saving + ' PC')
            
            '''Plot Parameters by Archive'''
            for d1 in dirs:
                d1_data.append(sep_data['All'][i][d1])
                temp,_ = pearsonr(sep_data['All'][i][d1][param],sep_data['All'][i][d1]['DIAMETER'])
                rsquared = round(temp**2,4)
                label = d1 + ' $r^2$ = ' + str(rsquared)
                label_list.append(label)
            main_plot(d1_data, [param,'DIAMETER'], title, ax_titles = [plabel,'Diameter'], labels = label_list, where = 'upper right', save_as = saving + ' A')

            
'''Extended Plots with Relation to Soma'''
regress = {}                      #setup regression dictionary
for i in comp_types:
    regress[i] = {}
    for connect in comp_list:     #setup correlation data for each compartment type
        param_comb = []; regress[i][connect] = []
        select_data = {}
        selection = initial_params if connect == 'Initial' else select_params
        for p1 in selection:
            
            '''Multiple Regression for Diameter''' 
            if p1 != 'DIAMETER':  #runs statsmodels OLS (multiple) regression for single and paired parameters
                model,predictions = ols_fit(comb_data[connect][i],p1,'DIAMETER')
                vals = regress_val(model)
                temp = [p1]
                temp.extend(vals)
                regress[i][connect].append(temp)
                for p2 in selection:          #checks for repeat combinations of parameters for easier analysis
                    if p2 != 'DIAMETER' and p2 != p1 and str(p2 + ' + ' + p1) not in param_comb:
                        param_comb.append(str(p1 + ' + ' + p2))
                        model,predictions = ols_fit(comb_data[connect][i],[p1,p2],'DIAMETER')
                        vals = regress_val(model)
                        temp = [str(p1 + ' + ' + p2)]
                        temp.extend(vals)
                        regress[i][connect].append(temp)
                        
            '''Plot Correlations to Diameter''' 
            select_data[p1] = comb_data[connect][i][p1]
        title = i + ' ' + connect + ' Feature Correlation' 
        corr(select_data, selection, title, save_as = title)  

#print regression data and sort by r-squared
for i in regress:
    for connect in regress[i]:
        print(i,connect)
        print('Parameter','F-Stat','R-Squared','Condition Number','AIC','BIC')
        regress[i][connect].sort(key=lambda x: x[2])
        for test in regress[i][connect][::-1]:
            print(test)



'''Additional Tools for Analysis'''
#connect = 'Continuing'; comp_type = 'Basal'
'''Transform Values after Visualization of Non-Linear Trends'''
#if Non-Linear trends appear in parameter plots to Radius/Diameter, can add transformed values as new parameter
'''
old_param = 'PATH_DIS'; new_param = 'Log_' + old_param
comb_data[connect][comp_type][new_param] = [np.log(i) for i in comb_data[connect][comp_type][old_param]]
for param in [old_param,new_param]:
    model,predictions = ols_fit(comb_data[connect][comp_type],param,'DIAMETER',extended = True)
    main_plot(comb_data[connect][comp_type],[param,'DIAMETER'], param + ' to ' + 'DIAMETER') #see if trend becomes linear
    if param == new_param:
        comb_data[connect][comp_type][new_param + '_Res'] = [yval - est for yval,est in zip(comb_data[connect][comp_type][new_param],predictions)]
        main_plot(comb_data[connect][comp_type],[param + '_Res','DIAMETER'], param + '_Res' + ' to ' + 'DIAMETER') #plot residuals
    #can extend residuals to remaining parameters
'''

'''Plot Two Variables with Diameter with Colorbar'''
'''
var = ['PARENT_DIA','DIAMETER','PATH_DIS'] #var as x,y,z
main_plot(comb_data[connect][comp_type], var, title, ax_titles = ['Parent Diameter','Diameter'], labels = var[2], plot_type = 'Color3d')
'''

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
'''
file_copies = list(file_list) 
random.Random(14).shuffle(file_copies) 
split_files = split_seq(file_copies,2)
training_files = split_files[0]
testing_files = split_files[1]

training_set = {} #will contain only needed parameters to train model
coeffs = {} #will contain coefficients from model training
training_all = {} #will contain all training values to verify fit
testing_all = {}  #will contain all testing values to verify fit
#testing to see if we can select which parameters to find and then iterate through regression code
test_params = {}
test_params['Apical'] = {'Initial':['PARENT_RAD','PATH_DIS'],'BP_Child':['PARENT_RAD','PATH_DIS'],'Continuing':['PARENT_RAD']}
test_params['Basal'] = {'Initial':['PARENT_RAD','PATH_DIS'],'BP_Child':['PARENT_RAD','PATH_DIS'],'Continuing':['PARENT_RAD']}
'''

#how to iterate between different comp. designations as some (1-SC)
#Testing complete_dict and compartment designation as multiplier
'''
for param in complete_dict:
    training_all[param] = []
    testing_all[param] = []

for num,fname in enumerate(complete_dict['FNAME']):
    for param in complete_dict:
        if fname in training_files:
            training_all[param].append(complete_dict[param][num])
        elif fname in testing_files:
            testing_all[param].append(complete_dict[param][num])
'''
'''
for comp_type in test_params:
    for connect in test_params[comp_type]:
        if connect == 'Initial': #(1-SC) for SC in complete_dict['DIRECT_SOMA']
        elif connect == 'BP_Child': #(BPC) for BPC in complete_dict['BP_CHILD']
        elif connect == 'Continuing': #(1-BPC) for BPC in complete_dict['BP_CHILD']
            
        model,predictions = ols_fit(training_set[comp_type][connect],test_params[comp_type][connect],'RADIUS',extended = True)

        #setup model for each compartment designation --> check if using complete_dict gives same results as comb_data
'''
'''
predict = {}
for comp_type in test_params:
    training_set[comp_type] = {}
    coeffs[comp_type] = {}
    predict[comp_type] = {'Train':[],'Test':[]}
    #for val in complete_dict:
    for connect in test_params[comp_type]:
        training_set[comp_type][connect] = {}
        coeffs[comp_type][connect] = {}
        for param in params:
            training_set[comp_type][connect][param] = []
            for num,i in enumerate(comb_data[connect][comp_type]['FNAME']):
                if i in training_files:      #for now, all parameter values available in training or testing sets
                    training_set[comp_type][connect][param].append(comb_data[connect][comp_type][param][num])
        #start multiple regression
        print('**************** ' + comp_type + ' ' + connect + ' ****************')
        model,_ = ols_fit(training_set[comp_type][connect],test_params[comp_type][connect],'RADIUS',extended = True)
        coeffs[comp_type][connect] = [coeff for coeff in model.params]
        for training_all, testing_all#...
        for num,val in enumerate(complete_dict):
            predict[comp_type]['Train'][num] = []
            if connect == 'Initial':
                newval = sum[newlist]
            elif connect == 'BP_Child':
            elif connect == 'Continuing':
'''
        

#somehow combine the values in complete_dict with coefficients and parameters from test_params
#hmm also have to combine with the training and testing sets...
#testing_set[comp_type][connect]['train'] = predictions



'''        
temp_list = []; num_lists = []
for num,test_param in enumerate(test_params[comp_type][connect]):
    temp = [coeffs[num] * val for val in testing_set[comp_type][connect][test_param]]
    temp_list.append(temp)
    num_lists.append(temp_list[num])
    #temp_list[num].append(temp)
    
total_list = list(chain(num_lists))  #<--- so close.....
testing_set[comp_type][connect]['test'] = [coeff[i] * test_param[i] for coeff,test_param in zip(coeffs,test_params)]
'''
#(1-SC) - Initial
#(BPC) - Branch Child
#(1-BPC) - Continuing
#model,predictions = ols_fit(comb_data[comb_type]['All'],[test1,test2],'RADIUS',extended = True)
#will need to add with test and train sets


#so I think it's working properly --> check test set with coefficients from train set
#somehow automate and check to see predictions to Radius values for fit within dictionary structure
#train_set['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_set['DIRECT_SOMA'],train_set['PARENT_RAD'])]
#train_set['Init_ND'] = [(1-SC)*ND for SC,ND in zip(train_set['DIRECT_SOMA'],train_set['NODE_DEGREE'])]
#model,predictions = ols_fit(train_set['comp_type'],['PARENT_RAD','NODE_DEGREE'],'RADIUS')                    

'''            
for connect in comb_data:
    for num,i in enumerate(comb_data[connect][comp_type]['FNAME']):
        if i in training_files:
            for param in params:
                training_set[connect][param].append(comb_data[connect][comp_type][param][num])
        else:
            for param in params:
                testing_set[connect][param].append(comb_data[connect][comp_type][param][num])
'''
'''
train_set = training_set 
test_set = testing_set
'''
'''Find Best Models to Predict Radius''' #compare across f-value, r2, and condition #
'''
#cell_type = 'Apical'
cell_type = 'Basal'
#cell_type = 'Striatal'
#cell_type = 'Purkinje'
comp_list = ['Direct','BP_Child','Non_BP']
for comp_type in comp_list: 
    #comp_type = 'Direct'
    save_name = cell_type + ' ' + comp_type + ' ' + 'Models'
    for i in cont_params:
        temp_list = []
        model,predictions = ols_fit(train_set[comp_type],['PARENT_RAD',i],'RADIUS')
        temp_dict = dict(model.params)
        for i in temp_dict:
            temp_list.append(i)
            #temp_list.append(temp_dict[i])
        temp_list.extend((model.fvalue, model.rsquared_adj,model.condition_number))
        #print(new_name,temp_list)
        save_file(temp_list, save_name)
'''
'''Save Coefficients to File'''

'''Striatal Dendrites'''
'''
temp_list = []
comp_type = 'Direct'
#train_set['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_set['DIRECT_SOMA'],train_set['PARENT_RAD'])]
#train_set['Init_ND'] = [(1-SC)*ND for SC,ND in zip(train_set['DIRECT_SOMA'],train_set['NODE_DEGREE'])]
model,predictions = ols_fit(train_set['comp_type'],['PARENT_RAD','NODE_DEGREE'],'RADIUS')
'''
'''
temp_dict = dict(model.params)
for i in temp_dict:
    add_name = comp_type + '_' + i
    temp_list.append((add_name, temp_dict[i]))
'''
'''
temp_list = add_coeff(temp_list,dict(model.params),comp_type)   #for some reason, more accurate to keep separate OLS per compartment type than all together
for i in temp_list:
    save_file(i,'SAVING_TEST')
'''
#save_file(model.params,'Striatal_Test')

#train_set['BP_Child_PR'] = [(1-BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['BP_Child_NO'] = [(1-BP)*NO for BP,NO in zip(train_set['BP_CHILD'],train_set['NODE_ORDER'])]
#model,predictions = ols_fit(train_set['BP_Child'],['PARENT_RAD','NODE_ORDER'],'RADIUS')
#save_file(model.params,'Striatal_Test')

#train_set['Non_BP_PR'] = [(BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['Non_BP_PD'] = [(BP)*PD for BP,PD in zip(train_set['BP_CHILD'],train_set['PATH_DIS'])]
#model,predictions = ols_fit(train_set['Non_BP'],['PARENT_RAD','PATH_DIS'],'RADIUS')
#save_file(model.params,'Striatal_Test')

#model,predictions = ols_fit(train_set,['Init_PR','Init_ND','BP_Child_PR','BP_Child_NO','Non_BP_PR','Non_BP_PD'],'RADIUS')

#train_coef = model.params
#train_set['predictions'] = predictions
#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*ND + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*NO + train_coef[4]*(BP)*PR + train_coef[5]*(BP)*PD for SC,BP,PR,ND,NO,PD in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['NODE_DEGREE'],test_set['NODE_ORDER'],test_set['PATH_DIS'])]


'''Purkinje Dendrites'''
'''
#train_set['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_set['DIRECT_SOMA'],train_set['PARENT_RAD'])]
#train_set['Init_PE'] = [(1-SC)*PE for SC,PE in zip(train_set['DIRECT_SOMA'],train_set['PATH_TO_END'])]
model,predictions = ols_fit(train_set['Direct'],['PARENT_RAD','PATH_TO_END'],'RADIUS')
save_file(model.params,'Purkinje_Test')

#train_set['BP_Child_PR'] = [(1-BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['BP_Child_PD'] = [(1-BP)*PD for BP,PD in zip(train_set['BP_CHILD'],train_set['PATH_DIS'])]
model,predictions = ols_fit(train_set['BP_Child'],['PARENT_RAD','PATH_DIS'],'RADIUS')
save_file(model.params,'Purkinje_Test')

#train_set['Non_BP_PR'] = [(BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['Non_BP_NO'] = [(BP)*NO for BP,NO in zip(train_set['BP_CHILD'],train_set['NODE_ORDER'])]
model,predictions = ols_fit(train_set['Non_BP'],['PARENT_RAD','NODE_ORDER'],'RADIUS')
save_file(model.params,'Purkinje_Test')


#model,predictions = ols_fit(train_set,['Init_PR','Init_PE','BP_Child_PR','BP_Child_PD','Non_BP_PR','Non_BP_NO'],'RADIUS')

#train_coef = model.params
#train_set['predictions'] = predictions
#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*PE + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PD + train_coef[4]*(BP)*PR + train_coef[5]*(BP)*NO for SC,BP,PR,PE,PD,NO in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['PATH_TO_END'],test_set['PATH_DIS'],test_set['NODE_ORDER'])]
'''

'''Basal Dendrites'''
'''
#train_set['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_set['DIRECT_SOMA'],train_set['PARENT_RAD'])]
#train_set['Init_BL'] = [(1-SC)*BL for SC,BL in zip(train_set['DIRECT_SOMA'],train_set['BRANCH_LEN'])]
model,predictions = ols_fit(train_set['Direct'],['PARENT_RAD','BRANCH_LEN'],'RADIUS')
#coefs = model.params
#train_set['Direct']['predictions'] = predictions
#test_set['Direct']['predictions'] = [model.params[0]*(1-SC)*PR + model.params[1]*(1-SC)*PE for SC,PR,PE in zip(test_set['Direct']['DIRECT_SOMA'],test_set['Direct']['PARENT_RAD'],test_set['Direct']['PATH_TO_END'])]
save_file(model.params,'Basal_Test')
save_file(model.params,'Pyramidal_Test')

#train_set['BP_Child_PR'] = [(1-BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['BP_Child_PD'] = [(1-BP)*PD for BP,PD in zip(train_set['BP_CHILD'],train_set['PATH_DIS'])]
model,predictions = ols_fit(train_set['BP_Child'],['PARENT_RAD','PATH_DIS'],'RADIUS')
#coefs = model.params
save_file(model.params,'Basal_Test')
save_file(model.params,'Pyramidal_Test')

#train_set['Non_BP_PR'] = [(BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
model,predictions = ols_fit(train_set['Non_BP'],['PARENT_RAD'],'RADIUS')
save_file(model.params,'Basal_Test')
save_file(model.params,'Pyramidal_Test')


#model,predictions = ols_fit(train_set,['Init_PR','Init_BL','BP_Child_PR','BP_Child_PD','Non_BP_PR'],'RADIUS')

#train_coef = model.params
#train_set['predictions'] = predictions

#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*BL + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PD + train_coef[4]*(BP)*PR for SC,BP,PR,BL,PD in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['BRANCH_LEN'],test_set['PATH_DIS'])]

#save_file(train_coef,'Basal_Model')
'''

'''Apical Dendrites'''
'''
training_set = {}
testing_set = {}
for connect in comb_data:
    training_set[connect] = {}
    testing_set[connect] = {}
    for param in params:
        training_set[connect][param] = []
        testing_set[connect][param] = []

#comp_type = apical_dict
comp_type = 'Apical'
for connect in comb_data:
    for num,i in enumerate(comb_data[connect][comp_type]['FNAME']):
        if i in training_files:
            for param in params:
                training_set[connect][param].append(comb_data[connect][comp_type][param][num])
        else:
            for param in params:
                testing_set[connect][param].append(comb_data[connect][comp_type][param][num])
        
train_set = training_set 
test_set = testing_set

#train_set['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_set['DIRECT_SOMA'],train_set['PARENT_RAD'])]
#train_set['Init_PE'] = [(1-SC)*BL for SC,BL in zip(train_set['DIRECT_SOMA'],train_set['PATH_TO_END'])]
model,predictions = ols_fit(train_set['Direct'],['PARENT_RAD','PATH_TO_END'],'RADIUS')
#coefs = model.params
#train_set['Direct']['predictions'] = predictions
#test_set['Direct']['predictions'] = [model.params[0]*(1-SC)*PR + model.params[1]*(1-SC)*PE for SC,PR,PE in zip(test_set['Direct']['DIRECT_SOMA'],test_set['Direct']['PARENT_RAD'],test_set['Direct']['PATH_TO_END'])]
save_file(model.params,'Apical_Test')
save_file(model.params,'Pyramidal_Test')

#train_set['BP_Child_PR'] = [(1-BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
#train_set['BP_Child_PE'] = [(1-BP)*PE for BP,PE in zip(train_set['BP_CHILD'],train_set['PATH_TO_END'])]
model,predictions = ols_fit(train_set['BP_Child'],['PARENT_RAD','PATH_TO_END'],'RADIUS')
#coefs = model.params
save_file(model.params,'Apical_Test')
save_file(model.params,'Pyramidal_Test')

#train_set['Non_BP_PR'] = [(BP)*PR for BP,PR in zip(train_set['BP_CHILD'],train_set['PARENT_RAD'])]
model,predictions = ols_fit(train_set['Non_BP'],['PARENT_RAD'],'RADIUS')
save_file(model.params,'Apical_Test')
save_file(model.params,'Pyramidal_Test')
'''
'''
#model,predictions = ols_fit(train_set,['Init_PR','Init_BL','BP_Child_PR','BP_Child_PE','Non_BP_PR'],'RADIUS')
#model,predictions = ols_fit(train_set,['Init_PR','Init_PD','BP_Child_PR','BP_Child_PE','Non_BP_PR'],'RADIUS')

#model,predictions = ols_fit(train_set,['Init_PR','Init_ND','BP_Child_PR','BP_Child_PE','Non_BP_PR'],'RADIUS')
#model,predictions = ols_fit(train_set,['Init_PR','Init_PE','BP_Child_PR','BP_Child_PE','Non_BP_PR'],'RADIUS')
#train_coef = model.params
#train_set['predictions'] = predictions

#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*BL + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PE + train_coef[4]*(BP)*PR for SC,BP,PR,BL,PE in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['BRANCH_LEN'],test_set['PATH_TO_END'])]
#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*PD + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PE + train_coef[4]*(BP)*PR for SC,BP,PR,PE,PD in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['PATH_TO_END'],test_set['PATH_DIS'])]

#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*ND + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PE + train_coef[4]*(BP)*PR for SC,BP,PR,PE,ND in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['PATH_TO_END'],test_set['NODE_DEGREE'])]
#test_set['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*(1-SC)*PE + train_coef[2]*(1-BP)*PR + train_coef[3]*(1-BP)*PE + train_coef[4]*(BP)*PR for SC,BP,PR,PE in zip(test_set['DIRECT_SOMA'],test_set['BP_CHILD'],test_set['PARENT_RAD'],test_set['PATH_TO_END'])]
'''
'''
r_list = []#; archive_r = {}
#for i in dirs:
    #archive_r[i] = []
    
for i in training_files:
    temp_pred = []; temp_rad = []
    for num, j in enumerate(train_set['FNAME']):
        if j == i:
            temp_pred.append(train_set['predictions'][num])
            temp_rad.append(train_set['RADIUS'][num])
    corr,_ = pearsonr(temp_pred,temp_rad)
    r_list.append(corr)
    
for i in testing_files:
    temp_pred = []; temp_rad = []
    for num, j in enumerate(test_set['FNAME']):
        if j == i:
            temp_pred.append(test_set['predictions'][num])
            temp_rad.append(test_set['RADIUS'][num])
    corr,_ = pearsonr(temp_pred,temp_rad)
    r_list.append(corr)

ave_r = np.average(r_list)
sd = np.std(r_list)
print(ave_r,sd)

f_label = 'F-Statistic : ' + '{:.2e}'.format(model.fvalue)
r_label = ' Average R : ' + '{:.4}'.format(ave_r)
additions = f_label +  r_label

simple_plot([train_set,test_set],['RADIUS','predictions'],[1,0.3], 'Apical Dendrites', legend = ['Training','Testing'], add = additions,labels = ['Original Radius','Predicted Radius'], where = 'upper left')

save_file(train_coef,'Apical_Model')
'''

'''
[SC*BP*PR + (1-SC)*train_coef[1]*PTE + (1-SC)*train_coef[2]*PR + (1-BP)*train_coef[3]*PTE + (1-BP)*train_coef[4]*PR for SC,BP,PR,PTE in zip(test_basal['DIRECT_SOMA'],test_basal['BRANCH_POINT'],test_basal['PARENT_RAD'],test_basal['PATH_TO_END'])]

#train_basal['Cont_PR'] = [(SC*BP)*PR for SC,BP,PR in zip(train_basal['DIRECT_SOMA'],train_basal['BRANCH_POINT'],train_basal['PARENT_RAD'])]
#train_basal['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PARENT_RAD'])]
#train_basal['Init_PTE'] = [(1-SC)*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PATH_TO_END'])]
#train_basal['Init_PD'] = [(1-SC)*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PATH_DIS'])]
#train_basal['Branch_PR'] = [(1-SC)*PR for SC,PR in zip(train_basal['BRANCH_POINT'],train_basal['PARENT_RAD'])]
#train_basal['Branch_PTE'] = [(1-SC)*PR for SC,PR in zip(train_basal['BRANCH_POINT'],train_basal['PATH_TO_END'])] #can run model regression fine

model,predictions = ols_fit(train_basal,['Cont_PR','Init_PTE','Init_PR','Branch_PTE','Branch_PR'],'RADIUS')      #get errors in predictions --> no plot yet
train_coef = model.params
train_basal['predictions'] = predictions
test_basal = [SC*BP*PR + (1-SC)*train_coef[1]*PTE + (1-SC)*train_coef[2]*PR + (1-BP)*train_coef[3]*PTE + (1-BP)*train_coef[4]*PR for SC,BP,PR,PTE in zip(test_basal['DIRECT_SOMA'],test_basal['BRANCH_POINT'],test_basal['PARENT_RAD'],test_basal['PATH_TO_END'])]
f_stat = '{:.2e}'.format(model.fvalue)
#r_val = '{:.4}'.format(model.rsquared_adj)
corr,_ = pearsonr(test_basal['predictions'],test_basal['RADIUS'])
r_val = corr*2
r_label = ' Adjusted $r^2$ = ' + r_val
additions = 'F-Statistic : ' + f_stat + r_label #' Adjusted R-Squared : ' + r_val
simple_plot([train_basal,test_basal],['RADIUS','predictions'],[1,0.3], 'Apical Hippocampus', legend = ['Training','Testing'], add = additions,labels = ['Radius','Predictions'])
'''

#####Extra Code for Possible Feature Analysis#####



'''Initial Plot Function to compare between different archives'''
'''
def merge_plot(data, to, extras):
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
'''

'''Initiate Fit and Residuals for Best Feature to Estimate Radius'''
#DA_PR, DA_PR_fit = initial_fit(comb_data['Direct']['Apical'], ['PARENT_RAD'],['RADIUS'], [func0,'y = mx + b'])
#IA_PR, IA_PR_fit = initial_fit(comb_data['Indirect']['Apical'], ['PARENT_RAD'], ['RADIUS'], [func0,'y = mx + b'])
#for param in IA_PR.keys():
    #simple_plot(IA_PR,[param,'PARENT_RAD_Res'], 1.0, 'Residual Plot : ' + param + ' to Parent Radius Residual')



#num_comps = len(basal_dict['CHILD'])



'''
training_basal = {}
testing_basal = {}
new_params = [i for i in params.keys() if i == 'NODE_DEGREE' or i == 'PATH_DIS' or i == 'BRANCH_LEN' or i == 'PATH_TO_END' or i == 'PARENT_RAD' or i == 'RADIUS']
for param in new_params:#header.keys():
    training_basal[param] = []
    testing_basal[param] = []
for num,i in enumerate(new_basal['Direct']['FNAME']):
    #(basal_dict['FNAME']):
    #for num in copy_list:
    if i in training_files:
        for param in new_params:#header.keys():
            training_basal[param].append(new_basal['Direct'][param][num])
            #(basal_dict[param][num])
    else:
        for param in new_params:#header.keys():
            testing_basal[param].append(new_basal['Direct'][param][num])
            #basal_dict[param][num])
        
train_basal = training_basal
test_basal = testing_basal
train_basal['Cont_PR'] = [SC*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PARENT_RAD'])]
train_basal['Init_PR'] = [(1-SC)*PR for SC,PR in zip(train_basal['DIRECT_SOMA'],train_basal['PARENT_RAD'])]
train_basal['Init_PD'] = [(1-SC)*PD for SC,PD in zip(train_basal['DIRECT_SOMA'],train_basal['PATH_DIS'])]
train_basal['Init_PTE'] = [(1-SC)*PTE for SC,PTE in zip(train_basal['DIRECT_SOMA'],train_basal['PATH_TO_END'])]
train_basal['Init_BL'] = [(1-SC)*BL for SC,BL in zip(train_basal['DIRECT_SOMA'],train_basal['BRANCH_LEN'])]

train_basal['Init_logPD'] = [(1-SC)*np.log(PD) for SC,PD in zip(train_basal['DIRECT_SOMA'],train_basal['PATH_DIS'])]
model,predictions = ols_fit(train_basal,['PARENT_RAD'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_PR'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_PTE'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_BL'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_PD'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_logPD'],'RADIUS')

model,predictions = ols_fit(train_basal,['Cont_PR','Init_BL','Init_PD'],'RADIUS')
model,predictions = ols_fit(train_basal,['Cont_PR','Init_PTE','Init_PD'],'RADIUS')

train_coef = model.params
train_basal['predictions'] = predictions
#test_apical['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*SC*PR + train_coef[2]*SC*np.log(PD) for SC,PR,PD in zip(test_apical['DIRECT_SOMA'],test_apical['PARENT_RAD'],test_apical['PATH_DIS'])]
#test_basal['predictions'] = [train_coef[0]*(1-SC)*PR + train_coef[1]*SC*PR for SC,PR in zip(test_basal['DIRECT_SOMA'],test_basal['PARENT_RAD'])]
test_basal['predictions'] = [train_coef[0]*SC*PR + train_coef[1]*(1-SC)*BL + train_coef[2]*(1-SC)*PD for SC,PR,BL,PD in zip(test_basal['DIRECT_SOMA'],test_basal['PARENT_RAD'],test_basal['BRANCH_LEN'],test_basal['PATH_DIS'])]
test_basal['predictions'] = [train_coef[0]*BL + train_coef[1]*PD for BL,PD in zip(test_basal['BRANCH_LEN'],test_basal['PATH_DIS'])]

f_stat = '{:.2e}'.format(model.fvalue)
r_val = '{:.4}'.format(model.rsquared_adj)
r_label = ' Adjusted $r^2$ = ' + r_val
additions = 'F-Statistic : ' + f_stat + r_label #' Adjusted R-Squared : ' + r_val
simple_plot([train_basal,test_basal],['RADIUS','predictions'],[1,1],'Initial Purkinje Dendrites', legend = ['Training Set','Testing Set'], add = additions,labels = ['Original Radius','Predicted Radius'])
simple_plot([train_basal,test_basal],['RADIUS','predictions'],[1,0.3], 'Apical Dendrite of Hippocampus: Initial Log Path Distance to Continuing Parent Radius', legend = ['Training','Testing'], add = additions,labels = ['Radius','Predictions'])
save_file(train_coef,'A_Coefs')
'''

#must change plot function to == 2:
#corr,_ = pearsonr(basal_dict['PARENT_RAD'],basal_dict['RADIUS'])
#new_basal['Direct']['log_PD'] = [np.log(PD) for PD in new_basal['Direct']['PATH_DIS']]
#new_apical['Direct']['log_PD'] = [np.log(PD) for PD in new_apical['Direct']['PATH_DIS']]

#corr,_ = pearsonr(new_basal['Direct']['log_PD'],new_basal['Direct']['RADIUS'])
#corr,_ = pearsonr(new_basal['Direct']['PARENT_RAD'],new_basal['Direct']['RADIUS'])
#corr,_ = pearsonr(new_apical['Direct']['PARENT_RAD'],new_apical['Direct']['RADIUS'])
#corr = round(corr**2,4)
#testing = ' R Squared = ' + str(corr)
#simple_plot(basal_dict,['PARENT_RAD','RADIUS'],.3,'Parent Radius to Radius: Hippocampal Proximal Dendrites',legend = testing, where = 'upper right', labels = ['Parent Radius','Radius'])
#simple_plot(new_basal['Direct'],['log_PD','RADIUS'],1,'Log Path Distance to Radius: Hippocampal Initial Proximal Dendrites',legend = testing, where = 'upper left', labels = ['Log Path Distance','Radius'])
#simple_plot(new_basal['Direct'],['log_PD','RADIUS'],1,'Log Path Distance to Radius: Cerebellar Initial Dendrites',legend = testing, where = 'upper left', labels = ['Log Path Distance','Radius'])
