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
#October 9, 2020

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
def main_plot(data, var, title, ax_titles = None, labels = None, save_as = None, plot_type = None, fit_line = None, add = None, where = None, size = None):
    plt.rc('font', size = 28) #default plot and font sizes
    if size:
        fig, ax = plt.subplots(figsize = (20,12))
    else:
        fig, ax = plt.subplots(figsize = (20,10))          
    if add:                                        #additional information to plot
        at = AnchoredText(add,
                          prop=dict(size=28), frameon=True,
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
    for i in dcorr:
        dcorr[i] = dcorr[i]**2                              #test to see if we can match r2 values from plots in correlation matrices
    fig,ax = plt.subplots(figsize = (20,15))
    sns.set(font_scale = 2.8)
    plt.rc('font', size = 26)
    rsquared = '$r^2$'
    #ax = sns.heatmap(dcorr,vmin = -1,vmax = 1, center = 0,cmap = sns.diverging_palette(20, 220, n=256),square = True,annot = True,cbar_kws={'shrink':1,'label':rsquared})#"Pearsons R"}
    ax = sns.heatmap(dcorr,vmin = 0,vmax = 1, center = 0.5,cmap = 'Blues',square = True,annot = True,cbar_kws={'shrink':1,'label':rsquared})#"Pearsons R"})
    ax.set_yticklabels([sub['short'] for sub in selection.values()],rotation = 0,fontsize = 30)
    ax.set_xticklabels([sub['short'] for sub in selection.values()],rotation = 0,horizontalalignment = 'right',fontsize = 30,ha = 'center')
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

def regress_val(model, split = None, r = None): 
    f_stat = '{:.2e}'.format(model.fvalue)
    r_squared = '{:.4}'.format(model.rsquared_adj)
    r_val = '{:.4}'.format(np.sqrt(float(r_squared))) 
    aic = '{:.4}'.format(model.aic)
    bic = '{:.4}'.format(model.bic)
    cond_num = '{:.4}'.format(model.condition_number)
    if split:
        return f_stat,r_val if r == None else f_stat,r_squared
    else:
        vals = [f_stat,r_squared,cond_num,aic,bic]
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

for connect in comp_list:
    sep_data[connect] = {}; comb_data[connect] = {}
    for i in comp_types:
        sep_data[connect][i] = {}; comb_data[connect][i] = {}
        for param in header.keys():
            comb_data[connect][i][param] = []
            if len(complete_dict.keys()) != len(header.keys()):
                complete_dict[param] = []
            for d1 in dirs:
                if connect == comp_list[0]:
                    complete_dict[param].extend(archive_dict[i][d1][param])
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
select_params = {'BRANCH_LEN': {'short':'TD','long':'Total Dendritic Length'},
                 'NODE_DEGREE': {'short':'TB','long':'Terminal Branch Order'},
                 'NODE_ORDER': {'short':'IB','long':'Initial Branch Order'},
                 'PARENT_DIA': {'short':'PD','long':'Parent Diameter'},
                 'PATH_DIS': {'short':'PS','long':'Path To Soma'},
                 'PATH_TO_END': {'short':'LP','long':'Longest Path to Terminal'},
                 'DIAMETER': {'short':'D','long':'Diameter'}}

initial_params = {key:value for (key,value) in select_params.items() if key != 'NODE_ORDER'}

'''Initial Plots and Parameter Analysis''' #save back in folder with fullpath os.cwd
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
        param_comb = []; regress[i][connect] = []; select_data = {}
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
        
        title = 'Correlation Dendrites'
        saving = i + ' ' + connect + ' Feature Correlation'
        corr(select_data, selection, 'Striatal Dendrites', save_as = saving)


#print regression data and sort by r-squared
for i in regress:
    for connect in regress[i]:
        print(i,connect)
        print('Parameter','F-Stat','R-Squared','Condition Number','AIC','BIC')
        regress[i][connect].sort(key=lambda x: x[2])
        for test in regress[i][connect][::-1]:
            print(test)

'''Separate Training and Testing Data'''
seed = 14 
file_copies = list(file_list)             #randomly separate files into training and testing sets
random.Random(seed).shuffle(file_copies)
split_files = split_seq(file_copies,2)
training_files = split_files[0]
testing_files = split_files[1]

'''Choose which Features to Estimate Radius for Equations'''  #estimating radius as found in .swc morphologies
test_params = {}      #'Apical' and 'Basal' are classification of compartment types within .swc morphologies
test_params['Apical'] = {'Initial':['PARENT_RAD','PATH_TO_END'],'BP_Child':['PARENT_RAD','PATH_TO_END'],'Continuing':['PARENT_RAD']}
test_params['Basal'] = {'Initial':['PARENT_RAD','BRANCH_LEN'],'BP_Child':['PARENT_RAD','PATH_DIS'],'Continuing':['PARENT_RAD']}
#test_params['Basal'] = {'Initial':['PARENT_RAD','NODE_ORDER'],'BP_Child':['PARENT_RAD','NODE_DEGREE'],'Continuing':['PARENT_RAD']}
#test_params['Basal'] = {'Initial':['PARENT_RAD','PATH_TO_END'],'BP_Child':['PARENT_RAD','PATH_DIS'],'Continuing':['PARENT_RAD']}

'''Create Equations to Predict Radius'''
pdict = {0:'Initial',1:'Continuing',2:'BP_Child'}   #original data (#'s) to dictionary organization (terms)
cdict = {'3':'Basal','4':'Apical'}
newrads = {fname:{} for fname in file_list}
models = {}; train_data = {}; train_anal = {}; test_anal = {}; train_plot = {}; test_plot = {}#; merge_data = {}

for i in comp_types: 
    models[i] = {'Initial':[],'Continuing':[],'BP_Child':[]} #saves each model for initial, continuing, and branch point children
    train_anal[i] = {'pred':{fname:[] for fname in training_files}, 'org':{fname:[] for fname in training_files}}
    test_anal[i] = {'pred':{fname:[] for fname in testing_files}, 'org':{fname:[] for fname in testing_files}}
    train_plot[i] = {'pred':[],'org':[]}; test_plot[i] = {'pred':[],'org':[]}
    for connect in pdict.values():  
        train_data = {p:[] for p in header.keys()}#,'Test':{p:[] for p in header.keys()}} #will temporarily hold data values in model formation
        for num,fname in enumerate(comb_data[connect][i]['FNAME']): 
            for param in header.keys():
                if fname in training_files:       #split data into training and testing sets to validate equations
                    train_data[param].append(comb_data[connect][i][param][num])
        
        model,predictions = ols_fit(train_data,test_params[i][connect],'RADIUS',extended = True)
        models[i][connect] = {feature:model.params[num] for num,feature in enumerate(test_params[i][connect])}

cd = complete_dict
for num,(comp,parent,fname,ctype,pcon) in enumerate(zip(cd['CHILD'],cd['PARENT'],cd['FNAME'],cd['TYPE'],cd['PAR_CONNECT'])):
    ctype = cdict[ctype]; pcon = pdict[pcon]
    newrad = 0
    for feature in test_params[ctype][pcon]:        
        '''new predictions with updating parent radius'''
        if feature == 'PARENT_RAD':                       #start with initial comps, with saved parent radius as soma val
            if pcon == 'Initial':                   
                newrad = newrad + models[ctype][pcon][feature] * complete_dict['PARENT_RAD'][num]
            else:                                         #if parent radius of any other comp, use updated radius
                newrad = newrad + models[ctype][pcon][feature] * newrads[fname][parent]
                #newrad = newrad + newrads[fname][parent]

        else:                                             #if other feature value, find and update radius
            newrad = newrad + models[ctype][pcon][feature] * complete_dict[feature][num]
        '''old predictions with non-updating parent radius'''
        #newrad = newrad + models[ctype][pcon][feature] * complete_dict[feature][num]

    newrads[fname][comp] = newrad   
    if fname in testing_files:  #save new predicted radius and original radius for R analysis and later plots
        #[x.append(y) for x,y in zip([test_anal[ctype]['pred'][fname],test_plot[ctype]['pred']])]
        #[x.append(complete_dict['RADIUS'][num]) for x in zip([test_anal[ctype]['org'][fname],test_plot[ctype]['org']])]
        test_anal[ctype]['pred'][fname].append(newrad)
        test_anal[ctype]['org'][fname].append(complete_dict['RADIUS'][num])
        test_plot[ctype]['pred'].append(newrad)
        test_plot[ctype]['org'].append(complete_dict['RADIUS'][num])
    elif fname in training_files:
        #[x.append(newrad) for x in zip([train_anal[ctype]['pred'][fname],train_plot[ctype]['pred']])]
        #[x.append(complete_dict['RADIUS'][num]) for x in zip([train_anal[ctype]['org'][fname],train_plot[ctype]['org']])]
        train_anal[ctype]['pred'][fname].append(newrad)
        train_anal[ctype]['org'][fname].append(complete_dict['RADIUS'][num])
        train_plot[ctype]['pred'].append(newrad)
        train_plot[ctype]['org'].append(complete_dict['RADIUS'][num])

for i in comp_types:
    test_r = []; train_r = []
    for fname in file_list:
        if fname in testing_files:          #calculate rsquared values for each file in training and testing sets
            r_val,_ = pearsonr(test_anal[i]['pred'][fname],test_anal[i]['org'][fname])
            test_r.append(r_val**2)
        elif fname in training_files:
            r_val,_ = pearsonr(train_anal[i]['pred'][fname],train_anal[i]['org'][fname])
            train_r.append(r_val**2)
    print(i, ' average',' stdev')
    print('Train',np.average(train_r),np.std(train_r))
    print('Test',np.average(test_r),np.std(test_r))  #calculate averaged rsquared for training, testing, and all files
    train_label = 'Average $R^2$ = ' + '{:.4}'.format(np.average(train_r)) + ' StDev = ' '{:.4}'.format(np.std(train_r))
    test_label = 'Average $R^2$ = ' + '{:.4}'.format(np.average(test_r)) + ' StDev = ' '{:.4}'.format(np.std(test_r))
    ave_r = '{:.4}'.format(np.average([np.average(train_r),np.average(test_r)]))
    title = i + ' Dendrites'
    additions = 'Average $R^2$ = ' + ave_r
    saving = i + ' Predicted Radii'                  #plot predicted vs. original radii
    main_plot([train_plot[i],test_plot[i]], ['org','pred'], title, ax_titles = ['Original Radius','Predicted Radius'], add = additions, labels = ['Training Set: ' + train_label,'Testing Set: ' + test_label], where = 'lower right', save_as=saving) 
