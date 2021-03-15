'''<Extracts data values from .swc files>
    Copyright (C) <2021>  <Jonathan Reed>

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
#        uses visual and statistical tools to analyze features predictive of node diameter
#        will output feature plots as .png and model equations as model.txt
#        where path is name of path to folder containing archive folder(s) with original .swc morphologies and .CNG_extract.txt files
#        .CNG_extract.txt files created from morph_feature_extract.py

#George Mason University
#Jonathan Reed
#March 14, 2021

import numpy as np
from scipy import optimize
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import statsmodels.api as sm
import random
import seaborn as sns
import argparse
import os
import glob
import re

'''Save plot as png'''            
def save_png(png, title):                       #default to save plots instead of plot to window
    png.savefig(title, bbox_inches = 'tight')   #remove empty space around plots
    print('File Created : ' + str(title))
    png.close()
    
'''Main Plot Function for Basic, Fit or 3d Color Plots'''
def main_plot(data, var, title = None, ax_titles = None, labels = None, save_as = None, plot_type = None, fit_line = None, add = None, where = None, marker_size = None):
    #plt.rc('font', size = 38)                       
    plt.rc('font', size = 34)                    #default plot and font sizes
    fig, ax = plt.subplots(figsize = (20,10))
    #fig, ax = plt.subplots(figsize = (10,10))
    if add:                                      #additional information to plot
        at = AnchoredText(add, prop=dict(size=28), frameon=True, loc='upper center')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") 
        ax.add_artist(at)
       
    if plot_type == 'Color3d':                    #plot 3 variables(x,y,z) with variable 'z' using colorbar
        scat = ax.scatter(data[var[0]], data[var[1]], c=data[var[2]], s=100, marker='o', label = labels)
        fig.colorbar(scat)
        
    elif plot_type == 'Fit':                      #plot initial fitted equation to data values
        plt.plot(data[var[0]], data[var[1]], 'o')
        plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) 
    
    else:
        for num in range(len(data)):
            if labels:                            #can modify markersize or alpha (transparency) of plotted points
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', label = labels[num], markersize = 15)
            else:
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', markersize = 15)

    if labels or fit_line:
        plt.legend(loc = where) if where else plt.legend()
    plt.xlabel(ax_titles[0]) if ax_titles else plt.xlabel(var[0])
    plt.ylabel(ax_titles[1]) if ax_titles else plt.ylabel(var[1])
    if title:
        plt.title(title)                           #will require 'title' if nothing passed in 'save_as'
    save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')

'''Plot Parameter Correlation to Diameter with Pearson's R'''
def corr(data,selection,title,save_as = None):
    dframe = pd.DataFrame(data,columns = selection.keys())  #dataframe reorders alphabetically feature columns by default 
    dcorr = dframe.corr()                                   #set feature order by select_params --> passed as 'selection'
    for i in dcorr:
        dcorr[i] = dcorr[i]**2                              #dataframe correlation as R2
    plt.rc('font', size = 30)            
    fig,ax = plt.subplots(figsize = (20,15))
    #sns.set(font_scale = 4)
    rsquared = '$R^2$'                                      #can modify heatmap for either (+/-) R or (+) R2 for Feature Correlations
    #ax = sns.heatmap(dcorr,vmin = -1,vmax = 1, center = 0,cmap = sns.diverging_palette(20, 220, n=256),square = True,annot = True,cbar_kws={'shrink':1,'label':rsquared})#"Pearsons R"}
    ax = sns.heatmap(dcorr,vmin = 0,vmax = 1, center = 0.5,cmap = 'Blues',square = True,annot = True,cbar_kws={'label':rsquared})
    ax.set_yticklabels([sub['short'] for sub in selection.values()],rotation = 0)
    ax.set_xticklabels([sub['short'] for sub in selection.values()],rotation = 0,horizontalalignment = 'right',ha = 'center')
    cbar = ax.collections[0].colorbar                        #selection will contain short and long versions of feature names
    plt.title(title, loc = 'left')   
    save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')
    
'''Function Definitions for Fit Equations'''
#fit selected variables to equations; can input more functions for other feature relationships (to diameter)
def func0(x,m,b):                                     
    return m*x+b

'''Fit Data to Selected Equation and Obtain Residuals'''
#ALTERNATIVE Method to plot predicted values and find residuals
def initial_fit(data, var, func):          #data as dictionary with feature data and var is list of selected features
    x = var[0]                             #assumes only fit to single x- and y- variable
    y = var[1]                             #will require modification if multiple 'x' or transformed 'x' features
    popt, pcov = optimize.curve_fit(func[0],data[x],data[y])    #fits equation to selected function if possible
    print('This is fit estimate for : ', func[1])
    print('popt',popt)
    print('pcov',pcov)
    print('Where x = ', x, ' and y = ', y) 
 
    '''Plot Data to Function and Find Residuals'''
    for i in data[x]:
        temp_max = round(max(data[i])) 
        temp_min = round(min(data[i]))
        x_range = np.arange(temp_min, temp_max + 1, 1)              #plot original feature values to line equation
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
#PREFERRED Method to predict feature variable values
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

'''Equation Values to compare OLS models'''
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
                        for point, val in enumerate(temp_line):                  #Child, Type, Parent are DESCRETE values defined by .swc morphology
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
        if line[header['TYPE']] == '4':                     #currently does not keep axon data if present in .swc morphology
            apical_list.append(line)                        
        elif line[header['TYPE']] == '3':
            basal_list.append(line)                        

    if not 'Basal' in archive_dict.keys():                  #assumes basal compartment 'designation' in all .swc iles
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
                    complete_dict[param].extend(archive_dict[i][d1][param])       #no longer keeps order as files were read
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
#if additional features to plot, add plot labels below
select_params = {'BRANCH_LEN': {'short':'TD','long':'Total Dendritic Length'},
                 'NODE_DEGREE': {'short':'TB','long':'Terminal Branch Order'},
                 'NODE_ORDER': {'short':'IB','long':'Initial Branch Order'},
                 'PARENT_DIA': {'short':'PD','long':'Parent Diameter'},
                 'PATH_DIS': {'short':'PS','long':'Path To Soma'},
                 'PATH_TO_END': {'short':'LP','long':'Longest Path to Terminal'},
                 'DIAMETER': {'short':'D','long':'Diameter'}}

initial_params = {key:value for (key,value) in select_params.items() if key != 'NODE_ORDER'}

'''Initial Plots and Parameter Analysis''' 
for i in comp_types:
    for param in select_params:
        if param != 'DIAMETER':
            d1_data = []; label_list = []
            plabel = select_params[param]['long']
            #title = i + ' ' + plabel + ' to Diameter'
            title = i + ' Dendrites'
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
                label = d1 + ' $R^2$ = ' + str(rsquared)
                label_list.append(label)
            main_plot(d1_data, [param,'DIAMETER'], title, ax_titles = [plabel,'Diameter'], labels = label_list, where = 'upper right', save_as = saving + ' A')

'''Plot Feature Correlations'''
regress = {}                      #setup regression dictionary
for i in comp_types:
    regress[i] = {}
    
    #Plot Correlation between 2 Features
    for connect in comp_list:     #setup correlation data for each compartment type
        param_comb = []; regress[i][connect] = []; select_data = {}
        selection = initial_params if connect == 'Initial' else select_params
        for p1 in selection:
            
            '''Multiple Regression for Diameter'''
            if p1 != 'DIAMETER':  #runs statsmodels OLS (multiple) regression for single and paired parameters
                model,predictions = ols_fit(comb_data[connect][i],p1,'DIAMETER') #be aware that single variable will still fit model
                vals = regress_val(model)                                        
                temp = [p1]
                temp.extend(vals)
                regress[i][connect].append(temp)
                for p2 in selection:          #checks for duplicates to print only unique feature (2) combinations
                    if p2 != 'DIAMETER' and p2 != p1 and str(p2 + ' + ' + p1) not in param_comb:
                        param_comb.append(str(p1 + ' + ' + p2))
                        model,predictions = ols_fit(comb_data[connect][i],[p1,p2],'DIAMETER')
                        vals = regress_val(model)
                        temp = [str(p1 + ' + ' + p2)]
                        temp.extend(vals)
                        regress[i][connect].append(temp)
                 
            select_data[p1] = comb_data[connect][i][p1]

        title = i + ' Dendrites'
        saving = i + ' ' + connect + ' Feature Correlation'
        corr(select_data, selection, title, save_as = saving)

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

print('train files',len(training_files))
print('test files',len(testing_files))

'''Choose which Features to Estimate Radius for Equations'''  #estimating radius as found in .swc morphologies
test_params = {}      #'Apical' and 'Basal' are classification of compartment types within .swc morphologies
test_params['Apical'] = {'Initial':['PARENT_DIA','PATH_TO_END'],'BP_Child':['PARENT_DIA','PATH_TO_END'],'Continuing':['PARENT_DIA']}
#test_params['Basal'] = {'Initial':['PARENT_DIA','BRANCH_LEN'],'BP_Child':['PARENT_DIA','PATH_DIS'],'Continuing':['PARENT_DIA']}
test_params['Basal'] = {'Initial':['PARENT_DIA','PATH_TO_END'],'BP_Child':['PARENT_DIA','PATH_DIS'],'Continuing':['PARENT_DIA','NODE_ORDER']}
#test_params['Basal'] = {'Initial':['PARENT_DIA','NODE_DEGREE'],'BP_Child':['PARENT_DIA','NODE_ORDER'],'Continuing':['PARENT_DIA']}

'''Create Equations (OLS Model) and Save Models to File'''
pdict = {0:'Initial',1:'BP_Child',2:'Continuing'}   #original data (#'s) to dictionary organization (terms)
cdict = {'3':'Basal','4':'Apical'}
newrads = {fname:{} for fname in file_list}; newrads_i = {fname:{} for fname in file_list}

model_file = open('model.txt','w')
models = {}; model_data = {}

for i in comp_types:
    models[i] = {'Initial':[],'Continuing':[],'BP_Child':[]} #saves model for initial, continuing, and branch point children
    model_data[i] = {ts:{'fname':{},'node_type':{nd:{'pred':[],'org':[]} for nd in pdict.values()}} for ts in ['Train','Test','Test_i']}

    #separate training and testing files in dictionary by file name
    #testing set split into test set including original diameters (Test_i)
    for ts in ['Train','Test','Test_i']:
        model_data[i][ts]['fname'] = {fn:{'pred':[],'org':[]} for fn in training_files} if ts == 'Train' else {fn:{'pred':[],'org':[]} for fn in testing_files}
        
    for connect in pdict.values():  
        train_data = {p:[] for p in header.keys()} #temporarily hold feature values
        for num,fname in enumerate(comb_data[connect][i]['FNAME']): 
            for param in header.keys():
                if fname in training_files:        #split data into training and testing sets to validate equations
                    train_data[param].append(comb_data[connect][i][param][num])
        
        model,_ = ols_fit(train_data,test_params[i][connect],'DIAMETER')#,extended = True)
        print(i, connect, model.params, 'R-Squared = ', model.rsquared_adj)
        models[i][connect] = {feature:model.params[num] for num,feature in enumerate(test_params[i][connect])}

        model_file.write(i + '\n')
        model_file.write(connect + '\n')
        for feature in models[i][connect]:
            line = ' '.join([feature, str(models[i][connect][feature])])
            model_file.write(line + '\n')

model_file.close()

'''Calculate New Predicted Values from Equations''' #includes updating Parent Diameters
for line in complete_list:
    #if line[header['FNAME']] == 'WT-1201MSN03.CNG_extract':
    comp = line[header['CHILD']]; parent = line[header['PARENT']]; fname = line[header['FNAME']]
    ctype = cdict[line[header['TYPE']]]; pcon = pdict[line[header['PAR_CONNECT']]]
    newrad = 0; newrad_i = 0
    for feature in test_params[ctype][pcon]:        
        '''new predictions with updating parent radius'''
        if feature == 'PARENT_DIA':                       #start with initial comps, with saved parent radius as soma val
            if pcon == 'Initial':                   
                newrad = newrad + models[ctype][pcon][feature] * line[header['PARENT_DIA']]
                newrad_i = newrad_i + line[header['DIAMETER']]
            else:                                         #if parent radius of any other comp, use updated radius
                newrad = newrad + models[ctype][pcon][feature] * newrads[fname][parent]
                newrad_i = newrad_i + models[ctype][pcon][feature] * newrads_i[fname][parent]
                #newrad = newrad + newrads[fname][parent]

        else:                                             #if other feature value, find and update radius
            newrad = newrad + models[ctype][pcon][feature] * line[header[feature]]
            if parent != '1': #only applies additional features to predicted initial radii
                newrad_i = newrad_i + models[ctype][pcon][feature] * line[header[feature]]
        '''old predictions with non-updating parent radius'''
        #newrad = newrad + models[ctype][pcon][feature] * complete_dict[feature][num]

    newrads[fname][comp] = newrad
    newrads_i[fname][comp] = newrad_i

    '''Save Predicted Values by File and by Node Type'''
    tlist = []
    tlist.append('Train') if fname in training_files else tlist.extend(['Test','Test_i'])

    #iterate through Train, Test, and Test_i dictionaries and save by file name and node type
    for ts in tlist:
        model_data[ctype][ts]['fname'][fname]['org'].append(line[header['DIAMETER']])   #save original width 
        model_data[ctype][ts]['node_type'][pcon]['org'].append(line[header['DIAMETER']])
        if ts != 'Test_i':
            model_data[ctype][ts]['fname'][fname]['pred'].append(newrad)                #save new width w/o original initial
            model_data[ctype][ts]['node_type'][pcon]['pred'].append(newrad)
        else:
            model_data[ctype][ts]['fname'][fname]['pred'].append(newrad_i)              #save new width w/ original initial
            model_data[ctype][ts]['node_type'][pcon]['pred'].append(newrad_i)

'''Calculate Correlation of Predicted to Original Diameter by File and by Node Type'''
corr_data = {} #correlation data from values found in model_data
for i in comp_types:
    corr_data[i] = {ts:{'fname':{},'node_type':{nd:[] for nd in pdict.values()}} for ts in ['Train','Test','Test_i']}
    
    for tt in ['Train','Test','Test_i']: #calculate correlation within node type across files
        corr_data[i][tt]['fname'] = {fname:[] for fname in training_files} if tt == 'Train' else {fname:[] for fname in testing_files}
        for nt in pdict.values():
            r_val,_ = pearsonr(model_data[i][tt]['node_type'][nt]['pred'],model_data[i][tt]['node_type'][nt]['org'])
            corr_data[i][tt]['node_type'][nt].append(r_val**2)

    for fname in file_list:  #calculate correlation within individual files 
        tlist = []
        tlist.append('Train') if fname in training_files else tlist.extend(['Test','Test_i'])
        for ts in tlist:
            r_val,_ = pearsonr(model_data[i][ts]['fname'][fname]['pred'],model_data[i][ts]['fname'][fname]['org'])
            corr_data[i][ts]['fname'][fname].append(r_val**2)

    '''Plot Predicted to Original Diameter Values and Correlations'''
    print(i, ' average',' stdev')
    t_labels = {}; to_plot = {}
    for pt in ['Train','Test','Test_i']:
        to_plot[pt] = {'org':[],'pred':[]}
        for fname in model_data[i][pt]['fname']:             #send single averaged R2 (and stdev) across all files to plot
            to_plot[pt]['org'].extend(model_data[i][pt]['fname'][fname]['org'])
            to_plot[pt]['pred'].extend(model_data[i][pt]['fname'][fname]['pred'])
        ave = np.average(list(flatten(corr_data[i][pt]['fname'].values()))); std = np.std(list(flatten(corr_data[i][pt]['fname'].values())))
        t_labels[pt] = '$R^2$ = ' + str(round(ave,2)) + ' +/- ' + str(round(std,2))
        
        print(pt,'By File',ave,std)
        print(corr_data[i][pt]['node_type'])
    title = i + ' Dendrites'
    #additions = 'Average $R^2$ = ' + ave_r
    saving = i + ' Predicted Diameter'
    main_plot([to_plot['Train'],to_plot['Test'],to_plot['Test_i']], ['org','pred'], title, ax_titles = ['Original Diameter ($\mu$m)','Predicted Diameter ($\mu$m)'], labels = ['Train         ' + t_labels['Train'],'Test          ' + t_labels['Test'],'Test + In. ' + t_labels['Test_i']], where = 'lower right', save_as=saving)# , add = additions)
