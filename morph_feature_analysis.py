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
#        optional parameter: --seed default=14 - for reproducing prior results on unix OS
#                            --train_test  npy file with training and testing files - for reproducing prior results on any OS

#George Mason University
#Jonathan Reed
#March 14, 2021

#Updated 2021 April 10
#Avrama Blackwell

import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import random
import argparse
import os
import morph_feature_utils as mfu
import morph_feature_plots as mfp

pd.set_option("display.max_columns", 20)
pd.set_option('display.width',240)

#global parameter, matches the swc designation
dend_types={'3':'Basal','4':'Apical'}

'''Organize Data by Parent Connection - Initial, Branch Children, Continuing Compartments'''
comp_list = ['Initial','BP_Child','Continuing','All']
connect_dict = {2:'Continuing',1:'BC', 0:'Initial'}

'''Select parameters to analyze with Compartment Diameter, also give abbreviations and long names'''
#if additional features to plot, add plot labels below
features = {'DIAMETER': {'short':'D','long':'Diameter'},
                 'PARENT_DIA': {'short':'PD','long':'Parent Diameter'},
                 'BRANCH_LEN': {'short':'TL','long':'Total Dendritic Length'},
                 'NODE_DEGREE': {'short':'TD','long':'Terminal Degree'},
                 'NODE_ORDER': {'short':'IB','long':'Initial Branch Order'},
                 'PATH_DIS': {'short':'PS','long':'Path To Soma'},
                 'PATH_TO_END': {'short':'LP','long':'Longest Path to Terminal'} }

###control output  - should be made additional parameters   
df_corr=False
parameter_plots=False
hist_features=[]#['PATH_DIS','PATH_TO_END','BRANCH_LEN','DIAMETER'] #set to empty to turn off these plots
corr_matrices=False
final_predictions=True
xcorr=False
rall_test=False

parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)
parser.add_argument("--seed",type=int, default=14)
parser.add_argument("--train_test",type=str)
args = parser.parse_args()                      
path = args.path

import glob
if not len(glob.glob(path)):  #if the path doesn't exist
    dirs=glob.glob(path+'*')  #use the wildcard character
    root=os.path.dirname(dirs[0])+'/'
    dirs=[os.path.basename(d1) for d1 in dirs]
else:
    '''Locates Archive Data'''
    root, dirs, files = list(os.walk(path))[0]           #all CNG_extract.txt files require identical parameters and length data
if len(dirs) == 0:
    if os.path.basename(path) == '':  #if single repo specified with no subdirectories
        path = root.rstrip('/')                      #fix for trailing '/' in path
    dirs = [os.path.basename(path)]   #then, the single repo is the same as the set of subdirectories
    root = os.path.dirname(path) + '/'  #and the root is the full path of the repo
subtitle=root.rstrip('/')+'_'
df,complete_list,header=mfu.read_morphologies(root,dirs)
nodes={ar:[] for ar in np.unique(df.ARCHIVE)}
for ar in np.unique(df.ARCHIVE):
    ardf=df[df.ARCHIVE==ar]
    for fname in np.unique(ardf.FNAME):
        fdf=ardf[ardf.FNAME==fname]
        nodes[ar].append(len(fdf))
    print(subtitle, ar, 'mean nodes', np.mean(nodes[ar]))

#Some parameters for the graphs
if 'Hippocampus' in path.split('/'):
    ax_max=10
    legend_loc='lower right'
elif 'Cerebellum' in path.split('/'):
    ax_max=8
    legend_loc='upper right'
elif 'Basal_Ganglia' in path.split('/'):
    ax_max=6
    legend_loc='upper right'
else:
    print('new archive')

#mfu.compare_EM(df,path)
################# correlations of various types using dictionaries
if df_corr:
    feat_corr={dtnum:{f:[] for f in hist_features+['PARENT_DIA']} for dtnum in np.unique(df.TYPE)}
    for dtnum in feat_corr.keys():
        ctdf=df[df.TYPE==dtnum] #added on 3/29
        for feat in hist_features+['PARENT_DIA']:
            for ar in np.unique(df['ARCHIVE']):
                feat_corr[dtnum][feat].append(round(ctdf[ctdf.ARCHIVE==ar][feat].corr(ctdf[ctdf.ARCHIVE==ar]['DIAMETER']),4))
        print(dend_types[dtnum],[ar for ar in np.unique(df['ARCHIVE'])])
        for feat in feat_corr[dtnum].keys():
            print(dtnum,feat,', corr by archive:',feat_corr[dtnum][feat])
        for ar in np.unique(df['ARCHIVE']):
             dfsubset=ctdf[(ctdf.ARCHIVE==ar) & (ctdf.PAR_CONNECT==1)]
             print(ar, dend_types[dtnum], 'corr: 3/2',pearsonr(dfsubset['DIAMETER']**(1.5),dfsubset['PARENT_DIA']**(1.5))[0],
                   'linear',pearsonr(dfsubset['DIAMETER'],dfsubset['PARENT_DIA'])[0])
        
    print(df.groupby(['TYPE','ARCHIVE'])[['PARENT_DIA','DIAMETER']].corr())
    print(df.groupby(['TYPE','PAR_CONNECT'])[['PARENT_DIA','DIAMETER']].corr())
    #print(df.groupby(['TYPE','ARCHIVE','PAR_CONNECT'])[['PARENT_DIA','DIAMETER']].corr())
    
    print(' ********** neuron-wise correlations')
    neuron_corr={a:{} for a in np.unique(df['ARCHIVE'])}
    for ar in neuron_corr.keys():
        for f in np.unique(df[df.ARCHIVE==ar]['FNAME']):
            neuron_corr[ar][f]=df[df.FNAME==f]['PARENT_DIA'].corr(df[df.FNAME==f]['DIAMETER'])
        print('     archive',ar,round(np.mean(list(neuron_corr[ar].values())),3))
initial_params = {key:value for (key,value) in features.items() if key != 'NODE_ORDER'}

####################################################3
'''Initial Plots and Parameter Analysis'''
for dtype in np.unique(df.TYPE): #apical vs basal
    for param in features:
        if param != 'DIAMETER':
            compdf=df[(df.TYPE==dtype)]
            plabel = features[param]['long']
            #title = i + ' ' + plabel + ' to Diameter'
            title = dend_types[dtype] + ' Dendrites'
            saving = subtitle+dend_types[dtype] + param+'_vs_Diam'
            
            '''Plot Parameters by Archive'''
            arch_labels={}
            drs_len=np.max([len(d1) for d1 in dirs])
            for d1 in dirs:
                archdf=compdf[(compdf['ARCHIVE']==d1)]
                temp,_ = pearsonr(archdf[param],archdf['DIAMETER'])
                arch_labels[d1]=label = d1.ljust(drs_len) + ' $R^2$ = ' + str(round(temp**2,4) )
            ct_labels={}
            ct_len=np.max([len(ctnm) for ctnm in connect_dict.values()])
            for ctnum in connect_dict.keys():
                ctdf=compdf[(compdf['PAR_CONNECT']==ctnum)]
                temp,_ = pearsonr(ctdf[param],ctdf['DIAMETER'])
                ct_labels[ctnum]=connect_dict[ctnum].ljust(ct_len) + ' $R^2$ = ' + str(round(temp**2,4))                
            if parameter_plots:
                print('fig title',saving)
                #Plot Parameters to Diameter
                mfp.df_plot(compdf, [param,'DIAMETER'], title=title, ax_titles = [plabel,'Diameter'], save_as = saving)
                
                #Plot Parameters by Parent Connection
                mfp.df_plot(compdf, [param,'DIAMETER'], labels=ct_labels, select_param='PAR_CONNECT',title=title, ax_titles = [plabel+' ($\mu$m)','Diameter'+' ($\mu$m)'], where = 'upper right', save_as = saving+'by_connect',ax_max=ax_max)
                #plot parameters by archive
                mfp.df_plot(compdf, [param,'DIAMETER'],labels=arch_labels, select_param='ARCHIVE',title=title, ax_titles = [plabel+' ($\mu$m)','Diameter'+' ($\mu$m)'],where='upper right', save_as = saving+'by_archive',ax_max=ax_max)

####################################################3
#plot histograms and test whether Rall's 3/2 power rule holds
if len(hist_features):
    histograms={};binset={}
    for feat in hist_features:
        if feat=='BRANCH_LEN':
            log=True
        else:
            log=False
        histograms[feat],binset[feat]=mfu.archive_histogram(df,dirs,dend_types,feat,features[feat]['long'],log=log)
    for dtype in histograms[feat]:
        saving = 'Histogram_'+subtitle+dend_types[dtype]
        mfp.plot_hist_set (histograms,binset,dtype,features,title=None,save_as=saving)

if rall_test:
    saving = 'Rall_test_'+subtitle#+dend_types[ctnum]
    parent_dia,sum_child_dia=mfp.Rall_plot(df,dirs,comp_list,save_as=saving)
    child_dia32={ct:{ar:[] for ar in dirs} for ct in sum_child_dia.keys()}
    parent_dia32={ct:{ar:[] for ar in dirs} for ct in parent_dia.keys()}
    for ctnum in child_dia32.keys():
        for ar in child_dia32[ctnum].keys():
            child_dia32[ctnum][ar]=[cd**1.5 for clist in sum_child_dia[ctnum][ar] for cd in clist]
            parent_dia32[ctnum][ar]=[cd**1.5 for clist in parent_dia[ctnum][ar] for cd in clist]
            print(ctnum,ar,round((pearsonr(child_dia32[ctnum][ar],parent_dia32[ctnum][ar])[0])**2,3))
####################################################3
''' multiple linear regression 
One method to improve this: loop through par2 in order of adjR - i.e., sort the dictionary '''
improvement=0.02 #this is arbitrary.  Use one parameter model unless 2 parameter model increases adjusted R2 by this amount
critical_CN =20 #condition numbers greater than this are worrisome
regress = {dend_types[ctnum]:{ct:[] for ct in comp_list[0:-1]} for ctnum in np.unique(df.TYPE)}
adjr={}
for ctnum in np.unique(df.TYPE):
    ct=dend_types[ctnum]
    #Plot Correlation between 2 Features
    for cn,connect in enumerate(comp_list[0:-1]):     #setup correlation data for each compartment type
        comb_df=df[(df.PAR_CONNECT==cn) & (df.TYPE==ctnum)]
        adjR={}
        if connect == 'Initial':
            indep_var={s:v for s,v in features.items() if s != 'DIAMETER' and s != 'NODE_ORDER'}
        else:
            indep_var={s:v for s,v in features.items() if s != 'DIAMETER'}
        for pnum,par1 in enumerate(indep_var.keys()):
            '''Multiple Regression for Diameter, runs statsmodels OLS (multiple) regression for parameters'''
            model,predictions = mfu.ols_fit(comb_df,[par1],'DIAMETER',connect=connect,extended=False,add_const=True) #be aware that single variable will still fit model
            vals = mfu.regress_val(model)
            regress[ct][connect].append([par1]+vals)
            adjR[par1]=model.rsquared_adj
        adjR=dict(sorted(adjR.items(),key=lambda item: item[1],reverse=True))
        print('FINISHED 1 par loop for', connect, adjR)
        par1='PARENT_DIA'; adjR1=adjR[par1]; bestR2=adjR1
        del adjR['PARENT_DIA']
        for par2,adjR2 in adjR.items():
            #regression for pairs of parameters - use 2nd parameter if pearonsr to diameter is higher than pearsonsR to parent_dia
            model,predictions = mfu.ols_fit(comb_df,[par1,par2],'DIAMETER',connect=connect,add_const=True)
            vals = mfu.regress_val(model)
            if model.rsquared_adj-max(adjR1,adjR2)>improvement:
                regress[ct][connect].append([str(par1 + ' + ' + par2)]+vals)
                bestR2=model.rsquared_adj
            else:
                print('PARAMETER PAIRS: not adding',connect,par1,'=',round(adjR1,4), 
                      par2,'=',round(adjR2,4),'vs combined:',round(model.rsquared_adj,4),
                      'condition',round(model.condition_number))
        # Additional loop to see if parameter combinatinos not using parent diameter are better
        for kk,(par1,adjR1) in enumerate(adjR.items()):
            for par2,adjR2 in list(adjR.items())[kk+1:]:
                model,predictions = mfu.ols_fit(comb_df,[par1,par2],'DIAMETER',connect=connect,add_const=True)
                vals = mfu.regress_val(model)
                if model.rsquared_adj-max(adjR1,adjR2)>improvement and model.rsquared_adj> bestR2:
                    regress[ct][connect].append([str(par1 + ' + ' + par2)]+vals)
                    bestR2=model.rsquared_adj
                else:
                    print(connect,'Not adding addtl PARAMETER PAIRS:', par1, 'R2=',round(adjR1,3),
                          par2,'R2=',round(adjR2,3), 'modelR2',round(model.rsquared_adj,4))
        print('FINISHED 2 param loop for', connect) 
        if corr_matrices:                
            title = ''#ct + ' Dendrites'
            saving = 'Corr_Matrix'+subtitle+dend_types[ctnum] +connect
            #saving = ct + ' ' + connect + ' Feature Correlation'
            mfp.corr(comb_df, indep_var, title, save_as = saving)

########### plot correlation by filename

fname_set={ar:np.unique(df[df.ARCHIVE==ar].FNAME) for ar in np.unique(df.ARCHIVE)}
if xcorr:
    saving = 'Cross_corr'+subtitle
    mean_xcorr,lags,decay=mfu.cross_corr(df,fname_set,dend_types)
    mfp.plot_xcorr(mean_xcorr,lags,decay,xmax=600,save_as=saving)
####################################################3
#print regression data and sort by r-squared
best_reg_params={i:{} for i in dend_types.values()}
for i in regress:
    for connect in regress[i]:
        if len(regress[i][connect]):
            print('FINAL RESULTS OF MULTI-REGRESION FOR', i,connect)
            print('Parameter','F-Stat','R-Squared','Condition Number','AIC','BIC','R2 of model with Intercept')
            regress[i][connect].sort(key=lambda x: x[2])
            for test in regress[i][connect][::-1]:
                print(test)
                best_reg_params[i][connect]=sorted(regress[i][connect], key=lambda x: x[2],reverse=True)[0][0].split(' + ')
        else:
            print('NOTHING SIGNIFICANT FOR', i,connect, 'using PARENT_DIA')
            best_reg_params[i][connect]=['PARENT_DIA']
      
####################################################3
'''Separate Training and Testing Data'''

#to reproduce model fit results, need the original split_files
#This loads the original training and testing files from the specified .npy file
if args.train_test:
    split_files=np.load(args.train_test, allow_pickle=True)
    training_files = split_files[0]
    testing_files = split_files[1]
else:
    file_copies = list(np.unique(df['FNAME']))             #randomly separate files into training and testing sets
    random.Random(args.seed).shuffle(file_copies)
    split_files = mfu.split_seq(file_copies,2)
    training_files = split_files[0]
    testing_files = split_files[1]

print('train files',len(training_files))
print('test files',len(testing_files))

######### neuron wise correlations and cross-correlations   ###########
fname_set={'train':training_files,'test':testing_files}
if xcorr:
   train_test_xcorr,train_test_lags,train_test_decay=mfu.cross_corr(df,fname_set,dend_types)
   mfp.plot_xcorr(train_test_xcorr,train_test_lags,train_test_decay,xmax=None,save_as='')

test_params=best_reg_params

add_const=False
'''Create Equations (OLS Model) and Save Models to File
Ideally this should go into a function, input as two sets of filenames: train and test
Then, could call this in a loop over archive, or do all archives'''
pdict = {0:'Initial',1:'BP_Child',2:'Continuing'}   #original data (#'s) to dictionary organization (terms)
newrads = {fname:{} for fname in np.unique(df['FNAME'])}; newrads_i = {fname:{} for fname in np.unique(df['FNAME'])}

model_file = open(subtitle+'model.txt','w')
models = {}
for ctnum in np.unique(df.TYPE): #apical vs basal
    i=dend_types[ctnum]
    models[i] = {tp: [] for tp in comp_list[0:-1]} #saves model for initial, continuing, and branch point children 
    #separate training and testing files in dictionary by file name
    #testing set split into test set including original diameters (Test_i)
    for cn,connect in pdict.items():  
        comb_df=df[(df.PAR_CONNECT==cn) & (df.TYPE==ctnum)]
        fname_dfs=[]
        for fname in training_files:
            fname_dfs.append(comb_df[comb_df.FNAME==fname])
        train_df=pd.concat(fname_dfs, ignore_index=True)
        model,_ = mfu.ols_fit(train_df,test_params[i][connect],'DIAMETER',connect=connect,add_const=add_const)#,extended = True)
        if add_const:
            test_params[i][connect]=test_params[i][connect]+['const']
        print(' *****for Table ',i, connect, model.params, model.pvalues, 'R-Squared = ', model.rsquared_adj)
        models[i][connect] = {feature:model.params[num] for num,feature in enumerate(test_params[i][connect])}

        model_file.write(i + '\n')
        model_file.write(connect + '\n')
        for feature in models[i][connect]:
            line = ' '.join([feature, str(models[i][connect][feature])])
            model_file.write(line + '\n')

model_file.close()

'''Calculate New Predicted Values from Equations
Calculations done recursively, starting from the initial segment
Thus, only the radius of the soma is used, or
radius of initial segment for that set of predictions
If this is part of shape_shifter, then shouldn't be done here.
Ideally, create a function to use model_file, OR 
if that function already exists in shape_shifter, import the function here
possibly loop over that function using the different archive models 
Then, could call this in a loop over archive, or do all archives''' #includes updating Parent Diameters

######## set up dictionaries to hold predictions
if final_predictions:
    model_data = {}
    for ctnum in np.unique(df.TYPE): #apical vs basal
        i=dend_types[ctnum]
        model_data[i] = {ts:{'fname':{},'node_type':{nd:{'pred':[],'org':[]} for nd in pdict.values()}} for ts in ['Train','Test','Test_i']}
        for ts in ['Train','Test','Test_i']:
            model_data[i][ts]['fname'] = {fn:{'pred':[],'org':[]} for fn in training_files} if ts == 'Train' else {fn:{'pred':[],'org':[]} for fn in testing_files}
    
    for line in complete_list:
        #if line[header['FNAME']] == 'WT-1201MSN03.CNG_extract':
        comp = line[header['CHILD']]; parent = line[header['PARENT']]; fname = line[header['FNAME']]
        ctype = dend_types[line[header['TYPE']]]; pcon = pdict[line[header['PAR_CONNECT']]]
        newrad = 0; newrad_i = 0
        for feature in test_params[ctype][pcon]:        
            '''new predictions with updating parent radius'''
            if feature == 'PARENT_DIA':                       #start with initial comps, with saved parent radius as soma val
                if pcon == 'Initial':                   
                    newrad = newrad + models[ctype][pcon][feature] * line[header['PARENT_DIA']]
                    newrad_i = line[header['DIAMETER']]
                else:                                         #if parent radius of any other comp, use updated radius
                    newrad = newrad + models[ctype][pcon][feature] * newrads[fname][parent]
                    newrad_i = newrad_i + models[ctype][pcon][feature] * newrads_i[fname][parent]
            elif feature=='const':
                 newrad = newrad+models[ctype][pcon][feature]
                 if pcon != 'Initial':
                     newrad_i = newrad+models[ctype][pcon][feature]
            else:                                             #if other feature value, find and update radius
                newrad = newrad + models[ctype][pcon][feature] * line[header[feature]]
                if pcon=='Initial':
                    newrad_i=line[header['DIAMETER']]
                if parent != '1': #only applies additional features to predicted initial radii
                    newrad_i = newrad_i + models[ctype][pcon][feature] * line[header[feature]]
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
    for i in model_data.keys():
        corr_data[i] = {ts:{'fname':{},'node_type':{nd:[] for nd in pdict.values()}} for ts in ['Train','Test','Test_i']}
        
        for tt in ['Train','Test','Test_i']: #calculate correlation within node type across files
            corr_data[i][tt]['fname'] = {fname:[] for fname in training_files} if tt == 'Train' else {fname:[] for fname in testing_files}
            for nt in pdict.values():
                r_val,_ = pearsonr(model_data[i][tt]['node_type'][nt]['pred'],model_data[i][tt]['node_type'][nt]['org'])
                corr_data[i][tt]['node_type'][nt].append(r_val**2)
    
        for fname in np.unique(df['FNAME']):  #calculate correlation within individual files 
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
            ave = np.average(list(mfu.flatten(corr_data[i][pt]['fname'].values()))); std = np.std(list(mfu.flatten(corr_data[i][pt]['fname'].values())))
            t_labels[pt] = '$R^2$ = ' + str(round(ave,2)) + ' +/- ' + str(round(std,2))
            
            print(pt,'By File',ave,std)
            print(corr_data[i][pt]['node_type'])
        title = ''#i + ' Dendrites'
        #additions = 'Average $R^2$ = ' + ave_r
        saving = 'Predict'+subtitle+i
        mfp.main_plot([to_plot['Train'],to_plot['Test'],to_plot['Test_i']], ['org','pred'], title, 
                  ax_titles = ['Original Diameter ($\mu$m)', 'Predicted Diameter ($\mu$m)'], 
                  labels = ['Train         ' + t_labels['Train'],'Test          ' + t_labels['Test'],'Test + In. ' + t_labels['Test_i']],
                  where = legend_loc, save_as=saving,ax_max=ax_max)# , add = additions)
