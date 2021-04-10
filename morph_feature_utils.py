'''<utilities for analyzing morphology features extracted from .swc files>
    Copyright (C) <2021>  <Jonathan Reed>,<Avrama Blackwell>

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

import numpy as np
from scipy import optimize
import pandas as pd
import os
import glob
import re
import statsmodels.api as sm

'''Function Definitions for Fit Equations'''
#fit selected variables to equations; can input more functions for other feature relationships (to diameter)
def func0(x,m,b):                                     
    return m*x+b

def exp1(t,tau1,A,B):
    return A*np.exp(-t/tau1)+B 

def exp2(t,tau1,tau2,A,B,C):
    return A*np.exp(-t/tau1)+B*np.exp(-t/tau2)+C 

def archive_histogram(df,dirs,dend_types,feature,xaxis_title,log=False):
    hist={}
    numsteps=50
    for dtype in np.unique(df.TYPE): #apical vs basal
        hist[dtype]={}
        minval=np.percentile(df[feature],1)
        maxval=np.percentile(df[feature],99)
        if log and minval>0:
            logbins=np.arange(np.log10(minval),np.log10(maxval),(np.log10(maxval)-np.log10(minval))/numsteps)
            bins=10**logbins
        else:
            bins=np.arange(minval,maxval,(maxval-minval)/numsteps)
        for jj,ar in enumerate(dirs):
            dia=df[(df['ARCHIVE']==ar) & (df['TYPE']==dtype)][feature] #edited on 3/29 to select by dtype
            hist[dtype][ar],_bins=np.histogram(dia,bins=bins,range=(minval,maxval))
    return hist,bins
    
'''OLS Regression to fit Variables'''
def ols_fit(data,xlabel,ylabel,extended = False,connect='',add_const=False):
    temp_df = pd.DataFrame(data)   #xlabel can be multiple names of parameters to fit towards Diameter
    X = temp_df[xlabel]
    Y = temp_df[ylabel]
    model = sm.OLS(Y,X).fit()
    if connect=='BP_Child' and 'PARENT_DIA' in xlabel: #1. FIXME, 2. compare automatic with Reed models
        diam_index=xlabel.index('PARENT_DIA')
        X1=temp_df[xlabel[diam_index]]**(3/2)
        if len(xlabel)>1:
            X2=temp_df[xlabel[len(xlabel)-diam_index-1]]
            X32=pd.concat([X1,X2],axis=1)
        else:
            X32=X1
        Y32=Y**(3/2)
        model32=sm.OLS(Y32,X32).fit()
        print('$$$$$$$$$$$$$$$ using 3/2 power for',connect,xlabel,
                  round(model32.rsquared_adj,3),'VERSUS',round(model.rsquared_adj,3))            
    XC = sm.add_constant(X)
    modelC = sm.OLS(Y,XC).fit()
    if len(xlabel)==1:
        if modelC.f_pvalue > 0.05:
            print('@@@@@@@@@@@@@',connect, xlabel,'NOT SIGNIFICANT !!!!!!!!!!!!!!')
            if extended:
                print(model.summary())     #by default no extended output printed within function
                print('THIS ONE MORE RELEVANT:::', modelC.summary())
            return model,[]
        else:
            print ('@@@ Connect',connect,'Indep', xlabel,'no Const',model.f_pvalue, model.pvalues.iloc[0],model.rsquared_adj,
                   '   \nwith Const',modelC.f_pvalue,modelC.pvalues.iloc[0],modelC.rsquared_adj)
            if extended:
                print(model.summary())     #by default no extended output printed within function
            return model,model.predict(X)
    elif add_const:
        return modelC,modelC.predict(XC)
    else:
        return model,model.predict(X)

def cross_corr(df,fname_lists,dend_types):
    lags={dend_types[ct]:{} for ct in np.unique(df.TYPE)}
    decay={dend_types[ct]:{} for ct in np.unique(df.TYPE)}
    mean_xcorr={dend_types[ct]:{} for ct in np.unique(df.TYPE)}
    estimate={dend_types[ct]:{} for ct in np.unique(df.TYPE)}
    for ctnum in np.unique(df.TYPE):
        ct=dend_types[ctnum]
        for ts in fname_lists.keys():
            xcorr={}
            for fn in fname_lists[ts]:
                df_subset=df[(df.TYPE==ctnum)&(df.FNAME==fn)]
                norm_data=df_subset.DIAMETER-np.mean(df_subset.DIAMETER)
                xcorr[fn]=np.correlate(norm_data,norm_data,'full')/len(norm_data)/np.var(norm_data)
            minlength=np.min([len(xc) for xc in xcorr.values()])//2
            mean_xcorr[ct][ts]=np.mean([xc[len(xc)//2+1:len(xc)//2+1+minlength] for xc in xcorr.values()],axis=0)
            lags[ct][ts]=np.arange(0,minlength)
            param_bounds2 = ([0,0,-100,-100,-100],[lags[ct][ts][-1],lags[ct][ts][-1],100, 100,100])
            param_bounds1 = ([0,-100,-100],[lags[ct][ts][-1],100,100])
            #FIXME send in either exp1 or exp2 function, and but then save either 1 or 2 taus, and param_bounds 1 or 2
            popt, pcov = optimize.curve_fit(exp2,lags[ct][ts],mean_xcorr[ct][ts],bounds=param_bounds2,maxfev=5000)
            print(ct,ts,'XCORR DECAY',sorted(popt[0:2]))
            #decay[ct][ts]=round(popt[0],3)
            decay[ct][ts]=(round(np.min(popt[0:2]),3),round(np.max(popt[0:2]),3))
            estimate[ct][ts]=exp2(lags[ct][ts],*popt)
    return mean_xcorr,lags,decay

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
    return [cd for clist in container for cd in clist]

'''Split N sequences into 'equal' sizes'''
def split_seq(seq, size):               #splits file_list into equal amounts for testing-training datasets
        newseq = []
        splitsize = 1.0/size*len(seq)
        for i in range(size):
            newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
        return newseq

def read_morphologies(root,dirs):
    header = {}                                
    complete_list = []
    
    '''Read in Morphology Data in Folder with Single/Multiple Archive(s)'''
    for d1 in dirs:                                      #can accept individual archives or multiple archives within single directory           
        fullpath = root + d1 + '/*CNG_extract.txt'                      
        print('Working on Directory : ', str(d1))
        for fname in glob.glob(fullpath):                               #locates and loops through _extract files in path
            temp_name =  os.path.basename(fname).split('.txt')[0]
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
                                if point != header['CHILD'] and point != header['TYPE'] and point != header['PARENT'] and point != header['TYPE'] :
                                    temp_line[point] = float(temp_line[point])
                                if point==header['PAR_CONNECT']:
                                    temp_line[point]=int(temp_line[point])
                            temp_line.extend([temp_line[header['RADIUS']]*2,temp_line[header['PARENT_RAD']]*2,temp_name,d1]) 
                            complete_list.append(temp_line)          #complete_list will be used to test equations changing radius
                                                                                                                                        
        '''Initial Data Separation by .swc Compartment Types (Basal/Apical if present) and Archive'''
        if len(temp_line) > len(header):                             
            for i in ['DIAMETER','PARENT_DIA','FNAME','ARCHIVE']:    #header with new added parameters
                header[i] = len(header)
    morph_df=pd.DataFrame(complete_list,columns=header.keys())
    return morph_df,complete_list,header

def compare_EM(df,path):
    ###################### Compare to EM data ################### 
    #proximal thick: from 0 to 100 microns, range 1.8-2.5
    #medial: from100 to 250 microns, range 1.6-2.2
    #distal: from 250 - 450, range 1.0-1.5
    #radiatum thin: 0.15-0.4
    #str LM - > 450 um, range 0.8-1.2; 0.3-0.8; 0.15-0.4
    #basal: proximal: 0.50-0.9; distal: 0.25-0.45

    import scipy
    from scipy.stats import percentileofscore
    megias={'ap':{'dist':[0,100],'diam':[0.45,2.5]},
            'am':{'dist':[100,250],'diam':[0.45,2.2]},
            'ad':{'dist':[250,450],'diam':[0.45,2.1]},
            'lm':{'dist':[450,1000],'diam':[0.15,1.2]}}
    wilson={'prim':{'nd':0,'diam':[1.25,2.5]}, #wilson reports 2.25
            'sec':{'nd':1,'diam':[0.75,1.25]}, #Wilson reports 1.0
            'tert':{'nd':2,'diam':[0.29,0.63]}} #Wilson reports 0.29 - 0.63

    if 'Hippocampus' in path.split('/'):
        for cat,vals in megias.items():
            ca1=df[(df.ARCHIVE=='Groen') & (df.TYPE=='4')]
            mintile=percentileofscore(ca1[(ca1.PATH_DIS>vals['dist'][0]) & (ca1.PATH_DIS<=vals['dist'][1])]['DIAMETER'],score=vals['diam'][0])
            maxtile=percentileofscore(ca1[(ca1.PATH_DIS>vals['dist'][0]) & (ca1.PATH_DIS<=vals['dist'][1])]['DIAMETER'],score=vals['diam'][1])
            print(cat, vals['dist'], 'min', round(mintile,3),'max',round(maxtile,3))
    elif 'Basal_Ganglia' in path.split('/'):
        for cat,vals in wilson.items():
            mintile=percentileofscore(df[(df.NODE_ORDER>=vals['nd'])]['DIAMETER'],score=vals['diam'][0])
            maxtile=percentileofscore(df[(df.NODE_ORDER>=vals['nd'])]['DIAMETER'],score=vals['diam'][1])
            print(cat, vals['nd'], 'min', round(mintile,1),'max',round(maxtile,1))

