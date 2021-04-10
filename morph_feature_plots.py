'''<utilities for graphing morphology features extracted from .swc files>
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
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
from scipy.stats import pearsonr
import morph_feature_utils as mfu

'''Save plot as png'''            
def save_png(png, title):                       #default to save plots instead of plot to window
    png.savefig(title, bbox_inches = 'tight')   #remove empty space around plots
    print('File Created : ' + str(title))
    png.close()
    
def main_plot(data, var, title = None, ax_titles = None, labels = None, save_as = None, plot_type = None, fit_line = None, add = None, where = None, ax_max = None):
    if save_as:
        scale=2
    else:
        scale=1
    plt.rc('font', size = scale*10)                       
    colors=['tab:blue','tab:orange','tab:purple']
    fig, ax = plt.subplots(figsize = (scale*5,scale*2.9))
    if add:                                      #additional information to plot
        at = AnchoredText(add, prop=dict(size=scale*14), frameon=True, loc='upper center')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") 
        ax.add_artist(at)
       
    if plot_type == 'Color3d':                    #plot 3 variables(x,y,z) with variable 'z' using colorbar
        scat = ax.scatter(data[var[0]], data[var[1]], c=data[var[2]], s=2, marker='o', label = labels)
        fig.colorbar(scat)
        
    elif plot_type == 'Fit':                      #plot initial fitted equation to data values
        plt.plot(data[var[0]], data[var[1]], 'o')
        plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) 
    else:
        for num in range(len(data)):
            if labels:                            #can modify markersize or alpha (transparency) of plotted points
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', color=colors[num],label = labels[num], markersize = scale*2)
            else:
                plt.plot(data[num][var[0]], data[num][var[1]], 'o', color=colors[num],markersize = scale*2)

    if labels or fit_line:
        plt.legend(loc = where) if where else plt.legend()
    plt.xlabel(ax_titles[0]) if ax_titles else plt.xlabel(var[0])
    plt.ylabel(ax_titles[1]) if ax_titles else plt.ylabel(var[1])
    if ax_max:
        plt.ylim([0,ax_max])
        plt.xlim([0,ax_max])
    plt.tight_layout()
    if save_as:
        fig.savefig(save_as+'.tiff')
        save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')

'''Main Plot Function for Basic, Fit or 3d Color Plots'''
def df_plot(data, var, title = None, select_param=None, ax_titles = None, labels = None, save_as = None, plot_type = None, fit_line = None, add = None, where = None, ax_max = None):
    if save_as:
        scale=2
    else:
        scale=1
    colors=['tab:blue','tab:orange','tab:purple']
    plt.rc('font', size = scale*10)  #34                     
    fig, ax = plt.subplots(figsize = (scale*5,scale*2.5)) #(20,10)
    if add:                                      #additional information to plot
        at = AnchoredText(add, prop=dict(size=28), frameon=True, loc='upper center')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") 
        ax.add_artist(at)
       
    if plot_type == 'Color3d':                    #plot 3 variables(x,y,z) with variable 'z' using colorbar
        scat = ax.scatter(data[var[0]], data[var[1]], c=data[var[2]], s=2, marker='o', label = labels)
        fig.colorbar(scat)
        
    elif plot_type == 'Fit':                      #plot initial fitted equation to data values
        plt.plot(data[var[0]], data[var[1]], 'o')
        plt.plot(fit_line[0], fit_line[1], color = 'orange', label = fit_line[2]) 
    else:
        if labels:                            #can modify markersize or alpha (transparency) of plotted points
            for i,(kind,label) in enumerate(labels.items()):
                subset=data[data[select_param]==kind]
                if len(subset):
                    plt.plot(subset[var[0]], subset[var[1]], 'o', color=colors[i],label = label, markersize = scale*2) #make 4 for inset
        else:
                plt.plot(data[var[0]], data[var[1]], '.', markersize = scale*2)

    ####### make axis scales an option
    '''
    if var[0]=='PARENT_DIA':
        plt.xscale('log')
        plt.yscale('log')
    '''
    if title:
        plt.title(title)                           #will require 'title' if nothing passed in 'save_as'

    if labels or fit_line:
        plt.legend(loc = where) if where else plt.legend() #turn off for inset
    plt.xlabel(ax_titles[0]) if ax_titles else plt.xlabel(var[0])
    plt.ylabel(ax_titles[1]) if ax_titles else plt.ylabel(var[1])
    if ax_max:
        plt.ylim([0,ax_max])
    plt.tight_layout()
    if save_as:
        fig.savefig(save_as+'.tiff')
        save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')

'''Plot Parameter Correlation to Diameter with Pearson's R'''
def corr(data,selection,title,save_as = None):
    if save_as:
        scale=2
    else:
        scale=1
    dframe = pd.DataFrame(data,columns = selection.keys())  #dataframe reorders alphabetically feature columns by default 
    dcorr = dframe.corr()                                   #set feature order by features --> passed as 'selection'
    for i in dcorr:
        dcorr[i] = round(dcorr[i]**2,4)                              #dataframe correlation as R2
    plt.rc('font', size = scale*10/1.2)            
    fig,ax = plt.subplots(figsize = (scale*5,scale*4))
    #sns.set(font_scale = 4)
    rsquared = '$R^2$'                                      #can modify heatmap for either (+/-) R or (+) R2 for Feature Correlations
    #ax = sns.heatmap(dcorr,vmin = -1,vmax = 1, center = 0,cmap = sns.diverging_palette(20, 220, n=256),square = True,annot = True,cbar_kws={'shrink':1,'label':rsquared})#"Pearsons R"}
    ax = sns.heatmap(dcorr,vmin = 0,vmax = 1, center = 0.5,cmap = 'Blues',square = True,annot = True,cbar_kws={'label':rsquared})
    ax.set_yticklabels([sub['short'] for sub in selection.values()],rotation = 0)
    ax.set_xticklabels([sub['short'] for sub in selection.values()],rotation = 0,horizontalalignment = 'right',ha = 'center')
    cbar = ax.collections[0].colorbar                        #selection will contain short and long versions of feature names
    plt.title(title, loc = 'left')
    plt.tight_layout()
    if save_as:
        fig.savefig(save_as+'.tiff')
        save_png(plt,save_as + '.png') if save_as else save_png(plt,title + '.png')
    
def plot_hist_set (histograms,binset,dtype,features,title=None,save_as=None):
    if save_as:
        scale=2
    else:
        scale=1
    rows=int(round(np.sqrt(len(histograms))))
    cols=int(np.ceil(len(histograms)/rows))
    plt.rc('font', size = 10)            
    fig,axis=plt.subplots(nrows=rows,ncols=cols,figsize = (scale*4,scale*3))
    ax=fig.axes
    if title:
        plt.title(title+ dend_types[dtype])
    for i,feat in enumerate(histograms.keys()):
        for ar in histograms[feat][dtype].keys():
            ax[i].plot(binset[feat][0:-1],histograms[feat][dtype][ar]/np.sum(histograms[feat][dtype][ar]),label=ar) #dividing by sum converts counts to frequency
        ax[i].set_xlabel(features[feat]['long'])
        ax[i].set_ylabel('fraction')
    ax[1].legend()
    plt.tight_layout()
    if len(save_as):
        fig.savefig(save_as+'.tiff')

def Rall_plot(df,dirs,comp_list,save_as=''):
    if save_as:
        scale=1.7
    else:
        scale=1    
    parent_dia32={ct:{ar:[] for ar in dirs} for ct in np.unique(df.TYPE)}
    child_dia32={ct:{ar:[] for ar in dirs} for ct in np.unique(df.TYPE)}
    for ct in np.unique(df.TYPE): #loop over apical or basal
        for ar in dirs:  #loop over archive
            bp_df=df[(df.PAR_CONNECT==float(comp_list.index('BP_Child'))) & (df.TYPE==ct) & (df.ARCHIVE==ar)] #extract branch point children
            for fn in np.unique(bp_df['FNAME']):
                rall_dict={}
                for parent in np.unique(bp_df[(bp_df.FNAME==fn) & (bp_df.TYPE==ct)]['PARENT']): 
                    child_df=bp_df[(bp_df.FNAME==fn) & (bp_df.PARENT==parent)] #find all children
                    rall_dict[parent]={'dia':child_df.iloc[0]['DIAMETER'],'child':list(child_df['DIAMETER'])}
                parent_dia32[ct][ar].append([float(a['dia']) for a in rall_dict.values()])
                child_dia=[a['child'] for a in rall_dict.values()]
                #transfer **1.5 for parent, to np.sum(child**1.5)**(2/3) to units of x & y axes are um
                child_dia32[ct][ar].append([np.sum([float(b)**1.5 for b in a])**(2/3) for a in child_dia ])
    ####### Plotting part ##########  Ideally put in separate function
    plt.rc('font', size = scale*10)  #34                     
    colors=['tab:blue','tab:orange','tab:purple']
    for ct in parent_dia32.keys():
        xymax=0
        fig=plt.figure()
        if not save_as:
            plt.title('Rall test, dend type '+ ct)
        for i,ar in enumerate(parent_dia32[ct].keys()):
            x=list(mfu.flatten(parent_dia32[ct][ar]))
            y=list(mfu.flatten(child_dia32[ct][ar]))
            corr=pearsonr(x,y)[0]
            xymax=np.max([np.max(x),np.max(y),xymax])
            print('corr of parent to (sum child^1.5)^(2/3)',ar, round(corr*corr,3))
            plt.plot(x,y,'.',color=colors[i],label=ar.ljust(10)+' $R^2$='+str(round(corr*corr,3)),markersize = scale*4)
 #dividing by sum converts counts to frequency
        plt.plot([0,xymax],[0,xymax],'gray')
        plt.xlabel('parent diameter ($\mu$m)')
        plt.ylabel('\u03A3 child diameter ($\mu$m)')
        plt.legend(loc='lower right')
        plt.tight_layout()
        if save_as:
            fig.savefig(save_as+ct+'.tiff')
    return parent_dia32,child_dia32

def plot_xcorr(mean_xcorr,lags,decay,xmax=None,save_as=''):
    if save_as:
        scale=2
    else:
        scale=1
    for ct in mean_xcorr.keys():
        plt.rc('font', size = scale*10)  #34                     
        fig, ax = plt.subplots(figsize = (scale*5,scale*2)) #(20,10)
        #plt.title('cross-correlation '+ct)
        min_xcorr=0
        min_xcorr=min(0,np.min([np.min(xc) for k in mean_xcorr.keys() for xc in mean_xcorr[k].values() ]))
        for ts in mean_xcorr[ct].keys():
            label=ts
            for i in range(len(decay[ct][ts])):
                label=label+', $\lambda$'+str(i+1)+'='+str(round(decay[ct][ts][i],1))
            plt.plot(lags[ct][ts],mean_xcorr[ct][ts],label=label,linewidth=scale*1)
            #plt.plot(lags[ct][ts],estimate[ct][ts],'gray')
        if xmax:
            plt.xlim([0,xmax])
        plt.ylim([min_xcorr,1])
        plt.xlabel('Lag (number of nodes)')
        plt.ylabel ('Correlation')
        plt.legend()
        plt.tight_layout()
        if len(save_as):
            fig.savefig(save_as+ct+'.tiff')

 
