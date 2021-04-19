'''
Jonathon Reed, Avrama Blackwell
Analysis of simulation of neuron morphologies
comparing original morphology with predicted radii
Using two types of current injection and 1 location of synaptic input
2021 March 12
'''
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.ticker import StrMethodFormatter
import numpy as np
from scipy import optimize
import os
import glob

convert2milli=1000

def save_png(png, title):               #default to save plots instead of plot to window
    png.savefig(title, bbox_inches = 'tight')
    print('File Created : ' + str(title))
    png.close()

def exp2(t,tau1,tau2,A,B,C):
    return A*np.exp(-t/tau1)+B*np.exp(-t/tau2)+C #0th tau --> single soma as RM*CM
                                                 #tau1 = 0th tau(1+(n*pi/elec_len)**2) where n = tau #
                                                 #tau --> short sims beginning at peak
                                                 #change in voltage --> long sims @ steady reach
def exp1(t,tau1,A,B):
    return A*np.exp(-t/tau1)+B #0th tau --> single soma as RM*CM

def main_plot(data, var, colors,title = '', ax_titles = None, labels = None, save_as = None, where = 'best', size = None,xmax=None):
    #plt.rc('font', size = 36) #default plot and font sizes
    if save_as:
        scale=2
    else:
        scale=1
    plt.rc('font', size = scale*9)
    fig=plt.figure(figsize = (scale*4,scale*2))     
    for key,traces in data.items():
        if labels:
            plt.plot(data[key][var[0]], data[key][var[1]], label = labels[key], linewidth = scale*1,color=colors[key])
        else:
            plt.plot(data[key][var[0]], data[key][var[1]], linewidth = scale*1,color=colors[key])
            
    if labels:
        plt.legend(loc = where)
    if xmax:
        plt.xlim([0,xmax])
    plt.xlabel(ax_titles[0]) if ax_titles else plt.xlabel(var[0])
    plt.ylabel(ax_titles[1]) if ax_titles else plt.ylabel(var[1])
    #plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    #plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    plt.tight_layout()
    if save_as:
        #save_png(plt,save_as + '.png')# if save_as else save_png(plt,title + '.png')
        fig.savefig(save_as+'.tiff')
    else:
        plt.title(title)
        plt.show()
    return plt
        
def show_curvefit(func,sim_data,popt,rkey,start,comp='somaVm',title=None):
    estimate = func(np.array(sim_data[rkey]['time'][start:]),*popt) #send in all curve_fit params for estimate
    #show result of curve fit
    plt.figure()
    if title:
        plt.title(title)
    else:
        plt.title(rkey)
    plt.plot(sim_data[rkey]['time'][start:],sim_data[rkey][comp][start:],'o')
    plt.plot(sim_data[rkey]['time'][start:], estimate)
    #save_png(plt,'check_eq.png')
    #plt.close()

_FUNCTIONS={
    'exp1': exp1,
    'exp2': exp2}
####### Uncomment the root file name (pattern) and synapse compartment of one of these #########
#only saving syntau and dVsyn in soma for one of the synapses
pattern = 'Purkinje-slice-ageP35-1.CNG_*' #matches quite well to pred+in.  Looks very similar to original plot (shape, not values)
syn=['12_3']  #synaptic input ,,'1342_3',
#pattern = 'WT-0201MSN03.CNG_*' #Lai striatum - matches quite well, including synaptic input
#syn=['19_3']#,'35_3',
#pattern = 'WT-P270-20-14ak.CNG_*' #Lindroos striatum - 
#syn=['35_3']#'194_3',
#pattern = 'ri06_mod.CNG_*' #Golding hippocampus - good fit, good synapse.  
#syn=['31_4']#,'3400_3'] #
#pattern = 'Adol-20100419cell1.CNG_*' #Groen - hippocampus - good fit, good synapse. 
#syn=['32_3']#,'112_3'] #response is ~0.1 mV 
### to fit to single exponential, change to 'exp1'
fit='exp2'
param_bounds = ([0,0,-100,-100,-100],[1000,1000, 100,100,100])#msec,msec,mV,mV,mV
showfit=True

#parameters to control plotting, specific to the set of simulations, but same for all file names
inj_set={'short' : '0.001', 'long': '0.8'}
syn_set={'syn'+str(i):syn[i] for i in range(len(syn))}
stim_set=dict(syn_set,**inj_set)
labels={'org':'Original','pred':'Predict','pred_i':'Predict+In.','zdia2':'2.0 Diameter'}
colors={'org':'k','pred':'b','pred_i':'g','zdia2':'gray'}
ax_titles = ['Time (ms)','Vm (mV)']
#dictionary to store summary measures
measures={'tau':{},'tau2':{},'deltaV':{},'syntau':{},'syntau2':{},'dVsyn':{}}
fit_cov={}
for stimtype,stim in stim_set.items():
    fnames=sorted(glob.glob(pattern+stim+'*.txt'))
    sim_data = {}
    for fname in fnames:
        rkey=fname.split(stim)[0].split('CNG_')[-1]
        header=open(fname).readline()
        data_keys=header.split()[1:]
        sim_data[rkey] = {dk:[] for dk in data_keys}
        print('processing',fname,' .p file type=',rkey,' data columns',data_keys)
        data=np.loadtxt(fname)
        for i, dk in enumerate(data_keys):
            sim_data[rkey][dk]=data[:,i]*convert2milli
        if len(data_keys)>2:
            if 'dVsyn'+'_'+data_keys[2] not in measures.keys():
                for kk in ['dVsyn','syntau','syntau2']:
                    measures[kk+'_'+data_keys[2]]={}
        if stimtype=='short': #determine time constants
            start = np.argmax(sim_data[rkey]['somaVm']) #point at which is highest voltage
            #print('exp start', start,sim_data[rkey]['time'][start])
            if fit=='exp2':
                #curve fit to double tau, needs parameter constraints
                popt,pcov = optimize.curve_fit(exp2,sim_data[rkey]['time'][start:],sim_data[rkey]['somaVm'][start:],bounds=param_bounds,maxfev=5000)
                measures['tau'][rkey]=max(popt[0],popt[1])
                measures['tau2'][rkey]=min(popt[0],popt[1])
            else:
                popt,pcov = optimize.curve_fit(exp1,sim_data[rkey]['time'][start:],sim_data[rkey]['somaVm'][start:])
                measures['tau'][rkey] = popt[0]
            fit_cov[rkey]=pcov[0,0]
            if showfit:
                show_curvefit(_FUNCTIONS[fit],sim_data,popt,rkey,start)
        elif stimtype=='long': #determine deltaV (proportional to input resistance)
            delay = 50e-3
            pt = np.max(np.where(np.array(sim_data[rkey]['time'])<delay))
            baseV=sim_data[rkey]['somaVm'][pt]
            peakpt=np.argmax(sim_data[rkey]['somaVm'])
            peakV=sim_data[rkey]['somaVm'][peakpt]
            measures['deltaV'][rkey] = peakV-baseV
        elif stimtype.startswith('syn'): #determine peak time and peak Vm
            delay=10e-3
            for comp in list(sim_data[rkey].keys())[1:]:
                peakpt=np.array(sim_data[rkey][comp]).argmax()
                print('   analysis of',stimtype,comp, 'for', rkey, 'peakpt=',peakpt)
                peakV=sim_data[rkey][comp][peakpt]
                basept = np.max(np.where(np.array(sim_data[rkey]['time'])<delay))
                baseV=sim_data[rkey][comp][basept]
                if comp=='somaVm':
                    add_to_key=''
                else:
                    add_to_key='_'+comp
                if fit=='exp2':
                    popt,pcov = optimize.curve_fit(exp2,sim_data[rkey]['time'][peakpt:],sim_data[rkey][comp][peakpt:],bounds=param_bounds,maxfev=5000)
                    measures['syntau2'+add_to_key][rkey]=min(popt[0],popt[1])
                    measures['syntau'+add_to_key][rkey]=max(popt[0],popt[1])
                else:
                    popt,pcov = optimize.curve_fit(exp1,sim_data[rkey]['time'][peakpt:],sim_data[rkey][comp][peakpt:])
                    measures['syntau'+add_to_key][rkey]=popt[0]
                measures['dVsyn'+add_to_key][rkey]=peakV-baseV
                if showfit:
                    show_curvefit(_FUNCTIONS[fit],sim_data,popt,rkey,peakpt,comp=comp,title=rkey+add_to_key)

    #plot all morphologies for one type of stimulation
    title=fname.split('.')[0]
    if title == 'ri06_mod':
        title = 'ri06'
    print(title)
    savename = title+'_'+stim_set[stimtype]
    main_plot(sim_data, data_keys, colors, title=title+' '+stim_set[stimtype], ax_titles=ax_titles, labels = labels, where = 'upper right',save_as=savename)
    if stimtype.startswith('syn'):
        savename = title+'_'+'dend_'+stim_set[stimtype]        
        main_plot(sim_data, np.array(data_keys)[[0,2]], colors, title='dend '+stim_set[stimtype], ax_titles=ax_titles, labels = labels, where = 'upper right',xmax=200, save_as=savename)

#calculate the difference as ratio abs((orig-new)/orig)
diff={k:{} for k in measures.keys()}
for meastype,measure in measures.items():
    if len(measure):
        for key in measure.keys():
            diff[meastype][key]=np.abs(measure[key]/measure['org']-1)
        print ('*******', meastype,'******')
        print([(k,round(m,5)) for k,m in measure.items()],'\ndiff = ',[(k,round(d,3)) for k,d in diff[meastype].items()])
        if meastype=='tau':
            print('fit', fit_cov[key])

##### Save measures to a file #####3            
outfname=pattern.split('.')[0]+'.measures'
f=open(outfname,'w')
for key in diff:
    f.write(key+'='+str(measures[key])+'\n')
    f.write('diff='+str(diff[key])+'\n')
f.close()

####### to save figure as tif:
#import PIL
#fig.savefig('fname.tif')
  




