#   module add python/2.7.15
#   module add moose
#   python
#>>> execfile('nameoffile.py')

import moose
import numpy as np
from numpy import linspace
from pprint import pprint
import moose_my_func as f
import moose_startparams as p
import sys


#filename, sim_duration, sim amplitute, simtime, optionally compartment name for synaptic stimulation
args = sys.argv[1:]
fname = args[0]
delay = float(args[1])
dur = float(args[2])
amp = float(args[3])
stime = float(args[4])
if fname.startswith('WT-P270'):
    print('Using Lindroos parameters')
    pas=p.params['Lindroos']
else:
    pas=p.params['default']
if len(args)>5:
    syncompname=args[5]
    stim='syn'
    amp=0
    outfname=fname+syncompname+'_RM'+str(round(pas['RM']))+'.txt'
    if len(args)>6:
        p.synparams['gbar']=float(args[6]) #overwrite default value
else:
    outfname=fname + str(dur)+'_RM'+str(round(pas['RM'])) + '.txt'
    stim='inject'
cname = 'cell'
cell = moose.loadModel(fname + '.p',cname)
comps = []
for comp in moose.wildcardFind('%s/#[TYPE=Compartment]'%(cell.name)):
    comps.append(f.adjust_comp(comp, pas['RM'], pas['CM'], pas['RA']))
    #comps.append(f.adjust_comp(comp, RM, CM, RA))
somaname = comps[0].path
injectcomp = moose.element(somaname)
soma = moose.element(somaname)

if stim=='syn':
    syncomp=moose.element(cell.path+'/'+syncompname)
    synproto=f.make_synchan(p.synparams)
    synchan,sh=f.add_synchan(p.synparams['channame'],syncomp.path,p.synparams['gbar'])
    spiketime=[delay]
    stimtab=f.create_tt(spiketime)
    f.connect_synchan(sh,stimtab)
if stim=='inject':
    soma_pulse = f.create_pulse('somapulse', soma, delay, dur, amp, p.delay1) #setup for current injection
    currenttab, vmtabs = f.output([soma],10*p.simdt,soma_pulse) #<--- plot/out dt
   
else:
    currenttab, vmtabs = f.output([soma,syncomp],10*p.simdt) #<--- plot/out dt
print('CURRENTTAB: ', currenttab,'VMTABS: ', vmtabs,' simdt (before hsolve)', soma.dt)

#################### set up hsolve
hsolve=moose.HSolve(cell.path+'/hsolve')
hsolve.dt=p.simdt
hsolve.target=somaname
print('after hsolve, simdt=', hsolve.dt)
moose.reinit()
moose.start(stime)
'''
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
plt.ion()
for vmtab in vmtabs:
    plt.plot(vmtab.vector*1000,label=vmtab.name)
plt.legend()
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
'''
header = 'time somaVm   '
if stim=='syn':
    header=header+syncompname+'Vm'
f.outfile(outfname, stime, header, vmtabs) 

'''
if swcfile: 
    cell = moose.loadModel(swcfile,cellname)
    #print('cell = ', cell, 'cellname = ', cellname, 'cell.name = ', cell.name)
    for comp in moose.wildcardFind('/' + cell.name + '/#[TYPE=Compartment]'):
        comp = f.adjust_comp(comp, pas['RM'], pas['CM'], pas['RA'])
    #for here soma from swcfile has to be changed
    somaname = cellname + '/soma[0]'  #this is the same path for either morphology file or create neuron

'''
