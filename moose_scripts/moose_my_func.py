import moose
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

def make_synchan(params):
    if not moose.exists('/library'):
        lib=moose.Neutral('/library')
    else:
        lib=moose.element('/library')
        
    synchan=moose.SynChan('/library/'+params['channame'])
    synchan.tau1=params['tau1']
    synchan.tau2=params['tau2']
    synchan.Ek=params['Erev']
    return synchan

def add_synchan(chanpath,syncomp,gbar):
    proto=moose.element('/library/'+chanpath)
    print('add_synchan',proto,syncomp,chanpath)
    synchan=moose.copy(proto,syncomp,chanpath)[0]
    synchan.Gbar=gbar
    m=moose.connect(syncomp,'channel',synchan,'channel')
    sh=moose.SimpleSynHandler(synchan.path+'/SH')
    moose.connect(sh,'activationOut',synchan,'activation')
    return synchan,sh

def create_tt(spiketime):
    stimtab=moose.TimeTable('presyn')
    stimtab.vector=spiketime
    return stimtab

def connect_synchan(sh,presyn):
    numsyn=sh.synapse.num
    sh.synapse.num=numsyn+1
    sh.synapse[numsyn].delay=1e-3
    m=moose.connect(presyn,'eventOut',sh.synapse[numsyn],'addSpike')
    print('connection',m)
    return

def create_comp(cellname, compname, RM, CM, RA, initVm, comp_len = None, comp_dia = None,  Vm = 70e-3):
    '''Create morphology compartments with parameter values'''
    '''Create new compartments'''
    if ~moose.exists(compname):
        comp = moose.Compartment(cellname + '/' + compname)
        print('comp created :' + comp.name)
        if comp.length != None and comp.diameter != None: #this is for reading in morphology file
            comp.diameter = comp_dia                      #already has diameter and length values 
            comp.length = comp_len                        #so new values only for created compartments
    comp.initVm = initVm
    SA = np.pi * 2 * (comp.diameter/2) * comp.length
    comp.Rm = RM/SA
    comp.Cm = CM*SA
    x_area = np.pi * comp.diameter * comp.diameter/4
    comp.Ra = RA*comp.length/x_area
    return comp

def adjust_comp(comp, RM, CM, RA):#, Em, initVm):
    #if moose.exists('/' + cellname + '/' + compname):
        #comp = moose.element('/' + cellname + '/' + compname) #Element
    #print('old Rm',comp.path,comp.Rm/1e6,'old Cm,Ra',comp.Cm*1e6,comp.Ra)
    SA = np.pi * comp.diameter * comp.length
    if SA == 0:
            print('comp',comp.name, 'surface area',SA,'comp length', comp.length, 'old Rm',comp.path,comp.Rm/1e6,'old Cm,Ra',comp.Cm*1e6,comp.Ra)
    xa = np.pi * comp.diameter * comp.diameter/4
    comp.Rm = RM/SA
    comp.Cm = CM*SA
    comp.Ra = RA*comp.length/xa
    return comp
    #print('new Rm',comp.path,comp.Rm/1e6,'new Cm,Ra',comp.Cm*1e6,comp.Ra)
    #comp.initVm = Vrest
    #comp.Em = Wrest
    #if cellname == '1_1':
        #comp.length = comp.diameter
    #return comp_name #create_comp returns comp, so comp_name is python var for comp

#creates pulse generator
#should I put creation of output directory in create_pulse func?
#since it needs 'output' to exist?
def create_pulse(pulse_name, comp, delay0, width, level, delay1):
    pulse = moose.PulseGen(pulse_name)
    print(pulse)
    pulse.delay[0] = delay0
    pulse.width[0] = width
    pulse.level[0] = level
    pulse.delay[1] = delay1
    moose.connect(pulse, 'output', comp, 'injectMsg')
    return pulse #do I need to return pulse?

#creates output as separate moose.Table
def output(comp,plot_dt,pulse=None):  
    if ~moose.exists('/data' ):
        data = moose.Neutral('/data' )
    if pulse:
        current_tab = moose.Table('/data' +  '/current' + comp[0].name)
        moose.setClock(current_tab.tick,plot_dt)
        moose.connect(current_tab, 'requestOut', pulse, 'getOutputValue')
    else:
        current_tab=[]
    vmtabs=[]
    for c in comp:
        vm_tab = moose.Table('/data' + '/Vm' + c.name)
        vmtabs.append(vm_tab)
        moose.setClock(vm_tab.tick,plot_dt)
        moose.connect(vm_tab, 'requestOut', c, 'getVm')
    return current_tab, vmtabs

def outfile(name, endtime, headtxt, vmtabs):
    ts = np.linspace(0, endtime, len(vmtabs[0].vector))
    data = [ts]+[vmtab.vector for vmtab in vmtabs]
    #print(data)
    np.savetxt(name, np.column_stack(data), header = headtxt)
    #return data

def save_png(png, title):               #default for all plots instead of plotting to window
    png.savefig(title + '.png')
    print('File Created : ' + str(title))
    png.close()

def make_plot(data,title,labels):
    fig,ax = plt.subplots(figsize = (20,10))
    plt.rc('font', size = 28)
    plt.plot(data[0],data[1])
    plt.title = title
    plt.xlabel = labels[0]
    plt.ylabel = labels[1]
    save_png(plt,title + '.png')

def compare_plot(data, title, legend, labels, add = None, marker_size = None, where = None):
    fig, ax = plt.subplots(figsize = (20,10))
    plt.rc('font', size = 28)
    
    if add:
        at = AnchoredText(add,
                  prop=dict(size=20), frameon=True,
                  loc='upper center',
                  )
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2") #maybe look for ways to remove border around annotation
        ax.add_artist(at)
    
    for num,i in enumerate(data):
        if num == 2:
            plt.plot(data[i][0],data[i][1],color='green',label = legend[num])
        elif num == 3:
            plt.plot(data[i][0],data[i][1],color='limegreen',label = legend[num])
        else:
            plt.plot(data[i][0],data[i][1], label = legend[num])
        #plt.plot(data[i][0],data[i][1], label = legend[num])
    #plt.plot(data[0][0],data[0][1], color = 'b', label = legend[0])#, markersize = 10)
    #plt.plot(data[1][0],data[1][1], color = 'orange', label = legend[1])#, markersize = 10)
    plt.legend()
    plt.xlabel(labels[0], fontsize = 28)
    ax.tick_params(axis='x', labelsize=28)
    ax.tick_params(axis='y', labelsize=28)
    plt.ylabel(labels[1], fontsize = 28)
    plt.title(title, fontsize = 28)
    save_png(plt,title)
