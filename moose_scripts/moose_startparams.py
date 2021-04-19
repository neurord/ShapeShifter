#this will be the location of all initial parameter values (Vm, RA, RM etc.)

#use available parameters from our morph comparisons
#convert to SI
#for Lindroos cell #ohm-m2 #Lindroos used 8, but could go as low as 1 to account for Kir
params={'Lindroos': {'RM':8.0,'RA':0.15,'CM':0.01}, #F/m^2 = 1 uF/cm^2 
        'default':{'RM': 1.6,'RA': 1.98,'CM' :0.0186}} #F/m2

synparams={'channame':'ampa','tau1':1e-3,'tau2':5e-3,'Erev':-5e-3,'gbar':1e-9}
#may need up to 5e-9
#sending in sim_params in line
#units of time in sec.

#short injection
#short_width = 1e-3
#short_simtime = 100e-3

#long injection
#long_width = 400e-3
#long_simtime = 800e-3

delay0 = 1e-3
#width = long_width#100e-3 #change from 100e-3 to see if spiking 
level = 1e-9 
delay1 = 1e9

#default simdt is 5e-5
#simdt = 1e-3
simdt = 1e-5

#adjusted values from swc/p file
#aVm = 70e-3
#aCM = .03
#aRM = 2.8
#aRA = 4.0
#initVm = -70e-3 #seems equivalent to Eleak 

