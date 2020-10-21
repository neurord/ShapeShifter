'''<Modifies morphology files (.p) for comp. simulations from types below>
   <Condenses .p files of neuron traces to maintain specific lambda (ac or dc) and removes 0 length compartments.>
    Copyright (C) <2016>  <Saivardhan Mada>
    Copyright (C) <2017>  <Avrama Blackwell>
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

#Usage: python shape_shifter.py --file 'filename.p' --type 'condense'
#Types:
#  '0'        just remove compartments of size 0
#  'condense' combine compartments of similar radius (specify --rad_diff 0 to only combine identical radii),
#               electrotonic length not to exceed max_len* lambda.
#  'expand'   to change single, long compartment (e.g. a Neuron software segment) into multiple smaller compartments
#               electrotonic length of subdivided compartments do not exceed max_len* lambda
#  'radii'    to change the radius value depending on the distance to the end of the longest branch
# can specify alternative values to default rm [4], cm [0.01], ri [2.5] units are SI
# can specify frequency (--f) for ac lambda calculation, current default is 0.1 hz
# can specify maximum electrotonic length (--max_len) current default is 0.1, specify 0 to use dc lambda
# can specify maximum difference (rad_diff) in radii for combining adjacent compartments, default is 0.1 (=10%)
# change info or debug (below) to get additional information while running program


#Ideal Usage Scenario:
#neuromorpho.org morphology files often inaccurate compartment (node) diameters

#run shape_shifter.py through type radii for more accurate diameters
#can currently accept either .swc or .p morphology files
#requires model.txt from morph_feature_analysis.py to load in MLR model

# Input:   python shape_shifter.py --file morph_file.swc --type radii --model model.txt
# Ouput:   morph_file_org.p
 
#run shape_shifter.py again through type condense to combine similar compartments together for future simulation
#can also pass additional parameters i.e. f (frequency) and radius difference from optional arguments above
# Input:   python shape_shifter.py --file filenameconvert.p --type condense --f 100 --rad_diff 0
# Ouput:   filenameconvertoutout.p

#TO DO:
# test remove_zeros and condense (--type 'condense') using additional p files
#  - validate voltage response with simulations 
# read (and use) rm, cm and ri from p file
# write rm, cm and ri to p file if (a) not already specified and (b) given as input parameters

#George Mason University
#Saivardhan Mada
#Jul 20th, 2016
#Avrama Blackwell
#Apr 28, 2017
#Jonathan Reed
#Oct 21, 2020

from __future__ import print_function, division
import sys
import os
import glob
import math
import argparse
import collections
import numpy as np
from scipy.stats import pearsonr
from extract_dynamic import create_dict
import matplotlib.pyplot as plt
#from graph_extract import simple_plot
#from graph_extract import save_png

#testing recursion
#sys.setrecursionlimit(1500)
debug=0   #default is 0
info=0    #default is 1
#Globals for compartment values stored in list structure
CHILD = 0; PARENT = 1; X = 2; Y = 3; Z = 4; DIA = 5; COMPLEN = 6; END_DIS = 7
class morph:
        
        def __init__(self,input_file):
                self.linelist = self.read_file(input_file)
                self.soma_count()
                self.remove_zeros()
               
        def read_file(self,input_file):                 #read either .swc or .p morphology file and convert to .p if needed for morph-changes
                copy_lines = []
                lines = open(input_file, 'r').readlines()
                for line in lines:                      #avoid comment lines
                    if line[0] !='*' and line[0] !='/' and line[0] != '\n' and line[0] != '\r' and line[0] != '#':
                        copy_line = line.split()
                        for x in range(2,6):            #change certain values to floats for calculations
                                copy_line[x] = float(copy_line[x])
                        copy_lines.append(copy_line)
                        
                _,extension = os.path.splitext(input_file)

                if extension == '.p': #default file format for morph-changes
                        return copy_lines
                        
                elif extension == '.swc': #change to .p file format
                        newlines = []; parlist = []
                        for x in range(5):
                                parlist.append(list(zip(*copy_lines))[x]) #list of parent values
                        for num,line in enumerate(copy_lines):            #.p file slightly different order/values to .swc
                            newline = line[:6]
                            newline[0] = (str(line[0]) + '_' + str(line[1]))
                            if num == 0:                                  #soma is always first line, with unique 'none' in .swc
                                newline[1] = 'none'
                            else:
                                parent_index = parlist[0].index(line[6])  #parlist[0] refers to parent location of comp at line[6]
                                parent_type = parlist[1][parent_index]
                                newline[1] = str(line[6]) + '_' + str(parent_type)
                            newline[5] = round(newline[5]*2,4)     #.p file format takes diameter values instead of radius
                            newlines.append(newline)
                        return newlines

        def soma_count(self):                            #count number of soma nodes, usually 3-pt or 1-pt
                soma_start = self.linelist[0]            #by default, the first soma node should connect to all other nodes
                soma_nodes = 0
                for line in self.linelist:               #check for all soma nodes with _1 ending
                        if '_1' in line[0]:
                                soma_nodes = soma_nodes + 1
                if soma_nodes == 1:
                        print('***1-pt soma detected***')
                if soma_nodes == 3:
                        print('***3-pt soma detected***')
                        comp1 = self.linelist[1]; comp2 = self.linelist[2]
                        check_if_parent = []
                        for num,comp in enumerate(self.linelist):   #remove extra soma nodes which do not connect to others (only visuallizes width)
                                check_if_parent.append('true')
                        if len(check_if_parent) == 0 and soma_start[DIA] == comp1[DIA] == comp2[DIA] and soma_start[CHILD] == comp1[PARENT] == comp2[PARENT]:
                                self_linelist = [line for line in self.linelist if line != comp1 and line != comp2]
                        else:
                                ('***Multiple soma node connections; possibly check morphology file***') #usually extra 3-pt soma nodes exact diameter equal to soma
                elif soma_nodes == 2 or soma_nodes > 3:
                        print('***More than 3-pt soma; possibly check morphology file***') #some morphology files contain multiple soma nodes
                        
        def remove_zeros(self): 
                newlines=[]; parent_dict = {}
                #Calculate compartment length from XYZ values
                for line in self.linelist:
                        complen = np.sqrt((line[X])**2 + (line[Y])**2 + (line[Z])**2)
                        line.append(complen) 
                        
                #Remove Zeroes and Replace oldparent id with newparent id
                for line in self.linelist:
                        #removes any compartments with 0 distance or diameter OTHER THAN soma
                        if((line[DIA] == 0  or line[COMPLEN] == 0) and line[PARENT] != 'none'):
                                #relabel remaining children to previous parent (or compartment that is deleted)
                                # if                        4,3 0 0 0 <--- compartment needs to be removed since size 0
                                # then                      5,4 x y z <--- original (CHILD) compartment needs PARENT change to existing compartment
                                # should be changed to      5,3 x y z <--- new parent connection to remaining children  
                                parent_dict[line[CHILD]] = line[PARENT]
                                print ("######## Zero size comp removal: orig num lines=", line[CHILD], ", new num lines=", line[PARENT])
                        else:
                                newlines.append(line)
                                
                #Call to change parent if removing 0_compartments
                if parent_dict: #checks if parent_dict has any items, THEN does recursion
                        print('Zero point compartments removed; remaining comparments require parent-update for reconnection')
                        print(parent_dict)
                        parent_len = len(parent_dict.items())
                        if parent_len < 100: #arbitrary number; 100 zero-point compartments
                                for num, line in enumerate(self.linelist):
                                        if line[PARENT] != 'none':
                                                counter = 0 
                                                line = self.replace_par(parent_dict, counter, line)
                                print(len(newlines))
                                self.linelist = newlines
                        else:
                                print('Many zero-point compartments --> check morphology file')

        #replace parent values of compartments whose original parents were zero-point compartments and removed 
        def replace_par(self,parent_dict, counter, line):
                print('Original line', line)
                if line[PARENT] in parent_dict:
                        print('Line child', line[CHILD], 'Line parent', line[PARENT])
                        counter = counter + 1
                        line[PARENT] = parent_dict[line[PARENT]]
                        self.replace_par(parent_dict, counter, line)
                else:
                        if counter>0:
                                print('fixed line', counter ,line)
                        return line



def calc_lambda (type1, RM, RI, CM, F):
        tm = RM * CM #time constant
        dc_factor = math.sqrt(RM/(RI*4.0))*1000 #1000 is for unit conversion to microns 
        partial_ac= 2.0*math.pi*F*tm 
        ac_factor = dc_factor * math.sqrt(2.0/(1.0+math.sqrt(1.0+math.pow(partial_ac,2.0))))
        return ac_factor 

def calc_electrotonic_len(line,factor):
        lamb = factor * math.sqrt(float(line[DIA])) #calculates ac or dc lambda from compartment radius
        if debug:
                print('    calc_len comp %s lambda=%.1f len=%.3f L=%.5f' % (line[CHILD], lamb, line[COMPLEN], line[COMPLEN]/lamb))
        return line[COMPLEN]/lamb

def write_file(output, line):
        write_line = [str(val) for val in line[CHILD:COMPLEN]]
        write_line = ' '.join(write_line) #converts from list of strings to single string
        output.write(write_line + '\n')   #writes every node as new line

def calc_newcomp(condense,surface_tot,Ltot,lamb_factor):
        diameter=math.pow(surface_tot/(Ltot*np.pi*lamb_factor),(2/3.))
        l=surface_tot/(np.pi*diameter)
        if debug:
                print('TSA: dia %.5f, l=%.2f tsa=%.1f Ltot=%.5f' % (diameter, l, np.pi*diameter*l,l/(lamb_factor*math.sqrt(diameter))))
        x = 0.0
        y = 0.0
        z = 0.0
        #total distance in x, y and z.  
        for comp in condense:
                x = np.sqrt(x)
                y = np.sqrt(y)
                z = np.sqrt(z)
        if debug:
                print ('xtot %.3f ytot %.3f ztot %.3f'%(x, y, z))
        if (x > 0) :
                theta = np.arctan(y/x)
        else:
                theta = 0
        if (x>0 or y>0 or z>0):
                phi = np.arccos(z/math.sqrt(x*x+y*y+z*z))
        else:
                phi = 0
        if debug:
                print ('theta %.4f phi %.4f ' %(theta,phi))
        x = np.round(l*math.cos(theta)*math.sin(phi),3)
        y = np.round(l*math.sin(theta)*math.sin(phi),3)
        z = np.round(l*math.cos(phi),3)
        return str(x),str(y),str(z),str(np.round(diameter,3))

def subdivide_comp(newline,segs):  #only for type expand --> non-working version
        newlen = 0; lineset = []
        print ("subdivide", line[CHILD], "x,y,z,len", line[X], line[Y], line[Z], line[COMPLEN], "into", segs, "segments")
        for j in [X,Y,Z]: 
                newline[j]=np.round(line[j]/segs,5) 
        seg_length = np.sqrt((newline[X])**2 + (newline[Y])**2 + (newline[Z])**2) 
        for i in range(segs):
                lineset.append(line)
                if i > 0:
                        lineset[i][PARENT] = lineset[i-1][CHILD]
                        lineset[i][CHILD] = lineset[0] + str(i)
                newlen=newlen+seg_length
                #save new parent_dict to connect subsequent compartment to new segmented one either here or in type = 'expand'
        if info:
                for i in range(segs):
                        print ("new seg", i, myname[i],PARENT[i], newcomp[i])
                print ("total length", newlen, "seg", seg_length )
        return lineset

def to_condense(condense,surface_tot,Ltot,lambda_factor,rad_diff,delta_rad,line,comp2,outfile):
        if len(condense):
                #cannot add any additional compartments.  Condense the set
                x,y,z,diameter=calc_newcomp(condense,surface_tot,Ltot,lambda_factor)
                if info:
                        print ('#######condense', condense, 'stop before', comp2[CHILD])
                        print ('new: x y z',x,y,z, diameter )
                newcomp = [condense[-1][CHILD], condense[0][PARENT], x, y, z, diameter]
                write_file(outfile, newcomp)
                condense = []  #reinitializing to condense next set of compartments
                Ltot = 0
                surface_tot = 0
        else:
                #cannot condense this comp and nothing in condense set.  Just print out line
                write_file(outfile, line)
                if info:
                        print ('not condensing',len(condense),Ltot, line[CHILD], "rad", rad_diff,delta_rad)
        return condense, Ltot, surface_tot

def branch_comp(linelist):
        parents = [line[PARENT] for line in linelist] #create a list of just parents
        parent_count = collections.Counter(parents)   #create dictionary giving number of times each parent occurs
        branch_points = [line[CHILD] for line in linelist if parent_count[line[CHILD]]>1] #set of parents that appear more than once
        #print('branch points', branch_points)
        return branch_points, parents

def end_comp(linelist, parents):
        unique_parents=set(parents) #set() eliminates duplicate parents
        children=[line[CHILD] for line in linelist] #create a list of just children
        end_points=list(set(children)-unique_parents)
        print(end_points)
        return end_points

def calc_enddis(dis_to_end, par_comp, line, ending = 0):
        if ending == 0:                                     #if endpoint
                dis_to_end = dis_to_end + line[COMPLEN]/2
        else:                                               #all other comps. except soma
                dis_to_end = dis_to_end + line[COMPLEN]/2 + par_comp
                par_comp = line[COMPLEN]/2
        end_dis = dis_to_end
        return dis_to_end, par_comp, end_dis


def condenser(m, type1, max_len, lambda_factor, rad_diff):
        #num_comps=0
        ######## type = '0' removes 0 length compartments, this is done in the class, so just write file
        filename,_ = os.path.splitext(h.file)
        if(type1 == "0"):
                removed = open(filename + '_removed.p','w')
                removed.write('*relative' + '\n')
                for line in m.linelist:
                        write_file(removed,line)
                removed.close()
                print('Removed Zero point Compartments : Created ' + filename + '_removed.p')

        ####### type = "expand" takes long compartments and subdivides into multiple "segments"  --> in progress non-working version
        if (type1 == "expand"):
                #Expands linelist to include new created segments
                
                new_par = {} #holds all compartments whose parent needs to change with included segments
                newlinelist = []
                for num,line in enumerate(m.linelist):
                        L_comp = calc_electrotonic_len(line,lambda_factor)
                        print ("max_len", max_len, "L", L_comp)
                        if L_comp > max_len:
                                print('Adding segments to: ', line[CHILD])
                                segs=int(math.ceil(L_comp/max_len))
                                newlines=subdivide_comp(line,segs)
                                for seg in newlines:
                                        print('Segment added: ', seg[CHILD])
                                new_par[m.linelist[num][CHILD]] = newlines[-1][CHILD]
                                newlinelist.append(newlines)
                        else:
                                newlinelist.append(line)
                        if info:
                                print ("OK", line)

                #Replaces parent-child of any compartment whose parent was expanded
                print(len(new_par),new_par)
                print(len(newlinelist), newlinelist)
                
                #for line in newlinelist:                #change parent-child here
                #m.replace_parent[comp1[0]]=myname[-1] maybe NOT call replace par as it would change ALL subsequent parent_child connections
                #ONLY want to change the final newlines comp to original subsequent one
                #must replace parent of sudsequent compartment with last segmented comp.
                #maybe do a search here for the added comp. maybe outside this for loop too?
                #maybe do the re.search here for ANY subsequent compartment that would possibly connect
                #maybe add the expand_list to correct position into linelist outside for loop
                #write statement for entire line_list
        branch_points, parents = branch_comp(m.linelist)
        ######## type = condense condenses branches with similar radius and combined electronic length < 0.1 lambda
        if (type1 == "condense"):    #if rad_diff = 0, only condenses branches with same radius
                #cvapp converts .swc to .p file format AND from absolute to relative coordinates
                #print(m.linelist)
                #print('Original Compartments : ' + len(m.linelist[0]))
                condensed = open(filename + '_condensed.p','w')
                condensed.write('*relative' + '\n')
                Ltot = 0; surface_tot = 0; condense = [] #list to hold compartments to be condensed
                for num,line in enumerate(m.linelist):
                        comp1 = line
                        if debug:
                                print ('**** begin', line[CHILD])
                        if num <= (len(m.linelist)-2):
                                comp2 = m.linelist[num+1]
                                len_comp1=calc_electrotonic_len(comp1,lambda_factor)
                                len_comp2=calc_electrotonic_len(comp2,lambda_factor)
                        if len(condense): #assumes condense has values
                                delta_rad=abs((condense[0][DIA] - comp2[DIA])/condense[0][DIA])
                                tot_len=len_comp2+Ltot
                        else:
                                delta_rad=abs((comp1[DIA] - comp2[DIA])/comp1[DIA])
                                tot_len=len_comp1 + len_comp2
                        #rad_diff default is 0.1, max_len default is 0.1, see parameters at bottom
                        if info:
                                print ('comp1', comp1[CHILD], comp1[PARENT], 'comp2', comp2[CHILD], comp2[PARENT])
                        if(delta_rad <= rad_diff and comp2[PARENT]==comp1[CHILD] and tot_len < max_len and comp1[CHILD] not in branch_points):
                                print('***comps to be condensed: ', comp1[CHILD],comp2[CHILD], '***')  #branch point list has to be created before this point
                                if info:
                                        print('delta_rad <= rad_diff for :', 'comp1', comp1, 'comp2', comp2)
                                length1 = math.sqrt(comp1[X]**2 + comp1[Y]**2 + comp1[Z]**2)
                                length2 = math.sqrt(comp2[X]**2 + comp2[Y]**2 + comp2[Z]**2)
                                #if the radii are almost the same and comp2 is attached to comp1
                                #if this is last line of file, comp2=comp1, thus comp2[1] != comp1[0], skip down to condensing or writing file
                                if len(condense) == 0:
                                        #if this is 1st compartment of a set, add it
                                        condense.append(comp1)
                                        Ltot = len_comp1
                                        surface_tot = math.pi * comp1[DIA] * length1
                                #always add the 2nd compartment to set of comartments to be condensed
                                condense.append(comp2)
                                Ltot=Ltot+len_comp2
                                surface_tot=surface_tot+(math.pi * comp2[DIA] * length2)    
                        else:
                                condense, Ltot, surface_tot = to_condense(condense,surface_tot,Ltot,lambda_factor,rad_diff,delta_rad,line,comp2,condensed)
                                
                #temp_name = m.filename.split('.swc')[0] + '_condense.swc'
                #print('Created : ' + m.filename  )#+ ' with ' + len(m.linelist) + ' Condensed Compartments')
                #print(len(m.outfile))
                #print ("finished,", len(m.linelist), "output compartments")
                print('Condensed Compartments : Created ' + filename + '_condensed.p')
        ######## type = radii changes radius with distance dependent diameter equation from Lindroos et al 2018 'Basal Ganglia Neuromodulation'
        if (type1 =="radii"):
                #read in model, as morph file already read in as m.linelist
                read_model = open(h.model,'r').readlines()
                models = {}
                for line in read_model:
                        line = line.split()
                        if line[0] == 'Apical' or line[0] == 'Basal':
                                comp = line[0]
                                if not line[0] in models.keys():
                                        models[comp] = {}
                        elif line[0] == 'Initial' or line[0] == 'BP_Child' or line[0] == 'Continuing':
                                connect = line[0]
                                if not line[0] in models[comp].keys():
                                        models[comp][connect] = {}
                        elif len(line) == 2: 
                                models[comp][connect][line[0]] = float(line[1])

                
                feature_file = new_file + '_extract.txt'
                features = create_dict(feature_file)
                for feature in features.keys():                    #remove trailing '\n' from final feature name
                        if '\n' in feature:
                                new_key = feature.strip('\n')
                                features[new_key] = features.pop(feature)
                                
                soma_count = 0
                newrads = {}; newrads_i = {}
                pdict = {0:'Initial',1:'BP_Child',2:'Continuing'}   #original data (#'s) to dictionary organization (terms)
                cdict = {3:'Basal', 4:'Apical'}
                
                for num,line in enumerate(m.linelist):
                        newrad = 0; newrad_i = 0
                        if '_1' in line[CHILD]:
                                soma_count = soma_count + 1
                        else:
                                new_num = num - soma_count   #features do not have soma extended values
                                ctype = cdict[features['TYPE'][new_num]]; pcon = pdict[features['PAR_CONNECT'][new_num]]
                                for feature in models[ctype][pcon]:
                                        if feature == 'PARENT_DIA':              
                                                if pcon == 'Initial':
                                                        newrad = newrad + models[ctype][pcon][feature] * (features['RADIUS'][new_num]*2)
                                                        newrad_i = newrad_i + (features['RADIUS'][new_num]*2)
                                                else:
                                                        newrad = newrad + models[ctype][pcon][feature] * newrads[line[PARENT]]
                                                        newrad_i = newrad_i + models[ctype][pcon][feature] * newrads_i[line[PARENT]]
                                        else:                              
                                                newrad = newrad + models[ctype][pcon][feature] * features[feature][num]
                                                if not '_1' in line[PARENT]:
                                                        newrad_i = newrad_i + models[ctype][pcon][feature] * features[feature][new_num]
                                newrads[line[CHILD]] = newrad
                                newrads_i[line[CHILD]] = newrad_i
                
                org_file = open(filename + '_org.p','w')
                pred_file = open(filename + '_pred.p','w')
                pred_i_file = open(filename + '_pred_i.p','w')
                zdia1_file = open(filename + '_zdia1.p','w')
                zdia2_file = open(filename + '_zdia2.p','w')
                
                for fname in [org_file,pred_file,pred_i_file,zdia1_file,zdia2_file]:
                        fname.write('*relative' + '\n')

                for num,line in enumerate(m.linelist):
                        if num >= soma_count:             
                                write_file(org_file,line)
                                pred_line = line
                                pred_line[DIA] = newrads[line[CHILD]]
                                write_file(pred_file,pred_line)
                                pred_i_line = line
                                pred_i_line[DIA] = newrads_i[line[CHILD]]
                                write_file(pred_i_file,pred_i_line)
                                zdia1_line = line
                                zdia1_line[DIA] = 1.0
                                write_file(zdia1_file,zdia1_line)
                                zdia2_line = line
                                zdia2_line[DIA] = 2.0
                                write_file(zdia2_file,zdia2_line)
                        else:
                                for fname in [org_file,pred_file,pred_i_file,zdia1_file,zdia2_file]:
                                        write_file(fname,line)   
                org_file.close()
                pred_file.close()
                pred_i_file.close()
                zdia1_file.close()
                zdia2_file.close()
                
 
if __name__ == '__main__':
        #set default parameters
        rm = 4.00 #ohms-m^2 
        ri = 2.50 #ohms-m
        cm = 0.010 #F/m^2
        max_len = 0.1 #electrotonic length
        f = .1 #Hz
        rad_diff=0.1

        #set up argument parsers
        parser = argparse.ArgumentParser()
        parser.add_argument('--file')
        parser.add_argument('--type', choices={'0','condense','expand','radii'}, default='condense') #should default be 0
        #parser.add_argument('--path', type = str) 
        parser.add_argument('--model')
        parser.add_argument('--save_as') #will need to reformat file if .swc output desired
        
        parser.add_argument('--rad_diff', default=rad_diff, type=float)
        parser.add_argument('--rm', default=rm, type=float)
        parser.add_argument('--ri', default=ri, type=float)
        parser.add_argument('--cm', default=cm, type=float)
        parser.add_argument('--f', default=f, type=float)
        parser.add_argument('--max_len', default=max_len, type=float)
        
        try:
	        args = ARGS.split() #in python: define space-separated ARGS string
	        do_exit = False
        except NameError: #if you are not in python, read in filename and other parameters
	        args = sys.argv[1:]
	        do_exit = True
        h = parser.parse_args(args)
        newmorph=morph(h.file)
        
        #calculate lambda
        lambd_factor=calc_lambda(h.type, h.rm, h.ri, h.cm, h.f)
        #print('params', h, 'lambda', lambd_factor)
        #Optionally, can condense multiple comps into one, expand large comp into multiple, or assign radii when there are none
        condenser(newmorph, h.type, h.max_len, lambd_factor, h.rad_diff)
