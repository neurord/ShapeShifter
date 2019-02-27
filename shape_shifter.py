'''<Condenses .p files of neuron traces to maintain specific lambda (ac or dc) and removes 0 length compartments.>
    Copyright (C) <2016>  <Saivardhan Mada>
    Copyright (C) <2017>  <Avrama Blackwell>
    Copyright (C) <2019> <Jonathan Reed>
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
#type can be:
#  '0' just remove compartments of size 0
#  'condense' combine compartments of similar radius (specify 0 to only combine identical radii),
#          electrotonic length not to exceed max_len* lambda.
#  'expand' to change single, long compartment (e.g. a Neuron software segment) into multiple smaller compartments
#           electrotonic length of subdivided compartments do not exceed max_len* lambda
#  'radii'    to change the diameter value depeneding on the distance to the end of the longest branch
# can specify alternative values to default rm [4], cm [0.01], ri [2.5] units are SI
# can specify frequency (--f) for ac lambda calculation, current default is 0.1 hz
# can specify maximum electrotonic length (--max_len) current default is 0.1, specify 0 to use dc lambda
# can specify maximum difference (rad_diff) in radii for combining adjacent compartments, default is 0.1 (=10%)
# change info or debug (below) to get additional information while running program

#TO DO:
# test remove_zeros and condense (--type 'condense') using additional p files
#  - validate voltage response with simulations 
# read (and use) rm, cm and ri from p file
# write rm, cm and ri to p file if (a) not already specified and (b) given as input parameters
# test for whether absolute or relative coordinates.  If absolute - print error message (until length calculation can be updated)

#George Mason University
#Saivardhan Mada
#Jul 20th, 2016
#Avrama Blackwell
#Apr 28, 2017

from __future__ import print_function, division
import sys
import math
import argparse
import collections
import numpy as np

debug=0   #default is 0
info=0    #default is 1
#Globals for compartment values stored in linelist
child = 0; parent = 1; X = 2; Y = 3; Z = 4; dia = 5; complen = 6; end_dis = 7 #globals make all caps

class morph:
        def __init__(self,pfile):
                self.readfile(pfile)
                self.remove_zeros()
                
        def readfile(self,pfile):
                '''Read from .p file and set to list'''
                self.pfile=pfile
                self.linelist = []
                out_name=pfile.split('.p')[0]+'out.p'
                self.outfile = open(out_name, 'w')
                lines=open(pfile, 'r').readlines()
                self.lines = lines
                print('***Number of compartments in .p file: ', len(self.lines), '***')

                for line in lines:
                        if line[0] !='*' and line[0] !='/' and line[0] != '\n' and line[0] != '\r':
                                if len(line) > 5:                           #alternative code instead of list comprehension below
                                        self.linelist.append(line.split())  #the empty lines between text actually have len = 2, so they are not caught in len = None
                                else:
                                        self.outfile.write('\n')
                                #print('a',len(line),line)
                        else:
                                self.outfile.write(line) #writes commented portion of original file to created file
                                #print('b',len(line),line)
                #self.linelist = [x for x in self.linelist if x]  #fixes prior error which added empty lists to linelist, see comments above
                
                '''Conversion from strings --> floats takes place here (XYZ and diameter)'''
                for num,line in enumerate(self.linelist):
                        for id, val in enumerate(line):
                                if id >= X and id <= dia:
                                        line[id] = float(line[id])
                                else:
                                        line[id] = str(line[id])
                                        
                '''Checks for 1-pt or 3-pt soma'''
                line = self.linelist[0]; comp1 = self.linelist[1]; comp2 = self.linelist[2]
                parcomplist = []
                if line[X] == 0 and line[Y] == 0 and line[Z] == 0:
                        print('***Soma comp. detected with XYZ coordinates 0,0,0***')
                        for num, comp in enumerate(self.linelist): #checks if 2 comps. connected to soma do not connect to other comps.
                                if comp[parent] == comp1[child] or comp[parent] == comp2[child]:
                                        parcomplist.append('true')
                        if (comp1[dia] == line[dia] and comp2[dia] == line[dia]) and (comp1[parent] == line[child] and comp2[parent] == line[child]):
                                #check if soma comps. have same diameter and check if soma comps. connected to soma[0]
                                if not parcomplist: 
                                        for val in [X,Y,Z]:
                                                line[val] = comp1[val]-comp2[val]   #Change XYZ values of soma; no longer 0,0,0
                                        del self.linelist[2]
                                        del self.linelist[1]
                                        print('***Soma comps adjusted --> compare original and output files***')
                else:
                        print('***Soma comp. detected w/o XYZ coordinates 0,0,0***')
                        
        def remove_zeros(self): #change distance to comp_len
                newlines=[]; parent_dict = {}
                '''Calculate line distance'''
                for line in self.linelist:
                        distance = np.sqrt((line[X])**2 + (line[Y])**2 + (line[Z])**2)
                        line.append(distance) #may move down to type = radii since not used elsewhere
                        
                '''Remove Zeroes and Replace oldparent id with newparent id'''
                for line in self.linelist:
                        #removes any compartments with 0 distance or diameter OTHER THAN soma
                        if((line[dia] == 0  or line[complen] == 0) and line[parent] != 'none'):
                                #relabel remaining children to previous parent (or compartment that is deleted)
                                #i.e. if                    4,3 0 0 0 <--- needs to be removed
                                # then                      5,4 x y z
                                # should be changed to      5,3 x y z <--- new parent connection to remaining children  
                                parent_dict[line[child]] = line[parent]
                                print ("######## Zero size comp removal: orig num lines=", line[child], ", new num lines=", line[parent])
                        else:
                                newlines.append(line)
                                
                '''Call to change parent if removing 0_compartments'''
                if parent_dict: #checks if parent_dict has any items, THEN does recursion
                        print('Zero point compartments removed; remaining comparments require parent-update for reconnection')
                        print(parent_dict)
                        parent_len = len(parent_dict.items())
                        if parent_len < 100: #just a placeholder for now; i.e. 100 zero-point compartments
                                for num, line in enumerate(self.linelist):
                                        if line[parent] != 'none':
                                                counter = 0 
                                                line = self.replace_par(parent_dict, counter, line)
                                print(len(newlines))
                                self.linelist = newlines
                        else:
                                print('Many zero-point compartments --> check morphology file')

        #replace parent values of compartments whose original parents were zero-point compartments and removed 
        def replace_par(self,parent_dict, counter, line):
                print('Original line', line)
                if line[parent] in parent_dict:
                        print('Line child', line[child], 'Line parent', line[parent])
                        counter = counter + 1
                        print(counter)
                        line[parent] = parent_dict[line[parent]]
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
        lamb = factor * math.sqrt(float(line[dia])) #calculates ac or dc lambda from compartment radius
        if debug:
                print('    calc_len comp %s lambda=%.1f len=%.3f L=%.5f' % (line[child], lamb, line[complen], line[complen]/lamb))
        return line[complen]/lamb

def write_file(output, line):
        for x in range(X,complen):
                line[x] = round(float(line[x]),4)
        write_line = [str(val) for val in line[child:complen]]
        write_line = ' '.join(write_line) #converts from list of strings to single string
        output.write(write_line + '\n')   #writes every compartment as new line

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

def subdivide_comp(newline,segs):
        #newlines is var for remove_zeros function
        newlen = 0; lineset = []
        print ("subdivide", line[child], "x,y,z,len", line[X], line[Y], line[Z], line[complen], "into", segs, "segments")
        for j in [X,Y,Z]: 
                newline[j]=np.round(line[j]/segs,5) 
        seg_length = np.sqrt((newline[X])**2 + (newline[Y])**2 + (newline[Z])**2) 
        #must add lines to extent of segmentation
        for i in range(segs):
                lineset.append(line)
                if i > 0:
                        lineset[i][parent] = lineset[i-1][child]
                        lineset[i][child] = lineset[0] + str(i)
                newlen=newlen+seg_length
                #save new parent_dict to connect subsequent compartment to new segmented one either here or in type = 'expand'
        if info:
                for i in range(segs):
                        print ("new seg", i, myname[i],parent[i], newcomp[i])
                print ("total length", newlen, "seg", seg_length )
        return lineset

def to_condense(condense,surface_tot,Ltot,lambda_factor,rad_diff,delta_rad,line,comp2,outfile):
        if len(condense):
                #cannot add any additional compartments.  Condense the set
                x,y,z,diameter=calc_newcomp(condense,surface_tot,Ltot,lambda_factor)
                if info:
                        print ('#######condense', condense, 'stop before', comp2[child])
                        print ('new: x y z',x,y,z, diameter )
                newcomp = [condense[-1][child], condense[0][parent], x, y, z, diameter]
                #print out newcomp if above works but not in file
                write_file(outfile, newcomp)
                condense = []  #reinitializing to condense next set of compartments
                Ltot = 0
                surface_tot = 0
        else:
                #cannot condense this comp and nothing in condense set.  Just print out line
                write_file(outfile, line)
                if info:
                        print ('not condensing',len(condense),Ltot, line[child], "rad", rad_diff,delta_rad)
        return condense, Ltot, surface_tot

def branch_comp(linelist):
        parents = [line[parent] for line in linelist] #create a list of just parents
        parent_count = collections.Counter(parents) #create dictionary giving number of times each parent occurs
        branch_points = [line[child] for line in linelist if parent_count[line[child]]>1] #set of parents that appear more than once
        print('branch points', branch_points)
        return branch_points, parents

def end_comp(linelist, parents):
        unique_parents=set(parents) #set() eliminates duplicate parents
        children=[line[child] for line in linelist] #create a list of just children
        end_points=list(set(children)-unique_parents)
        print(end_points)
        return end_points

def condenser(m, type1, max_len, lambda_factor, rad_diff):
        num_comps=0
        ######## type = '0' removes 0 length compartments, this is done in the class, so just write file
        if(type1 == "0"): 
                for line in m.linelist:
                        write_file(m.outfile, line)

        ####### type = "expand" takes long compartments and subdivides into multiple "segments"
        if (type1 == "expand"):
                '''Expands linelist to include new created segments'''
                new_par = {} #holds all compartments whose parent needs to change with included segments
                newlinelist = []
                for num,line in enumerate(m.linelist):
                        L_comp = calc_electrotonic_len(line,lambda_factor)
                        print ("max_len", max_len, "L", L_comp)
                        if L_comp > max_len:
                                print('Adding segments to: ', line[child])
                                segs=int(math.ceil(L_comp/max_len))
                                newlines=subdivide_comp(line,segs)
                                for seg in newlines:
                                        print('Segment added: ', seg[child])
                                new_par[m.linelist[num][child]] = newlines[-1][child]
                                newlinelist.append(newlines)
                        else:
                                newlinelist.append(line)
                        if info:
                                print ("OK", line)

                '''Replaces parent-child of any compartment whose parent was expanded'''
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

        branch_points, parents = branch_comp(m.linelist) #used in both type condense and type radii
        
        ######## type = condense condenses branches with similar radius and combined electronic length < 0.1 lambda
        if (type1 == "condense"):    #if rad_diff = 0, only condenses branches with same radius
                #cvapp converts .swc to .p file format AND from absolute to relative coordinates
                Ltot = 0; surface_tot = 0; condense = [] #list to hold compartments to be condensed
                for num,line in enumerate(m.linelist):
                        comp1 = line
                        if debug:
                                print ('**** begin', line[child])
                        if num <= (len(m.linelist)-2):
                                comp2 = m.linelist[num+1]
                                len_comp1=calc_electrotonic_len(comp1,lambda_factor)
                                len_comp2=calc_electrotonic_len(comp2,lambda_factor)
                        if len(condense): #this already assumes condense has values
                                delta_rad=abs((condense[0][dia] - comp2[dia])/condense[0][dia])
                                tot_len=len_comp2+Ltot
                        else:
                                delta_rad=abs((comp1[dia] - comp2[dia])/comp1[dia])
                                tot_len=len_comp1 + len_comp2
                        #rad_diff default is 0.1, max_len default is 0.1
                        if info:
                                print ('comp1', comp1[child], comp1[parent], 'comp2', comp2[child], comp2[parent])
                        if(delta_rad <= rad_diff and comp2[parent]==comp1[child] and tot_len < max_len and comp1[child] not in branch_points):
                                print('***comps to be condensed: ', comp1[child],comp2[child], '***')  #branch point list has to be created before this point
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
                                        surface_tot = math.pi * comp1[dia] * length1
                                #always add the 2nd compartment to set of comartments to be condensed
                                condense.append(comp2)
                                Ltot=Ltot+len_comp2
                                surface_tot=surface_tot+(math.pi * comp2[dia] * length2)    
                        else:
                                condense, Ltot, surface_tot = to_condense(condense,surface_tot,Ltot,lambda_factor,rad_diff,delta_rad,line,comp2,m.outfile)
        print ("finished,", len(m.linelist), "output compartments")

        ######## type = radii changes radius with distance dependent diameter equation from Lindroos et al 2018 'Basal Ganglia Neuromodulation'
        if (type1 =="radii"):
                end_points = end_comp(m.linelist, parents)
                print('branch points', branch_points)
                print('end points', end_points)
                for line in m.linelist:
                        line.append(0) #will hold temporary value of distance to end of longest branch at given compartment as line[end_dis] is filled
                for line in reversed(m.linelist):
                        if line[child] in end_points:
                                #initialize dis to end
                                #update par_comp
                                #call func.
                                dis_to_end = 0 #will be temporary value to added to end_dis as loop moves towards soma compartment
                                line[end_dis] = line[complen]
                                par_comp = line[complen]/2
                                dis_to_end = dis_to_end + line[complen]/2
                                print('endpoint found at', line[child], line[end_dis])                 #if line is endpoint --> begins end branch length
                        elif line[child] in branch_points:                                             #if line is branchpoint either:
                                if line[end_dis] == 0 or dis_to_end > line[end_dis]:                                                            #add to end branch length
                                        line[end_dis] = dis_to_end + line[complen]/2 + par_comp                                #if already has branch length (another branch)
                                        print('branchpoint found at', line[child], line[end_dis])
                                        dis_to_end = dis_to_end + line[complen]/2 + par_comp
                                        par_comp = line[complen]/2                                     #compare which branch length is larger
 
                        else:                                                                           #if line is normal compartment --> add end branch length 
                                dis_to_end = dis_to_end + line[complen]/2 + par_comp
                                par_comp = line[complen]/2
                                line[end_dis] = dis_to_end 

                for num, line in enumerate(m.linelist[1:]):
                        #assuming soma diameter does not need to be changed
                        d1 = 0.001 * line[end_dis]  + 0.87
                        dp = [x[dia] for x in m.linelist if x[child] == line[parent]] 
                        dp = dp[0]             #a compartment can be a parent for multiple children, but each will only have one parent itself
                        l = abs(line[complen])
                        d2 = dp * np.exp(-l * 0.08) #proximal diameters not reproduced accurately <-- put in commit message
                        line[dia] = max(d1,d2)
                        
                if info:
                        print('self', '| parent', '| compartment diameter', '| compartment distance to end')
                for line in m.linelist:
                        if info:
                                print(line[child],line[parent],line[dia],line[end_dis])
                        write_file(m.outfile, line)                       #will remove temporary values of complen and end branch distance as file is created
                #where: Lindroos et al. 2018 Basal Ganglia Neuromodulation
                #L = dendritic length of branch rooted from given node i.e. distance to the end of longest branch from the given compartmentm, end_dis
                #dp = diameter of parent node i.e. dia
                #l = distance to parent node i.e. complen
                
if __name__ == '__main__':
        #set default parameters
        rm = 4.00 #ohms-m^2 
        ri = 2.50 #ohms-m
        cm = 0.010 #F/m^2
        max_len = 0.1 #electrotonic length
        f = .1 #Hz
        rad_diff=0.1
        type1 = "condense"
        #
        #set up argument parsers
        parser = argparse.ArgumentParser()
        parser.add_argument('--file')
        parser.add_argument('--type', choices={'0','condense','expand','radii'}, default='condense')
        parser.add_argument('--rad_diff', default=rad_diff, type=float)
        parser.add_argument('--rm', default=rm, type=float)
        parser.add_argument('--ri', default=ri, type=float)
        parser.add_argument('--cm', default=cm, type=float)
        parser.add_argument('--f', default=f, type=float)
        parser.add_argument('--max_len', default=max_len, type=float)
        #
        try:
	        args = ARGS.split() #in python: define space-separated ARGS string
	        do_exit = False
        except NameError: #if you are not in python, read in filename and other parameters
	        args = sys.argv[1:]
	        do_exit = True
        h = parser.parse_args(args)
        #
        #reads p file, delete zero size compartments
        newmorph=morph(h.file)
        #
        #calculate lambda
        lambd_factor=calc_lambda(h.type, h.rm, h.ri, h.cm, h.f)
        print('params', h, 'lambda', lambd_factor)
        #Optionally, can condense multiple comps into one, expand large comp into multiple, or assign radii when there are none
        condenser(newmorph, h.type, h.max_len, lambd_factor, h.rad_diff)
