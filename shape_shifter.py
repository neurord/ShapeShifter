'''<Condenses .p files of neuron traces to maintain specific lambda (ac or dc) and removes 0 length compartments.>
    Copyright (C) <2016>  <Saivardhan Mada>
    Copyright (C) <2017>  <Avrama Blackwell>

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

#Usage: python shape_shifter.py --file 'filename.p' --type 'radii'
#type can be:
#  '0' just remove compartments of size 0
#  'condense' combine compartments of similar radius (specify 0 to only combine identical radii),
#          electrotonic length not to exceed max_len* lambda.
#  'expand' to change single, long compartment (e.g. a Neuron software segment) into multiple smaller compartments
#           electrotonic length of subdivided compartments do not exceed max_len* lambda
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

#Saivardhan Mada
#Jul 20th, 2016
#Avrama Blackwell
#George Mason University
#Apr 28, 2017

from __future__ import print_function, division
import sys
import math
import argparse
import numpy as np

debug=0   #default is 0
info=0    #default is 1

class morph:
        def __init__(self,pfile):
                #dictionary to hold number of branches of each compartment
                self.children={}
                #dictionary to hold new parent compartments for deleted voxels
                self.replace_parent = {}
                self.readfile(pfile)
                self.remove_zeros()
                self.ID_branches()
                #is this where to put 3-pt --> 1-pt soma code? and save as .p file?
                #could read in data of first couple lines (need to change soma anyway)
                
        def readfile(self,pfile):
                '''
                Reads from pfile as STRING values --> will convert later + include swc file specifics
                '''
                self.pfile=pfile
                self.linelist = []
                lines=open(pfile, 'r').readlines()
                #skip parameters (*), comments (/) and blank lines (\n or \r)
                self.lines=[line for line in lines if (line[0] !='*' and line[0] !='/' and line[0] != '\n' and line[0] != '\r')]
                out_name=pfile.split('.p')[0]+'out.p'
                self.outfile = open(out_name, 'w')
                '''
                Reads in line by line than Converts to list of list
                '''
                for num, line in enumerate(self.lines):
                        self.linelist.append(line)
                        self.linelist[num] = self.linelist[num].split()
                        self.outfile.write(line)
                '''
                Conversion from strings --> floats takes place here
                '''
                #.p file child and parent are STRINGS --> keep as strings for now
                #convert rest str. values into float (xyz, dia and dis)
                for num,line in enumerate(self.linelist):
                        for id, val in enumerate(line):
                                if id > 1:
                                        line[id] = float(line[id])
                #print(self.linelist)
                #.swc files have child/parent saved as int instead of names, maybe read as FLOATS??
        
        def remove_zeros(self):
                #set globals here for dictionary-like retrieval of split data in linelist
                child = 0; parent = 1; x = 2; y = 3; z = 4; dia = 5; dis = 6
                newlines=[]; parent_dict = {}; parent_list = [] #may want to convert parent dict to list
                print(len(self.linelist))
                for line in self.linelist:
                        print(line)
                for num,line in enumerate(self.linelist):
                        #maybe we could calculate distance HERE to then help with REMOVING ZEROS later
                        distance = np.sqrt((line[x])**2 + (line[y])**2 + (line[z])**2)
                        line.append(distance)
                '''
                Remove Zeroes and Replace oldparent id with newparent id
                '''
                for num,line in enumerate(self.linelist):
                        #removes (does not update new-list) any compartments with 0 distance or diameter OTHER THAN soma
                        if((line[dia] == 0  or line[dis] == 0) and line[parent] != 'none'):
                                #relabel remaining children to previous parent (or compartment that is deleted)
                                #i.e. if                    4,3 0 0 0 <--- needs to be removed
                                # then                      5,4 x y z
                                # should be changed to      5,3 x y z <--- new parent connection to remaining children  
                                parent_dict[line[child]] = line[parent]
                                print ("######## Zero size comp removal: orig num lines=", line[child], ", new num lines=", line[parent])
                        else:
                                newlines.append(line)
                if parent_dict: #checks if parent_dict has any items, THEN does recursion
                        print('Zero point compartments removed; remaining comparments require parent-update for reconnection')
                        print(parent_dict)
                        #maybe put a 'limit' if too many zero-point compartments exists in a single morphology
                        parent_len = len(parent_dict.items())
                        
                        if parent_len < 100: #just a placeholder for now; i.e. 100 zero-point compartments
                                for num, line in enumerate(self.linelist):
                                        #this should only iterate for each line in linelist
                                        if line[parent] != 'none':
                                                counter = 0 #in this organization the counter comes back to 0 for each line
                                                line = self.replace_par(parent_dict, counter, line, parent, child)
                                print(len(newlines))
                                self.linelist = newlines
                        else:
                                print('Many zero-point compartments --> check morphology file')

        #new function to calculate diameter from distance

        def ID_branches(self):
                for line in self.lines:
                        #count the number of children of each compartment
                        parent=line.split()[1]
                        if parent in self.children.keys():
                                self.children[parent]=self.children[parent]+1
                        else:
                                self.children[parent]=1
                                #create list of compartments that shouldn't be deleted
                                #1st (soma) compartment
                self.do_not_delete=[x.split()[0] for x in self.lines if x.split()[1]=='none']
                #branch points
                self.do_not_delete=self.do_not_delete+[comp for comp in self.children.keys() if self.children[comp]>1]
                print ("do not delete these parent branch compartments:", self.do_not_delete)

        #replace parent values of compartments whose original parents were zero-point compartments and removed 
        def replace_par(self,parent_dict, counter, line, parent, child):
                #so this changes each line, not the parent dictionary itself
                print('Original line', line)
                if line[parent] in parent_dict:
                        print('Line child', line[child], 'Line parent', line[parent])
                        counter = counter + 1
                        print(counter)
                        line[parent] = parent_dict[line[parent]]
                        self.replace_par(parent_dict, counter, line, parent, child)
                        #so will var line be updated and rest of arguments remain the same as recursion ocurrs?
                else:
                        if counter>0:
                                print('fixed line', counter,line)
                        return line
        
def calc_lambda (type1, RM, RI, CM, F):
	tm = RM * CM #time constant
	dc_factor = math.sqrt(RM/(RI*4.0))*1000 #1000 is for unit conversion to microns 
	partial_ac= 2.0*math.pi*F*tm 
        ac_factor = dc_factor * math.sqrt(2.0/(1.0+math.sqrt(1.0+math.pow(partial_ac,2.0))))
        return ac_factor 

def calc_length(comp): #is this the same as the distance calculated above?
        x=float(comp[2])
        y=float(comp[3])
        z=float(comp[4])
        length = math.sqrt(x*x+y*y+z*z)
        return length

def calc_electrotonic_len(comp,factor):
        lamb = factor * math.sqrt(float(comp[5])) #calculates ac or dc lambda from compartment radius
        length=calc_length(comp) #electronic length of compartment
        if debug:
                print('    calc_len comp %s lambda=%.1f len=%.3f L=%.5f' % (comp[0],lamb, length, length/lamb))
        return length/lamb

def write_line(out_file,replace_parent,comp,numcomps): # do i need this as well?
	if(replace_parent.has_key(comp[1])): # if parent refers to deleted voxel, replace the parent with the deleted voxels parent
                newline=comp[0]+" "+replace_parent[comp[1]]+"   "+comp[2]+"  "+comp[3]+"  "+comp[4]+"  "+comp[5]+" \n"
		print('!!!!old:',comp,'deleted:',comp[1], 'new:',newline)
	else:
                newline=comp[0]+" "+comp[1]+"   "+comp[2]+"  "+comp[3]+"  "+comp[4]+"  "+comp[5]+" \n"
	out_file.write(newline)
        numcomps=numcomps+1
        return numcomps

def calc_newcomp(condense,surface_tot,Ltot,lamb_factor):
        dia=math.pow(surface_tot/(Ltot*np.pi*lamb_factor),(2/3.))
        l=surface_tot/(np.pi*dia)
        if debug:
                print('TSA: dia %.5f, l=%.2f tsa=%.1f Ltot=%.5f' % (dia, l, np.pi*dia*l,l/(lamb_factor*math.sqrt(dia))))
	x = 0.0
	y = 0.0
	z = 0.0
        #total distance in x, y and z.  
	for comp in condense:
		x = x + float(comp[2])
		y = y + float(comp[3])
		z = z + float(comp[4])
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
        return str(x),str(y),str(z),str(np.round(dia,3))

def subdivide_comp(comp,segs):
        myname=[]
        parent=[]
        newcomp=[]
        xyz=[float(x) for x in comp[2:5]]
        length=math.sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])
        newlen=0
        print ("subdivide", comp[0], "x,y,z,len", xyz, length, "into", segs, "segments")
        new_xyz=np.zeros(3)
        for j in range(3):
                new_xyz[j]=np.round(xyz[j]/segs,5)
                seg_length=math.sqrt(new_xyz[0]*new_xyz[0] + new_xyz[1]*new_xyz[1] + new_xyz[2]*new_xyz[2])
        for i in range(segs):
                newcomp.append(['0','0','0',comp[5]])
                myname.append(comp[0]+'_'+str(i))
                if i==0:
                   parent.append(comp[1])
                else:
                   parent.append(myname[i-1])
                for j in range(3):
                        newcomp[i][j]=str(new_xyz[j])
                newlen=newlen+seg_length
        if info:
                for i in range(segs):
                        print ("new seg", i, myname[i],parent[i], newcomp[i])
                print ("total length", newlen, "seg", seg_length )
        return myname,parent,newcomp
        
def condenser(m, type1, max_len, lambda_factor, rad_diff):
        num_comps=0
        ######## type = 1 removes 0 length compartments, this is done in the class, so just write file
	if(type1 == "0"):
	        for line in m.lines:
                        m.outfile.write(line)
                        num_comps=num_comps+1
        ####### type = "expand" takes long compartments and subdivides into multiple "segments"
        if (type1 == "expand"):
 	        for line in m.lines:
		        comp1 = line.split()
                        L_comp1=calc_electrotonic_len(comp1,lambda_factor)
                        print ("max_len", max_len, "L", L_comp1)
                        if L_comp1>max_len:
                                segs=int(math.ceil(L_comp1/max_len))
                                myname,par,xyzdiam=subdivide_comp(comp1,segs)
                                m.replace_parent[comp1[0]]=myname[-1]
                                if info:
                                        print (m.replace_parent, "attach distal branches to", m.replace_parent[comp1[0]], "instead of", comp1[0])
                                for n,p,xyzd in zip(myname,par,xyzdiam):
                                        newcomp=[n, p, xyzd[0], xyzd[1], xyzd[2], xyzd[3]]
                                        num_comps=write_line(m.outfile,m.replace_parent,newcomp,num_comps)
                        else:
                                num_comps=write_line(m.outfile,m.replace_parent,comp1,num_comps)
                                if info:
                                        print ("OK", comp1[0])
        ######## type = condense condenses branches with similar radius and combined electronic length < 0.1 lambda
	if (type1 == "condense"):    #if rad_diff = 0, only condenses branches with same radius
	        condense = []
	        Ltot = 0
	        surface_tot = 0
	        for line in m.lines:
		        comp1 = line.split()
                        if debug:
                                print ('**** begin', comp1[0])
 		        line_num = m.lines.index(line)
                        if line_num <= (len(m.lines)-2):
			        comp2 = m.lines[line_num+1].split()
                                len_comp1=calc_electrotonic_len(comp1,lambda_factor)
                                len_comp2=calc_electrotonic_len(comp2,lambda_factor)
                        if len(condense):
                                delta_rad=abs(float(condense[0][5]) - float(comp2[5]))/float(condense[0][5])
                                tot_len=len_comp2+Ltot
                        else:
                                delta_rad=abs(float(comp1[5]) - float(comp2[5]))/float(comp1[5])
                                tot_len=len_comp1 + len_comp2
			if(delta_rad <= rad_diff and comp2[1]==comp1[0] and tot_len < max_len and comp1[0] not in m.do_not_delete):
                                #if the radii are almost the same and comp2 is attached to comp1 and not in do_not_delete list
                                #if this is last line of file, comp2=comp1, thus comp2[1] != comp1[0], skip down to condensing or writing file
                                if len(condense)==0:
                                        #if this is 1st compartment of a set, add it
                                        condense.append(comp1)
					Ltot = len_comp1
                                        #surface_area=pi*diam*len, comp[5] stores diameter, not radius
					surface_tot = math.pi * float(comp1[5]) * calc_length(comp1)
                                #always add the 2nd compartment to set of comartments to be condensed
                                condense.append(comp2)
                                Ltot=Ltot+len_comp2
                                surface_tot=surface_tot+(math.pi * float(comp2[5]) * calc_length(comp2))
				if(comp1[1] in m.replace_parent): #m.replace_parent.has_key(comp[1]))
					m.replace_parent[comp1[0]] = m.replace_parent[comp1[1]]
				else:
					m.replace_parent[comp1[0]] = comp1[1]
			else:
                                if len(condense):
                                        #cannot add any additional compartments.  Condense the set
                                        x,y,z,dia=calc_newcomp(condense,surface_tot,Ltot,lambda_factor)
                                        if info:
                                                print ('#######condense', condense, 'stop before', comp2[0])
                                                print ('new: x y z',x,y,z )
					newcomp=[condense[-1][0], condense[0][1], x, y, z, dia]
                                        num_comps=write_line(m.outfile,m.replace_parent,newcomp,num_comps)
					condense = []
					Ltot = 0
					surface_tot = 0
				else:
                                        #cannot condense this comp and nothing in condense set.  Just print out line
                                        num_comps=write_line(m.outfile,m.replace_parent,comp1,num_comps)
                                        if info:
                                                print ('not condensing',len(condense),Ltot, comp1[0], "rad", rad_diff,delta_rad)
        print ("finished,", num_comps, "output compartments")
        #if type1==('radii'):
        '''
        def diameter_create(self.linelist):
                #I need to somehow put a check to see if all subsequent values of diameter are equal with each other (or zero)
                for num,line in enumerate(self.linelist):
                         for id, val in enumerate(line):
                                if id > 1:
                                        line[id] = float(line[id])
                d1 = 0.001 * L + 0.87
                d2 = dp * np.exp(-l * 0.08)
                d = np.max(d1,d2)

                #where:
                #L = dendritic length of branch
                #dp = diameter of parent node
                #l = distance to parent node
        '''
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
	parser.add_argument('--type', choices={'0', 'condense','expand'}, default='condense')
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


'''extra code'''
'''
     if any()value in parent_dict.values() == key in parent_dict.keys():
                        value = parent_dict[key]
                        
                for key, value in parent_dict.items():
                        if key == value:
                                #replace value with key-ed value
                                value = key[value]
                        #how to detect if one 'self' would equal another's parent???
                #ideally by changing the parent_dict
                #I can use the new dictionary and change the values in my 'old' code below
                print(parent_dict)
'''
'''
                #converts parent dict to list --> still trying to just replace parent from dictionary in list BEFORE changing lines
                for key, value in parent_dict.items():
                        temp = [key,value]
                        parent_list.append(temp) #as list of lists
                
                #it seems like I would need to do recursion here as well... huh
                #change dictionary values here
                if #there is a match between separate keys and values
                for i, item in enumerate(parent_list):
                for j, val in enumerate(item):
                print(j)
                
                for num, line in enumerate(self.linelist):
                        if line[parent] != 'none':
                                #change line values here...
                #what if we combined list-format of dictionary and replacing dictionary values instead?
'''
