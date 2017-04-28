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

#Usage: python shape_shifter.py --file 'filename.p' --type 'dc'
#Cannot have blank lines in p file
#type can be:
#  '0' just remove compartments of size 0
#  'radii' combine compartments of similar radius (specify 0 to only combine identical radii),
#          electrotonic length not to exceed max_len* lambda. (Specify 0 to use dc lambda)
# can also specify alternative values to default rm [4], cm [0.01], ri [2.5] units are SI
# can also specify frequency (--f) for ac lambda calculation, current default is 0.1 hz
# can also specify maximum electrotonic length (--max_len) current default is 0.1
# change info or debug (below) to get additional information while running program

#TO DO:
# fix subdivide_comp to work with branches at an angle
# read (and use) rm, cm and ri from p file
# write rm, cm and ri to p file if (a) not already specified and (b) given as input parameters
# test for whether absolute or relative coordinates.  If absolute - print error message (until length calculation can be updated)
# deal with blank lines (that begin with spaces) in .p file

#Saivardhan Mada
#Jul 20th, 2016
#Avrama Blackwell
#George Mason University
#Apr 28, 2017

import sys
import math
import argparse
import numpy as np

debug=0
info=1

class morph:
        def __init__(self,pfile):
                #dictionary to hold number of branches of each compartment
                self.children={}
                #dictionary to hold deleted voxels
                self.replace_parent = {}
                self.readfile(pfile)
                self.remove_zeros()
                self.ID_branches()
        
        def readfile(self,pfile):
                self.pfile=pfile
                linelist=open(pfile, 'r').readlines()
                self.linelist=[line for line in linelist if (line[0] !='*' and line[0] !='/' and line[0] != '\n')]
                numheader=len(linelist)-len(self.linelist)
                out_name=pfile.split('.p')[0]+'out.p'
                self.outfile = open(out_name, 'w')
                for line in linelist[0:numheader]:
                        self.outfile.write(line)

        def remove_zeros(self):
                newlines=[]
                for line_num,line in enumerate(self.linelist):
		        comp1 = line.split()
		        if(line_num <= len(self.linelist)-2):
			        comp2 = self.linelist[line_num+1].split()
                                #Test further, possibly modify "0" to floating point 0
			        if(comp1[2] == "0" and comp1[3] == "0" and comp1[4] == "0" and comp1[1] != 'none'):
                                        #if all coordinates are 0's and not soma compartment
				        self.replace_parent[comp1[0]] = comp1[1]# saves into dictionary to check future voxels
                                        if info:
				                print "!!!!!!!eliminate", comp1,"change child of", comp2, "to", comp1[1]
                                                print "replace_parent comp dictionary", self.replace_parent
                                else:
                                        xyzd=comp1[2]+"  "+comp1[3]+"  "+comp1[4]+"  "+comp1[5]
	                                if(self.replace_parent.has_key(comp1[1])):
                                        # if voxel refers to deleted voxel change the connection to the deleted voxels connection
                                                newlines.append(comp1[0]+" "+self.replace_parent[comp1[1]]+"   "+xyzd+" \n")
                                        else:
                                                newlines.append(comp1[0]+" "+comp1[1]+"   "+xyzd+" \n")
                        else:
                                #don't forget to add the last line
                                newlines.append(line)
                #print self.linelist[-1], newlines[-1]
                print "######## Zero size comp removal: orig num lines=", len(self.linelist), ", new num lines=", len(newlines)
                self.linelist=newlines

        def ID_branches(self):
                for line in self.linelist:
                        #count the number of children of each compartment
                        parent=line.split()[1]
                        if parent in self.children.keys():
                                self.children[parent]=self.children[parent]+1
                        else:
                                self.children[parent]=1
                #create list of compartments that shouldn't be deleted
                #1st (soma) compartment
                self.do_not_delete=[x.split()[0] for x in self.linelist if x.split()[1]=='none']
                #branch points
                self.do_not_delete=self.do_not_delete+[comp for comp in self.children.keys() if self.children[comp]>1]
                print "do not delete these parent branch compartments:", self.do_not_delete

def calc_lambda (type1, RM, RI, CM, F):
	tm = RM * CM #time constant
	dc_factor = math.sqrt(RM/(RI*4.0))*1000 #1000 is for unit conversion to microns 
	partial_ac= 2.0*math.pi*F*tm 
 	ac_factor = dc_factor * math.sqrt(2.0/(1.0+math.sqrt(1.0+math.pow(partial_ac,2.0))))
        return ac_factor 

def calc_length(comp):
        x=float(comp[2])
        y=float(comp[3])
        z=float(comp[4])
        length = math.sqrt(x*x+y*y+z*z)
        return length

def calc_electrotonic_len(comp,factor):
        lamb = factor * math.sqrt(float(comp[5])) #calculates ac or dc lambda from compartment radius
        #electronic length of compartment
        length=calc_length(comp)
        if debug:
                print('    calc_len comp %s lambda=%.1f len=%.3f L=%.5f' % (comp[0],lamb, length, length/lamb))
        return length/lamb

def write_line(out_file,replace_parent,comp,numcomps):
	if(replace_parent.has_key(comp[1])):# if voxel refers to deleted voxel change the connection to the deleted voxels connection
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
        #total distance in x, y and z.  Possibly this could be used for ac/dc type, with l = total length 
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
        print xyz
        for i in range(segs):
                #initialize x,y,z to 0 0 0
                newcomp.append(['0','0','0',comp[5]])
                myname.append(comp[0]+'_'+str(i))
                if i==0:
                   parent.append(comp[1])
                else:
                   parent.append(myname[i-1])
                if np.sum(np.abs(xyz))==np.max(np.abs(xyz)):
                        maxval=max(np.abs(xyz))
                        if maxval==-min(xyz):
                                maxindex=xyz.index(-maxval)
                        else:
                                maxindex=xyz.index(maxval)
                        comp_index=maxindex+2
                        newcomp[i][maxindex]=str(np.round(float(comp[comp_index])/segs,5))
                else:
                        print "new approach needed"
        for i in range(segs):
                print myname[i],parent[i], newcomp[i]
        return myname,parent,newcomp
        

def condenser(m, type1, max_len, lambda_factor, rad_diff):
        num_comps=0
        ######## type = 1 removes 0 length compartments, this is done in class, so just write file
	if(type1 == "0"):
	        for line in m.linelist:
                        m.outfile.write(line)
                        num_comps=num_comps+1
        ####### type = "expand" takes long compartments and subdivides into multiple "segments"
        if (type1 == "expand"):
 	        for line in m.linelist:
		        comp1 = line.split()
                        L_comp1=calc_electrotonic_len(comp1,lambda_factor)
                        print "max_len", max_len, "L", L_comp1
                        if L_comp1>max_len:
                                segs=int(math.ceil(L_comp1/max_len))
                                myname,par,xyzdiam=subdivide_comp(comp1,segs)
                                m.replace_parent[comp1[0]]=myname[-1]
                                if info:
                                        print "!!!!Branch to long:", comp1[0],L_comp1,segs
                                        print m.replace_parent, "attach distal branches to", m.replace_parent[comp1[0]], "instead of", comp1[0]
                                for n,p,xyzd in zip(myname,par,xyzdiam):
                                        newcomp=[n, p, xyzd[0], xyzd[1], xyzd[2], xyzd[3]]
                                        num_comps=write_line(m.outfile,m.replace_parent,newcomp,num_comps)
                        else:
                                num_comps=write_line(m.outfile,m.replace_parent,comp1,num_comps)
                                if info:
                                        print "OK", comp1[0],L_comp1
        ######## type = radii condenses branches with similar radius and combined electronic length < 0.1 lambda
	if (type1 == "radii"):    #if rad_diff = 0, only condenses branches with same radius
	        condense = []
	        Ltot = 0
	        surface_tot = 0
	        for line in m.linelist:
		        comp1 = line.split()
                        if debug:
                                print '**** begin', comp1[0]
 		        line_num = m.linelist.index(line)
                        if line_num <= (len(m.linelist)-2):
			        comp2 = m.linelist[line_num+1].split()
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
                                #if this is last line of file, comp2=comp1, thus comp2[1] != comp1[0], don't enter this loop
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
                                                print '#######condense', condense, 'stop before', comp2[0]
                                                print 'new: x y z',x,y,z 
					newcomp=[condense[-1][0], condense[0][1], x, y, z, dia]
                                        num_comps=write_line(m.outfile,m.replace_parent,newcomp,num_comps)
					condense = []
					Ltot = 0
					surface_tot = 0
				else:
                                        #cannot condense this comp and nothing in condense set.  Just print out line
                                        num_comps=write_line(m.outfile,m.replace_parent,comp1,num_comps)
                                        if info:
                                                print 'not condensing',len(condense),Ltot, comp1[0], "rad", rad_diff,delta_rad
        print "finished,", num_comps, "output compartments"

if __name__ == '__main__':
        #set default parameters
        rm = 4.00 #ohms-m^2 
        ri = 2.50 #ohms-m
        cm = 0.010 #F/m^2
        max_len = 0.1 #electrotonic length
        f = .1 #Hz
        rad_diff=0.1
        type1 = "ac"

	#array of variables 
	variables = [type1, rm, ri, cm, f, max_len, rad_diff]
        
	#sets up arg parser for all the variables 
	parser = argparse.ArgumentParser()
        parser.add_argument('--file')
	parser.add_argument('--type', choices={'0', 'radii','expand'}, default='radii')
	parser.add_argument('--rad_diff', default=rad_diff)
	parser.add_argument('--rm', default=rm)
	parser.add_argument('--ri', default=ri)
	parser.add_argument('--cm', default=cm)
	parser.add_argument('--f', default=f)
	parser.add_argument('--max_len', default=max_len)
        
	#takes values from arg parser and adds them to the variable array 
	h = parser.parse_args()
	variables[0] = h.type
	variables[1] = float(h.rm)
	variables[2] = float(h.ri)
	variables[3] = float(h.cm)
	variables[4] = float(h.f)
	variables[5] = float(h.max_len)
	variables[6] = float(h.rad_diff)

	#reads p file
        newmorph=morph(h.file)

        #calculate lambda
        lambd_factor=calc_lambda(variables[0], variables[1], variables[2], variables[3], variables[4])
        #may not need variables 1-4 once radii version fixed (or eliminated)
	condenser(newmorph, variables[0], variables[5], lambd_factor, variables[6])


