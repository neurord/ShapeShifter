'''<Condenses .p files of neuron traces to maintain specific lambda (ac or dc) and removes 0 length compartments.>
    Copyright (C) <2016>  <Saivardhan Mada>

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

#TO DO:
#0. verify that ac or dc works correctly - may need to check for relative vs absolute coordinates
#2. deal with blank lines
#1. correctly deal with parents when compartment is deleted using 'radii'

#Saivardhan Mada
#Jul 20th, 2016
#runs using python3
import sys
import math
import argparse
import numpy as np

#dictionary to hold deleted voxels
deleted = {}
needtowrite = []

def calc_lambda (type1, RM, RI, CM, F):
	tm = RM * CM #time constant
	dc_factor = math.sqrt(RM/(RI*4.0))*1000 #1000 is for unit conversion to microns 
	partial_ac= 2.0*math.pi*F*tm 
 	ac_factor = dc_factor * math.sqrt(2.0/(1.0+math.sqrt(1.0+math.pow(partial_ac,2.0))))
        if type1=='ac':
                return ac_factor
        elif type1=='dc':
                return dc_factor
        else:
                return ac_factor #apparently used by type=radius

def write_line(out_file,deleted,comp):
	if(deleted.has_key(comp[1])):# if voxel refers to deleted voxel change the connection to the deleted voxels connection
                newline=comp[0]+" "+deleted[comp[1]]+"   "+comp[2]+"  "+comp[3]+"  "+comp[4]+"  "+comp[5]+" \n"
		print('old:',comp,'deleted:',deleted[comp[1]], 'new:',newline)
		out_file.write(newline)
	else:
                line=comp[0]+" "+comp[1]+"   "+comp[2]+"  "+comp[3]+"  "+comp[4]+"  "+comp[5]+" \n"
		out_file.write(line)

def condenser(input_file, output_file, type1, RM, RI, CM, F, max_len, lambda_factor):
	comp2 = []#blank variable to hold previous compartment
	tobecondensed = []
	future = {}
	lamb_tot = 0
	surface_tot = 0
	for line in input_file:
                if line.startswith('*') or line.startswith('//') or line.startswith(' '):
                        #if comment or parameter, print the line, read in next line
			output_file.write(line)
                else:
                        #otherwise, assign 1st line as comp1 and next line as comp2
		        comp1 = line.split()
		        line_num = input_file.index(line)
                        #print 'text',comp1, 'num',line_num
		        if(line_num <= len(input_file)-2):
			        comp2 = input_file[line_num+1].split()
                        #
                        #compare the two compartments to determine if they can be merged
			if(type1 == "0"): # type which removes 0 length compartments
				if(comp1[2] == "0" and comp1[3] == "0" and comp1[4] == "0" and comp1[1] != 'none'):#if all coordinates are 0's and not soma compartment
				        print(comp2,comp1)
					deleted[comp1[0]] = comp1[1]# saves into dictionary to check future voxels
                                        print (deleted)
                                else:
                                        write_line(output_file,deleted, comp1)
			elif(type1 == "ac" or type1 == "dc"): #condenses branches with same radius and minimum electronic length of .01 lambda [AC or DC]
				lamb = lambda_factor * math.sqrt(float(comp2[5])) #calculates ac or dc lambda from compartment radius
				# calculate the lengths and check for total electronic length (len/lambda) being less then 0.1 (max_len)
				if(comp1[5] == comp2[5] and comp2[1]==comp1[0]): # comp2 is attached to comp1 and radii are the same (change this to radii within 10% of each other)
					len_comp1 = (math.pow(float(comp1[2]),2.0)+math.pow(float(comp1[3]),2.0)+math.pow(float(comp1[4]),2.0))/lamb #electronic length of compartment 1
					len_comp2 = (math.pow(float(comp2[2]),2.0)+math.pow(float(comp2[3]),2.0)+math.pow(float(comp2[4]),2.0))/lamb #electronic length of compartment 2
					#print(len_comp1+len_comp2)
					if((len_comp2+len_comp1) < max_len): #if it is less then the designated electrotonic length combine the two compartments
						if(comp1[1] in deleted):
							deleted[comp1[0]] = deleted[comp1[1]]
						else: #delete comp1, keep track of deletion, update comp2 x,y,z if relative, no need if absolute coordinates?
							deleted[comp1[0]] = comp1[1]
                                                #not being written in case it get's combined in next line
						#output_file.write(comp2[0]+" "+deleted[comp1[0]]+" "+comp2[2]+" "+comp2[3]+" "+comp2[4]+" "+comp2[5]+"\n")
                                        else:
                                                write_line(output_file,deleted,comp1)
                                else:
                                        write_line(output_file,deleted,comp1)
			elif(type1 == "radii"): #condenses branches with same radius and minimum electronic length of .01 lambda [AC], calculates new coordinates (i.e., if relative coordinates)
				ten = .2*float(comp1[5])
				#basically change this to 10% and calculate the lens again and check for electronic length being less then 
				if(abs(float(comp1[5]) - float(comp2[5])) <= ten):#if the radius is the same
					ac1 = lambda_factor * math.sqrt(float(comp1[5]))#calculates lambda for ac
					ac2 = lambda_factor * math.sqrt(float(comp2[5]))#calculates lambda for ac
					len_comp1 = (math.pow(float(comp1[2]),2.0)+math.pow(float(comp1[3]),2.0)+math.pow(float(comp1[4]),2.0))/ac1 #electronic length of compartment 1
					len_comp2 = (math.pow(float(comp2[2]),2.0)+math.pow(float(comp2[3]),2.0)+math.pow(float(comp2[4]),2.0))/ac2 #electronic length of compartment 2
					#print(len_comp1+len_comp2)
					#comp1[5] = str((float(comp1[5]) + float(comp2[5]))/2.0)
					if((len_comp2+len_comp1)+lamb_tot < max_len):#if it is less then the designated lambada value condense the branch
						#output_file.write(comp2[0]+" "+deleted[comp2[1]]+" "+comp2[2]+" "+comp2[3]+" "+comp2[4]+" "+comp2[5]+"\n")
						#print("hey")
						lamb_tot = lamb_tot + (len_comp2+len_comp1)
						surface_tot = surface_tot + (2* math.pi * float(comp1[5]) * (math.pow(float(comp1[2]),2.0)+math.pow(float(comp1[3]),2.0)+math.pow(float(comp1[4]),2.0)))
						if(comp1[1] in deleted):
							deleted[comp1[0]] = deleted[comp1[1]]
							tobecondensed.append(comp1)
						else:
							deleted[comp1[0]] = comp1[1]
								
						#output_file.write(comp2[0]+" "+comp1[1]+" "+comp2[2]+" "+comp2[3]+" "+comp2[4]+" "+comp2[5]+"\n")
					elif(comp1[1] in deleted and lamb_tot > 0):
						l = math.pow((surface_tot * math.pow((2*RM/RI), 1/2) / (2*math.pi*lamb_tot)), (2/3))
						r = lamb_tot*math.pow(l, 1/2)/ math.pow(2*RM/RI, .5)
						x = 0.0
						y = 0.0
						z = 0.0
						for comp in tobecondensed:
							x = x + float(comp[2])
							y = y + float(comp[3])
							z = z + float(comp[4])
						#print(x, y, z)
						if((x-float(comp1[2]) > 0)):
							theta = np.arctan((y-float(comp1[3]))/(x-float(comp1[2])))
						else:
							theta = 0
						if(math.pow((z-float(comp1[4])), 2) > 0):
							phi = np.arctan(math.pow(((x-float(comp1[2]))**2 + (y-float(comp1[3]))**2), 1/2)/ math.pow((z-float(comp1[4])), 2))
						else:
							phi = 0
						x = l*math.cos(theta)*math.sin(phi) + float(comp1[2])
						y = l*math.sin(theta)*math.sin(phi) + float(comp1[3])
						z = l*math.cos(phi) + float(comp1[4])
						output_file.write(comp1[0]+" "+deleted[comp1[1]]+" "+str(x)+" "+str(y)+" "+str(z)+" "+str(r)+"\n")
						tobecondensed = []
						lamb_tot = 0
						surface_tot = 0
					else:
						lamb_tot = 0
						surface_tot = 0
						output_file.write(line)
                                else:
                                        write_line(output_file,deleted,comp1)


def main():
        #set default parameters
        rm = 4.00 #ohms-m^2 
        ri = 2.50 #ohms-m
        cm = 0.010 #F/m^2
        max_len = 0.1 #electrotonic length
        f = .1 #Hz
        type1 = "ac"

	#array of variables 
	variables = [type1, rm, ri, cm, f, max_len]
        
	#sets up arg parser for all the variables 
	parser = argparse.ArgumentParser()
        parser.add_argument('--file')
	parser.add_argument('--type', choices={'ac', 'dc', '0'}, default='radii')
	parser.add_argument('--rm', default=rm)
	parser.add_argument('--ri', default=ri)
	parser.add_argument('--cm', default=cm)
	parser.add_argument('--f', default=f)
	parser.add_argument('--max_len', default=max_len)
        
	#takes values from arg parser and adds them to the variable array 
	h = parser.parse_args()
        in_name=h.file
        out_name=h.file.split('.p')[0]+'out.p'
	variables[0] = h.type
	variables[1] = h.rm
	variables[2] = h.ri
	variables[3] = h.cm
	variables[4] = h.f
	variables[5] = h.max_len
        
	#reads p file
	input_file = open(in_name, 'r').readlines()
	output_file = open(out_name, 'w')

        #calculate lambda
        if h.type=='ac' or h.type=='dc' or h.type=='radii':
                lambd_factor=calc_lambda(variables[0:5])
        else:
                lambd_factor=0
        #may not need variables 1-4 once radii version fixed (or eliminated)
	condenser(input_file, output_file, variables[0], variables[1], variables[2], variables[3], variables[4], variables[5], lambd_factor)

main()

