'''<Extracts data values from .swc files>
    Copyright (C) <2019>  <Jonathan Reed>

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

#Usage:  python graph_extract.py --path 'name_of_path' --choice 'type of choice'
#        where morph_list is separate file with list of names of .swc files
#Output: morphology_extract.py

#George Mason University
#Jonathan Reed
#Sep 3, 2019

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import argparse
import os
import glob
import re
'''
def graph(temp_list, param):
    for num, param in enumerate(parameter):
        plt.figure(num)
        plt.plot(temp_list[header[param]], temp_list['RADIUS'])
    #plt.show()
    #send in x,y, possibly z for either single or all files within or between certain archives
'''    
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   #may change, but for now is glob-command to look for 
parser.add_argument("--choice", type = str, choices = {'optimize','graph'})
parser.add_argument("--parameters", type = str, nargs = '+', choices = {'PARENT_RAD','NODE_DEGREE','NODE_ORDER','PATH_DIS','PATH_TO_END','NUM_ENDS','BRANCH_LEN','HS'}) 

args = parser.parse_args()                       #would need to add additional params here manually
path = args.path                                
choice = args.choice
parameters = args.parameters

root, dirs, files = list(os.walk(path))[0]  #dirs are subdirectories with archive morphology files
                                            #make sure all extracted files have same parameters before graphing
archive_dict = {}; header = {}
for d1 in dirs:
    fullpath = path + d1 + '/*CNG_extract.txt'
    data_list = []; apical_list = []; basal_list = []
    print('Working on Directory : ', str(d1))
    for fcount,fname in enumerate(glob.glob(fullpath)):
        with open(fname) as f:
            for line in f:
                if line.strip():                                    #removes any empty lines
                    if '*C' in line and not header:
                        if 'XYZ' in line:                           
                            line = re.sub('XYZ','X; Y; Z', line)    #find and fix header for each file correctly
                        line = line.split('; ')                     
                        for num, val in enumerate(line):            
                            if '\n' in val:
                                val = re.sub('\n','',val)
                            if val == '*CHILD':
                                header['CHILD'] = num
                            else:
                                header[val] = num                          #once this is complete dictionary is created
                    elif line[0] != '*' and line[0] != ' ' and line[0] != '/n':                #for rest of data line by line
                        temp_line = line.split()                           #for some reason empty first list for first file... 
                        if len(temp_line) == len(header):                  #seamless way to convert certain string values to numerical ones
                            for point, val in enumerate(temp_line):        #Child, Type, Parent use DESCRETE numerical labelling
                                if point != header['CHILD'] and point != header['TYPE'] and point != header['PARENT']:
                                    temp_line[point] = float(temp_line[point])
                            #temp_line.append(d1)                           #if wanted could append name of archive to each data line
                        data_list.append(temp_line)                    #or push the temp_list to the archive dictionary at this point...

    for line in data_list:
        if line[header['TYPE']] == '4':                     #splits by node type (Apical or Basal)
            apical_list.append(line)
        elif line[header['TYPE']] == '3':
            basal_list.append(line)
            
    archive_dict[d1] = {}                                   #sets archive name (in dirs) as nested dictionary within complete archive dictionary
    basal_list = zip(*basal_list)
    archive_dict[d1]['Basal'] = basal_list                  #keeps archive --> separate oranization of Apical and Basal
    if apical_list:
        apical_list = zip(*apical_list)                        
        archive_dict[d1]['Apical'] = apical_list
        
for archive in archive_dict.keys():                         #for every archive in archive dictionary
    for param in parameters:         #for eachparameter set as argument
        plt.figure(num)
        #for cell_type in archive_dict[archive]:
            #plt.plot(archive_dict[archive][cell_type][header[param]], archive_dict[archive][cell_type][header['RADIUS']], 'o', label = cell_type)
        #maybe plot the different archives together as Apical and Basal 
        if 'Apical' in archive_dict[archive].keys():
            plt.plot(archive_dict[archive]['Apical'][header[param]], archive_dict[archive]['Apical'][header['RADIUS']], 'o', label = 'Apical')
        plt.plot(archive_dict[archive]['Basal'][header[param]], archive_dict[archive]['Basal'][header['RADIUS']], 'o', label = 'Basal')
        title = str(archive) + ' ' + str(param) #str(cell_type) 
        plt.xlabel(str(param))
        plt.ylabel('Radius')
        plt.title(title)
        plt.legend()
        plt.show()
        #plt.savefig(title + '.png')
     

'''
if choice == 'optimize':
    print(choice)
if choice == 'graph':
    print(choice)
'''
'''
    
        if node == 'Basal':
            plt.plot(combined_list[header[param]], combined_list[header['RADIUS']],'o', label = 'basal')  
        plt.xlabel(param)                                                                               
        plt.ylabel('Radius')
        plt.legend()
    plt.show()
'''
'''
    if len(Radius_a):
        #plot both apical and basal data
        plt.plot(HoSt_a, Radius_a, 'o', label = 'apical dendrite')
        plt.plot(HoSt_b, Radius_b, 'o', label = 'basal dendrite')
    else:
        #plot only basal data as no apical data exists in particular archive
        plt.plot(HoSt_b, Radius_b, 'o', label = 'basal dendrite')
    #Axes3D.plot_surface()look for seqeuntial color map

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if len(Radius_a):
        ax.scatter(End_Distance_a, Radius_a, HoSt_a, c='r', marker = 'o')
        ax.scatter(End_Distance_b, Radius_b, HoSt_b, c='b', marker = 'o')
    #else:
    ax.scatter(Parent_Rad_b, Radius_b, End_Distance_b, c='g', marker = 'o')
    ax.set_xlabel('Parent_Rad')
    ax.set_ylabel('Radius')
    ax.set_zlabel('End_Distance')
    plt.show()

    #either send data to separate for loop(plot entire data set) or function 
    #plt.ylabel('Radius')
    #plt.xlabel('HS')
    
    #plt.legend()
    #plt.savefig(morph_file + 'HS' + ".png", format="PNG")
    #print(morph_file + 'HS' + ".png")
    #plt.show()
  
    #3D plotting

    #X,Y,Z = [][][]

    X = 
    Y = 
    Z = 
    ax.scatter(X, Y, Z, marker = 'o')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    plt.plot(Degree, Radius, 'o')
    plt.ylabel('Radius')
    plt.xlabel('Degree')

    print(len(Radius))
    print(len(Distance))
    print(len(Degree))
    print(len(Order))
    print(len(Parent_Rad))
    print(len(End_Distance))
    print(len(Ends))
    print(len(HoSt))

    plt.savefig(morph_file + 'degree' + ".png", format="PNG")
    print(morph_file + 'degree' + ".png")
    plt.show()
    
    
    if Ends == Degree:
        print('y')
    else:
        print('n')
    
'''
