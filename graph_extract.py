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

#Usage:  python graph_extract.py --path /name/of/path --parameters PATH_TO_END BRANCH_LEN
#        where path is name of path to directory containing archive directories
#        where parameters requried as capitalized in code
#Output: morphology_extract.py

#George Mason University
#Jonathan Reed
#Sep 3, 2019

import numpy as np
#import scipy
from scipy.stats import pearsonr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D

import argparse
import os
import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   #may change, but for now is glob-command to look for 
parser.add_argument("--choice", type = str, choices = {'optimize','graph'})
#parser.add_argument("--parameters", type = str, nargs = '+', choices = {'PARENT_RAD','NODE_DEGREE','NODE_ORDER','PATH_DIS','PATH_TO_END','NUM_ENDS','BRANCH_LEN','HS'}) 

args = parser.parse_args()                       #would need to add additional params here manually
path = args.path                                
choice = args.choice
#parameters = args.parameters

root, dirs, files = list(os.walk(path))[0]  #dirs are subdirectories with archive morphology files
                                            #make sure all extracted files have same parameters before graphing
archive_dict = {}; header = {}; params = {} #params will be filled once str_list params are removed from header.keys()
str_list = ['PARENT','X','Y','Z','CHILD','TYPE']  #identifiers of particular compartments not useful for Radius comparision, would need to update if other identifiers used in .swc 
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
                        temp_line = line.split()
                        if len(temp_line) == len(header):
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
            #i want to change bottom to covert entire list of child, type and parent to strings at once instead of line by line


                 
    archive_dict[d1] = {}                                   #sets archive name (in dirs) as nested dictionary within complete archive dictionary
    basal_list = zip(*basal_list)                           # is zip making the data a tuple organization (yes)
    archive_dict[d1]['Basal'] = basal_list                  #keeps archive --> separate oranization of Apical and Basal
    if apical_list:
        apical_list = zip(*apical_list)                        
        archive_dict[d1]['Apical'] = apical_list            #change organization to archive_dict['Apical'][d1]

params = {key:header[key] for key in header if key not in str_list}  #<-- if user input parameters send here
'''
for num, param in enumerate(parameters):                  #make y value for plot specified (input) possibly function to plot against radius and separate residuals
    for archive in archive_dict.keys():
        for cell_type in archive_dict[archive]:            #only potential problem here would be that basal and apical would be on the same plot
            label = archive + ' ' + cell_type
            plt.plot(archive_dict[archive][cell_type][header[param]], archive_dict[archive][cell_type][header['RADIUS']], 'o', label = label)
    title = str(archive) + ' ' + str(param) #str(cell_type) 
    plt.xlabel(str(param))
    plt.ylabel('Radius')
    plt.title(title)
    plt.legend()
    plt.show()
'''
# temporarily commented out to try Pearson's R below...  #rearrange with how I am changing the dictionary organization archive_dict --> compartement type --> archive
'''
for param in params.keys():
    #maybe initialize here axis = fig.axes
    for archive in dirs:          #one plot per archive
        #would have to initialize fig,ax = plt.subplot()
        for num,cell_type in enumerate(archive_dict[d1]):
            #ax[num].plot()        #if i wanted to plot basal and apical on same plot...
            #plt.figure(num)
            figure(num = num, figsize = (14,8))
            label = archive + ' ' + cell_type
            plt.plot(archive_dict[archive][cell_type][header[param]], archive_dict[archive][cell_type][header['RADIUS']], 'o', label = label)
            title = cell_type + ' ' + str(param) + ' vs. ' + 'RADIUS' #str(root) + ' ' + str(param) #str(cell_type) 
            plt.xlabel(str(param))
            plt.ylabel('RADIUS')
            plt.title(title)
            plt.legend()
            #plt.savefig(title + '.png')                       #works better but still missing the basal graphs if both basal and apical data exist in particular dictionary
    #plt.clf()                                         #more explicit way to refer to figures
    plt.show()                      #plots as intended, but some remaining issues with saving figures automatically as would appear in window


''' #non-working as of now. Will try and correlate across multiple archives of same region for more general trends...
#corr = np.zeros(shape=(len(params),len(params)), dtype = object)
#how to use compartment type variable instead of specific 'Apical'
#loop over compartment type
corr = pd.DataFrame(index = params, columns = params)                #easier to visualize comparison between parameter correlations
#for key in archive_dict.keys():                                     #this would print each corr separately for each archive
for r,rparam in enumerate(params.keys()):                            #how to concatenate each archive for single corr table
   for c,cparam in enumerate(params.keys()):
       xlist = []; ylist = []                                        #option to create concatenation dictionary
       for key in archive_dict.keys():                               #why iss archive_dict[key]['Apical'][header[rparam]] tuple?
           xlist = list(archive_dict[key]['Basal'][header[rparam]]) + xlist #concatenation <-- i like this..
           ylist = list(archive_dict[key]['Basal'][header[cparam]]) + ylist #can change to list comprehension if wanted
       corr.loc[rparam,cparam] = round(pearsonr(xlist, ylist)[0], 4)         #in or out of loop for temp or permanent x and y list
       #xlist = []; ylist = []                                       #first for initial correlation for params
       #for d1 in dirs:                                              #take all list code out of loop to make permanent list
           #xlist =                                                  #plot residuals of Radius and compare to other params
           #add_val = archive_dict[d1]['Apical'][header[rparam]]
           #add_x = [val in zip(archive_dict[d1]['Apical'][header[rparam]], ]
           #add_y = [val for val in archive_dict[d1]['Apical'][header[cparam]]]
           #xlist.append(add_x[0])
           #ylist.append(add_y[0])
       #zip(*archive_dict[dirs]['Apical'][header[rparam]]) #??
       #corr.loc[rparam,cparam] = round(pearsonr(xlist,ylist)[0], 4) 
print(corr)

#maybe if correlation (linear) is good for particular combination do something else with that in code
#i.e. if corr value is above certain acceptance level of good relationship
'''
if choice == 'optimize':
    print(choice)
if choice == 'graph':
    print(choice)
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
