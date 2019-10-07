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

#Usage:  python data_extract.py --path /path/to/folder
#        where folder contains all .swc files intended for data_extract
#Output: morphology.CNG_extract.txt

#George Mason University
#Jonathan Reed
#Sep 5, 2019

import numpy as np
import matplotlib.pyplot as plt
import btmorph

import argparse
import datetime
import os
import glob

'''fuction to locate parent compartment and copy end_point and branch_point data within node degree'''
def find_parent(swc_tree, node_data, stats, local_list, b_length = 0): 
    end_point = node_data[END_POINT]
    branch_point = node_data[BRANCH_POINT]
    branch_len = node_data[BRANCH_LEN]
    path_to_end = node_data[PATH_TO_END]
    parent = node_data[COMP].parent
    
    for item in local_list:                              #will complete the endpoint(s) for subsequent node for each end/branch point sent to function
        if item[INDEX] == parent.index and item[INDEX] not in soma and item[COMP] not in stats._bif_points:
            item.append(end_point)  
            item.append(branch_point)
            comp_len = stats.get_pathlength_to_root(node_data[COMP]) - stats.get_pathlength_to_root(parent)
            item.append(branch_len + comp_len)                              #maybe I should keep this as second to last append, with longest end distance as last append
            item.append(path_to_end + comp_len)
            local_list = find_parent(swc_tree, item, stats, local_list)
    return local_list                                                    

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)  
args = parser.parse_args()                                                                                      
path = args.path
if not os.path.exists(path):
    print('Error in Path statement')
else:
    path = path + '*CNG.swc'
    for filename in glob.glob(path):
        swc_tree = btmorph.STree2()                         
        swc_tree.read_SWC_tree_from_file(str(filename))     #reads in .swc file, each subsequent node within SNode2 designation
        stats = btmorph.BTStats(swc_tree)                   #obtains morphometric data from each compartment separate from SNode2
             
        soma = []                 #relavant soma indices 
        morph_list = []           #data values to send to file
        local_list = []           #local values to calculate end_distance for each node
        degree_list = []          #all possible degree values for particular node
        INDEX = 0; COMP = 1; DEGREE = 2; END_POINT = 3; BRANCH_POINT = 4; BRANCH_LEN = 5; PATH_TO_END = 6
        for node in swc_tree:
            local_list.append([node.index, node, stats.degree_of_node(node)])
            if not stats.degree_of_node(node) in degree_list:
                degree_list.append(stats.degree_of_node(node))
            if node.get_content()['p3d'].type == 1:           #soma compartments have type == 1
                soma.append(node.index)
        degree_list.sort(key=lambda x: x)
        local_list.sort(key=lambda x: x[DEGREE])              #sorts list by degree of node

        for item in stats._end_points:
            self_node = [c for c in local_list if c[COMP] == item]                
            self_node[0].append([self_node[0][INDEX]])                             #endpoint(s) now being filled
            self_node[0].append([0])                                               #lets try appending 0 as node index begins at 1 for placeholder
            self_node[0].append(0)                                                 #longest path distance to end
            self_node[0].append(0)                                                 #total (branching) path distance to end
            local_list = find_parent(swc_tree, self_node[0], stats, local_list)     #end_points sent to function find_parent for Degree 1 nodes
                
        for num in degree_list[1:]:
        #for num in degree_list[1:2]:
            for item in stats._bif_points:
                if item.index not in soma and stats.degree_of_node(item) == num:
                    #these are the bifurication points
                    self_node = [c for c in local_list if c[COMP] == item]
                    child_1 = item.children[0]
                    child_2 = item.children[1]
                    end_1 = [c[END_POINT] for c in local_list if c[COMP] == child_1] #finds child from local_list and saves correct end_point
                    end_2 = [c[END_POINT] for c in local_list if c[COMP] == child_2]
                    self_node[0].append([end_1, end_2])
                    self_node[0][END_POINT] = list(flatten(self_node[0][END_POINT]))
                    branch_1 = [c[BRANCH_POINT] for c in local_list if c[COMP] == child_1] # basically append itself + any branchpoints downstream
                    branch_2 = [c[BRANCH_POINT] for c in local_list if c[COMP] == child_2] # branch len would be the branch_len of all children as well.
                    self_node[0].append([branch_1,branch_2,self_node[0][INDEX]]) #adds previous branches and self to branch list
                    self_node[0][BRANCH_POINT] = list(flatten(self_node[0][BRANCH_POINT]))
                    #maybe also add in branch_len which is longer here... but we have the previous ends IN ITS CHILDREN>>>>
                    path_1 = [c[PATH_TO_END] for c in local_list if c[COMP] == child_1]
                    path_2 = [c[PATH_TO_END] for c in local_list if c[COMP] == child_2]
                    path_1 = path_1[0] + (stats.get_pathlength_to_root(child_1) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    path_2 = path_2[0] + (stats.get_pathlength_to_root(child_2) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    branch_len = path_1 + path_2
                    self_node[0].append(branch_len)
                    path = path_1 if path_1 > path_2 else path_2
                    self_node[0].append(path)
                    #here we want the btmorph pathlength funtionality
                    local_list = find_parent(swc_tree, self_node[0], stats, local_list)
                    
        for item in local_list:
            if item[INDEX] not in soma:
                item[BRANCH_POINT] = [x for x in item[BRANCH_POINT] if x != 0]             #how to convert 'nested' arrray into list for flatten...
                #print(item[INDEX], item[DEGREE], item[END_POINT], item[BRANCH_POINT], item[BRANCH_LEN], item[PATH_TO_END]) #item[PATH_TO_END])
        
        #currently branch_points of empty list is only for endpoints as they stem from terminal branch --> not sent to file****
        
        CHILD = 0; TYPE = 1; XYZ = 2; RADIUS = 3; NODE_DEGREE = 4; NODE_ORDER = 5; PARENT = 6; PARENT_RAD = 7; PATH_DIS = 8; END_DIS = 9; NUM_ENDS = 10; HS = 11
        #Node_Degree seems to be the number of 'leafs' downstream
        #Node_Order is amount of  bifurications upstream
        #Path_Dis is distance from node to soma
        #End_Dis is distance from node to end of the longest branch 

        for node in swc_tree.get_nodes():
            if node.index not in soma:
                temp = []
                temp.append(node.index)
                temp.append(node.get_content()['p3d'].type)                      
                for num in range(3):
                    temp.append(node.get_content()['p3d'].xyz[num])                     #['p3d'] holds node coordinates, radius, and compartment subtype
                temp.append(node.get_content()['p3d'].radius)                           #subtype as soma(1), axon(2), basal(3) or apical(4) dendrite
                temp.append(stats.degree_of_node(node))
                temp.append(stats.order_of_node(node))
                temp.extend([node.parent.index,node.parent.get_content()['p3d'].radius,stats.get_pathlength_to_root(node)])                     
                self_node = [c for c in local_list if c[COMP] == node]
                #include more features here...
                temp.append(len(self_node[0][END_POINT]))                           #if there are out of index errors, check original .swc files
                temp.append(stats.local_horton_strahler(node))                      #one instance where btmorph did not catch single end_point
                temp.append(self_node[0][BRANCH_LEN])
                temp.append(self_node[0][PATH_TO_END])
                morph_list.append(temp)

        dirname = os.path.dirname(filename)
        filename = os.path.basename(filename)
        filename = filename.split('.swc')[0] + '_extract.txt'
        completename = dirname + '/' + filename
        outfile = open(completename,'w')
        outfile.write('*Original .swc file : ')
        outfile.write(filename.split('_extract.txt')[0] + '.swc'  + '\n')
        outfile.write('*Extracted .swc data on : ')
        outfile.write(str(datetime.datetime.now()) + '\n')
        outfile.write('\n')
        outfile.write('*CHILD; TYPE; X; Y; Z; RADIUS; NODE_DEGREE; NODE_ORDER; PARENT; PARENT_RAD; PATH_DIS; NUM_ENDS; HS; BRANCH_LEN; PATH_TO_END')
        outfile.write('\n')

        for line in morph_list:
            write_line = [str(val) for val in line]
            write_line = ' '.join(write_line)
            outfile.write(str(write_line) + '\n')

        print('File Created  :  ')
        print(filename)
        outfile.close()
 
