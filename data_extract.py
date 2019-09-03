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
#        where folder contains all .swc files intended for data_extract and btmorph directory
#Output: morphology.CNG_extract.txt

#George Mason University
#Jonathan Reed
#Sep 3, 2019

import numpy as np
import matplotlib.pyplot as plt
import btmorph

import argparse
import datetime
import os
import glob

def find_end(degree, node_data, stats, local_list):
    end_point = node_data[END_POINT]
    parent = node_data[COMP].parent.index
    for item in local_list:                              #will complete the endpoint(s) for subsequent node for each end/branch point sent to function
        if item[INDEX] == parent and item[INDEX] not in soma and item[COMP] not in stats._bif_points:
            item.append(end_point)
            local_list = find_end(degree, item, stats, local_list)
    return local_list

'''
def find_branch(node_data, local_list):
    branch_point = node_data[BRANCH_POINT]
    parent = node_data[COMP].parent.index
    for item in local_list:                              #will complete the endpoint(s) for subsequent node for each end/branch point sent to function
        if item[INDEX] == parent and item[INDEX] not in soma and item[COMP] not in stats._bif_points:
            item.append(branch_point)
            local_list = find_branch(node_data, local_list)
    return local_list
'''

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)   #may change, but for now is glob-command to look for
args = parser.parse_args()                                                                                      
path = args.path
if not os.path.exists(path):
    print('Error in Path statement')
else:
    path = path + '*CNG.swc'
    for filename in glob.glob(path):
        print(os.path.basename(filename))
        #filename = os.path.basename(filename)              
        swc_tree = btmorph.STree2()                         
        swc_tree.read_SWC_tree_from_file(str(filename))     #reads in .swc file, each subsequent node within SNode2 designation
        stats = btmorph.BTStats(swc_tree)                   #obtains morphometric data from each compartment separate from SNode2
             
        soma = []                 #relavant soma indices 
        morph_list = []           #data values to send to file
        local_list = []           #local values to calculate end_distance for each node
        degree_list = []          #all possible degree values for particular node
        COMP = 1; DEGREE = 2; INDEX = 0; END_POINT = 3; PATH = 4
        for node in swc_tree:
            local_list.append([node.index, node, stats.degree_of_node(node)])
            if not stats.degree_of_node(node) in degree_list:
                degree_list.append(stats.degree_of_node(node))
            if node.index == 1:                               #soma has largest degree_of_node
                max_degree = stats.degree_of_node(node)
            if node.get_content()['p3d'].type == 1:
                soma.append(node.index)
        degree_list.sort(key=lambda x: x)
        local_list.sort(key=lambda x: x[DEGREE])              #sorts list by degree of node

        min_degree = 1 #degree is initial degree designation for all terminal compartments (i.e. no branching, only branches)
        
        for item in stats._end_points:
            #print(item.index)
            self_node = [c for c in local_list if c[COMP] == item]                
            self_node[0].append([self_node[0][INDEX]])                            #endpoint(s) now being filled
            local_list = find_end(min_degree, self_node[0], stats, local_list)    #end_points sent to function find_parent for Degree 1 nodes

        for num in degree_list[1:]:         
            #print(num)
            for item in stats._bif_points:
                if stats.degree_of_node(item) == num:
                    #these are the bifurication points
                    self_node = [c for c in local_list if c[COMP] == item] 
                    child_1 = [c[END_POINT] for c in local_list if c[INDEX] == item.children[0].index] #finds child from local_list and saves correct end_point
                    child_2 = [c[END_POINT] for c in local_list if c[INDEX] == item.children[1].index]
                    self_node[0].append([child_1, child_2])
                    local_list = find_end(num, self_node[0], stats, local_list)

        for num, item in enumerate(local_list): 
            path = 0
            if item[DEGREE] > 1 and item[INDEX] not in soma:
                item[END_POINT] = list(flatten(item[END_POINT]))                   #contains all possible end_points downstream of node
                for val in item[END_POINT]:
                    end = swc_tree.get_node_with_index(val)                        #tests each end_point to see longest path for end_distance
                    if (stats.get_pathlength_to_root(end) - stats.get_pathlength_to_root(item[COMP])) > path:
                        path = stats.get_pathlength_to_root(end) - stats.get_pathlength_to_root(item[COMP])
                item.append(path)

            elif item[DEGREE] == 1 and item[INDEX] not in soma:                    #all end_points have degree == 1, so end_distance to self as end_point = 0
                #end = swc_tree.get_node_with_index(item[END_POINT])
                #path = stats.get_pathlength_to_root(end) - stats.get_pathlength_to_root(item[COMP])
                item.append(0)
                #item.append(path)
        '''
        for item in local_list:
            if item[DEGREE] < 3:
                print(item)
        '''
        CHILD = 0; TYPE = 1; XYZ = 2; RADIUS = 3; NODE_DEGREE = 4; NODE_ORDER = 5; PARENT = 6; PARENT_RAD = 7; PATH_DIS = 8; END_DIS = 9; NUM_ENDS = 10; HS = 11
        #Node_Degree seems to be the number of 'leafs' downstream
        #Node_Order is amount of  bifurications upstream
        #Path_Dis is distance from node to soma
        #End_Dis is distance from node to end of the longest branch 

        for node in swc_tree.get_nodes():
            temp = []
            temp.append(node.index)
            temp.append(node.get_content()['p3d'].type)                      
            for num in range(3):
                temp.append(node.get_content()['p3d'].xyz[num])                     #['p3d'] holds node coordinates, radius, and compartment subtype
            temp.append(node.get_content()['p3d'].radius)                           #subtype as soma(1), axon(2), basal(3) or apical(4) dendrite
            temp.append(stats.degree_of_node(node))
            temp.append(stats.order_of_node(node))
            if node.index not in soma:
                temp.extend([node.parent.index,node.parent.get_content()['p3d'].radius,stats.get_pathlength_to_root(node)])                     
                self_node = [c for c in local_list if c[COMP] == node]
                #include more features here...
                temp.append(self_node[0][PATH])
                temp.append(len(self_node[0][END_POINT]))                           #if there are out of index errors, check original .swc files
                temp.append(stats.local_horton_strahler(node))                      #one instance where btmorph did not catch single end_point
            morph_list.append(temp)

            '''
            for item in local_list:
                if item[DEGREE] == 1:
                    item.append(0)
                if item[DEGREE] == 2:
                    item.append(1)
            for num in degree_list[3:]:
                for item in stats._bif_points:
                    #if stats.degree_of_node(item) == 0:
                        #item.append(0)
                    #if stats.degree_of_node(item)  == 1:
                        #initiate branch points
                        #self_node = [c for c in local_list if c[COMP] == item]
                        #self_node.append(1)
                        #local_list = find_branch(self_node, local_list)
                    if stats.degree_of_node(item) == num:# and num > 1:
                        child_1 = [c[BRANCH_POINT] for c in local_list if c[INDEX] == item.children[0].index] #finds child from local_list and saves correct end_point
                        child_2 = [c[BRANCH_POINT] for c in local_list if c[INDEX] == item.children[1].index]
                        child_1 = child_1[0]
                        child_2 = child_2[0]
                        self_node = [c for c in local_list if c[COMP] == item]
                        self_node.append(child_1 + child_2)
                        local_list = find_branch(self_node, local_list)

            for item in local_list:
                if item[DEGREE] == 4:
                    #print(item)
                    print(item[BRANCH_POINT])
            '''  
        filename = os.path.basename(filename) #places extract files in location of btmorph command
        out_name = filename.split('.swc')[0] + '_extract.txt'
        outfile = open(out_name,'w')
        outfile.write('*Original .swc file : ')
        outfile.write(filename + '\n')
        outfile.write('*Extracted .swc data on : ')
        outfile.write(str(datetime.datetime.now()) + '\n')
        outfile.write('\n')
        outfile.write('*CHILD; TYPE; XYZ; RADIUS; NODE_DEGREE; NODE_ORDER; PARENT; PARENT_RAD; PATH_DIS; END_DIS; NUM_ENDS; HS')
        outfile.write('\n')

        for line in morph_list:
            write_line = [str(val) for val in line]
            write_line = ' '.join(write_line)
            outfile.write(str(write_line) + '\n')

        print('File Created  :  ')
        print(out_name)
        outfile.close()

        #next step is to set up graphing of the values to diameter to see relationships


