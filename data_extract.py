'''<Extracts data values from .swc files>
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

#####btmorph requires python2.7; not updated to work with python3#####

#Usage:  python2 morph_feature_extract.py --path /path/to/folder
#        where folder contains all CNG.swc files intended to extract feature data 
#Output: morphology.CNG_extract.txt

#George Mason University
#Jonathan Reed
#May 25, 2020

import numpy as np
import matplotlib.pyplot as plt
import btmorph

import argparse
import datetime
import os
import glob

'''Copy end compartment and branching compartment data for continuing compartments'''
def find_paths(swc_tree, node_data, stats, local_list): 
    for item in local_list:                              #will complete the endpoint(s) for subsequent node for each end/branch point sent to function
        if item[INDEX] == node_data[COMP].parent.index and item[INDEX] not in soma and item[COMP] not in stats._bif_points:
            comp_len = stats.get_pathlength_to_root(node_data[COMP]) - stats.get_pathlength_to_root(node_data[COMP].parent)
            new_length = node_data[BRANCH_LEN] + comp_len
            new_endpath = node_data[PATH_TO_END] + comp_len 
            item.extend((node_data[END_POINT],node_data[BRANCH_POINT],new_length,new_endpath))
            local_list = find_paths(swc_tree, item, stats, local_list)
    return local_list                                                    

'''Flatten nested lists of values into single list of values'''
def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

'''Start of Working Code'''
parser = argparse.ArgumentParser()
parser.add_argument("--path", type = str)  
args = parser.parse_args()                                                                                      
path = args.path
if not os.path.exists(path):
    print('Error in Path')
else:
    path = path + '*CNG.swc'
    
    '''Load .swc morphologies from folder'''
    for filename in glob.glob(path):
        soma = []; morph_list = []; local_list = []; degree_list = []
        #soma compartments, feature list, feature calculation values, compartment order from terminals toward soma (Node Degree)

        #Load morphology with btmorph
        swc_tree = btmorph.STree2()                         
        swc_tree.read_SWC_tree_from_file(str(filename))    
        stats = btmorph.BTStats(swc_tree)                   #morphometric data from each compartment in morphology
             
        INDEX = 0; COMP = 1; DEGREE = 2; END_POINT = 3; BRANCH_POINT = 4; BRANCH_LEN = 5; PATH_TO_END = 6
        for node in swc_tree:
            local_list.append([node.index, node, stats.degree_of_node(node)])
            if not stats.degree_of_node(node) in degree_list:    #create list of available node degrees for separation
                degree_list.append(stats.degree_of_node(node))
            if node.get_content()['p3d'].type == 1:              #separate soma compartment(s) from dendritic compartments
                soma.append(node.index)                          #soma compartments designated Type = 1 from .swc file format
        degree_list.sort(key=lambda x: x)
        local_list.sort(key=lambda x: x[DEGREE])                 #sorted list by node degree

        '''Calculate Feature Values associated with Path Distances'''
        #Start calculations from terminal branches with node degree = 0
        for item in stats._end_points:
            self_node = [c for c in local_list if c[COMP] == item]
            self_node[0].extend(([self_node[0][INDEX]],[0],0,0))
            #end points, branch points, total dendritic length, longest path to terminal end 
            local_list = find_paths(swc_tree, self_node[0], stats, local_list)     #calculate compartments from terminal ends to distal branch points

        #Start calculations from branch compartments and continuing (in between) compartments; traversing towards soma
        for num in degree_list[1:]:
            for item in stats._bif_points:                                            #branch/bifurcation compartments
                if item.index not in soma and stats.degree_of_node(item) == num:
                    self_node = [c for c in local_list if c[COMP] == item]
                    child_1 = [c for c in local_list if c[COMP] == item.children[0]]  #find children (subsequent) compartments of branch compartments
                    child_2 = [c for c in local_list if c[COMP] == item.children[1]] 

                    end_1 = child_1[0][END_POINT]
                    end_2 = child_2[0][END_POINT]
                    self_node[0].append([end_1, end_2])
                    self_node[0][END_POINT] = list(flatten(self_node[0][END_POINT]))        #find all endpoints which stem from branch compartments

                    branch_1 = child_1[0][BRANCH_POINT]
                    branch_2 = child_2[0][BRANCH_POINT]
                    self_node[0].append([branch_1,branch_2,self_node[0][INDEX]]) 
                    self_node[0][BRANCH_POINT] = list(flatten(self_node[0][BRANCH_POINT]))  #find all branchpoints which stem from branch compartments 

                    child_1_path = child_1[0][BRANCH_LEN]
                    child_2_path = child_2[0][BRANCH_LEN]
                    branch_path_1 = child_1_path + (stats.get_pathlength_to_root(child_1[0][COMP]) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    branch_path_2 = child_2_path + (stats.get_pathlength_to_root(child_2[0][COMP]) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    branch_len = branch_path_1 + branch_path_2
                    self_node[0].append(branch_len)                         #calculate the total dendritic length rooted at branch compartment

                    end_path_1 = child_1[0][PATH_TO_END] + (stats.get_pathlength_to_root(child_1[0][COMP]) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    end_path_2 = child_2[0][PATH_TO_END] + (stats.get_pathlength_to_root(child_2[0][COMP]) - stats.get_pathlength_to_root(self_node[0][COMP]))
                    path = end_path_1 if end_path_1 > end_path_2 else end_path_2
                    self_node[0].append(path)                               #calculate the longest path to terminal end rooted at branch compartment

                    #calculate feature values by adding compartment length for compartments in-between branch compartments
                    local_list = find_paths(swc_tree, self_node[0], stats, local_list)
        
        '''Calculate Initial Compartment as Feature Value'''
        direct_to_soma = {}
        for node in swc_tree.get_nodes():
            if node.index not in soma:
                if node.parent.index in soma:
                    direct_to_soma[node.index] = 0
                else:
                    direct_to_soma[node.index] = 1
            else:
                direct_to_soma[node.index] = -1

        '''Organize all Feature Values by Compartment'''
        for node in swc_tree.get_nodes():
            if node.index not in soma:                                     #Each compartment has the following identifiers/features:
                temp = []
                temp.append(node.index)                                    #'Child'           - Name
                temp.append(node.get_content()['p3d'].type)                #'Type'            - Designation; soma(1), axon(2), basal(3) or apical(4) dendrite
                for num in range(3):
                    temp.append(node.get_content()['p3d'].xyz[num])        #'XYZ'             - Coordinates
                temp.append(node.get_content()['p3d'].radius)              #'Radius'          - Radius             
                temp.append(stats.degree_of_node(node))                    #'Node Degree'     - Terminal Branch Order
                temp.append(stats.order_of_node(node))                     #'Node Order'      - Initial Branch Order
                temp.append(node.parent.index)                             #'Parent'          - Name of Parent Compartment
                temp.append(node.parent.get_content()['p3d'].radius)       #'Parent Radius'   - Radius of Parent Compartment
                temp.append(stats.get_pathlength_to_root(node))            #'Path Distance'   - Path from soma to Compartment

                self_node = [c for c in local_list if c[COMP] == node]     #Feature values from local list calculations
                temp.append(stats.local_horton_strahler(node))             #'Horton Strahler' - Horton Strahler Branch Order
                temp.append(self_node[0][BRANCH_LEN])                      #'Branch Length    - Total Dendritic Length rooted at Compartment
                temp.append(self_node[0][PATH_TO_END])                     #'Path to End'     - Longest Path to Terminal End
                temp.append(direct_to_soma[node.index])                    #'Direct Soma'     - Initial Compartments directly stemming from soma
                if node in stats._bif_points:                              #'Branch Point'    - Branching/Bifurcation Compartments
                    temp.append(0)  #if bifurcation or branching point
                else:
                    temp.append(1)  #if continuing segment
                if node.parent in stats._bif_points:                       #'Branch Children' - Compartments stemming after branching/bifurcation
                     temp.append(0) #if child to bifurcation or branching point
                else:
                    temp.append(1)  #if remaining segment
                morph_list.append(temp)

        '''Save to file'''
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
        outfile.write('*CHILD; TYPE; X; Y; Z; RADIUS; NODE_DEGREE; NODE_ORDER; PARENT; PARENT_RAD; PATH_DIS; HS; BRANCH_LEN; PATH_TO_END; DIRECT_SOMA; BRANCH_POINT; BP_CHILD')
        outfile.write('\n')

        for line in morph_list:
            write_line = [str(val) for val in line]
            write_line = ' '.join(write_line)
            outfile.write(str(write_line) + '\n')

        print('File Created  :  ')
        print(filename)
        outfile.close()
        
