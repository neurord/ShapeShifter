'''<Reformats .p data file format into .swc file format.>
   <Output .swc file>
    Copyright (C) <2019>  <Zhi Cheng Wu>

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

#Usage:  
#Output: filenameconvert.p

#George Mason University
#Jonathan Reed and Zhi Cheng Wu
#March 1, 2019

import argparse
import datetime
from convert_swc_pfile import read_data_comments
PARENT=1
SELF=0
DIAM=5

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--file")
    args = parser.parse_args()
    filename = args.file

    lines = open(filename, 'r').readlines()

    data_lines,comment_lines=read_data_comments(lines)
    
    comment_lines.append('# Original .p file : '+filename + '\n')
    comment_lines.append('# Modified to .swc file on : '+str(datetime.datetime.now()) + '\n')

    all_new_lines = []
    #rearranges values for .p format
    for num, line in enumerate(data_lines):
        newline = [line[x] for x in range(6)]
        newline[5] = round(float(line[5])/2,4) #replace diameter with radius
        #very first line has no parent, and is of type soma
        if num == 0:
            newline.append(-1) #parent
            newline[PARENT] = 1  #type
        #splicing underscore out and choosing portion of the value desired
        else:
            newline.append( newline[PARENT].split('_')[0]) #parent
            newline[PARENT] = newline[0].split('_')[1] #type
        newline[SELF] = newline[SELF].split('_')[0] #self
        #creates list of newline that can be recorded to search for previous neurons
        all_new_lines.append(newline)
        parentneuron = int(newline[6])
        #print(num,'*** relative:', newline,'parent=',parentneuron)
        if int(parentneuron) > 0:
            for x in range(2,5):
                #convert from relative to absolute coordinates, but not for first line
                #print(newline[x], 'add to', all_new_lines[parentneuron-1][x])
                newline[x] = round(float(newline[x])+float(all_new_lines[parentneuron-1][x]),4)
            #print('   abs',newline,'\n','  parent',parentneuron,' :',all_new_lines[parentneuron-1])

    out_name = filename.split('.p')[0] + '.swc'
    outfile = open(out_name,'w')
    for line in comment_lines:
        outfile.write(line)
    for line in all_new_lines:
        write_line = ' '.join([str(val) for val in line]) #converts from list of strings to single string
        outfile.write(write_line + '\n')
    outfile.close()

    print('Converted ' + str(len(data_lines)) + ' compartments from original file : ' + str(filename))
    print('Modified file created : ' + str(out_name))
