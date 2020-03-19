'''<Extracts data values from .txt .p .swc files>
    Copyright (C) <2020>  <Zhi Cheng Wu>

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

#Usage:  python extract_dynamic.py or import 
#        Have the data file being extracted within the same cwd, change file name in main class
#        Data file's first line of string must be the header of the data values of the columns separated by spaces
#        Data extracted are all changed to floats, only extract numerical values of data (no units, etc)
#        
#        Use as function through "from extract_dynamic import create_dict" in beginning line of main python code
#        create_dict function returns a dictionary so a variable has to be assigned when called: temp = create_dict('*.txt')
#        You can then manipulate the dictionary using the keys from the header of the data
#        Can work with *.swc, *.p, *.txt
#George Mason University
#Zhi Cheng Wu
#Feb. 27, 2020

import re

def create_dict(pythonfile):
    newFile = open(pythonfile,'r')
    content = newFile.readlines() #content now list of lines
    final_dict = {}
    fail = False
    temp_dict = {}
    report = 0
    header = None
    row = None
    for row in content:
        if row[0] == '*' and row[1] == 'C':
        #naming system of dict dynamic
            header = row.split('; ')
            #if '*CHILD' in header:
            header[0] = 'CHILD'
            for categories in header:
                temp_dict[categories] = []
            break
    final_dict = temp_dict #save point of category header keys
    for line in content:
        fail = False
        row = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        count = 0
        for categories in header:
            #first line of content has the header containing string to label the dict keys
            #using try in order to conserve processing power to determine whether
            try:
                temp_dict[categories].append(float(row[count]))
                count = count+1
            except:
                fail = True
        #if theres a point where a data value is missing/incorrectly named/extra data value in one line
        if fail == False:
            final_dict = temp_dict
        else:
            temp_dict = final_dict
            report += 1
    if report > len(header):
        print("Possible error in data set format.")
    return temp_dict
def main():
    temp = create_dict('Adol-20100419cell1.CNG_extract.txt')
    print(temp)
if __name__ == "__main__":
    main()
        
