Description
============
shape_shifter is a morphology converter which takes in Genesis (.p) files of neurons and simplifies the morphology by combining sets of compartments with the same (or similar) radius and with combined length less than specified value of electrotonic length, caculated using either AC or DC methods. Previous trials of this program has resulted in reduction of the orginal file by almost 90%.  shape_shifter also can be used to transform a simple morphology with extremely long compartments into smaller compartments with same topology (connectivity).

**input files**

Morphology files can be attained through http://neuromorpho.org/ and these .swc files on NeuroMorpho can be converted to .p files automatically using convert_swc_pfile.py.  We strongly recommend using this instead of cvapp as parsing issues are common causing a freeze of the program.  

**convert_swc_pfile usage**
 Input:  python convert_swc_pfile.py --file filename.swc
 Output: filenameconvert.p

**shape_shifter usage**
 Input:  
 ``` 
 python shape_shifter.py --file 'filename.p' --type 'radii' 
 ```
 Output: filenameout.p

type can be:
  - '0'        just remove compartments of size 0
  - 'condense' combine compartments of similar radius (specify --rad_diff 0 to only combine identical radii),
               electrotonic length not to exceed max_len* lambda.
  - 'expand'   to change single, long compartment (e.g. a Neuron software segment) into multiple smaller compartments
               electrotonic length of subdivided compartments do not exceed max_len* lambda
  - 'radii'    to change the diameter value depeneding on the distance to the end of the longest branch

+ optional arguments:

  - can specify alternative values to default rm [4], cm [0.01], ri [2.5] for calculating lambda, units are SI
  - can specify frequency (--f) for ac lambda calculation, default is 0.1 hz, specify 0 to use dc lambda
  - can specify maximum electrotonic length (--max_len), default is 0.1
  - can specify criteria maximum difference (rad_diff) in radii for combining adjacent compartments, default is 0.1 (=10%)
  
change info or debug parameter to get more or less information while running program
           

Ideal Usage Scenario
============
neuromorpho.org morphology files often inaccurate compartment (node) diameters

+ convert .swc file --> .p file using convert_swc_pfile.py

 - Input:  
 ``` 
 python convert_swc_pfile.py --file filename.swc 
 ```
 - Output:  filenameconvert.p

+ run shape_shifter.py through type radii for more accurate diameters

 - Input:   
 ``` 
 python shape_shifter.py --file filenameconvert.p --type radii
 ```
 - Ouput:   filenameconvertout.p
 
+ run shape_shifter.py again through type condense to combine similar compartments together for future simulation
 can also pass additional parameters i.e. f (frequency) and radius difference from optional arguments above
 
 - Input:   
 ``` 
 python shape_shifter.py --file filenameconvert.p --type condense --f 100 --rad_diff 0 
 ```
 - Ouput:   filenameconvertoutout.p

CVAPP - Description
============
Mentioned above, 
If using cvapp, provided in this repository. We strongly recommend that the 1st three soma lines in the cvapp output file be combined into a single coma compartment by hand prior to running shape_shifter, because cvapp output of neuromorpho files typically has 3 soma compartments, with the first, parent compartment (which cannot be eliminated by shape_shifter) having zero size. In addition, you may need to remove blank lines from the .p file to avoid error messages.

**CVAPP usage**

- in Unix:
  java -classpath cvapp.jar: $CLASSPATH cvapp
- on Windows:
  java -classpath cvapp.jar; $CLASSPATH cvapp

1. click on "file", then "open" in the dropdown, and then navigate to your swc file
2. click on "file", then "save as Genesis - flat", create filename
3. click on "file", then "quit" to exit cvapp
4. edit the output .p file as follows:
- add *cartesian, *asymmetrical, and all the resistance, voltage and capacitance values
- edit the lines specifying the soma. If the 1st compartment has 0,0,0 as coordinates, this will create a compartment of length 0 and volume 0 (not good. The next two compartments are likely specifying the left and right edges (or top and bottom edges) of the soma.  Change these three lines by
  - changing the 1st soma comp x, y, z to be equal the sum of the x, y, z values of the 2nd and 3d soma comps
  - delete the 2nd and 3d soma comps
5. run shape_shifter
