Description
============
shape_shifter is a morphology converter which takes in Genesis (.p) files of neurons and simplifies the morphology by combining sets of compartments with the same (or similar) radius and with combined length less than specified value of electrotonic length, caculated using either AC or DC methods. Previous trials of this program has resulted in reduction of the orginal file by almost 90%.  shape_shifter also can be used to transform a simple morphology with extremely long compartments into smaller compartments with same topology (connectivity).

**input files**

Morphology files can be attained through http://neuromorpho.org/ and these .swc files on NeuroMorpho can be converted to .p files using cvapp, provided in this repository. We strongly recommend that the 1st three soma lines in the cvapp output file be combined into a single coma compartment by hand prior to running shape_shifter, because cvapp output of neuromorpho files typically has 3 soma compartments, with the first, parent compartment (which cannot be eliminated by shape_shifter) having zero size. In addition, you may need to remove blank lines from the .p file to avoid error messages.

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
 
**shape_shifter usage**

python shape_shifter.py --file 'filename.p' --type 'radii'
+ type can be:

  - '0' - just remove compartments of size 0
  - 'radii' - combine compartments of similar radius (specify 0 for rad_diff to only combine identical radii),
          electrotonic length of condensed compartments do not exceed max_len* lambda 
  - 'expand' - to change single, long compartments (e.g. a Neuron software section) into multiple smaller compartments (e.g. Moose compartments or Neuron segments). Resulting electrotonic length of subdivided compartments does not exceed max_len* lambda

+ optional arguments:

  - can specify alternative values to default rm [4], cm [0.01], ri [2.5] for calculating lambda, units are SI
  - can specify frequency (--f) for ac lambda calculation, default is 0.1 hz, specify 0 to use dc lambda
  - can specify maximum electrotonic length (--max_len), default is 0.1
  - can specify criteria maximum difference (rad_diff) in radii for combining adjacent compartments, default is 0.1 (=10%)
  
change info or debug parameter to get more or less information while running program
           

