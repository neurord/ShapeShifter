## Description

Shape Shifter is a morphology converter which takes in Genesis (.p) files of neurons and simplifies the morphology by combining sets of compartments with the same (or similar) radius and with combined length less than specified value of electrotonic length, caculated using either AC or DC methods. Previous trials of this program has resulted in reduction of the orginal file by almost 90%.  Shape_shifter also can be used to transform a simple morphology with extremely long compartments into smaller compartments with same topology (connectivity).

Morphology files can be attained through http://neuromorpho.org/ and these .swc files on NeuroMorpho can be converted to .p files using cvapp (http://neuron.duke.edu/cells/usage.html). We strongly recommend that the 1st three soma lines in the cvapp output file be combined into a single coma compartment by hand prior to running shape_shifter, because cvapp output of neuromorpho files typically has 3 soma compartments, with the first, parent compartment (which cannot be eliminated by shape_shifter) having zero size.  

Usage: python shape_shifter.py --file 'filename.p' --type 'radii'

+ type can be:

  - '0' - just remove compartments of size 0
  - 'radii' - combine compartments of similar radius (specify 0 for rad_diff to only combine identical radii),
          electrotonic length of condensed compartments do not exceed max_len* lambda 
  - 'expand' - to change single, long compartments (e.g. a Neuron software section) into multiple smaller compartments (e.g. Moose compartments  or Neuron segments)
           electrotonic length of subdivided compartments do not exceed max_len* lambda
           
 + can specify alternative values to default rm [4], cm [0.01], ri [2.5] for calculating lambda, units are SI
 + can specify frequency (--f) for ac lambda calculation, default is 0.1 hz, specify 0 to use dc lambda
 + can specify maximum electrotonic length (--max_len), default is 0.1
 + can specify criteria maximum difference (rad_diff) in radii for combining adjacent compartments, default is 0.1 (=10%)
 
 change info or debug parameter to get more or less information while running program
