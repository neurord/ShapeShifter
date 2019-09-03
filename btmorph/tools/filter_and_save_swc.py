import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def filter_and_save_SWC(destination,filter,types=range(10),prefix="_filtered"):
    """
    Removes points from a SWC structure and saves the new SWC to a file.

    Can be used to remove unwanted structures that are identifiable by \
    the type-field in the SWC description. Specification of (standard) \
    SWC type fields can be found `here <http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html>`_.

    To select the basal dendrites only, use the argument `types=[1,3]`:\
    1 to select the soma and 3 for the basal dendrites themselves.

    Parameters
    -----------
    destination : string
        string with the location of where to find the SWC files.
    types : list of int
        types that are to be remained in the SWC file.
    """
    # change to a directory of choice for I/O
    pwd = os.getcwd()
    os.chdir(destination)    

    # load morphologies and initialize statistics
    all_f = glob.glob(filter)

    for f in all_f:
        print "processing file: ",f
        temp_tree = btmorph.STree2()
        temp_tree.read_SWC_tree_from_file(f,types=types)
        outN = f.split(".swc")[0]+prefix+".swc"
        temp_tree.write_SWC_tree_to_file(outN)

    os.chdir(pwd)
