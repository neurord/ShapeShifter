import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def _get_node_features(stats,node,term=False):
    O = stats.order_of_node(node)
    D = stats.order_of_node(node)
    ED = stats.get_Euclidean_length_to_root(node)
    PL = stats.get_pathlength_to_root(node)
    SPL = stats.get_segment_pathlength(node)
    SED = stats.get_segment_Euclidean_length(node)
    if term:
        return O,D,ED,PL,SPL,SED
    else:
        PA = stats.partition_asymmetry(node)
        AMP = stats.bifurcation_angle_vec(node)
        return O,D,ED,PL,SPL,SED,PA,AMP

def perform_2D_analysis(destination,filter="*.swc",max_n=None):
    """
    Wrapper function to perform an analysis of the vector features of one neuronal morphology (in the SWC format and with 3-point soma)

    For both the terminal points and the branching points the following features are recorded

    - Order of the node
    - degree of the node
    - Euclidean distance to the soma
    - path length to the soma
    - pathlength of the segment (coming in to a node)
    - Euclidean distance of the segment (coming in the a node)
    - branch angle amplitude [branch points only]
    
    Two text files are generated, for terminal and branching points. Each row corresponds to a node (point) and the six
    columns correspond to the features above.

    Parameters
    -----------
    destination : string
        string with the location of where to find the SWC files.
    """
    
    # change to a directory of choice for I/O
    os.chdir(destination)
    
    # load morphologies and initialize statistics
    all_f = glob.glob(filter)
    if not max_n == None:
        all_f = all_f[:max_n]
    swc_trees = {}
    individual_stats = {}
    for f in all_f:
        print "f: ", f
        cell_name = f.split(filter)[0]
        temp_tree = btmorph.STree2()
        temp_tree.read_SWC_tree_from_file(f)
        swc_trees[cell_name] = temp_tree
        temp_stats = btmorph.BTStats(temp_tree)
        individual_stats[cell_name] = temp_stats

    """ 2D features calculated per neuron and for the whole population.
        For each bifurcation and terminal point, the following features \
        at that location will be recorded: order, degree, asymmetry_index, \
        path length, euclidean distance,(bifurcation angle, local)"""

    for cell_name in individual_stats:
        print "analyzing cell %s" % cell_name
        bif_nodes = individual_stats[cell_name]._bif_points
        term_nodes = individual_stats[cell_name]._end_points

        term_to_write = ""
        for node in term_nodes:
            O,D,ED,PL,SPL,SED= \
              _get_node_features(individual_stats[cell_name],node,term=True)
            term_to_write += "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (O,D,ED,PL,SPL,SED)
            
        term_writer = open(cell_name+"_terms_multivariate.txt","w")
        term_writer.write(term_to_write)
        term_writer.flush()
        term_writer.close()

        bif_to_write = ""
        for node in bif_nodes:
            O,D,ED,PL,SPL,SED,PA,AMP= \
              _get_node_features(individual_stats[cell_name],node)
            bif_to_write += "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (O,D,ED,PL,SPL,SED,PA,AMP)
        bif_writer = open(cell_name+"_bifs_multivariate.txt","w")
        bif_writer.write(bif_to_write)
        bif_writer.flush()
        bif_writer.close()

if __name__=="__main__":
    destination = "."
    if len(sys.argv) == 2:
        destination = sys.argv[1]
    perform_2D_analysis(destination)
