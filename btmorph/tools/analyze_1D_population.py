import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def perform_1D_population_analysis(destination,filter="*.swc",depth="Y",bar=[200,200,200],post_name=None,max_n=None):
    """
    Wrapper function to perform a complete analysis of a population of neuronal morphologies stored in SWC format (and three-point soma).

    Computes the following features:

    - # bifurcations
    - # terminals
    - # stems
    - total length
    - total volume
    - total surface
    - max centrifugal order
    - inter bifurcation length
    - terminal segment length
    - euclidean distance between terminal tips and soma
    - path legth between terminal tips and soma

    For each feature a list is created with all values collected from all morphologies.
    
    These vectors are saved as python Pickle objects.

    At the same time a histogram is generated to display the data.

    Parameters
    -----------
    destination : string
        string with the location of where to find the SWC files.
    filter :  string
        filter to select SWC files. Default is "*.swc". See `glob
        <https://docs.python.org/2/library/glob.html>`_ documentation for more advanced filters
    depth : string
        Dimension that indicates "depth"/distance from top. Default is "Y"; NeuroMac generated files use "Z".
    bar : array of float
        Include a scale-bar with the specified dimensions
    max_n : int
        To perform the analysis on a subset of morphologies specify a \
        number. Default is None and all morphologies will be taken \
        into account.
    post_name : string
        string appended to the file name when saving. Default None
    """
    # change to a directory of choice for I/O
    pwd = os.getcwd()
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
        plt.clf()
        btmorph.plot_2D_SWC(f,outN=cell_name+"_2D.pdf",align=True,depth=depth,show_axis=True)
        #btmorph.true_2D_projections_equal(f,outN=f.split(".swc")[0]+"_projections.pdf",depth=depth,bar=bar)
        plt.close()

    """ 1D features, without dependencies on other quantities
        Both global (1 per neuron) and local (N per neuron) features"""

    one_d_stats = {}

    # bulk statistics of global features
    funcs = ["total_length","no_bifurcations","no_terminals","no_stems",\
             "total_volume","total_surface"]
    for cell_name in individual_stats:
        for func in funcs:
            method = getattr(individual_stats[cell_name], func)
            ret = method()
            if isinstance(ret,tuple):
                ret = ret[0]
            if func in one_d_stats:
                one_d_stats[func].append(ret)
            else:
                one_d_stats[func] = []
                one_d_stats[func].append(ret)

    # extract the max order            
    for cell_name in individual_stats:
        max_order = 0
        for node in swc_trees[cell_name].get_nodes():
            order = individual_stats[cell_name].order_of_node(node)
            if order > max_order:
                max_order = order
        if "max_order" in one_d_stats:
            one_d_stats["max_order"].append(max_order)
        else:
            one_d_stats["max_order"] = []
            one_d_stats["max_order"].append(max_order)

    # extract the "inter-bifurcation segment length"        
    for cell_name in individual_stats:
        bif_nodes = individual_stats[cell_name]._bif_points
        for node in bif_nodes:
            L = individual_stats[cell_name].get_segment_pathlength(node)
            if "inter_bif_L" in one_d_stats:
                one_d_stats["inter_bif_L"].append(L)
            else:
                one_d_stats["inter_bif_L"] = []
                one_d_stats["inter_bif_L"].append(L)

    # extract the terminal segment lengths
    # extract the path length / euclidean distance to terminal tips
    for cell_name in individual_stats:
        term_nodes = individual_stats[cell_name]._end_points
        for node in term_nodes:
            L = individual_stats[cell_name].get_segment_pathlength(node)
            path_L = individual_stats[cell_name].get_pathlength_to_root(node)
            eucl_L = individual_stats[cell_name].get_Euclidean_length_to_root(node)
            if "term_L" in one_d_stats:
                one_d_stats["term_L"].append(L)
            else:
                one_d_stats["term_L"] = []
                one_d_stats["term_L"].append(L)

            if "term_path_L" in one_d_stats:
                one_d_stats["term_path_L"].append(path_L)
            else:
                one_d_stats["term_path_L"] = []
                one_d_stats["term_path_L"].append(path_L)

            if "term_eucl_L" in one_d_stats:
                one_d_stats["term_eucl_L"].append(eucl_L)
            else:
                one_d_stats["term_eucl_L"] = []
                one_d_stats["term_eucl_L"].append(eucl_L)


    for func in one_d_stats:
        plt.figure(0)
        plt.clf()
        plt.hist(one_d_stats[func])
        plt.xlabel(func)
        plt.ylabel("#")
        if post_name == None:
            plt.savefig("pop_1D_"+func+"_hist.pdf")
            p_name = "pop_1D_"+func+"_raw.pkl"
        else:
            plt.savefig("pop_1D_"+func+"_hist"+post_name+".pdf")
            p_name = "pop_1D_"+func+"_raw"+post_name+".pkl"            
        pickle.dump(one_d_stats[func], open(p_name,"w"))

    # and change back to the original directory where the command was given
    os.chdir(pwd)

if __name__=="__main__":
    destination = "."
    if len(sys.argv) == 2:
        destination = sys.argv[1]
    perform_1D_population_analysis(destination)

    
