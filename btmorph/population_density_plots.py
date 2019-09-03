import sys,time
sys.setrecursionlimit(10000)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

import os,glob
import btmorph

def population_density_projection(destination=".",filter="*.swc",outN=None,depth="Z",precision=[10,10,10]):
    """
    Plot a 2D heat-map of all neurites in a population. A population can \
    be sepcified by setting the `destination` and `filter` arguments.\

    Parameters
    ----------
    destination : string
        Set the directory in which the population is located. Default "." sets the current working directory.
    filter : string
        Filter the file-names to constrain the population. Default "*.swc" \
        will select all SWC files in one directory. The filter uses the \
        `glob <https://docs.python.org/2/library/glob.html>`_ module \
        so make sure glob understands your filter!
    outN : string
        File name to save the figure to. Default is None and the figure will not be saved by default.
    depth : string
        Indicate the axis from superficial to deep. In most files this is "Y".\
        Default, however, is set to the mathematical standard "Z".
    precision : list of ints
        Set the bin size of the heatmap. Default [10,10,10] (in micron)
    """

    all_SWC = {}

    min_x,min_y,min_z = 1000000,1000000,1000000
    max_x, max_y, max_z = -1,-1,-1
    
    for file_name in glob.glob(destination+"/"+filter):
        print "Processing file: ", file_name
        x = open(file_name,'r')
        soma = None
        SWC = {}
        for line in x :
            if(not line.startswith('#')) :
                splits = line.split()
                index = int(splits[0])
                if index == 1:
                    soma_x = float(splits[2])
                    soma_y = float(splits[3])
                    soma_z = float(splits[4])

                n_type = int(splits[1])
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z
                r = float(splits[5])
                parent = int(splits[-1])
                SWC[index] = (x,y,z,r,parent,n_type)

                # keep track of the overall dimensions of the population                
                max_x = x if x > max_x else max_x
                max_y = y if y > max_y else max_y
                max_z = z if z > max_z else max_z
                min_x = x if x < min_x else min_x
                min_y = y if y < min_y else min_y
                min_z = z if z < min_z else min_z
        all_SWC[file_name] = SWC

    print "%f< X < %f" % (min_x,max_x)
    print "%f< Y < %f" % (min_y,max_y)
    print "%f< Z < %f" % (min_z,max_z)

    if depth == "Y":
        front_range_x = np.arange(min_x,max_x,precision[0])
        front_range_y = np.arange(min_y,max_y,precision[1])
    else:
        front_range_x = np.arange(min_x,max_x,precision[0])
        front_range_y = np.arange(min_z,max_z,precision[2])        
    front_dens = np.zeros((len(front_range_x),len(front_range_y)))
    
    for file_name in all_SWC.keys():
        SWC = all_SWC[file_name]
        for index in SWC.keys():
            C = SWC[index]
            # front: if depth=="Z": XZ plane, else XY plane
            if depth=="Y":
                X = C[0]
                Y = C[1]
            else:
                X = C[0]
                Y = C[2]    
            x_ind = np.nonzero( (X>=front_range_x) & (X<=front_range_x+precision[0]))
            y_ind = np.nonzero( (Y>=front_range_y) & (Y<=front_range_y+precision[1]))
            for x_i in x_ind:
                for y_i in y_ind:
                    front_dens[x_i,y_i] = front_dens[x_i,y_i] + 1
    plt.figure()
    if(depth == "Y"):
        plt.imshow(front_dens.T,origin="lower",aspect="auto",extent=[min_x,max_x,min_y,max_y])
    else:
        plt.imshow(front_dens.T,origin="lower",aspect="auto",extent=[min_x,max_x,min_z,max_z])
        
    if not outN == None:
        plt.savefig(outN)

def population_2D_density_projections(destination=".",filter="*.swc",\
                                      outN=None,depth="Z",\
                                      limits=None,\
                                      precision=[10,10,10]):
    """
    Plot a pseudo-3D heat-map of all neurites in a population by means \
    of three 2D plots. A population can be sepcified by setting the \
    `destination` and `filter` arguments.

    Parameters
    -----------
    destination : string
        Set the directory in which the population is located. Default "." sets the current working directory.
    filter : string
        Filter the file-names to constrain the population. Default "*.swc" \
        will select all SWC files in one directory. The filter uses the \
        `glob <https://docs.python.org/2/library/glob.html>`_ module \
        so make sure glob understands your filter!
    outN : string
        File name to save the figure to. Default is None and the figure will not be saved by default.
    depth : string
        Indicate the axis from superficial to deep. In most files this is "Y".\
        Default, however, is set to the mathematical standard "Z".
    precision : list of ints
        Set the bin size of the heatmap. Default [10,10,10] (in micron)
    """
                                      
    # quickly run theriough all data to find the limits and collect data
    all_SWC = {}
    min_x,min_y,min_z = 1000000,1000000,1000000
    max_x, max_y, max_z = -1,-1,-1
    for file_name in glob.glob(destination+"/"+filter):
        print "Processing file: ", file_name
        x = open(file_name,'r')
        soma = None
        SWC = {}
        for line in x :
            if(not line.startswith('#')) :
                splits = line.split()
                index = int(splits[0])
                if index == 1:
                    soma_x = float(splits[2])
                    soma_y = float(splits[3])
                    soma_z = float(splits[4])

                n_type = int(splits[1])
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z
                r = float(splits[5])
                parent = int(splits[-1])
                SWC[index] = (x,y,z,r,parent,n_type)

                # keep track of the overall dimensions of the population                
                max_x = x if x > max_x else max_x
                max_y = y if y > max_y else max_y
                max_z = z if z > max_z else max_z
                min_x = x if x < min_x else min_x
                min_y = y if y < min_y else min_y
                min_z = z if z < min_z else min_z
        all_SWC[file_name] = SWC

    print "%f< X < %f" % (min_x,max_x)
    print "%f< Y < %f" % (min_y,max_y)
    print "%f< Z < %f" % (min_z,max_z)

    if depth == "Y":
        front_range_x = np.arange(min_x,max_x,precision[0])
        front_range_y = np.arange(min_y,max_y,precision[1])
        side_range_x = np.arange(min_z,max_z,precision[2])
        side_range_y = np.arange(min_y,max_y,precision[1])
        top_range_x = np.arange(min_x,max_x,precision[0])
        top_range_y = np.arange(min_z,max_z,precision[2])
    else:
        front_range_x = np.arange(min_x,max_x,precision[0])
        front_range_y = np.arange(min_z,max_z,precision[2])
        side_range_x = np.arange(min_y,max_y,precision[1])
        side_range_y = np.arange(min_z,max_z,precision[2])
        top_range_x = np.arange(min_x,max_x,precision[0])
        top_range_y = np.arange(min_y,max_y,precision[1])
        
    front_dens = np.zeros((len(front_range_x),len(front_range_y)))
    side_dens = np.zeros((len(side_range_x),len(side_range_y)))
    top_dens = np.zeros((len(top_range_x),len(top_range_y)))

    for file_name in all_SWC.keys():
        SWC = all_SWC[file_name]
        for index in SWC.keys():
            C = SWC[index]
            # front: if depth=="Z": XZ plane, else XY plane
            if depth=="Y":
                X = C[0]#x
                Y = C[1]#y
                x_ind = np.nonzero( (X>=front_range_x) & (X<=front_range_x+precision[0]))
                y_ind = np.nonzero( (Y>=front_range_y) & (Y<=front_range_y+precision[1]))
            else:
                X = C[0]#x
                Y = C[2]#z
                x_ind = np.nonzero( (X>=front_range_x) & (X<=front_range_x+precision[0]))
                y_ind = np.nonzero( (Y>=front_range_y) & (Y<=front_range_y+precision[2])) # precision 2, z-axis
            
            for x_i in x_ind:
                for y_i in y_ind:
                    front_dens[x_i,y_i] = front_dens[x_i,y_i] + 1

            # side: if depth=="Z": YZ plane, else ZY plane
            if depth=="Y":
                X = C[2]#z
                Y = C[1]#y
                x_ind = np.nonzero( (X>=side_range_x) & (X<=side_range_x+precision[2]))
                y_ind = np.nonzero( (Y>=side_range_y) & (Y<=side_range_y+precision[1]))                
            else:
                X = C[1]#y
                Y = C[2]#z  
            x_ind = np.nonzero( (X>=side_range_x) & (X<=side_range_x+precision[1]))
            y_ind = np.nonzero( (Y>=side_range_y) & (Y<=side_range_y+precision[2]))
            for x_i in x_ind:
                for y_i in y_ind:
                    side_dens[x_i,y_i] = side_dens[x_i,y_i] + 1

            # top: if depth=="Z": XY plane, else XZ plane
            if depth=="Y":
                X = C[0]#x
                Y = C[2]#z
                x_ind = np.nonzero( (X>=top_range_x) & (X<=top_range_x+precision[0]))
                y_ind = np.nonzero( (Y>=top_range_y) & (Y<=top_range_y+precision[2]))                
            else:
                X = C[0]#x
                Y = C[1]#y  
            x_ind = np.nonzero( (X>=top_range_x) & (X<=top_range_x+precision[0]))
            y_ind = np.nonzero( (Y>=top_range_y) & (Y<=top_range_y+precision[1]))
            for x_i in x_ind:
                for y_i in y_ind:
                    top_dens[x_i,y_i] = top_dens[x_i,y_i] + 1                       
    plt.figure()
    ax = plt.subplot2grid((2,2),(0,0))
    if(depth == "Y"):
        extent = [min_x,max_x,min_y,max_y] if limits == None else limits
        plt.imshow(front_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("X")
        plt.ylabel("Y")
    else:
        extent = [min_x,max_x,min_z,max_z] if limits == None else limits
        plt.imshow(front_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("X")
        plt.ylabel("Z")

    ax = plt.subplot2grid((2,2),(0,1))
    if(depth == "Y"):
        extent = [min_z,max_z,min_y,max_y] if limits == None else limits
        plt.imshow(side_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("Z")
        plt.ylabel("Y")
    else:
        extent =  [min_y,max_y,min_z,max_z] if limits == None else limits
        plt.imshow(side_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("Y")
        plt.ylabel("Z")

    ax = plt.subplot2grid((2,2),(1,1))
    if(depth == "Y"):
        extent =  [min_x,max_x,min_z,max_z] if limits == None else limits
        plt.imshow(top_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("X")
        plt.ylabel("Z")
    else:
        extent =  [min_x,max_x,min_y,max_y] if limits == None else limits
        plt.imshow(top_dens.T,origin="lower",aspect="auto",extent=extent)
        plt.xlabel("X")
        plt.ylabel("Y")

    plt.tight_layout()
        
    if not outN == None:
        plt.savefig(outN)
            
