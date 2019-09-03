"""
Generate animation of morphologies using matplotlib.

Usage is restricted to morphologies in the sWC format with the
three-point soma `standard <http://neuromorpho.org/neuroMorpho/SomaFormat.html>`_
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

import btmorph
from btmorph import config

def animate_SWC_rotation(file_name,offset=None,align=True,\
                         color_scheme="neuromorpho",filter=range(10),
                         depth="Y",out_n="rotate_animation"):
    """
    Rotate illustration of a neuronal morphology over 360 degrees.

    Parameters
    ----------
    file_n : string
        Filename of SWC description. **Must be using the 3-point soma**
        description. If the file is not in this format, use
        :py:func:`btmorph.btstructs2.STree2.read_SWC_tree_from_file` followed by
        :py:func:`btmorph.btstructs2.STree2.write_SWC_tree_to_file`, which will convert the
        description to three-point soma.
    offset : array_like
        Move the structure by the specified [x,y,z]
    color_scheme : string
        Default or "neuromorpho"
    filter_range : array_like
        List of integers indicating the SWC types to be included
        (1:soma, 2:axon, 3:basal,4:apical,...). By default, all SWC
        types are included
    depth : string
        Specify which direction (axis) represents "depth". On NeuroMorpho.org
        this is the Y-axis (and is used here by default). In NeuroMaC,
        depth is on the Z-axis.
    out_n : string
        File name of the generated gif file. The ".gif" extension is
        added automatically.
    """
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    #ax._axis3don = False
    ax.set_axis_off()

    if color_scheme == 'default':
        my_color_list = config.c_scheme_default['neurite']
    elif color_scheme == 'neuromorpho':
        my_color_list = config.c_scheme_nm['neurite']
    print 'my_color_list: ', my_color_list    

    """
    TODO?
    Can btmorph.plot_3D_SWC be integrated here? Or, is is special
    code for neuromorpho.org?
    btmorph.plot_3D_SWC(file_name=file_n)
    """
    
    x = open(file_name,'r')
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
            if not offset == None :
                x = float(splits[2])-offset[0]
                y = float(splits[3])-offset[1]
                z = float(splits[4])-offset[2]
            elif not align:
                x = float(splits[2])
                y = float(splits[3])
                z = float(splits[4])
            else:
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z 
            r = float(splits[5])
            parent = int(splits[-1])
            if n_type in filter:
                SWC[index] = (x,y,z,r,parent,n_type)

    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]
                
        if(index <= 3) :
            print 'do not draw the soma and its CNG, !!! 2 !!! point descriptions'
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            if depth=="Y": # default in neuromorpho.org files
                pl = plt.plot([p_x,c_x],[p_z,c_z],[p_y,c_y],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
                ax.set_xlabel('X')
                ax.set_ylabel('Z')
                ax.set_zlabel('Y')                
            elif depth=="Z": # used in NeuroMaC files
                pl = plt.plot([p_x,c_x],[p_y,c_y],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                
    if color_scheme=="neuromorpho":
        ax.set_axis_bgcolor(config.c_scheme_nm['bg'])
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
    
    anim = animation.FuncAnimation(fig, _animate_rotation,fargs=(ax,), frames=60)
    anim.save(out_n+".gif", writer='imagemagick', fps=4);
    #anim.save('demoanimation.mp4', writer='ffmpeg', fps=4);    

"""
Private functions below
"""

def _animate_rotation(nframe,fargs):
    fargs.view_init(elev=0, azim=nframe*6)

    
