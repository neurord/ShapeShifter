"""
File contains:

    - :class:`P3D2`
    - :class:`SNode2`
    - :class:`STree2`
    
B. Torben-Nielsen (from legacy code). Daniele Linaro contributed the
iterators in  :class:`STree2`.
"""
import numpy as np
import matplotlib.pyplot as plt
import math

class P3D2(object) :
    """
    Basic container to represent and store 3D information
    """
    def __init__(self,xyz,radius,type=7) :
        """ Constructor.

        Parameters
        -----------

        xyz : numpy.array
            3D location
        radius : float
        type : int
            Type asscoiated with the segment according to SWC standards
        """
        self.xyz = xyz
        self.radius = radius
        self.type = type

    def __str__(self) :
        return "P3D2 [%.2f %.2f %.2f], R=%.2f" % (self.xyz[0],self.xyz[1],self.xyz[2],self.radius)

class SNode2(object) :
    """
    Simple Node for use with a simple Tree (STree)
    
    By design, the "content" should be a dictionary. (2013-03-08)
    """
    
    def __init__(self,index) :
        """
        Constructor.

        Parameters
        -----------
        index : int
           Index, unique name of the :class:`SNode2`
        """
        self.parent = None
        self.index = index
        self.children = []
        self.content = {}

    def get_parent(self) :
        """
        Return the parent node of this one.

        Returns
        -------
        parent : :class:`SNode2`
           In case of the root, None is returned.Otherwise a :class:`SNode2` is returned
        """          
        return self.__parent

    def set_parent(self,parent) :
        """
        Set the parent node of a given other node

        Parameters
        ----------
        node : :class:`SNode2`
        """
        self.__parent = parent

    parent = property(get_parent, set_parent)
    
    def get_index(self) :
        """
        Return the index of this node

        Returns
        -------
        index : int
        """
        return self.__index
        
    def set_index(self,index) :
        """
        Set the unqiue name of a node

        Parameters
        ----------

        index : int
        """
        self.__index = index
    
    index = property(get_index, set_index)

    def get_children(self) :
        """
        Return the children nodes of this one (if any)

        Returns
        -------
        children : list :class:`SNode2`
           In case of a leaf an empty list is returned
        """                  
        return self.__children
        
    def set_children(self, children) :
        """
        Set the children nodes of this one

        Parameters
        ----------

        children: list :class:`SNode2`
        """
        self.__children = children

    children = property(get_children, set_children)

    def get_content(self) :
        """
        Return the content dict of a :class:`SNode2`

        Returns
        -------
        parent : :class:`SNode2`
           In case of the root, None is returned.Otherwise a :class:`SNode2` is returned
        """                  
        return self.__content
    
    def set_content(self,content) :
        """
        Set the content of a node. The content must be a dict

        Parameters
        ----------
        content : dict
            dict with content. For use in btmorph at least a 'p3d' entry should be present
        """        
        if isinstance(content,dict) :
            self.__content = content 
        else :
            raise Exception("SNode2.set_content must receive a dict")    

    content = property(get_content, set_content)

    def add_child(self,child_node) :
        """
        add a child to the children list of a given node

        Parameters
        -----------
        node :  :class:`SNode2`
        """
        self.children.append(child_node)
            
    def make_empty(self):
        """
        Clear the node. Unclear why I ever implemented this. Probably to cover up some failed garbage collection
        """
        self.parent = None
        self.content = {}
        self.children = []
            
    def remove_child(self, child) :
        """
        Remove a child node from the list of children of a specific node

        Parameters
        -----------
        node :  :class:`SNode2`
            If the child doesn't exist, you get into problems.
        """
        self.children.remove(child)
        
    def __str__(self) :
        return 'SNode2 (ID: '+str(self.index)+')'

    def __lt__(self,other):
        if self.index < other.index :
            return True
    def __le__(self,other):
        if self.index <= other.index :
            return True
    def __gt__(self,other):
        if self.index > other.index :
            return True
    def __ge__(self,other):
        if self.index >= other.index :
            return True
    
    def __copy__(self) : # customization of copy.copy
        ret = SNode2(self.index)
        for child in self.children :
            ret.add_child(child)
        ret.content = self.content
        ret.parent = self.parent
        return ret

class STree2(object) :
    '''
    Simple tree for use with a simple Node (:class:`SNode2`).

    While the class is designed to contain binary trees (for neuronal morphologies) the number of children is not limited.
    As such, this is a generic implementation of a tree structure as a linked list.
    '''
    
    def __init__(self):
        """
        Default constructor. No arguments are passed.
        """
        self.root = None

    def __iter__(self):
        nodes = []
        self._gather_nodes(self.root,nodes)
        for n in nodes:
            yield n
    
    def __getitem__(self,index):
        return self._find_node(self.root,index)

    def set_root(self,node) :
        """
        Set the root node of the tree

        Parameters
        -----------
        node : :class:`SNode2`
            to-be-root node
        """
        if not node is None:
            node.parent = None
        self.__root = node
    def get_root(self) :
        """
        Obtain the root node

        Returns
        -------
        root : :class:`SNode2`
        """
        return self.__root
    root = property(get_root, set_root)
    
    def is_root(self,node) :
        """
        Check whether a node is the root node

        Returns
        --------
        is_root : boolean
            True is the queried node is the root, False otherwise
        """
        if node.parent is None:
            return True
        else :
            return False
            
    def is_leaf(self,node) :
        """
        Check whether a node is a leaf node, i.e., a node without children

        Returns
        --------
        is_leaf : boolean
            True is the queried node is a leaf, False otherwise
        """
        if len(node.children) == 0  :
            return True
        else :
            return False
            
    def add_node_with_parent(self,node,parent) :
        """
        Add a node to the tree under a specific parent node

        Parameters
        -----------
        node : :class:`SNode2`
            node to be added
        parent : :class:`SNode2`
            parent node of the newly added node
        """
        node.parent = parent
        if not parent is None:
            parent.add_child(node)
        
    def remove_node(self,node) :
        """
        Remove a node from the tree

        Parameters
        -----------
        node : :class:`SNode2`
            node to be removed
        """
        node.parent.remove_child(node)
        self._deep_remove(node)
        
            
    def _deep_remove(self,node) :
        children = node.children
        node.make_empty()
        for child in children :
            self._deep_remove(child)        

    def get_nodes(self) :
        """
        Obtain a list of all nodes int the tree

        Returns
        -------
        all_nodes : list of :class:`SNode2`
        """
        n = []
        self._gather_nodes(self.root,n) 
        return n

    def get_sub_tree(self,fake_root) :
        """
        Obtain the subtree starting from the given node

        Parameters
        -----------
        fake_root : :class:`SNode2`
            Node which becomes the new root of the subtree

        Returns
        -------
        sub_tree :  STree2
            New tree with the node from the first argument as root node
        """
        ret = STree2()
        cp = fake_root.__copy__()
        cp.parent = None
        ret.root = cp
        return ret

    def _gather_nodes(self,node,node_list) :
        if not node is None:
            node_list.append(node)
            for child in node.children :
                self._gather_nodes(child,node_list)
    
    def get_node_with_index(self, index) :
        """
        Get a node with a specific name. The name is always an integer

        Parameters
        ----------
        index : int
            Name of the node to be found

        Returns
        -------
        node : :class:`SNode2`
            Node with the specific index
        """
        return self._find_node(self.root,index)
        
    def get_node_in_subtree(self,index,fake_root) :
        """
        Get a node with a specific name in a the subtree rooted at fake_root. The name is always an integer

        Parameters
        ----------
        index : int
            Name of the node to be found
        fake_root: :class:`SNode2`
            Root node of the subtree in which the node with a given index is searched for

        Returns
        -------
        node : :class:`SNode2`
            Node with the specific index
        """
        return self._find_node(fake_root,index)
        
    def _find_node(self,node,index) :
        """
        Sweet breadth-first/stack iteration to replace the recursive call. 
        Traverses the tree until it finds the node you are looking for.

        Parameters
        -----------

        
        Returns
        -------
        node : :class:`SNode2`
            when found and None when not found
        """
        stack = []; 
        stack.append(node)
        while(len(stack) != 0) :
            for child in stack :
                if child.index == index  :
                    return child
                else :
                    stack.remove(child)
                    for cchild in child.children :
                        stack.append(cchild)
        return None # Not found!
        
    def degree_of_node(self,node) :
        """
        Get the degree of a given node. The degree is defined as the number of leaf nodes in the subtree rooted at this node.

        Parameters
        ----------
        node : :class:`SNode2`
            Node of which the degree is to be computed.

        Returns
        -------
        degree : int
        """
        sub_tree = self.get_sub_tree(node)
        st_nodes = sub_tree.get_nodes()
        leafs = 0
        for n in st_nodes :
            if sub_tree.is_leaf(n) :
                leafs = leafs +1
        return leafs
        
    def order_of_node(self,node) :
        """
        Get the order of a given node. The order or centrifugal order is defined as 0 for the root and increased with any bifurcation.
        Hence, a node with 2 branch points on the shortest path between that node and the root has order 2.

        Parameters
        ----------
        node : :class:`SNode2`
            Node of which the order is to be computed.

        Returns
        -------
        order : int
        """
        ptr = self.path_to_root(node)
        order = 0
        for n in ptr :
            if len(n.children) > 1  :
                order = order+1
        # order is on [0,max_order] thus subtract 1 from this calculation
        return order-1 
                
    def path_to_root(self,node) :
        """
        Find and return the path between a node and the root.

        Parameters
        ----------
        node : :class:`SNode2`
            Node at which the path starts

        Returns
        -------
        path : list of :class:`SNode2`
            list of :class:`SNode2` with the provided node and the root as first and last entry, respectively.
        """
        n = []
        self._go_up_from(node,n)            
        return n
        
    def _go_up_from(self,node,n):
        n.append(node)
        if not node.parent is None :
            self._go_up_from(node.parent,n)

    def path_between_nodes(self,from_node,to_node) :
        """
        Find the path between two nodes. The from_node needs to be of higher \
        order than the to_node. In case there is no path between the nodes, \
        the path from the from_node to the soma is given.

        Parameters
        -----------
        from_node : :class:`SNode2`
        to_node : :class:`SNode2`
        """
        n = []
        self._go_up_from_until(from_node,to_node,n)
        return n

    def _go_up_from_until(self,from_node,to_node,n) :
        n.append(from_node)
        if from_node == to_node :
            return
        if not from_node.parent is None :
            self._go_up_from_until(from_node.parent,to_node,n)

    def write_SWC_tree_to_file(self,file_n) :
        """
        Non-specific for a tree.

        Used to write an SWC file from a morphology stored in this
        :class:`STree2`. Output uses the 3-point soma standard.

        Parameters
        ----------
        file_n : str
            name of the file to open
        """
        writer = open(file_n,'w')
        nodes = self.get_nodes()
        nodes.sort()

        # 3 point soma representation (See Neuromoprho.org FAQ)
        s1p = nodes[0].content["p3d"]
        s1_xyz = s1p.xyz
        s2p = nodes[1].content["p3d"]
        s2_xyz = s2p.xyz
        s3p = nodes[2].content["p3d"]
        s3_xyz = s3p.xyz
        soma_str = "1 1 " +str(s1_xyz[0]) + " " + str(s1_xyz[1]) + \
          " " + str(s1_xyz[2]) + " " + str(s1p.radius) + " -1\n" + \
          "2 1 " +str(s2_xyz[0]) + " " + str(s2_xyz[1]) + \
          " " + str(s2_xyz[2]) + " " + str(s2p.radius) + " 1\n" + \
          "3 1 " +str(s3_xyz[0]) + " " + str(s3_xyz[1]) + \
          " " + str(s3_xyz[2]) + " " + str(s3p.radius) + " 1\n"
        writer.write(soma_str)
        writer.flush()
        
        # add the soma compartment, then enter the loop
        for node in nodes[3:] :
            p3d = node.content['p3d'] # update 2013-03-08
            xyz = p3d.xyz
            radius = p3d.radius
            tt = p3d.type
            p3d_string = str(node.index)+' '+str(tt) + ' ' + \
                str(xyz[0]) + ' ' + str(xyz[1])+ ' ' + str(xyz[2]) + \
                ' ' + str(radius) + ' ' \
                + str(node.parent.index)
            # print 'p3d_string: ', p3d_string
            writer.write( p3d_string + '\n' )
            writer.flush()
        writer.close()          
        #print 'STree::writeSWCTreeToFile -> finished. Tree in >',fileN,'<'
        
    def read_SWC_tree_from_file(self,file_n,types=range(1,10)) :
        """
        Non-specific for a "tree data structure"
        Read and load a morphology from an SWC file and parse it into
        an STree2 object. 

        On the NeuroMorpho.org website, 5 types of somadescriptions are 
        considered (http://neuromorpho.org/neuroMorpho/SomaFormat.html).
        The "3-point soma" is the standard and most files are converted
        to this format during a curation step. btmorph follows this default
        specificationand the *internal structure of btmorph implements
        the 3-point soma*.

        However, two other options to describe the soma
        are still allowed and available, namely:
        - soma absent: btmorph adds a 3-point soma in between of [TO DEFINE/TODO]
        - multiple cylinder: [TO DEFINE/TODO]
        
        Parameters
        -----------
        file_n : str
            name of the file to open
        """
        # check soma-representation: 3-point soma or a non-standard representation
        soma_type = self._determine_soma_type(file_n)
        #print "STree2::read_SWC_tree_from_file found soma_type=%i" % soma_type
        
        file = open(file_n,'r')
        all_nodes = dict()
        for line in file :
            if not line.startswith('#') :
                split = line.split()
                index = int(split[0].rstrip())
                swc_type = int(split[1].rstrip())
                x = float(split[2].rstrip())
                y = float(split[3].rstrip())
                z = float(split[4].rstrip())
                radius = float(split[5].rstrip())
                parent_index = int(split[6].rstrip())

                if swc_type in types:
                    tP3D = P3D2(np.array([x,y,z]),radius,swc_type)
                    t_node = SNode2(index)
                    t_node.content = {'p3d': tP3D}
                    all_nodes[index] = (swc_type,t_node,parent_index)
                else:
                    print type,index

        #print "len(all_nodes): ", len(all_nodes)

        # IF 3-point soma representation
        if soma_type == 1:
            for index,(swc_type,node,parent_index) in all_nodes.items() :
                if index == 1:
                    self.root = node
                elif index in (2,3):
                    # the 3-point soma representation (http://neuromorpho.org/neuroMorpho/SomaFormat.html)
                    self.add_node_with_parent(node,self.root)
                else:
                    parent_node = all_nodes[parent_index][1]
                    self.add_node_with_parent(node,parent_node)
        # IF multiple cylinder soma representation
        elif soma_type ==2:
            self.root = all_nodes[1][1]
            
            # get all some info
            soma_cylinders = []
            connected_to_root = []
            for index,(swc_type,node,parent_index) in all_nodes.items() :
                if swc_type == 1 and not index == 1:
                    soma_cylinders.append((node,parent_index))
                    if index > 1 :
                        connected_to_root.append(index)

            # make soma
            s_node_1, s_node_2 = self._make_soma_from_cylinders(soma_cylinders,all_nodes)
            
            # add soma
            self.root = all_nodes[1][1]
            self.root.content["p3d"].radius = s_node_1.content["p3d"].radius
            self.add_node_with_parent(s_node_1,self.root)
            self.add_node_with_parent(s_node_2,self.root)

            # add the other points            
            for index,(swc_type,node,parent_index) in all_nodes.items() :
                if swc_type == 1:
                    pass
                else:
                    parent_node = all_nodes[parent_index][1]
                    if parent_node.index in connected_to_root:
                        self.add_node_with_parent(node,self.root)
                    else:
                        self.add_node_with_parent(node,parent_node)
                    
        return self

    def _make_soma_from_cylinders(self,soma_cylinders,all_nodes):
        """Now construct 3-point soma
        step 1: calculate surface of all cylinders
        step 2: make 3-point representation with the same surface"""
        
        total_surf = 0
        for (node,parent_index) in soma_cylinders:
            n = node.content["p3d"]
            p = all_nodes[parent_index][1].content["p3d"]
            H = np.sqrt(np.sum((n.xyz-p.xyz)**2))
            surf = 2*np.pi*p.radius*H
            #print "(node %i) surf as cylinder:  %f (R=%f, H=%f), P=%s" % (node.index,surf,n.radius,H,p)
            total_surf = total_surf+surf
        print "found 'multiple cylinder soma' w/ total soma surface=", total_surf

        # define apropriate radius
        radius=np.sqrt(total_surf/(4*np.pi))
        #print "found radius: ", radius

        s_node_1 = SNode2(2)
        r = self.root.content["p3d"]
        rp = r.xyz
        s_p_1 = P3D2(np.array([rp[0],rp[1]-radius,rp[2]]),radius,1)
        s_node_1.content = {'p3d': s_p_1}
        s_node_2 = SNode2(3)
        s_p_2 = P3D2(np.array([rp[0],rp[1]+radius,rp[2]]),radius,1)
        s_node_2.content = {'p3d': s_p_2}
        
        return s_node_1,s_node_2
    
    def _determine_soma_type(self,file_n):
        """
        Costly method to determine the soma type used in the SWC file.
        This method searches the whole file for soma entries.  

        Parameters
        ----------
        file_n : string
            Name of the file containing the SWC description

        Returns
        -------
        soma_type : int
            Integer indicating one of the su[pported SWC soma formats.
            1: Default three-point soma, 2: multiple cylinder description,
            3: otherwise [not suported in btmorph]
        """
        file = open(file_n,"r")
        somas = 0
        for line in file:
            if not line.startswith('#') :
                split = line.split()
                index = int(split[0].rstrip())
                s_type = int(split[1].rstrip())
                if s_type == 1 :
                    somas = somas +1
        file.close()
        if somas == 3:
            return 1
        elif somas < 3:
            return 3
        else:
            return 2
        

    def __str__(self) :
        return "STree2 ("+str(len(self.get_nodes()))+" nodes)"
