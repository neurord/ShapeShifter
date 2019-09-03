"""
3D box counting method implementation
Contains
    BoxCounter
    
by Irina Reshodko
"""

from grids import VoxelGrid
import math

class BoxCounter:
    """
    3D box counting functionality
    """
    def __init__(self, vg):
        """
        Initialize
        
        Parameters
        -------
        vg : :class:'btmorph.btstructs2.VoxelGrid'
            Ready to use voxel grid
        """
        # Resolution is power of two
        self.vg = vg
        [dx, dy, dz] = vg.res
        if vg.res[2] == 0 or vg.res[2] == 1:
            startDim = min(vg.res[:-1])
        else:
            startDim = min(vg.res)
        nSizes = int(math.log(startDim, 2)) + 1
        self.countVals = [[]]*nSizes
        self.coverageVals = [0]*nSizes
        for i in range(0, len(self.countVals)):
            self.countVals[i] = []
            
    def grid_coverage(self, startDim, coords = None):
        """
        Box counting for coverage (standard grid method)
        Works faster than grid_count
        
        Parameters
        ----------
        startDim : int
        Dimension of the box in voxels
        coords : tuple of ints
        Grid position
        """
        if startDim < 1:
            return -1
        # Should be a power of two
        if not VoxelGrid.is_power_of_two(startDim):
            return -1
        # startDim should not be greater than smallest dimension
        if self.vg.res[2] == 0 or self.vg.res[2] == 1:
            if startDim > min(self.vg.res[:-1]):
                return -1
        elif startDim > min(self.vg.res):
            return -1
        if coords != None:
            if 1 == startDim:
                if True == self.vg[coords]:
                    return True
                else:
                    return False
            else:
                m = int(math.log(startDim, 2))
                newDim = startDim/2
                new_c = [coords[0]*2, coords[1]*2, coords[2]*2]
                if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                    new_c[2] = coords[2]
                tmp = False
                for i in [new_c[0], new_c[0] + 1]:
                    for j in [new_c[1], new_c[1] + 1]:
                        if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                            k = new_c[2]
                            tmp = self.grid_coverage(newDim, (i,j,k)) or tmp
                        else:
                            for k in [new_c[2], new_c[2] + 1]:
                                tmp = self.grid_coverage(newDim, (i,j,k)) or tmp
                self.coverageVals[m] += tmp
                return tmp
        else:
            dx = self.vg.res[0]/startDim
            dy = self.vg.res[1]/startDim            
            if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                    dz = 1
            else:
                dz = self.vg.res[2]/startDim
            for i in range(0, dx):
                for j in range(0, dy):
                    for k in range(0, dz):                        
                        self.grid_coverage(startDim, (i, j, k))

                
        
    def grid_count(self, startDim, coords = None):
        """
        Box counting for grid method (standard)
        
        Parameters
        ----------
        startDim : int
        Dimension of the box in voxels
        coords : tuple of ints
        Grid position
        """
        if startDim < 1:
            return -1
        # Should be a power of two
        if not VoxelGrid.is_power_of_two(startDim):
            return -1
        # startDim should not be greater than smallest dimension
        if self.vg.res[2] == 0 or self.vg.res[2] == 1:
            if startDim > min(self.vg.res[:-1]):
                return -1
        elif startDim > min(self.vg.res):
            return -1
        encBX = self.vg.encompassingBox[0]
        encBY = self.vg.encompassingBox[1]
        encBZ = self.vg.encompassingBox[2]
        if coords != None:
            if 1 == startDim:
                if True == self.vg[coords]:
                    return 1
                else:
                    return 0
            else:
                m = int(math.log(startDim, 2))
                newDim = startDim/2
                new_c = [coords[0]*2, coords[1]*2, coords[2]*2]
                rangeX = [math.trunc(encBX[0]/float(newDim)), int(math.ceil(encBX[1]/float(newDim)))]
                rangeY = [math.trunc(encBY[0]/float(newDim)), int(math.ceil(encBY[1]/float(newDim)))]
                rangeZ = [math.trunc(encBZ[0]/float(newDim)), int(math.ceil(encBZ[1]/float(newDim)))]
                if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                    new_c[2] = coords[2]
                    rangeZ = [0, 1]
                s = 0
                for i in [new_c[0], new_c[0]+1]:#[max(new_c[0], rangeX[0]), min(new_c[0] + 1, rangeX[1])]:
                    for j in [new_c[1], new_c[1]+1]:#[max(new_c[1], rangeY[0]), min(new_c[1] + 1, rangeY[1])]:
                        if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                            k = new_c[2]
                            s += self.grid_count(newDim, (i,j,k))
                        else:
                            for k in [new_c[2], new_c[2]+1]:#[max(new_c[2], rangeZ[0]), min(new_c[2] + 1, rangeZ[1])]:
                                s += self.grid_count(newDim, (i,j,k))
                if s > 0 :
                    self.countVals[m].append(s)
                return s
        else:            
            #dx = self.vg.res[0]/startDim
            #dy = self.vg.res[1]/startDim
            rangeX = [math.trunc(encBX[0]/float(startDim)), int(math.ceil(encBX[1]/float(startDim)))]
            rangeY = [math.trunc(encBY[0]/float(startDim)), int(math.ceil(encBY[1]/float(startDim)))]
            if self.vg.res[2] == 0 or self.vg.res[2] == 1:
                    #dz = 1
                    rangeZ = [0, 1]
            else:
                #dz = self.vg.res[2]/startDim
                rangeZ = [math.trunc(encBZ[0]/float(startDim)), int(math.ceil(encBZ[1]/float(startDim)))]
            for i in range(rangeX[0], rangeX[1]):
                for j in range(rangeY[0], rangeY[1]):
                    for k in range(rangeZ[0], rangeZ[1]):                        
                        self.grid_count(startDim, (i, j, k))
        
