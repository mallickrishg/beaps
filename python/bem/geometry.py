""" Constructing unit vectors and other geoemtric parameters for a sequence of points
Written by Rishav Mallick, Caltech Seismolab, 2023
"""
import numpy as np

class geometry:
    def __init__(self, xnode, znode, G = 30e3, nu = 0.25):
        # xnode is provided with east/right side being positive
        # znode is provided with up being positive
        self.G = G
        self.nu = nu

        self.xn = xnode
        self.zn = znode
        self.N = len(xnode)-1
        self.W = np.sqrt((self.xn[:-1]-self.xn[1:])**2 + (self.zn[:-1]-self.zn[1:])**2)

        self.x,self.z = self.compute_midpoints()
        # compute dip angles
        self.dip = -np.degrees(np.arctan2(-znode[:-1]+znode[1:], -xnode[:-1]+xnode[1:]))
        # compute unit vectors
        self.dv,self.nv = self.compute_unitvectors()

    def compute_midpoints(self):
        x = (self.xn[:-1]+self.xn[1:])/2
        z = (self.zn[:-1]+self.zn[1:])/2
        return x, z
    
    def compute_unitvectors(self):
        # Compute the unit slip vector dv
        dv = np.array([-np.cos(np.radians(self.dip)), np.sin(np.radians(self.dip))])
        # Compute unit normal vector nv
        nv = np.array([np.sin(np.radians(self.dip)), np.cos(np.radians(self.dip))])
        return dv.T,nv.T
    
class observations:
    def __init__(self,x,z):
        self.x = x
        self.z = z
        self.N = len(x)
        self.nv = np.zeros((self.N,2))