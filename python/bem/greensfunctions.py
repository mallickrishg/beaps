""" Definition of displacement and stress greens functions for 2-d problems
Written by Rishav Mallick, Caltech Seismolab, 2023
"""
import numpy as np

# Half space greens functions for anti-plane geometry
# Stress for fault slip
def sigmahs_antiplane(x2, x3, y2, y3, Wf, dip, G):
    s12 = G * ( \
        -(x3-y3) / ((x2-y2)**2 + (x3-y3)**2) + (x3+y3) / ((x2-y2)**2 + (x3+y3)**2) \
        + (x3-y3-Wf*np.sin(np.deg2rad(dip))) / ((x2-y2-Wf*np.cos(np.deg2rad(dip)))**2 + (x3-y3-Wf*np.sin(np.deg2rad(dip)))**2) - (x3+y3+Wf*np.sin(np.deg2rad(dip))) / ((x2-y2-Wf*np.cos(np.deg2rad(dip)))**2 + (x3+y3+Wf*np.sin(np.deg2rad(dip)))**2) \
        ) / (2 * np.pi)
    
    s13 = G * ( \
        (x2-y2) / ((x2-y2)**2 + (x3-y3)**2) - (x2-y2) / ((x2-y2)**2 + (x3+y3)**2) \
        - (x2-y2-Wf*np.cos(np.deg2rad(dip))) / ((x2-y2-Wf*np.cos(np.deg2rad(dip)))**2 + (x3-y3-Wf*np.sin(np.deg2rad(dip)))**2) + (x2-y2-Wf*np.cos(np.deg2rad(dip))) / ((x2-y2-Wf*np.cos(np.deg2rad(dip)))**2 + (x3+y3+Wf*np.sin(np.deg2rad(dip)))**2) \
        ) / (2 * np.pi)
    
    return s12, s13

# Displacements for fault slip
def uhs_antiplane(x2, x3, y2, y3, W, dip):
    u1_1 = (np.arctan2((x3-y3), (x2-y2)) - np.arctan2((x3-y3-W*np.sin(np.deg2rad(dip))), (x2-y2-W*np.cos(np.deg2rad(dip))))) / (2 * np.pi)
    u1_2 = (-np.arctan2((x3+y3), (x2-y2)) + np.arctan2((x3+y3+W*np.sin(np.deg2rad(dip))), (x2-y2-W*np.cos(np.deg2rad(dip))))) / (2 * np.pi)
    
    u1_1[u1_1 > 0.5] = u1_1[u1_1 > 0.5] - 1
    
    u1 = u1_1 + u1_2
    return u1

# Full space greens functions for plane strain geometry
def stressfs_plane(x, y, xe, ye, a, dip, Ds, Dn, nu, E):
    # Arguments: (input)
    #   x & y  - The observation points locations in real Cartesian coords.  
    #  xe & ye - The element midpoint location in real coords. 
    #    a     - The elements half length
    #    dip   - dip angle in degrees
    #  Dn,Ds   - The defined displacement of each element.(normal and shear)
            #    Dn+ is opening, Ds+ is left lateral shearing. 
    #    nu    - The Poisson's ratio
    #    E     - The Young's modulus E = 2*G*nu/(1-nu)
    #  Arguments: (output)
    #  Stress - Is the stress caused by the movement of the dislocation at the observataion points. [Sxx,Syy,Sxy].
    
    Beta = -np.deg2rad(dip)
    sm = E / (2 * (1 + nu))
    con = 1 / (4 * np.pi * (1 - nu))
    cons = 2 * sm
    H = a
    Dxb = Ds
    Dyb = -Dn
    sb = np.sin(Beta)
    cb = np.cos(Beta)
    s2b = np.sin(2 * Beta)
    c2b = np.cos(2 * Beta)
    XB = (x - xe) * cb + (y - ye) * sb
    YB = -(x - xe) * sb + (y - ye) * cb

    Y2 = YB ** 2
    XMH = XB - H
    XPH = XB + H
    XMH2 = XMH ** 2
    XPH2 = XPH ** 2
    R1S = XMH2 + Y2
    R1S2 = R1S ** 2
    R2S = XPH2 + Y2
    R2S2 = R2S ** 2

    FF4 = con * (YB / R1S - YB / R2S)
    FF5 = con * (XMH / R1S - XPH / R2S)
    FF6 = con * ((XMH2 - Y2) / R1S2 - (XPH2 - Y2) / R2S2)
    FF7 = 2 * con * YB * (XMH / R1S2 - XPH / R2S2)

    Sxx = cons * Dxb * (2 * (cb * cb) * FF4 + s2b * FF5 + YB * (c2b * FF6 - s2b * FF7)) + cons * Dyb * (-FF5 + YB * (s2b * FF6 + c2b * FF7))
    Syy = cons * Dxb * (2 * (sb * sb) * FF4 - s2b * FF5 - YB * (c2b * FF6 - s2b * FF7)) + cons * Dyb * (-FF5 - YB * (s2b * FF6 + c2b * FF7))
    Sxy = cons * Dxb * (s2b * FF4 - c2b * FF5 + YB * (s2b * FF6 + c2b * FF7)) + cons * Dyb * (-YB * (c2b * FF6 - s2b * FF7))
    # Stress = np.column_stack((Sxx.ravel(), Syy.ravel(), Sxy.ravel()))
    return Sxx,Syy,Sxy

# Displacements for fault slip
def ufs_plane(x, y, xe, ye, a, dip, Ds, Dn, nu):
    # Arguments: (input)
    #   x & y  - The observation points locations in real Cartesian coords.  
    #  xe & ye - The element midpoint location in real coords. 
    #    a     - The elements half length
    #    dip   - dip angle in degrees
    #  Dn,Ds   - The defined displacement of each element.(normal and shear)
            #    Dn+ is opening, Ds+ is left lateral shearing. 
    #    nu    - The Poisson's ratio
    
    length = x.shape
    lengthrow = length[0]
    if x.ndim==1:
        lengthcol = 1
    else:
        lengthrow = length[0]
        lengthcol = length[1]

    Beta = -np.deg2rad(dip)

    # Define material constant used in calculating influence coefficients.
    con = 1 / (4 * np.pi * (1 - nu))
    H = a
    Dxb = Ds
    Dyb = -Dn
    sb = np.sin(Beta)
    cb = np.cos(Beta)
    # Define array of local coordinates for the observation grid relative to the midpoint and orientation of the ith element.
    XB = (x - xe) * cb + (y - ye) * sb
    YB = -(x - xe) * sb + (y - ye) * cb

    # Flag all observation points for which Yb is not 0  
    i1 = np.where(YB)
    i2 = np.where((YB == 0) & (np.abs(XB) < H))
    i3 = np.where((YB == 0) & (np.abs(XB) > H))

    # Calculate derivatives of the function f(x,y)
    Y2 = np.power(YB, 2)
    XMH = XB - H
    XPH = XB + H
    XMH2 = np.power(XMH, 2)
    XPH2 = np.power(XPH, 2)
    R1S = XMH2 + Y2
    R2S = XPH2 + Y2
    
    FF2 = con * (np.log(np.sqrt(R1S)) - np.log(np.sqrt(R2S)))
    # Steve Martels Solution to elements lying on same plane
    # FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
    # FB3 = difference of arc tangents for all other pts.
    FF3 = np.zeros_like(YB)
    FF3[i1] = np.arctan2(YB[i1], XMH[i1]) - np.arctan2(YB[i1], XPH[i1])
    FF3[i2] = np.pi#*np.ones_like(i2)
    FF3[i3] = 0#np.zeros_like(i3)
    
    FF3 = -con*FF3
    # FF3 = FF3.reshape(lengthrow, lengthcol)

    FF4 = con*(YB/R1S - YB/R2S)
    FF5 = con*(XMH/R1S - XPH/R2S)

    pr1 = 1 - 2 * nu
    pr2 = 2 - 2 * nu
    Ux = Dxb * (-pr1 * sb * FF2 + pr2 * cb * FF3 + YB * (sb * FF4 - cb * FF5)) + \
         Dyb * (-pr1 * cb * FF2 - pr2 * sb * FF3 - YB * (cb * FF4 + sb * FF5))
    Uy = Dxb * (+pr1 * cb * FF2 + pr2 * sb * FF3 - YB * (cb * FF4 + sb * FF5)) + \
         Dyb * (-pr1 * sb * FF2 + pr2 * cb * FF3 - YB * (sb * FF4 - cb * FF5))
    # Disp = np.hstack((Ux.reshape(-1, 1), Uy.reshape(-1, 1)))
    return Ux,Uy  

def compute_tractionkernel(src,rcv):
    # Compute Stress components for each source at the locations provided and then project it to shear and normal direction for the surface
    # src - structure that contains source geometry and parameters
    #       x,z
    #       W,dip
    #       dv,nv (unit vectors)
    #       G,nu - elastic material parameters
    # rcv - structure that contains location and geometry of observation points
    #       x,z,dv,nv,N
    # returns 4 traction kernels

    
    # kernels from shear source
    Kts_s = np.zeros((rcv.N, src.N))
    Ktn_s = np.zeros((rcv.N, src.N))

    # kernels from normal source
    Kts_n = np.zeros((rcv.N, src.N))
    Ktn_n = np.zeros((rcv.N, src.N))

    for i in range(src.N):        
        # compute tractions for shear source
        Sxx, Szz, Sxz = stressfs_plane(rcv.x, rcv.z, src.x[i], src.z[i], src.W[i]/2, src.dip[i], 1, 0, src.nu, 2*src.G*(1+src.nu))

        t = np.vstack((Sxx * rcv.nv[:, 0] + Sxz * rcv.nv[:, 1],
                    Sxz * rcv.nv[:, 0] + Szz * rcv.nv[:, 1])).T

        Kts_s[:, i] = rcv.dv[:, 0] * t[:, 0] + rcv.dv[:, 1] * t[:, 1]
        Ktn_s[:, i] = rcv.nv[:, 0] * t[:, 0] + rcv.nv[:, 1] * t[:, 1]

        # compute tractions for normal source
        Sxx, Szz, Sxz = stressfs_plane(rcv.x, rcv.z, src.x[i], src.z[i], src.W[i]/2, src.dip[i], 0, 1, src.nu, 2*src.G*(1+src.nu))

        t = np.vstack((Sxx * rcv.nv[:, 0] + Sxz * rcv.nv[:, 1],
                    Sxz * rcv.nv[:, 0] + Szz * rcv.nv[:, 1])).T

        Kts_n[:, i] = rcv.dv[:, 0] * t[:, 0] + rcv.dv[:, 1] * t[:, 1]
        Ktn_n[:, i] = rcv.nv[:, 0] * t[:, 0] + rcv.nv[:, 1] * t[:, 1]

    return Kts_s,Ktn_s,Kts_n,Ktn_n

def compute_displacementkernel(src,rcv):
    # Compute Displacement components for each source at the locations provided
    # src - structure that contains source geometry and parameters
    #       x,z
    #       W,dip
    #       dv,nv (unit vectors)
    #       G,nu - elastic material parameters
    # rcv - structure that contains location and geometry of observation points
    #       x,z,N
    # returns 4 displacement kernels

    delta = 1e-6# small shift for self-displacement in m

    # kernels from shear source
    Gx_s = np.zeros((rcv.N, src.N))
    Gz_s = np.zeros((rcv.N, src.N))

    # kernels from normal source
    Gx_n = np.zeros((rcv.N, src.N))
    Gz_n = np.zeros((rcv.N, src.N))

    for i in range(src.N):        
        # compute tractions for shear source
        x = rcv.x + 0
        z = rcv.z + 0
        r = np.sqrt((x-src.x[i])**2 + (z-src.z[i])**2)
        x[r==0] = x[r==0] - rcv.nv[r==0,0]*delta
        z[r==0] = z[r==0] - rcv.nv[r==0,1]*delta
        
        Ux,Uz = ufs_plane(x, z, src.x[i], src.z[i], src.W[i]/2, src.dip[i], 1, 0, src.nu)
        # need to add a check for self-displacements. Move the (x,z) by an infinitesimal amount in the nv direction
        Gx_s[:, i] = Ux
        Gz_s[:, i] = Uz

        # compute tractions for normal source
        Ux,Uz = ufs_plane(x, z, src.x[i], src.z[i], src.W[i]/2, src.dip[i], 0, 1, src.nu)

        Gx_n[:, i] = Ux
        Gz_n[:, i] = Uz

    return Gx_s,Gz_s,Gx_n,Gz_n