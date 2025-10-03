import numpy as np


def grid(N:int, chord:float = 1, epsilon:float = 0.1):
    epsilon = epsilon*chord
    dx = chord/N
    I = np.arange(1, N+1)
    
    # colocation point
    xc  = chord/N *(I - 0.25)    
    zc  = 4* epsilon* xc / chord *(1 - xc/chord)
    
    # Vortex point
    x  = chord/N *(I - 0.75)    
    z  = 4* epsilon* x / chord* (1 - x/chord)
    
    # Normal at colocation point
    detadx = 4* epsilon/ chord* (1 - 2*xc/chord)
    
    enx = -detadx/np.sqrt(1+detadx**2) 
    enz = 1/np.sqrt(1+detadx**2)
    
    cp = np.array([xc, zc]).T
    vp =  np.array([x, z]).T
    n = np.array([enx, enz]).T    
    return cp, vp, n



def vor2D(p:np.ndarray, p_j:np.ndarray, Gamma:float, tol:float = 1e-3):
    x, z = p
    xj, zj = p_j
    rx = (x - xj)
    rz = (z - zj)
    r = np.sqrt(rx**2  + (z - zj)**2)
    U, W = 0, 0
    if r > tol:
        V = 0.5/np.pi*Gamma/r
        U = V*(rz/r)
        W = V*(-rx/r)
    return U, W

def influCoeff(cp:np.ndarray, vp:np.ndarray, n:np.ndarray):
    N = cp.shape(0)
    U = np.zeros((N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            u, w = vor2D(cp[i], vp[j], 1)
            U[i, j] = u
            W[i, j] = w
    