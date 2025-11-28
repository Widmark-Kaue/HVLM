import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

def plot_mesh(X, Y, Z, vcp):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    
    for i in range(X.shape[0]):
        ax.plot(X[i,:], Y[i,:], Z[i,:], 'k-')
    for j in range(X.shape[1]):
        ax.plot(X[:,j], Y[:,j], Z[:,j], 'k-')
        
    for k in range(vcp.shape[0]):
        # Plot center, quarter chord and three quarter chord points
        ax.plot(vcp[k, 0, 0], vcp[k, 0, 1], vcp[k, 0, 2], 'ko', markersize = 4)
        ax.plot(vcp[k, 1, 0], vcp[k, 1, 1], vcp[k, 1, 2], 'ro', markersize = 4)
        # ax.plot(vcp[k, 2, 0], vcp[k, 2, 1], vcp[k, 2, 2], 'bo', markersize = 4)
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect([1, 1, 0.1])
    plt.axis('equal')
    plt.show(block = False)

def meshPlanar(nspan:int, nchord:int, geometry:np.ndarray):
    # geometry: (N_stations, 3): [x_le, y, chord]
    
    assert geometry.shape[1] == 3, "Geometry must have columns [x_le, y, chord]"
    
    y_geo = geometry[:, 1]
    y = np.linspace(y_geo[0], y_geo[-1], nspan+1)
    
    # Interpolações
    x_le = np.interp(y, y_geo, geometry[:,0])
    chord  = np.interp(y, y_geo, geometry[:,2])
    
    # Malha (nchord+1, nspan+1)
    X = np.zeros((nchord+1, nspan+1))
    Y = np.zeros((nchord+1, nspan+1))
    Z = np.zeros((nchord+1, nspan+1))
    
    for j in range(nspan+1):
        for i in range(nchord+1):
            eta = i / nchord   # 0 → 1
            X[i,j] = x_le[j] + chord[j] * eta
            Y[i,j] = y[j]
            Z[i,j] = 0.0
    
    # Painéis
    n_panels = nchord * nspan
    panels_id = np.arange(n_panels).reshape(nchord, nspan)
    panels = np.zeros((n_panels, 4, 3))
    k = 0
    
    for j in range(nspan): 
        for i in range(nchord):  # percorre primeiro a corda depois o span
            P1 = [X[i,   j],   Y[i,   j],   Z[i,   j]]
            P2 = [X[i+1, j],   Y[i+1, j],   Z[i+1, j]]
            P3 = [X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]]
            P4 = [X[i,   j+1], Y[i,   j+1], Z[i,   j+1]]
            
            panels[k,:,:] = np.array([P1, P2, P3, P4])
            k += 1

    return X, Y, Z, panels

def vortexAndControlPoint(panels:np.ndarray):
    
    npanels = panels.shape[0]
    
    vcp = np.zeros((npanels, 3, 3)) # (panel x [center, quarter chord, three quarter chord])
    
    for k in range(npanels):
        c_2 = np.mean(panels[k, :,:],axis=0)
        c_4 = c_2.copy()
        c_34 = c_2.copy()
        
        c_4[0] = c_2[0] - c_2[0]/2
        c_34[0] = c_2[0] + c_2[0]/2
        
        vcp[k, 0, :] = c_2 
        vcp[k, 1, :] = c_4
        vcp[k, 2, :] = c_34 
    
    return vcp
    