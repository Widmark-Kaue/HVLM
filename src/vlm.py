import numpy as np
import matplotlib.pyplot as plt

def plot_mesh(mesh:dict,title:str = ''):
    
    X = mesh['X']
    Y = mesh['Y']
    Z = mesh['Z']
    
    normals = mesh['normals']
    vcp = mesh['vcp']
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    scale = 0.05*np.max(X)
    for i in range(X.shape[0]):
        ax.plot(X[i,:], Y[i,:], Z[i,:], 'k-')
    for j in range(X.shape[1]):
        ax.plot(X[:,j], Y[:,j], Z[:,j], 'k-')
        
    for k in range(vcp.shape[0]):
        nx, ny, nz = normals[k]
        cx, cy, cz = vcp[k, 2,:]
        ax.plot(vcp[k, 0, 0], vcp[k, 0, 1], vcp[k, 0, 2], 'ko', markersize = 4) # center point
        ax.plot(vcp[k, 1, 0], vcp[k, 1, 1], vcp[k, 1, 2], 'ro', markersize = 4) # quarter point
        ax.plot(vcp[k, 2, 0], vcp[k, 2, 1], vcp[k, 2, 2], 'bo', markersize = 4) # three quarter point

        ax.quiver(cx, cy, cz, nx, ny, nz, color = 'gray', length = scale)
        
        
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z") # type: ignore
    ax.set_box_aspect([1, 1, 0.1]) # type: ignore
    fig.suptitle(title)
    plt.axis('equal')
    plt.show(block = False)

def geometry(span:float, AR:float, sections:list=[], taper:float = 1, sweep:float = 0, twist:float = 0, symmetry:bool = True):
    
    
    S = span**2/AR
    cmean = S/span
    cr = 2*cmean/(1+taper)
    ct = cr*taper
    
    if sections == []:
        nsection = 2
    else:
        nsection = len(sections) + 1
    if symmetry:
        y = np.linspace(0, span/2, nsection)
        chord = cr + 2*np.abs(y/span)*(ct - cr)
    else:
        y = np.linspace(-span/2, span/2, nsection)
        chord = cr + np.abs(y/span)*(ct - cr)
    
    sweep = np.deg2rad(sweep)
    twist = np.deg2rad(twist)
    x_le = np.tan(sweep)*y
    z_le = np.tan(twist)*y
    geometry = dict(
        leading_edge = np.column_stack([x_le, y, z_le]),
        chord = chord,
        camber = np.zeros_like(chord),
        span = span,
        AR = AR,
        S = S,
        cmean = cmean,
        sections = sections,
        sweep = np.rad2deg(sweep),
        twist = np.rad2deg(twist),
        symmetry = symmetry
    )
    
    return geometry

def meshPlanar(nspan:int, nchord:int, geometry:dict):
    # leading_edge: (N_stations, 3): [x_le, y, z_le]
    
    leading_edge = geometry['leading_edge']
    
    
    assert leading_edge.shape[1] == 3, "leading_edge must have columns [x_le, y, z_le]"
    
    y_geo = leading_edge[:, 1]
    # span = y_geo[-1] - y
    y = np.linspace(y_geo[0], y_geo[-1], nspan+1)
    
    # Interpolações
    x_le = np.interp(y, y_geo, leading_edge[:,0])
    z_le = np.interp(y, y_geo, leading_edge[:, 2])
    chord  = np.interp(y, y_geo, geometry['chord'])
    
    # Malha (nchord+1, nspan+1)
    X = np.zeros((nchord+1, nspan+1))
    Y = np.zeros((nchord+1, nspan+1))
    Z = np.zeros((nchord+1, nspan+1))
    
    for j in range(nspan+1):
        for i in range(nchord+1):
            eta = i / nchord   # 0 → 1
            X[i,j] = x_le[j] + chord[j] * eta
            Y[i,j] = y[j]
            Z[i,j] = z_le[j]
    
    # Painéis
    n_panels = nchord * nspan
    panels_id = np.arange(n_panels).reshape(nchord, nspan)
    panels = np.zeros((n_panels, 4, 3))
    panels_span = np.zeros((n_panels, 1))
    normals = np.zeros((n_panels, 3))
    area = np.zeros((n_panels,1))
    
    
    for j in range(nspan): 
        for i in range(nchord):  # percorre primeiro a corda depois o span
            
            # Pontos do painel
            P1 = np.array([X[i,   j],   Y[i,   j],   Z[i,   j]])
            P2 = np.array([X[i+1, j],   Y[i+1, j],   Z[i+1, j]])
            P3 = np.array([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            P4 = np.array([X[i,   j+1], Y[i,   j+1], Z[i,   j+1]])
            
            k = panels_id[i, j]
            panels[k,:,:] = np.array([P1, P2, P3, P4])
            
            # Normais
            Ak = P3 - P1
            Bk = P4 - P2          
            normal = np.cross(Ak, Bk)
            area[k] = np.linalg.norm(normal)
            normal_unit = normal/area[k]
            normals[k, :] = normal_unit
            panels_span[k] = P4[1] - P1[1] # talvez não seja o suficiente para asas enflechadas
            
    vcp = vortexAndControlPoint(panels)
    MESH = dict(
        X = X, Y = Y, Z = Z,
        panels = panels, 
        panels_span = panels_span,
        panels_id = panels_id,
        normals = normals,
        area = area,
        vcp = vcp,
        geometry = geometry
    )     

    return MESH

def vortexAndControlPoint(panels:np.ndarray):
    
    npanels = panels.shape[0]
    
    vcp = np.zeros((npanels, 3, 3)) # (panel x [center, quarter chord, three quarter chord])
    
    for k in range(npanels):
        center_point = np.mean(panels[k, :,:],axis=0)
        quarter_point = center_point.copy()
        three_quarter_point = center_point.copy()
        
        half_chord = center_point[0] - panels[k,0,0] # final - inicial
        
        quarter_point[0] = center_point[0] - half_chord/2
        three_quarter_point[0] = center_point[0] + half_chord/2
        vcp[k, 0, :] = center_point 
        vcp[k, 1, :] = quarter_point
        vcp[k, 2, :] = three_quarter_point 
    
    return vcp

def vortexl(p:np.ndarray, p1:np.ndarray, p2:np.ndarray, Gamma:float = 1, tol:float = 1e-6) -> np.ndarray:
    """_summary_

    Parameters
    ----------
    p : np.ndarray
        Point where the velocity will be compute
    p1 : np.ndarray
        Initial point of the vortex line
    p2 : np.ndarray
        Final point of the vortex line
    Gamma : float, optional
        Vorticity, by default 1
    tol : float, optional
        Tolerance to avoid singularities, by default 1e-3

    Returns
    -------
    np.ndarray
        Induced velocity
    """
    # Distances vectors
    r1 = p - p1
    r2 = p - p2
    r0 = p2 - p1
    
    # Unit vectos
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    
    assert r1_norm > tol, "vortexl: Singular condition, r1 is close enough of the core"
    assert r2_norm > tol, "vortexl: Singular condition, r2 is close enough of the core"
    
    r1_unit = r1/r1_norm
    r2_unit = r2/r2_norm
    
    # Cross product
    r1xr2 = np.cross(r1, r2)
    r1xr2_norm_square = np.linalg.norm(r1xr2)**2
    
    # distances = np.array([r1_norm, r2_norm, r1xr2_norm_square])
    # msg = ''
    # if np.any(distances < tol):
    #     msg = 'vortexl: Singular condition'
    #     print(msg)
    #     K = Gamma/(4*np.pi)/r1xr2_norm_square *(r0 @ r1_unit - r0 @ r2_unit)
    # else:
    #     K = Gamma/(4*np.pi)/r1xr2_norm_square *(r0 @ r1_unit - r0 @ r2_unit)
    
    # induced velocity
    K = Gamma/(4*np.pi)/r1xr2_norm_square *(r0 @ r1_unit - r0 @ r2_unit)
    q12 = K*r1xr2
    
    return q12
    
def hshoe(p:np.ndarray, pa:np.ndarray, pb:np.ndarray, pc:np.ndarray, pd:np.ndarray, Gamma:float = 1, tol:float = 1e-3)-> tuple:
    """_summary_

    Parameters
    ----------
    p : np.ndarray
        _description_
    pa : np.ndarray
        _description_
    pb : np.ndarray
        _description_
    pc : np.ndarray
        _description_
    pd : np.ndarray
        _description_
    Gamma : float, optional
        _description_, by default 1
    tol : float, optional
        _description_, by default 1e-3

    Returns
    -------
    tuple
        _description_
    """
    iv1 = vortexl(p, pa, pb, Gamma = Gamma, tol=tol)
    iv2 = vortexl(p, pb, pc, Gamma = Gamma, tol=tol)
    iv3 = vortexl(p, pc, pd, Gamma = Gamma, tol=tol)
    
    aij = iv1+iv2+iv3
    bij = iv1+iv3       # induced velocity by wake
    
    return aij, bij

def influence_coefficients(Vinf:np.ndarray, l_inf:float, mesh:dict, symmetry:bool = False)->tuple:
    """_summary_

    Parameters
    ----------
    Vinf : np.ndarray
        _description_
    l_inf : float
        _description_
    mesh: dict
        Dictionare with:
            panels : np.ndarray
                _description_
            panels_span:np.ndarray
                _description_
            vcp : np.ndarray
                _description_
            normals : np.ndarray
                _description_

    Returns
    -------
    tuple
        _description_
    """

    panels = mesh['panels']
    vcp = mesh['vcp']
    normals = mesh['normals']
    panels_span = mesh['panels_span']
    
    npanels = panels.shape[0]
    
    a = np.zeros((npanels, npanels))
    b = np.zeros((npanels, npanels))
    RHS = np.zeros(npanels)
    for ki in range(npanels):
        p = vcp[ki, 2, :]
        RHS[ki] = -Vinf @ normals[ki]
        for kj in range(npanels):
            # span = panels[kj, 3, 1] - panels[kj, 0, 1]
            span = panels_span[kj,0]
            
            # Define points of the horseshoe vortex
            pa = vcp[kj, 1, :] + np.array([0, -span/2, 0]) # point at infinite
            
            pa[0] = l_inf
            
            pb = vcp[kj, 1, :] + np.array([0, -span/2, 0])
            pc = vcp[kj, 1, :] + np.array([0, +span/2, 0])
            pd = vcp[kj, 1, :] + np.array([0, +span/2, 0]) # point at infinite
            
            pd[0] = l_inf
            
            aij, bij = hshoe(p, pa, pb, pc, pd)
            
            # Symmetry condition
            if symmetry:
                aux =np.array([1, -1, 1]) 
                ps = aux*p
                aij_sy, bij_sy = hshoe(ps, pa, pb, pc, pd)
                aij = aij + aux*aij_sy
                bij = bij + aux*bij_sy
                
            
            
            a[ki, kj] = aij @ normals[ki]
            b[ki, kj] = bij @ normals[ki]
    
    return a, b, RHS


def coefficients(Vinf:float, rho:float, Gamma:np.ndarray, wind:np.ndarray, mesh:dict):
    geometry = mesh['geometry']
    panels_span = mesh['panels_span'].reshape(Gamma.shape)
    
    aux = 2 if geometry['symmetry'] else 1
    
    # Ref parameters
    S = geometry['S']
    c = geometry['cmean']
    q = 0.5*rho*Vinf**2
    
    # Lift
    L = rho*Vinf*Gamma*panels_span
    L = np.sum(L)
    CL = aux*L/(q*S*c)
    
    # Induced drag
    D = -rho/2 * np.sum(Gamma*wind*panels_span)
    CDi = aux*D/(q*S*c)
    
    return CL, CDi

def run_polar(Vinf:float, rho:float, alpha:np.ndarray, mesh:dict):
    
    span = mesh['geometry']['span']
    sy = mesh['geometry']['symmetry']
    
    alpha_rad = np.deg2rad(alpha)
    CL = np.zeros(alpha.shape)
    CDi = np.zeros(alpha.shape)
    for i in range(len(alpha)):
        aoa = alpha_rad[i]
        V  = Vinf*np.array([np.cos(aoa), 0, np.sin(aoa)])
        A, B, RHS = influence_coefficients(V,20*span, mesh, symmetry=sy)
        Gamma = np.linalg.solve(A,RHS)
        wind = B @ Gamma
        
        cl, cdi = coefficients(Vinf, rho, Gamma, wind, mesh)
        CL[i] = cl
        CDi[i] = cdi
    return CL, CDi






    
    