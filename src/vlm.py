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
    panels_span = np.zeros((n_panels, 3))
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
            panels_span[k] = P4 - P1 # talvez não seja o suficiente para asas enflechadas
            
    vcp = vortexAndControlPoint(panels)
    MESH = dict(
        X = X, Y = Y, Z = Z,
        nspan = nspan,
        nchord = nchord,
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

def vortexl(p:np.ndarray, p1:np.ndarray, p2:np.ndarray, Gamma:float = 1, tol:float = 1e-10) -> np.ndarray:
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
    
    distances = np.array([r1_norm, r2_norm, r1xr2_norm_square])
    msg = ''
    if np.any(distances < tol):
        msg = 'vortexl: Singular condition'
        # print(msg)
        K = 0
    else:
        K = Gamma/(4*np.pi)/r1xr2_norm_square *(r0 @ r1_unit - r0 @ r2_unit)
    
    # induced velocity
    # K = Gamma/(4*np.pi)/r1xr2_norm_square *(r0 @ r1_unit - r0 @ r2_unit)
    q12 = K*r1xr2
    
    return q12

def vortex2(p,p1,p2,Gamma = 1.0):
    """
    Induced velocity (u, v, w) at point (x,y,z)
    due to a finite vortex segment from (x1,y1,z1) to (x2,y2,z2)
    with circulation gama (per unit length).
    """

    PI = np.pi
    RCUT = 1.0e-10
    x, y, z = p
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    

    # R1 x R2
    r1r2x = (y - y1)*(z - z2) - (z - z1)*(y - y2)
    r1r2y = -((x - x1)*(z - z2) - (z - z1)*(x - x2))
    r1r2z = (x - x1)*(y - y2) - (y - y1)*(x - x2)

    # |R1 x R2|^2
    square = r1r2x**2 + r1r2y**2 + r1r2z**2

    # |R1| and |R2|
    r1 = np.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
    r2 = np.sqrt((x - x2)**2 + (y - y2)**2 + (z - z2)**2)

    # Cutoff (same logic as GOTO 1)
    if (r1 < RCUT) or (r2 < RCUT) or (square < RCUT):
        return np.zeros(3)

    # R0 · R1 and R0 · R2
    r0r1 = (x2 - x1)*(x - x1) + (y2 - y1)*(y - y1) + (z2 - z1)*(z - z1)
    r0r2 = (x2 - x1)*(x - x2) + (y2 - y1)*(y - y2) + (z2 - z1)*(z - z2)

    coef = Gamma / (4.0 * PI * square) * (r0r1 / r1 - r0r2 / r2)

    u = r1r2x * coef
    v = r1r2y * coef
    w = r1r2z * coef
    q = np.array([u, v, w])

    return q
       
def hshoe(p:np.ndarray, pa:np.ndarray, pb:np.ndarray, pc:np.ndarray, pd:np.ndarray, Gamma:float = 1, tol:float = 1e-10)-> tuple:
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
    
    # iv1 = vortex2(p, pa, pb, Gamma)
    # iv2 = vortex2(p, pb, pc, Gamma)
    # iv3 = vortex2(p, pc, pd, Gamma)
    
    aij = iv1+iv2+iv3
    bij = iv1+iv3       # induced velocity by wake
    
    return aij, bij

def influence_coefficients(mesh:dict,l_inf:float)->tuple:
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
    symmetry = mesh['geometry']['symmetry']
    panels = mesh['panels']
    vcp = mesh['vcp']
    normals = mesh['normals']
    panels_span = mesh['panels_span']
    
    npanels = panels.shape[0]
    
    A = np.zeros((npanels, npanels))
    B = np.zeros((npanels, npanels))
    for ki in range(npanels):
        p = vcp[ki, 2, :]       # point where the induced velocitys will be evaluate
        for kj in range(npanels):
            # span = panels[kj, 3, 1] - panels[kj, 0, 1]
            span = panels_span[kj,1]
            
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
                
            
            
            A[ki, kj] = aij @ normals[ki]
            B[ki, kj] = bij @ normals[ki]
    
    return A, B

def RHS_matrix(Vinf:np.ndarray, normals:np.ndarray):
    npanels = normals.shape[0]
    RHS = np.zeros(npanels)
    for k in range(npanels):
        RHS[k] = -Vinf @ normals[k]  
    
    return RHS

def coefficients(V:np.ndarray, rho:float, Gamma:np.ndarray, wind:np.ndarray, mesh:dict):
    
    # Unit vectors
    ex = np.array([1, 0, 0])
    ey = np.array([0, 1, 0])
    ez = np.array([0, 0, 1])
    
    # Geometry and mesh parameters
    nspan = mesh['nspan']
    nchord = mesh['nchord']
    geometry = mesh['geometry']
    panels_span = mesh['panels_span']
    normals = mesh['normals']
    
    aux = 2 if geometry['symmetry'] else 1
    
    
    # # Lifting-surface convention
    # wind = wind.reshape(nchord, nspan) 
    # Gamma = Gamma.reshape(nchord, nspan).copy()
    
    # Ref parameters
    Vinf = np.linalg.norm(V)
    Vunit = V.reshape(1, 3)/Vinf
    aoa = np.arccos(V@ex /Vinf)
    S = geometry['S']
    c = geometry['cmean']
    q = 0.5*rho*Vinf**2
    
    # Correction panels_span
    deltaS = panels_span  - (panels_span @ Vunit.copy().reshape(3, 1))*Vunit.repeat(panels_span.shape[0], axis=0) 
    deltay = np.linalg.norm(deltaS,axis=1)
    # panels_span_proj = panels_span @ ey.reshape(3, 1)
    # panels_span_proj = panels_span_proj.reshape(Gamma.shape)
    
    # Correction of donwash velocity
    angle = normals @ ez.reshape(3, 1)
    angle = np.arccos(angle).reshape(wind.shape)
    wind = wind*np.cos(aoa + angle)
    wind = wind.reshape(Gamma.shape)
    
    #### Lift
    L = rho*Vinf*np.sum(Gamma*deltay)
    CL = aux*L/(q*S)
    
    ### Induced drag
    D = -rho * np.sum(Gamma*wind*deltay)
    CDi = aux*D/(q*S)
    
    return CL, CDi

def run_polar(Vinf:float, alpha:np.ndarray, mesh:dict, rho:float = 1.225):
    
    span = mesh['geometry']['span']
    normals = mesh['normals']
    npan = normals.shape[0]
    
    alpha_rad = np.deg2rad(alpha)
    CL = np.zeros(alpha.shape)
    CDi = np.zeros(alpha.shape)
    # CDi2 = np.zeros(alpha.shape)
    A, B = influence_coefficients(mesh,20*span)
    for i in range(len(alpha)):
        aoa = alpha_rad[i]
        V  = Vinf*np.array([np.cos(aoa), 0, np.sin(aoa)])
        RHS = RHS_matrix(V, normals)
        
        # Compute circulation
        Gamma = np.linalg.solve(A,RHS)
        
        # Compute induced velocity by wake
        wind = B @ Gamma
        
        # cl, cdi = coefficients(Vinf, rho, Gamma, wind, mesh)
        cl,cdi = coefficients(V, rho, Gamma, wind, mesh)
        CL[i] = cl
        CDi[i] = cdi
    
    polar = np.column_stack([alpha, CL, CDi])
    return polar





    
    