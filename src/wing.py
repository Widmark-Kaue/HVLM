import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from scipy.interpolate import interp1d
from pathlib import Path


@dataclass
class Geometry:
    span:float
    AR:float
    sections:list[str]
    taper:float = 1
    sweep:float = 0
    twist:float = 0 
    symmetry:bool = True
    camber_points:int = 50
    S:float = field(init=False, default=1.0)
    leading_edge:np.ndarray = field(init = False, repr = False, default_factory=lambda: np.empty(1))
    chord:np.ndarray = field(init = False, repr = False, default_factory=lambda: np.empty(1))
    cmean:float = field(init=False, repr=False)
    camber:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    
    def __post_init__(self):
        
        span = self.span
        taper = self.taper
        sweep  = self.sweep
        twist = self.twist
        AR = self.AR
        sections = self.sections
        symmetry = self.symmetry
        
        S = span**2/AR
        cmean = S/span
        cr = 2*cmean/(1+taper)
        ct = cr*taper
        
        # number of sections
        nsection = len(sections)
                
        # define spanwise and chord distribution  
        if symmetry:
            y = np.linspace(0, span/2, nsection)
            chord = cr + 2*np.abs(y/span)*(ct - cr)
        else:
            y = np.linspace(-span/2, span/2, nsection)
            chord = cr + np.abs(y/span)*(ct - cr)
        
        # Evaluate camber
        camber = np.zeros((self.camber_points, nsection))
        for i, section in enumerate(sections):
            if 'naca' in section.lower():
                number = section.split(' ')[-1]
                camber[:, i] = self.NacaAirfoil(number)
        
                
        # Modify x and z leading edge coord accord with sweep and twist
        sweep = np.deg2rad(sweep)
        twist = np.deg2rad(twist)
        x_le = np.tan(sweep)*y
        z_le = np.tan(twist)*y
        
        # add properties
        self.leading_edge = np.column_stack([x_le, y, z_le])
        self.chord = chord
        self.camber = camber
        self.cmean = cmean
        self.S = S
    
    def NacaAirfoil(self, number:str, spacing:str = 'Uniform'):
        m = float(number[0])/100
        p = float(number[1])/10
        
        x = np.linspace(0, 1, self.camber_points)    
        y = np.zeros_like(x)
        pos1 = x < p
        pos2 = (x >= p) * (x < 1)
        x1 = x[pos1]
        x2 = x[pos2] 
        
        if p != 0:
            y[pos1] = m/(p**2)   * (2*p*x1 - x1**2)
        y[pos2] = m/(1-p)**2 * (1 - 2*p + 2*p*x2 - x2**2) 
        return y

@dataclass
class Wing(Geometry):
    name:str = field(default='')
    nspan:int = field(init=False)
    nchord:int = field(init = False)
    X:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    Y:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    Z:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    panels:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    panels_span:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    panels_id:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    area:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    normals:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    vcp:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    horseshoe:np.ndarray = field(init=False, repr=False, default_factory=lambda: np.empty(1))
    
    
    def __vortexAndControlPoint__(self, panels:np.ndarray, panels_span:np.ndarray):
        npanels = panels.shape[0]
        
        vcp = np.zeros((npanels, 3, 3)) # (panel x [center, quarter chord, three quarter chord])
        horseshoe = np.zeros((npanels, 2, 3)) # panels x [pb, pc] 
        for k in range(npanels):
            span= panels_span[k,1]
            P1 = panels[k,0]
            P2 = panels[k,1]
            P3 = panels[k,2]
            P4 = panels[k,3]
            Lem = 0.5*(P1+P4)
            Tem = 0.5*(P2+P3)
            
            # Panel middle points
            center_point = np.mean(panels[k, :,:],axis=0)
            quarter_point = Lem + 0.25*(Tem - Lem)
            three_quarter_point = Lem + 0.75*(Tem - Lem)
            # center_point = np.mean(panels[k, :,:],axis=0)
            # quarter_point = center_point.copy()
            # three_quarter_point = center_point.copy()
            
            # half_chord = center_point[0] - panels[k,0,0] # final - inicial
            
            # quarter_point[0] = center_point[0] - half_chord/2
            # three_quarter_point[0] = center_point[0] + half_chord/2
            vcp[k, 0, :] = center_point 
            vcp[k, 1, :] = quarter_point
            vcp[k, 2, :] = three_quarter_point
            
            # Define points of the horseshoe vortex
            # pb = vcp[k, 1, :] + np.array([0, -span/2, 0])
            # pc = vcp[k, 1, :] + np.array([0, +span/2, 0])
            pa = P1 + 0.25*(P2 - P1)
            pb = P4 + 0.25*(P3 - P4)
            
            horseshoe[k, 0, :] = pa
            horseshoe[k, 1, :] = pb
            
        return vcp, horseshoe
    
    def mesh(self, nspan:int, nchord:int):
        leading_edge = self.leading_edge
        chord = self.chord
        camber = self.camber
        
        assert leading_edge.shape[1] == 3, "leading_edge must have columns [x_le, y, z_le]"
        
        # Geometry points
        y_geo = leading_edge[:, 1]
        x_c_geo = np.linspace(0, 1, self.camber_points)

        # Mesh points
        y = np.linspace(y_geo[0], y_geo[-1], nspan+1)
        x_c = np.linspace(0, 1, nchord+1)
        
        # Interpolations
        x_le = np.interp(y, y_geo, leading_edge[:,0])
        z_le = np.interp(y, y_geo, leading_edge[:, 2])
        chord  = np.interp(y, y_geo, chord)
        camber_fy = interp1d(y_geo, camber, axis=1)(y)
        camber_fx = interp1d(x_c_geo, camber_fy, axis=0)(x_c)
        
        # Malha (nchord+1, nspan+1)
        X = np.zeros((nchord+1, nspan+1))
        Y = np.zeros((nchord+1, nspan+1))
        Z = np.zeros((nchord+1, nspan+1))
        
        for j in range(nspan+1):
            for i in range(nchord+1):
                eta = i / nchord   # 0 → 1
                X[i,j] = x_le[j] + chord[j] * eta
                Y[i,j] = y[j]
                Z[i,j] = z_le[j] + camber_fx[i, j]*chord[j]
        
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
        
        vcp, hshoe = self.__vortexAndControlPoint__(panels, panels_span=panels_span)
        #add mesh properties
        self.X = X
        self.Y = Y
        self.Z = Z
        self.nspan = nspan
        self.nchord = nchord
        self.panels = panels 
        self.panels_span = panels_span
        self.panels_id = panels_id
        self.normals = normals
        self.area = area
        self.vcp = vcp
        self.horseshoe = hshoe
    
    def plot_mesh(self, title:bool = False, savefig:bool = False):
        #%
        path_images = Path('images')
        case = self.name
        
        X = self.X
        Y = self.Y
        Z = self.Z
        
        normals = self.normals
        vcp = self.vcp
        
        
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111, projection='3d')

        scale = 0.05*np.max(X)
        for i in range(X.shape[0]):
            ax.plot(X[i,:], Y[i,:], Z[i,:], 'k-')
        for j in range(X.shape[1]):
            ax.plot(X[:,j], Y[:,j], Z[:,j], 'k-')
        
        if self.symmetry:
            for i in range(X.shape[0]):
                ax.plot(X[i,:], -Y[i,:], Z[i,:], 'k-')
            for j in range(X.shape[1]):
                ax.plot(X[:,j], -Y[:,j], Z[:,j], 'k-')
            
        for k in range(vcp.shape[0]):
            nx, ny, nz = normals[k]
            cx, cy, cz = vcp[k, 2,:]
            ax.plot(vcp[k, 0, 0], vcp[k, 0, 1], vcp[k, 0, 2], 'ko', markersize = 4) # center point
            ax.plot(vcp[k, 1, 0], vcp[k, 1, 1], vcp[k, 1, 2], 'ro', markersize = 4) # quarter point
            ax.plot(vcp[k, 2, 0], vcp[k, 2, 1], vcp[k, 2, 2], 'bo', markersize = 4) # three quarter point

            ax.quiver(cx, cy, cz, nx, ny, nz, color = 'gray', length = scale)
            
        
        titleName = self.name if title else ''   
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z") # type: ignore
        ax.set_box_aspect([1, 1, 0.1]) # type: ignore
        plt.axis('equal')
        plt.tight_layout()
        if savefig:
            plt.savefig(path_images.joinpath(f'{case}3d.pdf'), dpi = 600, format = 'pdf')
        fig.suptitle(titleName)
        plt.show(block = False)
        