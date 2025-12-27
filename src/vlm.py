from src.wing import *

@dataclass
class VLM:
    wing: Wing
    rho:float = 1.225
    V:np.ndarray = field(init=True, repr=True, default_factory=lambda:np.zeros(3))
    
    
    # Vortex lattice parameters
    Gamma:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1))
    wind:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1))
    A:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1))
    B:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1))
    RHS:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1))
    Force:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1)) 
    polar:np.ndarray = field(init=False, repr=False, default_factory=lambda:np.empty(1)) 
    CL:float = field(init=False, default=0)
    CDi:float = field(init=False, default=0)
    
    
    # tolerance to induced velocity methods
    tol:float = 1e-10

    #################################################################
    # Properties
    #################################################################
    @property
    def ex(self): # Global unit vectors x
        return np.array([1, 0, 0])
    
    @property
    def ey(self): # Global unit vectors y
        return np.array([0, 1, 0])
    
    @property
    def ez(self): # Global unit vectors z
        return np.array([0, 0, 1])
    
    @property
    def Vinf(self):
        return np.linalg.norm(self.V)
    
    @property
    def alpha(self): # atan(Vz/Vx)
          return np.arctan(self.V[2]/self.V[0]) 
    
    @property
    def beta(self):# atan(-Vy/sqrt(Vx^2 + Vz^2))
          return np.arctan(-self.V[1]/np.sqrt(self.V[0]**2 + self.V[2]**2))
    @property
    def rotation_matrix(self):
        cosAlpha = np.cos(self.alpha)
        sinAlpha = np.sin(self.alpha)
        
        T1 = np.array([[ cosAlpha, 0, sinAlpha],
                       [        0, 1,       0 ],
                       [-sinAlpha, 0, cosAlpha]])
        
        return T1

    #################################################################
    # automatic functions
    #################################################################
    def run_polar(self, Vinf:float, alpha:np.ndarray):
    
        span = self.wing.span
        
        alpha_rad = np.deg2rad(alpha)
        CL = np.zeros(alpha.shape)
        CDi = np.zeros(alpha.shape)
        # CDi2 = np.zeros(alpha.shape)
        A, B = self.influence_coefficients_hshoe()
        for i in range(len(alpha)):
            aoa = alpha_rad[i]
            self.V  = Vinf*np.array([np.cos(aoa), 0, np.sin(aoa)])
            RHS = self.RHS_matrix()
            
            # Compute circulation
            self.Gamma = np.linalg.solve(A,RHS)

            # Compute induced velocity by wake
            self.wind = B @ self.Gamma
            
            cl, cdi = self.coefficients()
            
            CL[i] = cl
            CDi[i] = cdi
        
        polar = np.column_stack([alpha, CL, CDi])
        self.polar = polar
        return polar
    #################################################################
    # post-process functions 
    #################################################################
    
    def coefficients(self):
        # unit vectos
        ex = self.ex
        ez = self.ez
        
        # Geometry and mesh parameters
        nspan = self.wing.nspan
        nchord = self.wing.nchord
        panels_span = self.wing.panels_span
        vcp = self.wing.vcp
        npanels = panels_span.shape[0]
        aux = 2 if self.wing.symmetry else 1
        
        #VLM parameters
        Gamma = self.Gamma
        
        # Ref parameters
        Vinf = self.Vinf
        rho = self.rho
        S = self.wing.S
        q = 0.5*rho*Vinf**2
        
        Force = np.zeros_like(panels_span)
        wind_normal = np.zeros_like(Gamma)
        for k in range(npanels):
            p = vcp[k,1,:]          # middle point in bounded vortex line
            
            l = self.wing.horseshoe[k,1] - self.wing.horseshoe[k,0] # rb -ra 
            Vj, Vwj = self.induced_velocity(p, skip = k)
            Vi = Vj + self.V                # \sum{r=1-N} Gamma * Vj(r) + Vinf
            Vwi = Vwj
            ViXl = np.cross(Vi, l)
            VwiXl = np.cross(Vwi, l)
            Force[k] = self.rho*ViXl*Gamma[k]
            # wind_normal[k] = Vwi[2]* (l@self.ey)
            wind_normal[k] = -np.linalg.norm(VwiXl)
                    
        Force = aux*Force/(q*S)
        self.Force = Force
        TotalForce = np.sum(Force, axis = 0)
        
        Coef = self.rotation_matrix @ TotalForce.T
        self.CL = Coef[2]
        self.CDi = Coef[0]
        
        
        #Compute drag in Treffzt plane
        Gamma_re = Gamma.copy().reshape(nchord, nspan)
        wind_normal = wind_normal.reshape(nchord, nspan)
        D = -0.5*self.rho*np.sum(Gamma_re*wind_normal)
        self.CDi = aux*D/(q*S)
        return self.CL, self.CDi
        
        
    #################################################################
    # Vortex Lattice matrix construction functions 
    #################################################################
    def influence_coefficients_hshoe(self)->tuple:
        
        symmetry = self.wing.symmetry
        panels = self.wing.panels
        vcp = self.wing.vcp
        normals = self.wing.normals
        horseshoePoints = self.wing.horseshoe
            
        npanels = panels.shape[0]
        
        A = np.zeros((npanels, npanels))
        B = np.zeros((npanels, npanels))
        for ki in range(npanels):
            p = vcp[ki, 2, :]       # point where the induced velocitys will be evaluate
            for kj in range(npanels):
                pb = horseshoePoints[kj, 0]
                pc = horseshoePoints[kj, 1]
                
                aij, bij = self.hshoe(p,pb, pc)     
                
                if symmetry:
                    aux =np.array([1, -1, 1]) 
                    ps = aux*p
                    aij_sy, bij_sy = self.hshoe(ps, pb, pc)     
                    
                    aij = aij + aux*aij_sy
                    bij = bij + aux*bij_sy
                    
                
                
                A[ki, kj] = aij @ normals[ki]
                B[ki, kj] = bij @ normals[ki]
        self.A = A
        self.B = B
        return A, B

    def RHS_matrix(self):
        normals = self.wing.normals
        npanels = normals.shape[0]
        RHS = np.zeros(npanels)
        Vinf = self.V
        for k in range(npanels):
            RHS[k] = -Vinf @ normals[k]  
        
        self.RHS = RHS
        return RHS
   
    def induced_velocity(self, p:np.ndarray, skip:int = -1):
        npanels = self.wing.normals.shape[0]
        
        Vi = np.zeros(3)
        Vw = np.zeros(3)
        for k in range(npanels):
            pa = self.wing.horseshoe[k, 0, :]    
            pb = self.wing.horseshoe[k, 1, :]
            # if k == skip:
            #     _, vw = self.hshoeDrela(p, pa, pb)
            #     Vi+=vw*self.Gamma[k]
            #     Vw+=vw*self.Gamma[k]
            #     continue
            vi,vw = self.hshoe(p, pa, pb)
            Vi+=vi*self.Gamma[k]
            Vw+=vw*self.Gamma[k]
        return Vi, Vw
       
   
    #################################################################
    # Induced velocity by a horseshoe element functions 
    #################################################################
   
    
    def hshoe(self, p:np.ndarray, pa:np.ndarray, pb:np.ndarray, **args):
        ex = self.ex
        a = p - pa
        b = p - pb
        anorm = np.linalg.norm(a)
        bnorm = np.linalg.norm(b)
        
        aXb = np.cross(a, b)
        aXx = np.cross(a, ex)
        bXx = np.cross(b, ex)
        
        if np.all(np.abs(0.5*(pa+pb) - p) < self.tol):
            bounded_vortex = np.zeros(3)
        else:
            bounded_vortex = aXb/(anorm*bnorm + a@b)*(1/anorm + 1/bnorm)
            
        a_point_trailing_leg = aXx/(anorm - a@ex) * 1/anorm
        b_point_trailing_leg = bXx/(bnorm - b@ex) * 1/bnorm
        Vi = 1/(4*np.pi)*(bounded_vortex+a_point_trailing_leg - b_point_trailing_leg)
        Vwake = 1/(4*np.pi)*(a_point_trailing_leg - b_point_trailing_leg)
        return Vi, Vwake
    
   