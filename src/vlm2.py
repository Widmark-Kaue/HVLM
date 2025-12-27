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
        
        cosBeta = np.cos(self.beta)
        sinBeta = np.sin(self.beta)
        
        T1 = np.array([[ cosAlpha, 0, sinAlpha],
                       [        0, 1,       0 ],
                       [-sinAlpha, 0, cosAlpha]])
        
        return T1

    #################################################################
    # automatic functions
    #################################################################
    def run_polar(self, Vinf:float, alpha:np.ndarray, ref:str = 'Katz', coef:str = 'drela'):
    
        span = self.wing.span
        
        alpha_rad = np.deg2rad(alpha)
        CL = np.zeros(alpha.shape)
        CDi = np.zeros(alpha.shape)
        # CDi2 = np.zeros(alpha.shape)
        A, B = self.influence_coefficients_hshoe(20*span, ref = ref)
        for i in range(len(alpha)):
            aoa = alpha_rad[i]
            self.V  = Vinf*np.array([np.cos(aoa), 0, np.sin(aoa)])
            RHS = self.RHS_matrix()
            
            # Compute circulation
            self.Gamma = np.linalg.solve(A,RHS)

            # Compute induced velocity by wake
            self.wind = B @ self.Gamma
            
            
            # cl, cdi = coefficients(Vinf, rho, Gamma, wind, mesh)
            if coef.lower() == 'drela':
                cl, cdi = self.coefficients2()
            else:
                cl,cdi = self.coefficients()
            
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
        normals = self.wing.normals
        
        aux = 2 if self.wing.symmetry else 1
        
        #VLM parameters
        Gamma = self.Gamma
        wind = self.wind
        rho = self.rho
        V = self.V
        
        # # Lifting-surface convention
        # wind = wind.reshape(nchord, nspan) 
        # Gamma = Gamma.reshape(nchord, nspan).copy()
        
        # Ref parameters
        Vinf = self.Vinf
        Vunit = V.reshape(1, 3)/Vinf
        aoa = self.alpha
        S = self.wing.S
        c = self.wing.cmean
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
        self.CL = CL
        ### Induced drag
        D = -rho * np.sum(Gamma*wind*deltay)
        CDi = aux*D/(q*S)
        self.CDi = CDi
        return CL, CDi
    
    def coefficients2(self):
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
        
        # self.CL = np.sum(TotalForce @ ez)
        # # self.CDi = -np.sum(TotalForce @ ex)
        #Compute drag in Treffzt plane
        Gamma_re = Gamma.copy().reshape(nchord, nspan)
        # Gamma_re = np.sum(Gamma_re, axis = 0)
        wind_normal = wind_normal.reshape(nchord, nspan)
        # wind_normal = np.sum(wind_normal, axis = 0)
        D = -0.5*self.rho*np.sum(Gamma_re*wind_normal)
        self.CDi = aux*D/(q*S)
        return self.CL, self.CDi
        
        
    #################################################################
    # Vortex Lattice matrix construction functions 
    #################################################################
    def influence_coefficients_hshoe(self, l_inf:float, ref:str = 'Katz')->tuple:
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
        symmetry = self.wing.symmetry
        panels = self.wing.panels
        vcp = self.wing.vcp
        normals = self.wing.normals
        horseshoePoints = self.wing.horseshoe
            
        npanels = panels.shape[0]
        
        if ref.lower() == 'katz':
            funcH = lambda p, pa, pb, pc, pd: self.hshoe(p, pa, pb, pc, pd)
        else:
            funcH = lambda p, pa, pb, pc, pd: self.hshoeDrela(p, pb, pc)         
        
        
        A = np.zeros((npanels, npanels))
        B = np.zeros((npanels, npanels))
        for ki in range(npanels):
            p = vcp[ki, 2, :]       # point where the induced velocitys will be evaluate
            for kj in range(npanels):
                pb = horseshoePoints[kj, 0]
                pc = horseshoePoints[kj, 1]
                
                pa = pb.copy()
                pd = pc.copy()
                pa[0] = l_inf
                pd[0] = l_inf
                
                aij, bij = funcH(p, pa, pb, pc, pd)     
                # if method == 1:
                #     aij, bij = hshoe(p, pa, pb, pc, pd)
                # else:
                #     aij, bij = hshoeDrela(p, pb, pc)
                # Symmetry condition
                if symmetry:
                    aux =np.array([1, -1, 1]) 
                    ps = aux*p
                    aij_sy, bij_sy = funcH(ps, pa, pb, pc, pd)     
                    # if method ==1:
                    #     aij_sy, bij_sy = hshoe(ps, pa, pb, pc, pd)
                    # else:
                    #     aij_sy, bij_sy = hshoeDrela(ps, pb, pc)
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
            vi,vw = self.hshoeDrela(p, pa, pb)
            Vi+=vi*self.Gamma[k]
            Vw+=vw*self.Gamma[k]
        return Vi, Vw
       
   
    #################################################################
    # Induced velocity by a horseshoe element functions 
    #################################################################
    def vortexl(self,
        p:np.ndarray, 
        p1:np.ndarray, 
        p2:np.ndarray, 
        Gamma:float = 1, 
        ) -> np.ndarray:
        
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
        tol = self.tol
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
    
    def vortex2(self, 
                p:np.ndarray,
                p1:np.ndarray,
                p2:np.ndarray,
                Gamma = 1.0):
        """
        Induced velocity (u, v, w) at point (x,y,z)
        due to a finite vortex segment from (x1,y1,z1) to (x2,y2,z2)
        with circulation gama (per unit length).
        """

        PI = np.pi
        RCUT = self.tol
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
    
    def hshoeDrela(self, p:np.ndarray, pa:np.ndarray, pb:np.ndarray, **args):
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
    
    def hshoe(self, 
              p:np.ndarray, 
              pa:np.ndarray, 
              pb:np.ndarray, 
              pc:np.ndarray, 
              pd:np.ndarray, 
              Gamma:float = 1)-> tuple:
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
        
        iv1 = self.vortexl(p, pa, pb, Gamma = Gamma)
        iv2 = self.vortexl(p, pb, pc, Gamma = Gamma)
        iv3 = self.vortexl(p, pc, pd, Gamma = Gamma)
        
        # iv1 = vortex2(p, pa, pb, Gamma)
        # iv2 = vortex2(p, pb, pc, Gamma)
        # iv3 = vortex2(p, pc, pd, Gamma)
        aij = iv1+iv2+iv3
        bij = iv1+iv3       # induced velocity by wake
            
        
        return aij, bij