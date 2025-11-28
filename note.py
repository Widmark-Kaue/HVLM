import numpy as np
from scipy.linalg import solve

def panel(P1:np.ndarray, P2:np.ndarray, P3:np.ndarray, P4:np.ndarray)-> tuple:
    
    # Vetores auxiliares
    Pa = P2 - P3
    Pb = P4 - P1
    
    A = np.cross(Pa, Pb)
    norm = np.linalg.norm(A)
    if norm == 0.0:
        C = np.zeros_like(A)
    else:
        C = A/norm
        
    # cálculo de áreas auxiliares para estimar área do painel
    E = P3 - P1
    F = P2 - P1
    
    s1 = np.cross(F, Pb)
    s2 = np.cross(Pb, e)

   

    area = 0.5*(np.linalg.norm(s1) + np.linalg.norm(s2))
    return C, area


def grid(Nchord, Nspan, B, X, SN1, CS1, CH, QF=None, QC=None, DS=None, DXW=None):
    """
    Tradução do SUBROUTINE GRID do Fortran.
    Parâmetros de entrada:
      IB, JB : inteiros (n_chordwise_boxes, n_spanwise_boxes)
      B      : semi-envergadura (ou span, conforme o código)
      X      : array-like com 4 elementos [X_root_LE, X_tip_LE, X_tip_TE, X_root_TE]
      SN1, CS1: sin(alpha), cos(alpha) (usados para inclinação no eixo z)
      CH     : altura sobre o solo (usada no cálculo de z)
      QF, QC, DS : arrays opcionais pre-alocadas (veja shapes abaixo)
      DXW    : extensão do wake (se None, será definido como 100.*B)

    Retorna:
      QF, QC, DS, S, C, AR, DXW
    Shapes esperadas (se não passadas, serão criadas automaticamente):
      QF: (IB+2, JB+1, 3)
      QC: (IB,   JB,   3)
      DS: (IB,   JB,   4)  # [nx, ny, nz, area]
    """
    X = np.asarray(X, dtype=float)
    if DXW is None:
        DXW = 100.0 * B

    IB1 = IB + 1
    IB2 = IB + 2
    JB1 = JB + 1

    # aloca arrays se necessário
    if QF is None:
        QF = np.zeros((IB+2, JB+1, 3), dtype=float)  # indices 0..IB+1, 0..JB
    if QC is None:
        QC = np.zeros((IB, JB, 3), dtype=float)      # collocation: 0..IB-1, 0..JB-1
    if DS is None:
        DS = np.zeros((IB, JB, 4), dtype=float)      # nx,ny,nz,area

    # coordenadas dos vértices e pontos de wake
    DY = B / float(JB)
    for j in range(JB1):               # j = 0..JB  <-> Fortran J=1..JB1
        YLE = DY * j
        XLE = X[0] + (X[1] - X[0]) * YLE / B
        XTE = X[3] + (X[2] - X[3]) * YLE / B
        DX = (XTE - XLE) / float(IB)

        # Fortran: DO I=1,IB1 ; QF(I,J,1)=(XLE+DX*(I-0.75))*CS1
        # Python index i = I-1 -> fator (I - 0.75) = (i + 1 - 0.75) = i + 0.25
        for i in range(IB1):          # i = 0..IB  <-> Fortran I=1..IB1
            x_coord = XLE + DX * (i + 0.25)
            QF[i, j, 0] = x_coord * CS1
            QF[i, j, 1] = YLE
            QF[i, j, 2] = - QF[i, j, 0] * SN1 + CH

        # wake far-field point (Fortran QF(IB2,J,*) = ...)
        # Fortran IB2 index -> python index IB+1
        QF[IB+1, j, 0] = XTE + DXW
        # QF(IB2,J,2) = QF(IB1,J,2) where IB1 -> python index IB
        QF[IB+1, j, 1] = QF[IB, j, 1]
        QF[IB+1, j, 2] = QF[IB, j, 2]

    # collocation points QC and cálculo de vetores normais / área via panel()
    for j in range(JB):                 # j = 0..JB-1  <-> Fortran J=1..JB
        for i in range(IB):             # i = 0..IB-1  <-> Fortran I=1..IB
            # médias dos 4 cantos do elemento
            qc_x = (QF[i, j, 0] + QF[i, j+1, 0] + QF[i+1, j+1, 0] + QF[i+1, j, 0]) / 4.0
            qc_y = (QF[i, j, 1] + QF[i, j+1, 1] + QF[i+1, j+1, 1] + QF[i+1, j, 1]) / 4.0
            qc_z = (QF[i, j, 2] + QF[i, j+1, 2] + QF[i+1, j+1, 2] + QF[i+1, j, 2]) / 4.0
            QC[i, j, 0] = qc_x
            QC[i, j, 1] = qc_y
            QC[i, j, 2] = qc_z

            # chama panel() para calcular normal e área
            c1, c2, c3, area = panel(
                QF[i, j, 0], QF[i, j, 1], QF[i, j, 2],
                QF[i+1, j, 0], QF[i+1, j, 1], QF[i+1, j, 2],
                QF[i, j+1, 0], QF[i, j+1, 1], QF[i, j+1, 2],
                QF[i+1, j+1, 0], QF[i+1, j+1, 1], QF[i+1, j+1, 2]
            )
            DS[i, j, 0] = c1
            DS[i, j, 1] = c2
            DS[i, j, 2] = c3
            DS[i, j, 3] = area

    # S - área total da semi-asas seg. Fortran: S=0.5*(X(3)-X(2)+X(4)-X(1))*B
    S = 0.5 * ((X[2] - X[1]) + (X[3] - X[0])) * B
    C = S / B
    AR = 2.0 * B * B / S

    return QF, QC, DS, S, C, AR, DXW
Uso rápido (exemplo)
python
Copiar código
import numpy as np
IB = 4
JB = 13
B = 13.0
X = [0.0, 0.0, 4.0, 4.0]   # conforme seu código Fortran
alpha = 5.0 * np.pi/180.0
SN1 = np.sin(alpha)
CS1 = np.cos(alpha)
CH = 1000.0

QF, QC, DS, S, C, AR, DXW = grid(IB, JB, B, X, SN1, CS1, CH)
print("S, C, AR:", S, C, AR)
    

