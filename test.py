#%% Imports
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import use
use('TkAgg')

#%% Parâmetros da asa
b = 4           # envergadura
cr = 0.5        # corda na raiz
taper = 0.5     # razão de afilamento
ct = taper * cr # corda na ponta
sweep = 1      # Enflechamento (em graus) em 1/4 de corda
diedral = 0     # Diedro
twist = 0       # Torção (pode ser fixa ou ao longo da envergadura)


#%% Distribuições ao longo da envergadura
n = 10
y = np.linspace(0, b/2, n)
chord = cr - 2*(cr - ct)*y/b # Distribuição trapeizodal de corda

x_qc = cr/4 + y*np.tan(np.deg2rad(sweep))
x_le = x_qc - 0.25*chord
x_te = x_le + chord  

xr = np.concatenate((x_le, x_te[::-1]))
yr = np.concatenate((y, y[::-1]))

 
xl = xr[::-1]
yl = -yr


X = np.concatenate((xl, xr))
Y = np.concatenate((yl, yr))
Z = np.zeros(len(X))

plt.figure()
plt.plot(X,Y, 'k')
plt.fill(X, Y, alpha = 0.15)
plt.axis('equal')
plt.show(block = False)

# ax = plt.figure().add_subplot(projection = '3d')
# ax.plot(X, Y, Z, 'k')
# plt.axis('equal')
# plt.show(block = False)
