import numpy as np
import matplotlib.pyplot as plt

Nx = 200
Ny = 200

x = np.linspace(0,1,Nx)
y = np.linspace(0,1,Ny)

its = 500

V = np.zeros((Nx,Ny))
V_old = V

for it in range(its):
    V_old = V
    for i in range(Nx):
        for j in range(Ny):
            if i == 0 or j == 0 or i == Nx-1 or j == Ny-1:
                V[i,j] = 0
            elif np.sqrt((x[i]-.2)**2 + (y[j]-.2)**2) <= .1:
                V[i,j] = 1
            else:
                V[j,i] = .25*(V_old[i+1,j]+V_old[i-1,j]+V_old[i,j+1]+V_old[i,j-1])

plt.imshow(V)

plt.show()

