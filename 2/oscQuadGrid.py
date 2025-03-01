import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from datetime import datetime as dt
from seeding import seed
import matplotlib.pyplot as plt
from plot import colorline
import matplotlib.colors as mcolors
import matplotlib.cm as cm

start = dt.now()

def f(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
    I = np.array([[a_I, 0, 0], 
              [0, -a_I, 0],
              [0, 0, 0]])
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1**2 - x_2**2)))

    E_1 = I
    E_2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            E_2[i,j] = I[0,i] * x_1 * r_V[j] + I[1,i] * x_2 * r_V[j] + I[0,j] * x_1 * r_V[i] + I[1,j] * x_2 * r_V[i]

    #E_2 = np.array([[x_1**2 * a_I, 1/2 * a_I * x_1 * x_2, 1/2 * a_I * x_1 * x_3], 
    #                    [1/2 * a_I * x_1 * x_2, -1 * x_2**2 * a_I, -1/2 * a_I * x_2 * x_3], 
    #                    [1/2 * a_I * x_1 * x_3, -1/2 * a_I * x_2 * x_3, 0]])
    

    E_3 = np.identity(3) * IpqXpXq
    E_4 = np.zeros((3,3))
    for i in range(len(I)):
        for j in range(len(I)):
            E_4[j,i] = (r_V[i] * r_V[j])

    E_4 *= IpqXpXq


    E = (-6 * E_1) + (30 * E_2 / (r**2)) + (15 * E_3 / (r**2)) + (-105 * E_4 / (r**4))
    
    if give:
        return E, E_1, E_2, E_3, E_4
    return E

    I = np.array([[0, a_I, 0], 
              [a_I, 0, 0],
              [0, 0, 0]])
    
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1 * x_2 + x_1 * x_2)))

    E_1 = I
    E_2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            E_2[i,j] = I[0,i] * x_1 * r_V[j] + I[1,i] * x_2 * r_V[j] + I[0,j] * x_1 * r_V[i] + I[1,j] * x_2 * r_V[i]

    #E_2 = np.array([[x_1**2 * a_I, 1/2 * a_I * x_1 * x_2, 1/2 * a_I * x_1 * x_3], 
    #                    [1/2 * a_I * x_1 * x_2, -1 * x_2**2 * a_I, -1/2 * a_I * x_2 * x_3], 
    #                    [1/2 * a_I * x_1 * x_3, -1/2 * a_I * x_2 * x_3, 0]])
    

    E_3 = np.identity(3) * IpqXpXq
    E_4 = np.zeros((3,3))
    for i in range(len(I)):
        for j in range(len(I)):
            E_4[j,i] = (r_V[i] * r_V[j])

    E_4 *= IpqXpXq


    E = (-6 * E_1) + (30 * E_2 / (r**2)) + (15 * E_3 / (r**2)) + (-105 * E_4 / (r**4))
    
    if give:
        return E, E_1, E_2, E_3, E_4
    return E

def e_solve(E, icity, val_return = False, s = None, r_vec = None): # r_past dotted with e vects
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    # Create array with False for values with some sign as icity, true otherwise
    e_valcheck = (e_vals * icity < 0)

    # Mask all evals with Trues, find the index
    index = np.argmax(ma(e_vals * icity, e_valcheck))    

    # Access correct e vec
    r_change = e_vecs[:, index]

    if val_return: # Needed for saving eigenvalue
        mag = e_vals[index]
        #with open("output.txt", "a") as file:
        #    file.write(str(e_vals) + ' ' + str(index) + ' ' + str(r_vec) + ' ' + str(e_valcheck) + ' ' + str(e_vals[index]) + ' ' + str(s) + '\n')
        return r_change, mag
    else:
        return r_change

p1 = np.array([10,10,0])
a_I = 1

x_cords = 200
y_cords = 200

pos_norm = mcolors.Normalize(vmin=0, vmax=100)
neg_norm = mcolors.Normalize(vmin=0, vmax=100)

x_list = np.linspace(-80,80,x_cords)
y_list = np.linspace(-80,80,y_cords)

pos_grid = np.meshgrid(x_list,y_list)[0]
neg_grid = np.meshgrid(x_list,y_list)[0]

for i, x in enumerate(x_list):
    for j, y in enumerate(y_list):
        r_V = np.array([x,y,0])

        E = f(r_V)

        _, pos_e_val = e_solve(E, 1, True)
        _, neg_e_val = e_solve(E, -1, True)

        pos_grid[j,i] = pos_e_val
        neg_grid[j,i] = neg_e_val * -1

pos_fig, pos_ax = plt.subplots(1)
neg_fig, neg_ax = plt.subplots(1)

pos_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = pos_norm), ax = pos_ax, label = 'eigenvalue')
pos_ax.imshow(pos_grid, 'jet')
pos_fig.savefig('plots/posGrid.png')
neg_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = neg_norm), ax = neg_ax, label = 'eigenvalue')
neg_ax.imshow(neg_grid, 'jet')
neg_fig.savefig('plots/negGrid.png')


print('Runtime: %s' %(dt.now() - start))
