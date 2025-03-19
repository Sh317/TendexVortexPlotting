from ctypes import *
from datetime import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np

from util.plot import colorline
from util.seeding import seed


#Runtime
start = dt.now()

class vect(Structure):
    _fields_ = [('x', c_double*10000), ('y', c_double*10000), ('z', c_double*10000), ('m', c_double*10000), ('its', c_int)]

#Value setup
title, x0, y0, z0 = seed(4) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random
curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'c':[],'lg':[]}
seeds = len(x0)
num_its = 50
delta_0 = 10e-2
h0 = 10e-2
safety = .9
ending_tolerance = .5
pos_color = 'red'
neg_color = 'blue'
icity = 1
R = 8
sigma = 1
vX = 1

# load C++
rka_iter_E = CDLL("./cppScripts/EFieldCalc").rka_iter
rka_iter_B = CDLL("./cppScripts/BFieldCalc").rka_iter
rka_iter_E.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double]
rka_iter_B.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double]
rka_iter_E.restype = vect
rka_iter_B.restype = vect

# Start plot
pos_fig, pos_ax = plt.subplots(2, sharex=True)
neg_fig, neg_ax = plt.subplots(2, sharex=True)

#fig, ax = plt.subplots((2,2), sharex=True)


#Color maps
pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('lightsteelblue'), ('navy')], N=100)
neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('mistyrose'), ('darkred')], N=100)
blank_pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('blue'), ('blue')], N=1)
blank_neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('red'), ('red')], N=1)
E_norm = mcolors.LogNorm(vmin=10e-16, vmax=10e-7)
B_norm = mcolors.LogNorm(vmin=10e-7, vmax=10e-1)
width = False
color = True


for i in range(seeds):
    for icity in [1,-1]:
        # Starting point of each field line
        x, y, z = x0[i], y0[i], z0[i]
    
        vect_c = rka_iter_E(R, sigma, vX, x, y, z, num_its, icity, ending_tolerance, delta_0, safety, h0)
        its = vect_c.its
        x = vect_c.x[0:its]
        y = vect_c.y[0:its]
        z = vect_c.z[0:its]
        m = np.abs(np.array(vect_c.m[0:its]))

        if icity == 1:
            if width:
                lc = colorline(pos_ax[0], x, y, m, norm = E_norm, widths = m, cmap=blank_pos_cmap)
            elif color:
                lc = colorline(pos_ax[0], x, y, m, norm = E_norm, cmap='jet')
        else:
            if width:
                lc = colorline(neg_ax[0], x, y, m, norm = E_norm, widths = m, cmap=blank_neg_cmap)
            elif color:
                lc = colorline(neg_ax[0], x, y, m, norm = E_norm, cmap='jet')

for i in range(seeds):
    for icity in [1,-1]:
        # Starting point of each field line
        x, y, z = x0[i], y0[i], z0[i]
    
        vect_c = rka_iter_B(R, sigma, vX, x, y, z, num_its, icity, ending_tolerance, delta_0, safety, h0)
        its = vect_c.its
        x = vect_c.x[0:its]
        y = vect_c.y[0:its]
        z = vect_c.z[0:its]
        m = np.abs(np.array(vect_c.m[0:its]))
        #print(min(m))
        #print(max(m))
        #print()

        if icity == 1:
            if width:
                lc = colorline(pos_ax[1], x, y, m, norm = B_norm, widths = m, cmap=blank_pos_cmap)
            elif color:
                lc = colorline(pos_ax[1], x, y, m, norm = B_norm, cmap='jet')
        else:
            if width:
                lc = colorline(neg_ax[1], x, y, m, norm = B_norm, widths = m, cmap=blank_neg_cmap)
            elif color:
                lc = colorline(neg_ax[1], x, y, m, norm = B_norm, cmap='jet')

print(dt.now() - start) 

if color:
    pos_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = E_norm),ax=pos_ax[0])
    neg_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = E_norm),ax=neg_ax[0])

    pos_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = B_norm),ax=pos_ax[1])
    neg_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = B_norm),ax=neg_ax[1])

lim = 12

pos_ax[0].set_xlim(-lim,lim)
pos_ax[0].set_ylim(-lim,lim)
pos_ax[1].set_xlim(-lim,lim)
pos_ax[1].set_ylim(-lim,lim)
neg_ax[0].set_xlim(-lim,lim)
neg_ax[0].set_ylim(-lim,lim)
neg_ax[1].set_xlim(-lim,lim)
neg_ax[1].set_ylim(-lim,lim)

pos_ax[0].set_title('Positive Electric component')
pos_ax[1].set_title('Positive Magnetic component')
neg_ax[0].set_title('Negative Electric component')
neg_ax[1].set_title('Negative Magnetic component')




pos_fig.savefig('pos.png', dpi=300)
neg_fig.savefig('neg.png', dpi=300)