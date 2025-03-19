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

#x0 = [-40]
#y0 = [40]
#z0 = [0]

seeds = len(x0)
num_its = 2000
delta_0 = 10e-3
h0 = 10e-3
safety = .9
ending_tolerance = 1.0
pos_color = 'red'
neg_color = 'blue'
icity = 0
sepX = 0
sepY = .1

# load C++
rka_iter = CDLL("./cppScripts/rka_iter_90").rka_iter_double
rka_iter.argtypes = [c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double]
rka_iter.restype = vect

# Start plot
pos_fig, pos_ax = plt.subplots(1)
neg_fig, neg_ax = plt.subplots(1)

#Color maps
pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('lightsteelblue'), ('navy')], N=100)
neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('mistyrose'), ('darkred')], N=100)
blank_pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('blue'), ('blue')], N=1)
blank_neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('red'), ('red')], N=1)
pos_norm = mcolors.Normalize(vmin=0, vmax=80)
neg_norm = mcolors.Normalize(vmin=0, vmax=80)
width = False
color = True


for i in range(seeds):
    for icity in [1,-1]:
        # Starting point of each field line
        x, y, z = x0[i], y0[i], z0[i]
    
        vect_c = rka_iter(sepX, sepY, x, y, z, num_its, icity, ending_tolerance, delta_0, safety, h0)
        its = vect_c.its
        x = vect_c.x[0:its]
        y = vect_c.y[0:its]
        z = vect_c.z[0:its]
        m = np.abs(np.array(vect_c.m[0:its]))

        if icity == 1:
            if width:
                lc = colorline(pos_ax, x, y, m, widths = m, cmap=blank_pos_cmap, norm = pos_norm)
            elif color:
                lc = colorline(pos_ax, x, y, m, cmap='jet', norm = pos_norm)
            else:
                pos_ax.plot(x,y,color='blue')
        else:
            if width:
                lc = colorline(neg_ax, x, y, m, widths = m, cmap=blank_neg_cmap, norm = neg_norm)
            elif color:
                lc = colorline(neg_ax, x, y, m, cmap='jet', norm = pos_norm)
            else:
                neg_ax.plot(x,y,color='red')

print(dt.now() - start)

if color:
    pos_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = pos_norm), ax=pos_ax, label = 'pos eigenvalue')
    neg_fig.colorbar(cm.ScalarMappable(cmap='jet', norm = neg_norm), ax=neg_ax, label = 'neg eigenvalue * -1')
lim = 200
pos_ax.scatter([sepX,-sepX], [sepY,-sepY])
pos_ax.set_xlim(-lim,lim)
pos_ax.set_ylim(-lim,lim)
neg_ax.scatter([sepX,-sepX], [sepY,-sepY])
neg_ax.set_xlim(-lim,lim)
neg_ax.set_ylim(-lim,lim)
neg_ax.set_title('2 Dipole, x sep: %s y sep: %s, 180 degree difference' %(sepX, sepY))


pos_fig.savefig('pos.png', dpi=500)
neg_fig.savefig('neg.png', dpi=500)