from ctypes import *
from datetime import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import matplotlib.animation as animation
from functools import partial

from util.plot import colorline
from util.seeding import seed

#Runtime
start = dt.now()

class vects(Structure):
    _fields_ = [('x', (c_double*1000)*100), 
                ('y', (c_double*1000)*100), 
                ('z', (c_double*1000)*100), 
                ('m', (c_double*1000)*100), 
                ('its', c_int*100)]

class vect(Structure):
    _fields_ = [('x', c_double*10000), ('y', c_double*10000), ('z', c_double*10000), ('m', c_double*10000), ('its', c_int)]


# load C++
rka_iter = CDLL("./cppScripts/rka_iter_custT").rka_iter_double_seeded
rka_iter.argtypes = [c_double, c_double, c_double, c_int, c_int, c_int, c_int, c_double, c_double, c_double, c_double]
rka_iter.restype = vects

# Value setup
num_its = 100 # Iterations per seed point
delta_0 = 10e-2 # Min deviation from next step
h0 = 10e-3 #Starting step size
safety = .98 # Scaling factor for step size change
ending_tolerance = 1.0 # How far from poles to step simulation
icity = 1

# Color maps
cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('lightsteelblue'), ('navy')], N=100)
norm = mcolors.Normalize(vmin=0, vmax=80)

# Animation setup
seeds = 800
seed_max = 150
lim = 250 # x y bounds for plot
frame_num = 50 # Num frames for animation
artists = []
t = .78
print(t)
sep = 10.0
v = rka_iter(t, sep, sep, seeds, seed_max, num_its, icity, ending_tolerance, delta_0, safety, h0)
print(v.x[0][0])
'''

its = vect_c.its
x_list = vect_c.x
y_list = vect_c.y
z_list = vect_c.z
m_list = vect.c.m
m = np.abs(np.array(vect_c.m[0:its]))

lc = colorline(ax, x, y, m, cmap='jet', norm = norm)

ax.scatter([-sep,0], [-sep,0])
ax.set_xlim(-lim,lim)
ax.set_ylim(-lim,lim)

title = 'pi* ' + str(int(((t) / np.pi)*100)/100) +' S:' + str(int(sep))
fig.suptitle(title)

plt.show()
'''