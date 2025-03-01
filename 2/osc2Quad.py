import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from datetime import datetime as dt
from seeding import seed
import matplotlib.pyplot as plt
from plot import colorline
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#Runtime
start = dt.now()

#with open("output.txt", "w") as file:
#    file.write('start: ' + str(start) + '\n')



# RK Coeffecients
a = {'2':1/5,         '3':3/10,     '4':3/5,         '5':1,            '6': 7/8}
b = {'21':1/5, 
     '31':3/40,       '32':9/40, 
     '41':3/10,       '42':-9/10,   '43':6/5, 
     '51':-11/54,     '52':5/2,     '53':-70/27,    '54': 35/27, 
     '61':1631/55296, '62':175/512, '63':575/13824, '64':44275/110592, '65': 253/4096}
c_dif = {'1':37/378 - 2825/27648, '3':250/621 - 18575/48384, '4':125/621 - 13525/55296, '5':-277/14336, '6':512/1771 - 1/4}
c = {'1':37/378, '3':250/621, '4':125/621, '6':512/1771}

def f1(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
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

def f2(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
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
    
def f(r_V): #assuming I33 = 0
    E_1 = f1(r_V + p1)
    E_2 = f2(r_V - p1)

    return E_1 + E_2
    
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

# Model setup
h0 = 10e-6 #Scaling for E field movement
num_its = 500
ending_tolerance = 2 #How far from 0 to break
delta_0 = 10e-4
safety = .98
a_I = 2
p1 = np.array([10,10,0])

# Look setup
num_frames = 1
time_per_frame = 50
show_seeds = True
limit = 400
plot_x_range = [-limit,limit]
plot_y_range = [-limit,limit]
plot_z_range = [-limit,limit]
pos_color = 'red'
neg_color = 'blue'

#Seeding
title, x0, y0, z0 = seed(4) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random

title += '_2_poles_xSep:%s_ySep:%s' % (p1[0], p1[1])
#x0 = [20.]
#y0 = [20.]
#z0 = [0.]
seeds = len(x0)

#Time setup
time_steps = np.linspace(0, 2 * np.pi, num_frames)

frames = []
for t in time_steps:
    curve = {'x':[],'y':[],'z':[],'m':[], 'i':[]}

    for i in range(seeds):
        #with open("output.txt", "a") as file:
        #    file.write('______NEW SEED______ \n')
        for icity in [1,-1]:
            h = h0
            
            tot = 1 # Check to see if rk adaptive fails outside of seed points

            # Starting point of each field line
            x, y, z = x0[i], y0[i], z0[i]
            r_vec = np.array([x, y, z])
            r_start = np.array([x, y, z])

            # Accumulate segments for each line
            line_x, line_y, line_z, line_m = [], [], [], []
            
            r_past = np.zeros(3)

            for n in range(num_its // 2):
                r_mag = np.linalg.norm(r_vec)
                
                if True:                
                    #RK adaptive stepsize
                    step_sizing = True
                    while step_sizing:
                        #k1
                        E = f(r_vec)
                        k1 = h * e_solve(E, icity)

                        #k2
                        E_temp = f(r_vec + k1 * b['21'])
                        k2 = h * e_solve(E, icity)

                        #k3
                        E_temp = f(r_vec + k1 * b['31'] + k2 * b['32'])
                        k3 = h * e_solve(E, icity)

                        #k4
                        E_temp = f(r_vec + k1 * b['41'] + k2 * b['42'] + k3 * b['43'])
                        k4 = h * e_solve(E, icity)

                        #k5
                        E_temp = f(r_vec + k1 * b['51'] + k2 * b['52'] + k3 * b['53'] + k4 * b['54'])
                        k5 = h * e_solve(E, icity)

                        #k6
                        E_temp = f(r_vec + k1 * b['61'] + k2 * b['62'] + k3 * b['63'] + k4 * b['64'] + k5 * b['65'])
                        k6 = h * e_solve(E, icity)

                        delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
                        if np.linalg.norm(delta) == 0:
                            break
                            
                        if (np.abs(delta) > np.abs(delta_0)).any():
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
                        else:
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
                            r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
                            step_sizing = False
                else:
                    E = f(r_vec)

                    e_vec = e_solve(E, r_vec, r_past, icity)
                    r_change = e_vec
                
                #Wrap up
                r_mag = np.linalg.norm(r_change)

                # Fix np egien solver sign issue
                if np.dot(r_change/r_mag, r_past) < -.95:# and icity != -1:
                    r_change *= -1
                
                # Progress
                r_vec += r_change
                r_past = r_change/r_mag
                r_vec += r_change
                e_vec, e_val = e_solve(f(r_vec), icity, True, s = np.dot(r_change/r_mag, r_past), r_vec = r_vec)
                
                # Append point to line
                line_x.extend([r_vec[0]])
                line_y.extend([r_vec[1]])
                line_z.extend([0])
                
                line_m.extend([np.log10(np.abs(e_val)+1)])

                if (np.abs(np.linalg.norm(r_vec - p1)) < ending_tolerance) or (np.abs(np.linalg.norm(r_vec + p1)) < ending_tolerance): #Break since hit dipole
                    break
                elif r_mag == 0:
                    print('oopsie daisies')
                    break
                elif np.linalg.norm(r_vec) > limit:
                    break
        
            # Add segments for the current line
            curve['x'].append(line_x)
            curve['y'].append(line_y)
            curve['z'].append(line_z)
            curve['m'].append(np.array(line_m))
            curve['i'].append(icity)


print('Runtime, pre plot: %s' %(dt.now() - start))

# Normalize the width
if False:
    max_list = []
    min_list = []
    for i in curve['m']:
        i = np.array(i)
        max_list += [max(i)]
        min_list += [min(i)]

    tot_max = max(max_list)
    tot_min = min(min_list)

    for i, l in enumerate(curve['m']):
        l = (np.array(l) - tot_min) / (tot_max - tot_min)
        #urve['m'][i] = l

#for i,n in enumerate(curve['m']):
#    plt.cla()
#    plt.plot(n)
#    plt.savefig('testplots/' + str(i) + '.png')

# Set up the figure layout with animation controls
pos_fig, pos_ax = plt.subplots(1)
neg_fig, neg_ax = plt.subplots(1)

weight = True
color = False
pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('lightsteelblue'), ('navy')], N=100)
neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('mistyrose'), ('darkred')], N=100)

if weight:
    blank_pos_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('blue'), ('blue')], N=1)
    blank_neg_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('red'), ('red')], N=1)
    
    for x, y, m in zip(curve['x'], curve['y'], curve['m']):
        if curve['i'][i] == 1:
            lc = colorline(pos_ax, x, y, m, widths = m, cmap=blank_pos_cmap)
        if curve['i'][i] == -1:
            lc = colorline(neg_ax, x, y, m, widths = m, cmap=blank_neg_cmap)
elif color:    
    for i, (x, y, m) in enumerate(zip(curve['x'], curve['y'], curve['m'])):
        if curve['i'][i] == 1:
            lc = colorline(pos_ax, x, y, m, cmap=pos_cmap)
        else:
            lc = colorline(neg_ax, x, y, m, cmap=neg_cmap)
        
    pos_norm = mcolors.Normalize(vmin=0, vmax=100)
    neg_norm = mcolors.Normalize(vmin=0, vmax=-100)

    pos_fig.colorbar(cm.ScalarMappable(cmap=pos_cmap, norm = pos_norm), ax=pos_ax)
    neg_fig.colorbar(cm.ScalarMappable(cmap=neg_cmap, norm = neg_norm), ax=neg_ax)
    
elif False:
    x = np.linspace(0,10,10)
    y = np.linspace(0,10,10)
    m = np.linspace(0,1,10)
    custom_cmap = mcolors.LinearSegmentedColormap.from_list('pos', [('pink'), ('purple')], N=100)

    lc = colorline(x, y, m, cmap=custom_cmap)
else:
    
    for i, (line_x, line_y, line_m) in enumerate(zip(curve['x'], curve['y'], curve['m'])):
        if segment_colors[i] == 'red':
            ax.plot(line_x, line_y, color = 'red')

pos_ax.scatter([p1[0],-p1[0]], [p1[1],-p1[1]])
pos_ax.set_xlim(-120,120)
pos_ax.set_ylim(-120,120)
neg_ax.scatter([p1[0],-p1[0]], [p1[1],-p1[1]])
neg_ax.set_xlim(-120,120)
neg_ax.set_ylim(-120,120)
pos_fig.savefig('pos' + title + '.png', dpi=300)
neg_fig.savefig('neg' + title + '.png', dpi=300)

print('Runtime: %s' %(dt.now() - start))
