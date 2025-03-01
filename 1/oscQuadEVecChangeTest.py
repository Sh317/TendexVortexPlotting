import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from datetime import datetime as dt
from seeding import seed
from scipy.stats import linregress as lr
#from plot import plot
import matplotlib.pyplot as plt
from eigenvectorTest import generate_grid_from_function
import os


#Runtime
start = dt.now()

# RK Coeffecients
a = {'2':1/5,         '3':3/10,     '4':3/5,         '5':1,            '6': 7/8}
b = {'21':1/5, 
     '31':3/40,       '32':9/40, 
     '41':3/10,       '42':-9/10,   '43':6/5, 
     '51':-11/54,     '52':5/2,     '53':-70/27,    '54': 35/27, 
     '61':1631/55296, '62':175/512, '63':575/13824, '64':44275/110592, '65': 253/4096}
c_dif = {'1':37/378 - 2825/27648, '3':250/621 - 18575/48384, '4':125/621 - 13525/55296, '5':-277/14336, '6':512/1771 - 1/4}
c = {'1':37/378, '3':250/621, '4':125/621, '6':512/1771}

def dist1P(x,y):
    return ((2.2710024192177376) * x + 0.029466068242026466 - y) / np.sqrt(2.2710024192177376**2 + 1)
def dist2P(x,y):
    return (-.439579 * x - y) / np.sqrt(.439579)

def g1P(x):
    return (2.2749) * x
def g2P(x):
    return -.439579 * x
def g1N(x):
    return 0.44212777173202983 * x - 0.023953870557449797
def g2N(x):
    return -1/(0.44212777173202983) * x - 0.023953870557449797
def g1(x):
    return (0.441294996407) * x - 0.106966301356
def g2(x):
    return (0.45267421941) * x + 0.00496808575137
def g3(x):
    return -1/(0.441294996407) * x - 0.106966301356
def g4(x):
    return -1/(0.45267421941) * x + 0.00496808575137

def f(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
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
        return E
    return E

def e_solve(E, r_vec, r_past, icity, mag_return = False, give= False): # r_past dotted with e vects
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    #dots = np.zeros(3)
    #dots_check = np.zeros(3)
    

    #for i in range(3):
    #    e_vec = e_vecs[:, i]
    #    dots_check[i] = True

    #    dot = np.dot(r_past, e_vec)
    #    dot = int(dot * 10e9) / 10e9

    #    if e_vals[i] * icity > 0:
    #        dots[i] = np.abs(dot)
    #        dots_check[i] = False

    #dots = ma(dots, dots_check)
    
    #index = np.argmax(dots)
    #print(e_vals)
    e_valcheck = np.ma.greater(e_vals * -1 * icity, 0)

    index = np.argmin(e_valcheck)
    if e_vals[index] < 0:
        print('sfsfgasfga')

    r_change = e_vecs[:, index]

    if mag_return:
        #mag = e_vals[index]
        return r_change, e_vals, index, e_vecs
    elif give:
        return r_change, e_vals, index
    else:
        return r_change

# Model setup
h0 = 10e9 #Scaling for E field movement
num_its = 10000
ending_tolerance = 1 #How far from 0 to break
delta_0 = 10e-3
safety = .9
a_I = 1
I = np.array([[a_I, 0, 0], 
              [0, -a_I, 0],
              [0, 0, 0]])

# Look setup
num_frames = 1
time_per_frame = 50
show_seeds = True
limit = 500
plot_x_range = [-limit,limit]
plot_y_range = [-limit,limit]
plot_z_range = [-limit,limit]
pos_color = 'red'
neg_color = 'blue'

mags_on = []
mags_off = []
dist_on = []
dist_off = []


#Seeding
title, x0, y0, z0 = seed(0) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random

x_deg = []
y_deg = []
z_deg = []
#x0, y0, z0 = generate_grid_from_function(g1P, -30, 30, 8, .5)
#x0 = np.linspace(0,10, 100)
#x0 = [5, 20, 50]
y0 = [-2.1889632013843534]#, g2P(5)+.01]
#x0 = [7.5]
x0 = [4.979681851608]
#y0 = [-2.5]
#x0 = [-2.5]
#y0 = [7.5]
#x0 = [-7.5]
#y0 = [-12.5]
z0 = [0]
#y0 = g1(x0)
#z0 = np.zeros(100)
seeds = len(x0)
independent = False
#Time setup
time_steps = np.linspace(0, 2 * np.pi, num_frames)
seed_names = []
frames = []

dist = []
delt_e_vec = []
past_e_vec = []

for t in time_steps:
    
    curve_x, curve_y, curve_z = [], [], []
    curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'c':[],'lg':[]}
    segment_colors, segment_widths = [], []

    for i in range(seeds):
        for icity in [1]:
            for factor in [1]: #Moving forwards orbackwards
                h = h0
                tot = 1 # Check to see if rk adaptive fails outside of seed points

                # Starting point of each field line
                x, y, z = x0[i], y0[i], z0[i]
                r_vec = np.array([x, y, z])
                r_start = np.array([x, y, z])
                initial = True

                # Accumulate segments for each line
                line_x, line_y, line_z, line_c, line_m = [], [], [], [], []

                #legend_group = '(' + str(int(x*10)/10) + ', ' + str(int(y*10)/10) + ', ' + str(int(z*10)/10) + ')' + ' EVal Sign:' + str(icity)
                legend_group = 'Sign: '+ str(icity) + ' Movement Factor: ' + str(factor)
                seed_names += ['x: %s y: %s z: %s' % (x, y, z)]
                
                r_past = np.zeros(3)

                for n in range(num_its // 2):
                    r_mag = np.linalg.norm(r_vec)
                    
                    if r_mag > 0:
                                            
                        #RK adaptive stepsize
                        step_sizing = True
                        if False:
                            while step_sizing:
                                if independent:
                                    k1 = h * f(r_vec)
                                    k2 = h * f(r_vec + a['2'] * h)
                                    k3 = h * f(r_vec + a['3'] * h)
                                    k4 = h * f(r_vec + a['4'] * h)
                                    k5 = h * f(r_vec + a['5'] * h)
                                    k6 = h * f(r_vec + a['6'] * h)

                                else:
                                    #k1
                                    E = f(r_vec)
                                    k1 = h * e_solve(E, r_vec, r_past, icity)

                                    #k2
                                    E_temp = f(r_vec + k1 * b['21'])
                                    k2 = h * e_solve(E_temp, r_vec, r_past, icity)

                                    #k3
                                    E_temp = f(r_vec + k1 * b['31'] + k2 * b['32'])
                                    k3 = h * e_solve(E_temp, r_vec, r_past, icity)

                                    #k4
                                    E_temp = f(r_vec + k1 * b['41'] + k2 * b['42'] + k3 * b['43'])
                                    k4 = h * e_solve(E_temp, r_vec, r_past, icity)

                                    #k5
                                    E_temp = f(r_vec + k1 * b['51'] + k2 * b['52'] + k3 * b['53'] + k4 * b['54'])
                                    k5 = h * e_solve(E_temp, r_vec, r_past, icity)

                                    #k6
                                    E_temp = f(r_vec + k1 * b['61'] + k2 * b['62'] + k3 * b['63'] + k4 * b['64'] + k5 * b['65'])
                                    k6 = h * e_solve(E_temp, r_vec, r_past, icity)

                                delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
                                if np.linalg.norm(delta) == 0:
                                    tot += 1
                                    if tot > seeds:
                                        print('??', tot)
                                    break
                                    
                                if (np.abs(delta) > np.abs(delta_0)).any():
                                    h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
                                else:
                                    h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
                                    r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
                                    step_sizing = False

                        elif False:
                            E = f(r_vec)
                            k1 = h * e_solve(E, r_vec, r_past, icity)

                            E = f(r_vec + k1/2)
                            k2 = h * e_solve(E, r_vec, r_past, icity)

                            E = f(r_vec + k2/2)
                            k3 = h * e_solve(E, r_vec, r_past, icity)

                            E = f(r_vec + k3)
                            k4 = h * e_solve(E, r_vec, r_past, icity)

                            r_change =  k1/6 + k2/3 + k3/3 + k4/6
                        
                        else:
                            h = 10e-11#np.random.rand() * 10
                            E = f(r_vec)
                            #if np.si
                            e_vec = e_solve(E, r_vec, r_past, icity)
                            r_change = h * e_vec
                            _, e_vals, index, e_vecs = e_solve(E, r_vec, r_past, icity, True)
                            #os.system('cls' if os.name == 'nt' else 'clear')
                            print('e_vecs: %s' % e_vecs)
                            print('e_vals: %s' % e_vals)
                            print('x: %s y: %s' % (r_vec[0], r_vec[1]))
                            print('dist: %s' % dist2P(r_vec[0], r_vec[1]))

                            dist += [dist2P(r_vec[0], r_vec[1])]
                            if initial:
                                initial = False
                                delt_e_vec += [0]
                            else:
                                
                                delt_e_vec += [np.linalg.norm(e_vec - past_e_vec)]
                            past_e_vec = e_vec
                            
                        #Wrap up
                        if independent:
                            r_change = e_solve(r_change, r_vec, r_past, icity)
                        r_vec += r_change * factor
                        r_past = r_change
                        _, mags, index, vecs = e_solve(f(r_vec), r_vec, r_past, icity, True)
                        E_T = f(r_vec, True)
                        e_val = mags[index]
                        
                        # Append point to line
                        x_t = r_vec[0]
                        y_t = r_vec[1]
                        E_T = np.linalg.norm(E_T)
                        z_t = 0
                        line_x.extend([x_t])
                        line_y.extend([y_t])
                        line_z.extend([z_t])

                        if np.abs(g1P(x_t) - y_t) < 1.5 or np.abs(g2P(x_t) - y_t) < 1.5:
                            mags_on += [E_T]
                            #print(vecs[index], r_change)
                            
                            if dist1P(x_t, y_t) < dist2P(x_t, y_t):
                                dist_on += [dist1P(x_t, y_t)]
                            else:
                                dist_on += [dist2P(x_t, y_t)]
                            x_deg += [x_t]
                            y_deg += [y_t]
                            z_deg += [0]
                        else:
                            #print(vecs)
                            mags_off += [E_T]
                            if dist1P(x_t, y_t) < dist2P(x_t, y_t):
                                dist_off += [dist1P(x_t, y_t)]
                            else:
                                dist_off += [dist2P(x_t, y_t)]

                        #line_z.extend([np.log10(np.abs(e_val))])

                        if icity == 1:
                            line_c.extend(['b'])
                        else:
                            line_c.extend(['r'])
                        line_m.extend([1])

                        if (np.abs(np.linalg.norm(r_vec)) < ending_tolerance): #Break since hit dipole
                            break
            
                # Add segments for the current line
                curve['x'].append(line_x)
                curve['y'].append(line_y)
                curve['z'].append(line_z)
                curve['m'].append(np.array(line_m))
                curve['c'].append(line_c)
                curve['f'].append(factor)
                curve['lg'].append(legend_group)
                if icity == 1:
                    segment_colors += [pos_color]
                else:
                    segment_colors += [neg_color]

    # Add data for this frame
    frame_data = [
        go.Scatter3d(
            x=line_x,
            y=line_y,
            z=line_z,
            mode="lines",
            line=dict(color=segment_colors[i]),
            showlegend=True,
            legendgroup=curve['lg'][i],
            legendgrouptitle_text=curve['lg'][i],
            name = ('Fwd' if curve['f'][i] > 0 else 'Bwd')
        ) for i, (line_x, line_y, line_z) in enumerate(zip(curve['x'], curve['y'], curve['z']))
    ]

    frames.append(go.Frame(data=frame_data, name=f"frame{t}"))

print('Runtime, pre plot: %s' %(dt.now() - start))

# Normalize the width
#curve['m'] = np.array(curve['m'])
tot_min = 10e5
tot_max = 0
for i in curve['m']:
    if min(np.abs(i)) < tot_min:
        tot_min = min(np.abs(i))
    if max(np.abs(i)) > tot_max:
        tot_max = max(np.abs(i))

for i, l in enumerate(curve['m']):
    l = (l - tot_min) / (tot_max - tot_min)
    l = (l * 9) + 1
    curve['m'][i] = l

# Set up the figure layout with animation controls
if True:
    fig = go.Figure(
        data=frame_data,
        layout=go.Layout(
            title=title,
            scene=dict(
                xaxis_title="x",
                yaxis_title="y",
                zaxis_title="z",
                aspectratio=dict(x=1, y=1, z=1),
                xaxis_range=plot_x_range,
                yaxis_range=plot_y_range,
                zaxis_range=plot_z_range
            ),
            width=1000,
            height=800,
            updatemenus=[
                dict(
                    type="buttons",
                    showactive=False,
                    buttons=[
                        dict(
                            label="Play",
                            method="animate",
                            args=[
                                None,
                                dict(
                                    frame=dict(duration=200, redraw=True),
                                    fromcurrent=True,
                                    mode="immediate",
                                )
                            ],
                        ),
                        dict(
                            label="Pause",
                            method="animate",
                            args=[
                                [None],
                                dict(
                                    frame=dict(duration=0, redraw=False),
                                    mode="immediate",
                                    transition=dict(duration=0),
                                )
                            ],
                        ),
                        dict(
                            label="Reset",
                            method="animate",
                            args=[
                                [None],
                                dict(
                                    frame=dict(duration=0, redraw=True),
                                    mode="immediate",
                                    transition=dict(duration=0),
                                )
                            ],
                        )
                    ],
                )
            ],
            sliders=[{
                "active": 0,
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {
                    "font": {"size": 16},
                    "prefix": "Time: ",
                    "visible": True,
                    "xanchor": "right"
                },
                "transition": {"duration": time_per_frame},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [[f"frame{time}"], {"frame": {"duration": 100, "redraw": True}, "mode": "immediate"}],
                        "label": f"{time:.2f}",
                        "method": "animate"
                    } for time in time_steps
                ]
            }]
        ),
        frames=frames
    )

    if show_seeds:
        x = np.linspace(-100,100,200)
        z = np.zeros(1000)

        fig.add_trace(go.Scatter3d(x = x0, y = y0, z = z0, name='Seed points', line=dict(width=0)))
        fig.add_trace(go.Scatter3d(x = x_deg, y = y_deg, z = z_deg, name='Seed points', line=dict(width=0)))
        
        fig.add_trace(go.Scatter3d(x=x, y=g1P(x), z=z, name = 'g1P', mode="lines", line=dict(color='black'), showlegend=True))
        fig.add_trace(go.Scatter3d(x=x, y=g2P(x), z=z, name = 'g2P', mode="lines", line=dict(color='black'), showlegend=True))
        fig.add_trace(go.Scatter3d(x=x, y=g1N(x), z=z, name = 'g1N', mode="lines", line=dict(color='black'), showlegend=True))
        fig.add_trace(go.Scatter3d(x=x, y=g2N(x), z=z, name = 'g2N', mode="lines", line=dict(color='black'), showlegend=True))
    #fig.write_html("%s.html" %title)
    fig.show()

x_list = curve['x']
y_list = curve['y']
z_list = curve['z']
c_list = curve['c']

# Create the figure and axis
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#for x, y, z, c, name in zip(x_list, y_list, z_list, c_list, seed_names):
    #ax.plot(x,y,z, label = name)
#ax.plot(x_deg, y_deg, z_deg)
#fig.legend()
plt.plot(dist, delt_e_vec)
plt.xlabel('distance from the line')
plt.ylabel('magnitude of eigenvector change')
plt.show()

#plt.scatter(dist_off, mags_off, label = 'off', s = 1)
#plt.scatter(dist_on, mags_on, label = 'on', s = 1)
#plt.legend()
#plt.xlabel('dist')
#plt.ylabel('mag')

plt.show()


print('Runtime: %s' %(dt.now() - start))
