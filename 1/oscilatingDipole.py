import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from datetime import datetime as dt
from seeding import seed
from plot import plot

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

def g1(x):
    return (2.27135912155) * x + 0.24483709465
def g2(x):
    return (2.26291473703) * x - 0.179012299421
def g3(x):
    return -1/(2.27135912155) * x + 0.24483709465
def g4(x):
    return -1/(2.26291473703) * x - 0.179012299421
def rot(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

def f_old(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
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


    E = (-6 * E_1 / (r**5)) + (30 * E_2 / (r**7)) + (15 * E_3 / (r**7)) + (-105 * E_4 / (r**9))
    
    if give:
        return E, E_1, E_2, E_3, E_4
    return E

def f(r_V): #Assuming off diagonal = 0, assuming I33 = 0
    if False:
        p = 5
        p1 = np.array([p,0,0])
        p2 = np.array([-p,0,0])

        E_1 = f_old(r_V - p1)
        E_2 = f_old(r_V - p2)

        return E_1 + E_2

    else:
        return f_old(r_V)
    
x_deg = []
y_deg = []
z_deg = []

def e_solve(E, r_vec, r_past, icity, mag_return = False): # r_past dotted with e vects
    global x_deg
    global y_deg
    global z_deg
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    e0 = int(e_vals[0] * 10e3) / 10e3
    e1 = int(e_vals[1] * 10e3) / 10e3
    e2 = int(e_vals[2] * 10e5) / 10e5

    dots = np.zeros(3)
    dots_check = np.zeros(3)

    for i in range(3):
        e_vec = e_vecs[:, i]
        dots_check[i] = True

        dot = np.dot(r_past, e_vec)
        dot = int(dot * 10e9) / 10e9

        if e_vals[i] * icity > 0:
            dots[i] = np.abs(dot)
            dots_check[i] = False

    dots = ma(dots, dots_check)
    
    index = np.argmax(dots)

    x_t = r_vec[0]
    y_t = r_vec[1]
    #if (g2(x_t) < y_t and y_t < g1(x_t)) or (g4(x_t) < y_t and y_t < g3(x_t)):
    #    print('Next good')
    #else:
    #    for i in f_old(r_vec, True):
    #        print(i)#

    #    print(e_vals, index)
    #print()
#
    r_change = e_vecs[:, index]

    if mag_return:
        mag = e_vals[index]
        return r_change, mag
    else:
        return r_change


# Model setup
h0 = 10e-2 #Scaling for E field movement
num_its = 4000
ending_tolerance = 1 #How far from 0 to break
delta_0 = 10e-3
safety = .75
a_I = 200
I = np.array([[a_I, 0, 0], 
              [0, -1*a_I, 0],
              [0, 0, 0]])

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
title, x0, y0, z0 = seed(0) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random
#x0 = np.linspace(0,10, 20)
#y0 = np.column_stack((g1(x0), g2(x0)))
#y0 = np.mean(y0,1)
#y0 = g2(x0)
#z0 = np.zeros(20)
seeds = len(x0)

#Time setup
time_steps = np.linspace(0, 2 * np.pi, num_frames)

frames = []
for t in time_steps:
    
    curve_x, curve_y, curve_z = [], [], []
    curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'c':[],'lg':[]}
    segment_colors, segment_widths = [], []

    for i in range(seeds):
        for icity in [1,-1]:
            for factor in [1,-1]: #Moving forwards orbackwards
                h = h0
                tot = 1 # Check to see if rk adaptive fails outside of seed points

                # Starting point of each field line
                x, y, z = x0[i], y0[i], z0[i]
                r_vec = np.array([x, y, z])
                r_start = np.array([x, y, z])

                # Accumulate segments for each line
                line_x, line_y, line_z, line_c, line_m = [], [], [], [], []

                #legend_group = '(' + str(int(x*10)/10) + ', ' + str(int(y*10)/10) + ', ' + str(int(z*10)/10) + ')' + ' EVal Sign:' + str(icity)
                legend_group = 'Sign: '+ str(icity) + ' Movement Factor: ' + str(factor)
                
                r_past = np.zeros(3)

                for n in range(num_its // 2):
                    r_mag = np.linalg.norm(r_vec)
                    
                    if r_mag > 0:
                                            
                        #RK adaptive stepsize
                        step_sizing = True
                        while step_sizing:
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
                        
                        #Wrap up
                        r_vec += r_change * factor
                        r_past = r_change
                        _, e_val = e_solve(f(r_vec), r_vec, r_past, icity, True)
                        
                        # Append point to line
                        x_t = r_vec[0]
                        y_t = r_vec[1]
                        z_t = 0#np.log10(np.abs(e_val))
                        line_x.extend([x_t])
                        line_y.extend([y_t])
                        line_z.extend([z_t])

                        #if (g2(x_t) < y_t and y_t < g1(x_t)) or (g4(x_t) < y_t and y_t < g3(x_t)):
                        #    x_deg += [x_t]
                        #    y_deg += [y_t]
                        #    z_deg += [z_t]
                            
                        #line_z.extend([np.log10(np.abs(e_val))])

                        if icity == 1:
                            line_c.extend(['b'])
                        else:
                            line_c.extend(['r'])
                        line_m.extend([np.log10(np.abs(e_val)+1)])

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
    fig.add_trace(go.Scatter3d(x=x_deg, y=y_deg, z=z_deg, name = 'in zone', line=dict(color='black', width=0), showlegend=True))
    fig.add_trace(go.Scatter3d(x=x, y=g3(x), z=z, name = 'g1', mode="lines", line=dict(color='black'), showlegend=True))
    fig.add_trace(go.Scatter3d(x=x, y=g4(x), z=z, name = 'g2', mode="lines", line=dict(color='black'), showlegend=True))
    fig.add_trace(go.Scatter3d(x=(x + g3(x)) / np.sqrt(2), y=(g3(x) - x) / np.sqrt(2), z=z, name = 'g1`', mode="lines", line=dict(color='black'), showlegend=True))
    fig.add_trace(go.Scatter3d(x=(x + g4(x)) / np.sqrt(2), y=(g4(x) - x) / np.sqrt(2), z=z, name = 'g2`', mode="lines", line=dict(color='black'), showlegend=True))
        #return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    
    #fig.add_trace(go.Scatter3d(x=x, y=g3(x), z=z, name = 'g3', mode="lines", line=dict(color='black'), showlegend=True))
    #fig.add_trace(go.Scatter3d(x=x, y=g4(x), z=z, name = 'g4', mode="lines", line=dict(color='black'), showlegend=True))

#fig.write_html("%s.html" %title)
fig.show()

print('Runtime: %s' %(dt.now() - start))
