import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from datetime import datetime as dt
from seeding import seed#plane_seeds, circular_seeds, spherical_seeds, helical_seeds

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

def f(r_V): #Assuming off diagonal = 0, assuming I33 = 0
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1**2 - x_2**2)))

    E_1 = I
    E_2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            E_2[j,i] = I[0,i] * x_1 * r_V[j] + I[1,i] * x_2 * r_V[j] + I[0,j] * x_1 * r_V[i] + I[1,j] * x_2 * r_V[i]

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
    
    return E

def e_solve(E, past_e_vecs, past_e_vals):
    e_vals, e_vecs = np.linalg.eig(E)

    # Find index with smallest change in angle
    # 1. Create 3x3 array dots. Dots[i, j] is the dot between ith past e vector and jth current e vector
    # 2. Find argmax

    dots = np.zeros((3,3))
    dots_check = np.zeros((3,3))

    for i in range(3):
        for j in range(3):
            past_e_vec = past_e_vecs[:, i]
            e_vec = e_vecs[:, j]
            dots_check[i, j] = True

            dot = np.dot(past_e_vec, e_vec)
            dot = int(dot * 10e5) / 10e5
            #print(i, past_e_vec, j, e_vec, np.arccos(dot))
            if past_e_vals[i] * icity * dot > 0 and e_vals[j] * icity * dot > 0:
                dots[i, j] = np.abs(dot)
                dots_check[i, j] = False
                        

    dots = ma(dots, dots_check)
    
    index = np.argmax(dots) // 3

    r_change = past_e_vecs[:, index]# * past_e_vals[index]

    past_e_vals = e_vals
    past_e_vecs = e_vecs

    return r_change


# Model setup
h0 = 10e-2 #Scaling for E field movement
num_its = 4000
ending_tolerance = 1 #How far from 0 to break
delta_0 = 10e-3
safety = .75
a_I = .5
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
title, x0, y0, z0 = seed(0) # 0: Plane 1: Circular 2: Spherical 3: Helical
seeds = len(x0)

#Time setup
time_steps = np.linspace(0, 2 * np.pi, num_frames)

frames = []
for t in time_steps:
    
    curve_x, curve_y, curve_z = [], [], []
    curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'lg':[]}
    segment_colors, segment_widths = [], []

    for i in range(seeds):
        for icity in [1]:
            for factor in [1,-1]: #Moving forwards orbackwards
                h = h0
                tot = 1 # Check to see if rk adaptive fails outside of seed points

                # Starting point of each field line
                x, y, z = x0[i], y0[i], z0[i]
                r_vec = np.array([x, y, z])
                r_start = np.array([x, y, z])

                # Accumulate segments for each line
                line_x, line_y, line_z = [x], [y], [z]
                line_mag = []
                legend_group = '(' + str(int(x*10)/10) + ', ' + str(int(y*10)/10) + ', ' + str(int(z*10)/10) + ')'# + ' EVal Sign:' + str(icity)
                #legend_group = 'Sign: '+ str(icity) + ' Movement Factor: ' + str(factor)

                # Set index and icity
                past_e_vecs = np.zeros((3,3))
                past_e_vals = np.zeros(3) + icity
                r_change = 0
                initial = True
                
                for n in range(num_its // 2):
                    r_mag = np.linalg.norm(r_vec)
                    
                    if r_mag > 0:
                                            
                        #RK adaptive stepsize
                        step_sizing = True
                        while step_sizing:
                            #k1
                            E_temp = f(r_vec)
                            k1 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

                            #k2
                            E_temp = f(r_vec + k1 * b['21'])
                            k2 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

                            #k3
                            E_temp = f(r_vec + k1 * b['31'] + k2 * b['32'])
                            k3 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

                            #k4
                            E_temp = f(r_vec + k1 * b['41'] + k2 * b['42'] + k3 * b['43'])
                            k4 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

                            #k5
                            E_temp = f(r_vec + k1 * b['51'] + k2 * b['52'] + k3 * b['53'] + k4 * b['54'])
                            k5 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

                            #k6
                            E_temp = f(r_vec + k1 * b['61'] + k2 * b['62'] + k3 * b['63'] + k4 * b['64'] + k5 * b['65'])
                            k6 = h * e_solve(E_temp, past_e_vecs, past_e_vals)

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
                        E = f(r_vec)
                        past_e_vals, past_e_vecs = np.linalg.eig(E)

                        r_vec += r_change * factor
                        r_mag = np.linalg.norm(r_vec)
                        
                        # Append segment to line
                        line_x.extend([r_vec[0]])
                        line_y.extend([r_vec[1]])
                        line_z.extend([r_vec[2]])

                        # Color and thickness based on direction and magnitude
                        line_mag.append(r_mag)

                        if (np.abs(np.linalg.norm(r_vec)) < ending_tolerance): #Break since hit dipole
                            break
                        if r_vec[0] > 0:
                            m = r_vec[1] / r_vec[0]
                            if m - .43 < .2 or m - 2.2 < .2:
                                break
            
                # Add segments for the current line
                curve['x'].append(line_x)
                curve['y'].append(line_y)
                curve['z'].append(line_z)
                curve['m'].append(line_mag)
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

if True:
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
        fig.add_trace(go.Scatter3d(x = x0, y = y0, z = z0, name='Seed points', line=dict(width=0)))

    #fig.write_html("%s.html" %title)
    fig.show()

    print('Runtime: %s' %(dt.now() - start))