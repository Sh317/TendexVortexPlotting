import numpy as np
from numpy.ma import masked_array as ma
import plotly.graph_objects as go
from seeding import seed
from datetime import datetime as dt

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

def f(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1**2 - x_2**2)))

    E_1 = I
    #E_2 = np.zeros((3,3))
    #for i in range(3):
    #    for j in range(3):
    #        E_2[i,j] = I[0,i] * x_1 * r_V[j] + I[1,i] * x_2 * r_V[j] + I[0,j] * x_1 * r_V[i] + I[1,j] * x_2 * r_V[i]
    E_2 = np.array([[2 * a_I * x_1 * x_1, 0, 0], 
                    [0,-2 * a_I * x_2 * x_2, 0],
                    [0, 0, 0]])
    E_3 = np.identity(3) * IpqXpXq
    E_4 = np.zeros((3,3))
    for i in range(len(I)):
        for j in range(len(I)):
            E_4[j,i] = (r_V[i] * r_V[j])
    
    E_4 *= IpqXpXq
    E_4 = np.array([[IpqXpXq * x_1 * x_1, IpqXpXq * x_1 * x_2, 0], 
                    [IpqXpXq * x_1 * x_2, IpqXpXq * x_2 * x_2, 0],
                    [0, 0, 0]])

   


    E = (-6 * E_1) + (30 * E_2 / (r**2)) + (15 * E_3 / (r**2)) + (-105 * E_4 / (r**4))
    
    if give:
        return E
    return E

def e_solve(E, r_vec, r_past, icity, mag_return = False, give= False): # r_past dotted with e vects
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    e_valcheck = np.ma.greater(e_vals * -1 * icity, 0)

    index = np.argmin(e_valcheck)

    r_change = e_vecs[:, index]

    if mag_return:
        #mag = e_vals[index]
        return r_change, e_vals, index, e_vecs
    elif give:
        return r_change, e_vals, index
    else:
        return r_change

title, x0, y0, z0 = seed(4) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random
x0 = [5.0]
y0 = [5.0]
z0 = [0.0]
curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'c':[],'lg':[]}
seeds = len(x0)
num_its = 500
delta_0 = 10e-2
h0 = 10e-2
safety = .9
a_I = 1
I = np.array([[a_I, 0, 0], 
              [0, -a_I, 0],
              [0, 0, 0]])
time_steps = np.linspace(0, 2 * np.pi, 1)
pos_color = 'red'
neg_color = 'blue'

for i in range(seeds):
    # Starting point of each field line
    x, y, z = x0[i], y0[i], z0[i]
    r_vec = np.array([x, y, z])
    r_start = np.array([x, y, z])

    line_x, line_y, line_z, line_c, line_m = [], [], [], [], []

    r_past = np.zeros(3)
    r_change = np.zeros(3)
    for icity in [1]:
        for n in range(num_its // 2):
            step_sizing = True
            h = h0
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

            r_mag = np.linalg.norm(r_change)
            if np.dot(r_change/r_mag, r_past) < -.99:
                r_change *= -1
            print(r_change)
            r_vec += r_change
            r_past = r_change/r_mag

            line_x.extend([r_vec[0]])
            line_y.extend([r_vec[1]])
            line_z.extend([r_vec[2]])

            line_m.extend([1])

            if (np.abs(np.linalg.norm(r_vec)) < .5): #Break since hit dipole
                break

    curve['x'].append(line_x)
    curve['y'].append(line_y)
    curve['z'].append(line_z)
    curve['m'].append(line_m)
    #curve['lg'] += [icity]
    if icity == 1:
        curve['c'] += [pos_color]
    else:
        curve['c'] += [neg_color]

frames = []

frame_data = [
    go.Scatter3d(
        x=line_x,
        y=line_y,
        z=line_z,
        mode="lines",
        line=dict(color=curve['c'][i]),
        showlegend=False
        #legendgroup=curve['lg'][i]
        #legendgrouptitle_text=curve['lg'][i],
        #name = ('Fwd' if curve['f'][i] > 0 else 'Bwd')
    ) for i, (line_x, line_y, line_z) in enumerate(zip(curve['x'], curve['y'], curve['z']))
]

frames.append(go.Frame(data=frame_data, name=f"frame"))

fig = go.Figure(
    data=frame_data,
    layout=go.Layout(
        title=title,
        scene=dict(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title="z",
            aspectratio=dict(x=1, y=1, z=1),
            xaxis_range=[-400,400],
            yaxis_range=[-400,400],
            zaxis_range=[-400,400]
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
            "transition": {"duration": 30},
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

print(dt.now() - start)

fig.show()