import numpy as np
import plotly.graph_objects as go
from datetime import datetime as dt

#Runtime
start = dt.now()

#Constants
omega = .1 * 2 * np.pi
k = omega
lam = np.pi / omega

# RK Coeffecients
a = {'2':1/5,         '3':3/10,     '4':3/5,         '5':1,            '6': 7/8}
b = {'21':1/5, 
     '31':3/40,       '32':9/40, 
     '41':3/10,       '42':-9/10,   '43':6/5, 
     '51':-11/54,     '52':5/2,     '53':-70/27,    '54': 35/27, 
     '61':1631/55296, '62':175/512, '63':575/13824, '64':44275/110592, '65': 253/4096}
c_dif = {'1':37/378 - 2825/27648, '3':250/621 - 18575/48384, '4':125/621 - 13525/55296, '5':-277/14336, '6':512/1771 - 1/4}
c = {'1':37/378, '3':250/621, '4':125/621, '6':512/1771}

# The near (static) zone: d << r << lam
# The intermediate (induction) zone: d << r ~ lam
# The far (radiation) zone: d << lam << r

def random_seeds(seeds, distance):
    x_list = []
    y_list = []
    z_list = []

    for i in range(seeds//4):
        x = np.random.random() * distance
        y = np.random.random() * distance
        z = 0

        # Quadrant 1
        x_list += [x]
        y_list += [y]
        z_list += [z]#[np.random.random() * distance - distance / 2]

        # Quadrant 2
        x_list += [-1 * x]
        y_list += [y]
        z_list += [z]#[np.random.random() * distance - distance / 2]

        # Quadrant 3
        x_list += [-1 * x]
        y_list += [-1 * y]
        z_list += [z]#[np.random.random() * distance - distance / 2]

        # Quadrant 4
        x_list += [x]
        y_list += [-1 * y]
        z_list += [z]#[np.random.random() * distance - distance / 2]
        if False:
            # Quadrant 5
            x_list += [x]
            y_list += [y]
            z_list += [-1 * z]#[np.random.random() * distance - distance / 2]

            # Quadrant 6
            x_list += [-1 * x]
            y_list += [y]
            z_list += [-1 * z]#[np.random.random() * distance - distance / 2]

            # Quadrant 7
            x_list += [-1 * x]
            y_list += [-1 * y]
            z_list += [-1 * z]#[np.random.random() * distance - distance / 2]

            # Quadrant 8
            x_list += [x]
            y_list += [-1 * y]
            z_list += [-1 * z]#[np.random.random() * distance - distance / 2]

    return x_list, y_list, z_list

def f(r_V, p1, p2, factor):
    r = np.linalg.norm(r_V)
    n = r_V / r

    E_vec1 = k**2 * np.cross(np.cross(n, p1), n) * np.exp(1j * k * r) / r
    E_vec1 += (3 * n * np.dot(n, p1) - p1) * (1/ r**3 - (1j * k / r**2)) * np.exp(1j * k * r)

    E_vec2 = k**2 * np.cross(np.cross(n, p2), n) * np.exp(1j * k * r) / r
    E_vec2 += (3 * n * np.dot(n, p2) - p2) * ((1/ r**3) - (1j * k / r**2)) * np.exp(1j * k * r)

    E = (1 / (4 * np.pi)) * (E_vec1 + E_vec2)
    return E * factor

# Model setup
seed_type = 3 # 0: Plane 1: Circular 2: Spherical
h0 = 10e-3 #Scaling for E field movement
ending_tolerance = .05 #How far from 0 to break
delta_0 = 10e-2
safety = .95

# Look setup
num_frames = 5
time_per_frame = 50
show_seeds = True
rotate = True
limit = 120
plot_x_range = [-limit,limit]
plot_y_range = [-limit,limit]
plot_z_range = [-limit,limit]
pos_color = 'red'
neg_color = 'blue'


#Dipole setup
p1_initial = np.array([np.sqrt(2)/2, np.sqrt(2)/2, 0])
p2_initial = np.array([-np.sqrt(2)/2, np.sqrt(2)/2, 0])

seed_num = 300
distance = 10
num_its = 10
limit = distance + 5
plot_x_range = [-limit,limit]
plot_y_range = [-limit,limit]
plot_z_range = [-limit,limit]
if rotate:
    title = "Rotating dipole. Random seeding. Lambda: %s. RK Adaptive" % (lam)
else:
    title = "Two dipole. Random seeding. Lambda: %s. RK Adaptive" % (lam)

#Time setup
time_steps = np.linspace(0, np.pi, num_frames)

frames = []
x0, y0, z0 = random_seeds(seed_num,distance)
seeds = len(x0)
p1_past = 0
p2_past = 0
for t in time_steps:
    p1 = p1_initial * np.exp(1j * omega * t)
    p2 = p2_initial * np.exp(1j * omega * t)
    
    if rotate:
        p2 *= 1j
    
    curve_x, curve_y, curve_z = [], [], []
    curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'lg':[]}
    segment_colors, segment_widths = [], []

    for i in range(seeds):
        break_factor = False
        for factor in [1, -1]: #Moving forwards orbackwards
            if factor == -1 and break_factor:
                break
            h = h0

            # Starting point of each field line
            x, y, z = x0[i], y0[i], z0[i]
            x_start, y_start, z_start = x, y, z
            its_run = 0
            r_vec = np.array([x, y, z])
            r_start = np.array([x, y, z])

            # Accumulate segments for each line
            line_x, line_y, line_z = [x], [y], [z]
            line_mag = []
            legend_group = '(' + str(int(x*10)/10) + ', ' + str(int(y*10)/10) + ', ' + str(int(z*10)/10) + ')'

            for n in range(num_its // 2):
                r_mag = np.linalg.norm(r_vec)
                
                if r_mag > 0:
                    step_sizing = True
                    r_change = 0

                    #RK adaptive stepsize
                    while step_sizing:
                        k1 = h * f(r_vec, p1, p2, factor)
                        k2 = h * f(r_vec + b['21'] * k1, p1, p2, factor)
                        k3 = h * f(r_vec + b['31'] * k1 + b['32'] * k2, p1, p2, factor)
                        k4 = h * f(r_vec + b['41'] * k1 + b['42'] * k2 + b['43'] * k3, p1, p2, factor)
                        k5 = h * f(r_vec + b['51'] * k1 + b['52'] * k2 + b['53'] * k3 + b['54'] * k4, p1, p2, factor)
                        k6 = h * f(r_vec + b['61'] * k1 + b['62'] * k2 + b['63'] * k3 + b['64'] * k4 + b['65'] * k5, p1, p2, factor)

                        delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
                        if (np.abs(delta) > np.abs(delta_0)).any():
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
                        else:
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
                            r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
                            step_sizing = False

                    
                    #Wrap up
                    r_vec += np.real(r_change)
                    r_mag = np.linalg.norm(f(r_vec, p1, p2, 1))
                    
                    # Append segment to line
                    line_x.extend([r_vec[0]])
                    line_y.extend([r_vec[1]])
                    line_z.extend([r_vec[2]])

                    # Color and thickness based on direction and magnitude
                    line_mag.append(r_mag)

                    its_run += 1

                    if (np.abs(r_vec) < ending_tolerance).all(): #Break since hit dipole
                        break
                    if (np.abs(r_vec-r_start) < ending_tolerance).all() and n > 50: #Break since far zone
                        break_factor = True
                        break

            h = h0

            # Starting point of each field line
            fact_temp = -1 * factor
            r_vec_end = r_vec.copy()
            r_vec = r_vec_end

            r_dif = 10**10
            tot_away = 0
            base_dif = 0

            found_end = False

            while not found_end:
                r_mag = np.linalg.norm(r_vec)
                
                if r_mag > 0:
                    step_sizing = True
                    r_change = 0

                    #RK adaptive stepsize
                    while step_sizing:
                        k1 = h * f(r_vec, p1, p2, fact_temp)
                        k2 = h * f(r_vec + b['21'] * k1, p1, p2, fact_temp)
                        k3 = h * f(r_vec + b['31'] * k1 + b['32'] * k2, p1, p2, fact_temp)
                        k4 = h * f(r_vec + b['41'] * k1 + b['42'] * k2 + b['43'] * k3, p1, p2, fact_temp)
                        k5 = h * f(r_vec + b['51'] * k1 + b['52'] * k2 + b['53'] * k3 + b['54'] * k4, p1, p2, fact_temp)
                        k6 = h * f(r_vec + b['61'] * k1 + b['62'] * k2 + b['63'] * k3 + b['64'] * k4 + b['65'] * k5, p1, p2, fact_temp)

                        delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
                        if (np.abs(delta) > np.abs(delta_0)).any():
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
                        else:
                            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
                            r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
                            step_sizing = False

                    
                    #Wrap up
                    r_vec += np.real(r_change)

                    tot_away += 1
                    if np.linalg.norm(r_vec - r_start) < r_dif:
                        r_dif = np.linalg.norm(r_vec - r_start)
                        tot_away = 0
                    elif tot_away > 2:
                        base_dif = (np.linalg.norm((r_vec - r_start)))
                        break
            
            if t > 0:
                r_vec = r_vec_end

                r_dif = 10**10
                tot_away = 0

                found_end = False

                while not found_end:
                    r_mag = np.linalg.norm(r_vec)
                    
                    if r_mag > 0:
                        step_sizing = True
                        r_change = 0

                        #RK adaptive stepsize
                        while step_sizing:
                            k1 = h * f(r_vec, p1_past, p2_past, fact_temp)
                            k2 = h * f(r_vec + b['21'] * k1, p1_past, p2_past, fact_temp)
                            k3 = h * f(r_vec + b['31'] * k1 + b['32'] * k2, p1_past, p2_past, fact_temp)
                            k4 = h * f(r_vec + b['41'] * k1 + b['42'] * k2 + b['43'] * k3, p1_past, p2_past, fact_temp)
                            k5 = h * f(r_vec + b['51'] * k1 + b['52'] * k2 + b['53'] * k3 + b['54'] * k4, p1_past, p2_past, fact_temp)
                            k6 = h * f(r_vec + b['61'] * k1 + b['62'] * k2 + b['63'] * k3 + b['64'] * k4 + b['65'] * k5, p1_past, p2_past, fact_temp)

                            delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
                            if (np.abs(delta) > np.abs(delta_0)).any():
                                h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
                            else:
                                h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
                                r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
                                step_sizing = False

                        #Wrap up
                        r_vec += np.real(r_change)

                        tot_away += 1
                        
                        if np.linalg.norm(r_vec - r_start) < r_dif:
                            r_dif = np.linalg.norm(r_vec - r_start)
                            tot_away = 0
                        elif tot_away > 2:
                            #print()
                            #print((np.linalg.norm((r_vec - r_start)/r_vec)))
                            break

            # Add segments for the current line
            curve['x'].append(line_x)
            curve['y'].append(line_y)
            curve['z'].append(line_z)
            curve['m'].append(line_mag)
            curve['f'].append(factor)
            curve['lg'].append(legend_group)
    
    # Color everything
    for i, (line_x, line_y, line_z) in enumerate(zip(curve['x'], curve['y'], curve['z'])):
        # Calculate beginning of line vector
        try:
            l = np.array([line_x[5] - line_x[0], line_y[5] - line_y[0], line_z[5] - line_z[0]])
        except:
            segment_colors += ['grey']
            continue

        # Calculate dipole vector
        d = p1 + p2 

        # Angle says if we where plotting anti or towards dipole
        angle = np.real(np.dot(l, d) / np.linalg.norm(l) * np.linalg.norm(d)) * curve['f'][i]

        if angle > 0:
            color = pos_color
        else:
            color = neg_color
        
        segment_colors += [color]

    # Width norm
    w_Max = np.average(curve['m'][0])
    w_Min = np.average(curve['m'][0])
    widths = []

    # Get min and max mag for all lines (Cant np bc inhomogeneus)
    for i, line_mag in enumerate(curve['m']):
        w = np.average(line_mag)
        if w > w_Max:
            w_Max = w
        elif w < w_Min:
            w_Min = w

    # Get relative magnitude for each line
    for i, line_mag in enumerate(curve['m']):
        if w_Max == w_Min:
            widths += [0]
        else:
            w = np.average(line_mag)
            if np.isnan(w):
                widths += [0]
            else:
                widths += [(7* (w - w_Min) / (w_Max - w_Min)) + 3]

    # Add data for this frame
    frame_data = [
        go.Scatter3d(
            x=line_x,
            y=line_y,
            z=line_z,
            mode="lines",
            line=dict(width=widths[i], color=segment_colors[i]),
            showlegend=True,
            legendgroup=curve['lg'][i],
            legendgrouptitle_text=curve['lg'][i],
            name = ('Fwd' if curve['f'][i] > 0 else 'Bwd')
        ) for i, (line_x, line_y, line_z) in enumerate(zip(curve['x'], curve['y'], curve['z']))
    ]
    
    if show_seeds:
        p1T = np.real(p1)
        p2T = np.real(p2)
        frame_data += [
            go.Scatter3d(x = [p1T[0] / 2], y = [p1T[1] / 2], z = [p1T[2] / 2], line=dict(color=pos_color), legendgroup='1', legendgrouptitle_text='Dipole 1', name='pos'),
            go.Scatter3d(x = [p2T[0] / 2], y = [p2T[1] / 2], z = [p2T[2] / 2], line=dict(color=pos_color), legendgroup='2', legendgrouptitle_text='Dipole 2', name='pos'),
            go.Scatter3d(x = [-p1T[0] / 2], y = [-p1T[1] / 2], z = [-p1T[2] / 2], line=dict(color=neg_color), legendgroup='1', legendgrouptitle_text='Dipole 1', name='neg'),
            go.Scatter3d(x = [-p2T[0] / 2], y = [-p2T[1] / 2], z = [-p2T[2] / 2], line=dict(color=neg_color), legendgroup='2', legendgrouptitle_text='Dipole 2', name='neg')
        ]
    frames.append(go.Frame(data=frame_data, name=f"frame{t}"))

    p1_past = p1
    p2_past = p2

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
            #xaxis_range=plot_x_range,
            #yaxis_range=plot_y_range,
            #zaxis_range=plot_z_range
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

fig.show()

print('Runtime: %s' %(dt.now() - start))