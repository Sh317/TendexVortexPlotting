import numpy as np
import plotly.graph_objects as go
from datetime import datetime as dt

#Runtime
start = dt.now()

#Constants
omega = 2          
c = 3e8            
epsilon_0 = 8.85e-12
k = omega / c

# Funcs to generate seeds
def planeSeeds(xSeeds, ySeeds, zSeeds, seed_distance):
    totX_dist = (xSeeds - 1) / 2
    totY_dist = (ySeeds - 1) / 2
    totZ_dist = (zSeeds - 1) / 2

    x0 = []
    y0 = []
    z0 = []

    for x in range(xSeeds):
        for y in range(ySeeds):
            for z in range(zSeeds):
                x0 += [(x - totX_dist) * seed_distance]
                y0 += [(y - totY_dist) * seed_distance]
                z0 += [(z - totZ_dist) * seed_distance]

    return x0, y0, z0

def f(r_V, p1, k):
    r = np.linalg.norm(r_V)
    n = r_V / r

    E_vec1 = k**2 * np.cross(np.cross(n, p1), n) * np.exp(1j * k * r) / r
    E_vec1 += (3 * n * np.dot(n, p1) - p1) * (1/ r**3) * np.exp(1j * k * r)

    E_vec2 = k**2 * np.cross(np.cross(n, p2), n) * np.exp(1j * k * r) / r
    E_vec2 += (3 * n * np.dot(n, p2) - p2) * (1/ r**3) * np.exp(1j * k * r)

    E = (1 / (4 * np.pi * epsilon_0)) * (E_vec1 + E_vec2)
    return E


# Model setup
showSeeds = True
xSeeds = 10
zSeeds = 3
ySeeds = 2
seed_distance = 3 #Spacing between seeds
h = 10**-10 #Scaling for E field movement
num_frames = 1
num_its = 100
ending_tolerance = .5 #How far from 0 to break

#Dipole setup
q = 1
p0 = 1
p1x = np.sqrt(2)
p1y = np.sqrt(2)
p1z = 0
p2x = -np.sqrt(2)
p2y = np.sqrt(2)
p2z = 0

#Plot setup
plot_x_range = [-seed_distance * xSeeds,seed_distance * xSeeds]
plot_y_range = [-seed_distance * xSeeds,seed_distance * xSeeds]
plot_z_range = [-seed_distance * zSeeds,seed_distance * zSeeds]
title = "Alternating dipole rectangular seeding (%sx%s), 1 seed spacing. RK4" % (xSeeds, zSeeds)

#Seeding
x0, y0, z0 = planeSeeds(xSeeds, ySeeds, zSeeds, seed_distance)
seeds = len(x0)

#Time setup
time_steps = np.linspace(0, 10, num_frames)

frames = []
for t in time_steps:
    p1 = np.array([p1x, p1y, p1z]) * np.exp(1j * t)
    p2 = np.array([p2x, p2y, p2z]) * np.exp(1j * t)
    
    curve_x, curve_y, curve_z = [], [], []
    segment_colors, segment_widths = [], []

    for i in range(seeds):
        for factor in [1,-1]: #Moving forwards orbackwards
            # Starting point of each field line
            x, y, z = x0[i], y0[i], z0[i]

            #Set past variables
            E_vec = 0
            sign = 'a'

            # Accumulate segments for each line
            line_x, line_y, line_z = [x], [y], [z]

            for n in range(num_its // 2):
                r_vec = np.array([x, y, z])
                r_mag = np.linalg.norm(r_vec)
                
                if r_mag > 0:
                    #RK4
                    k1 = h * f(r_vec, p1, k)
                    k2 = h * f(r_vec + h/2, p1, k)
                    k3 = h * f(r_vec + h/2, p1, k)
                    k4 = h * f(r_vec + h, p1, k)

                    E_vec = E_vec + (k1 / 6) + (k2 / 3) + (k3 / 3) + (k4 / 6)
                    
                    #Wrap up
                    E_vec = np.real(E_vec)
                    E_mag = np.linalg.norm(E_vec)
                    
                    # Compute next point in field direction
                    x += factor * E_vec[0]
                    y += factor * E_vec[1]
                    z += factor * E_vec[2]
                    
                    # Append segment to line
                    line_x.extend([x])
                    line_y.extend([y])
                    line_z.extend([z])

                    # Color and thickness based on direction and magnitude
                    segment_widths.append(E_mag)


                    if sign != np.sign(E_mag) and sign != 'a':
                        print(E_mag)

                    sign = np.sign(E_mag)

                    if np.abs(x) < ending_tolerance and np.abs(y) < ending_tolerance and np.abs(z) < ending_tolerance:
                        break
        
            # Add segments for the current line
            curve_x.append(line_x)
            curve_y.append(line_y)
            curve_z.append(line_z)
    
    # Color everything
    for i, (line_x, line_y, line_z) in enumerate(zip(curve_x, curve_y, curve_z)):
        lY = len(line_y)
        if line_y[lY//4] > line_y[3 * lY//4]:
            if np.sign(p1[1]) > 0:
                color = "blue"
            else:
                color = "red"
        else:
            if np.sign(p1[1]) > 0:
                color = "red"
            else:
                color = "blue"
        segment_colors += [color]

    # Width norm
    w_Max = max(segment_widths)
    w_Min = min(segment_widths)

    width = []
    for w in segment_widths:
        width += [(7* (w - w_Min) / (w_Max - w_Min)) + 3]

    # Add data for this frame
    frame_data = [
        go.Scatter3d(
            x=line_x,
            y=line_y,
            z=line_z,
            mode="lines",
            line=dict(width=width[i], color=segment_colors[i]),
            showlegend=False
        ) for i, (line_x, line_y, line_z) in enumerate(zip(curve_x, curve_y, curve_z))
    ]
    frames.append(go.Frame(data=frame_data, name=f"frame{t}"))

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
        width=700,
        height=700,
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
                                frame=dict(duration=100, redraw=True),
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
            "transition": {"duration": 100},
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

if showSeeds:
    fig.add_trace(go.Scatter3d(x = x0, y = y0, z = z0))
    fig.add_trace(go.Scatter3d(x = [-p1x / 2,p1x / 2], y = [-p1y / 2,p1y / 2], z = [-p1z / 2,p1z / 2]))
    fig.add_trace(go.Scatter3d(x = [-p2x / 2,p2x / 2], y = [-p2y / 2,p2y / 2], z = [-p2z / 2,p2z / 2]))

#fig.write_html("%s.html" %title)
fig.show()

print('Runtime: %s' %(dt.now() - start))