import numpy as np
import plotly.graph_objects as go
from seeding import seed
from ctypes import *
from datetime import datetime as dt

#Runtime
start = dt.now()

class vect(Structure):
    _fields_ = [('x', c_double*10000), ('y', c_double*10000), ('z', c_double*10000), ('m', c_double*10000), ('its', c_int)]

title, x0, y0, z0 = seed(4) # 0: Plane 1: Circular 2: Spherical 3: Helical 4: Random


curve = {'x':[],'y':[],'z':[],'f':[],'m':[],'c':[],'lg':[]}
seeds = len(x0)
num_its = 1000
delta_0 = 10e-2
h0 = 10e-2
safety = .9
ending_tolerance = 1.0
time_steps = np.linspace(0, 2 * np.pi, 1)
pos_color = 'red'
neg_color = 'blue'

icity = 1

# load C++
rka_iter = CDLL("rka_iter").rka_iter
rka_iter.argtypes = [c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double]
rka_iter.restype = vect

for i in range(seeds):
    # Starting point of each field line
    x, y, z = x0[i], y0[i], z0[i]
    vect_c = rka_iter(x, y, z, num_its, icity, ending_tolerance, delta_0, safety, h0)
    its = vect_c.its
    curve['x'].append(vect_c.x[0:its])
    curve['y'].append(vect_c.y[0:its])
    curve['z'].append(vect_c.z[0:its])
    curve['m'].append(vect_c.m[0:its])

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