import numpy as np
import plotly.express as px
import math
import pandas as pd
from datetime import datetime as dt

#
#
#  +                -
#  |----------------|
#        100 units

start = dt.now()

d = 100
seeds = 50
seed_distance = d/100

def fibonacci_sphere(samples):

    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x* seed_distance - d/2, y * seed_distance, z * seed_distance))

    return points

seedPoints = pd.DataFrame(fibonacci_sphere(seeds), columns=['x','y','z'])

lines = {'i':[],'x':[],'y':[],'z':[]}
i = 0
for row in seedPoints.iterrows():
    row = row[1]

    done = False
    x = row['x']
    y = row['y']
    z = row['z']

    r = np.array([x,y,z])
    r1 = np.array([-d/2,0,0])
    r2 = np.array([d/2,0,0])

    while not done:
        x = r[0]
        y = r[1]
        z = r[2]

        E = ((r - r1) / (np.linalg.norm(r - r1))**3) - ((r - r2) / (np.linalg.norm(r - r2))**3)

        r += (E/np.linalg.norm(E))

        lines['i'] += [i]
        lines['x'] += [x]
        lines['y'] += [y]
        lines['z'] += [z]

        if np.abs(x-(d/2)) < 1 and np.abs(y) < 1 and np.abs(z) < 1:
            done = True

    i += 1

print('Processing Time: %s' %(dt.now() - start))


fig = px.line_3d(lines, x = 'x', y = 'y', z = 'z', line_group='i') # animation_frame = 



fig.update_layout(scene = dict(xaxis = dict(range=[-200,200],),yaxis = dict(range=[-200,200]),zaxis = dict(range=[-200,200])))
fig.update_layout(scene_aspectmode='cube')

fig.write_html("file.html")

fig.show()
