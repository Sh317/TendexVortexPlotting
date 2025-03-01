import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
import pandas as pd
from datetime import datetime as dt

start = dt.now()

def planeSeeds(seeds, seed_distance):
    seeds = int(seeds**.5)
    width = seed_distance * seeds

    l = []

    for x in range(seeds):
        for y in range(seeds):
            if x * seed_distance != 0 or y * seed_distance != 0:
                l += [[(x * seed_distance) - (width / 2), 0, (y * seed_distance) - (width / 2)]]

    return l


def fibonacci_sphere(samples, seed_distance,d):

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



#Config
d = 100
seeds = 10
seed_distance = 2



seedPoints = pd.DataFrame(planeSeeds(seeds, seed_distance), columns=['x','y','z'])
lines = {'lineNum':[],'x':[],'y':[],'z':[], 'c':[]}

i = 0

for row in seedPoints.iterrows():
    row = row[1]

    done = False
    reached = False
    
    x = row['x']
    y = row['y']
    z = row['z']

    r = np.array([x,y,z])
    
    
    eps = 8.8541878128*10**-12

    j = 0
    num = 0
    rho = 1.602*10**-19
    rho = 10**-7
    while not done:
        #E = ((r - r1) / (np.linalg.norm(r - r1))**3) - ((r - r2) / (np.linalg.norm(r - r2))**3)
        rMag = np.linalg.norm(r)
        
        p = np.array([0,(np.pi / 4) * rho * 1**4,0])
        n = r / rMag
        E = (3 * n * np.dot(n,p) - p)/ (rMag**3) * (1/(4 * np.pi * eps))
        print(r, E/10)
        lines['x'][-1]
        

        lines['lineNum'] += [i]
        lines['x'] += [r[0]]
        lines['y'] += [r[1]]
        lines['z'] += [r[2]]
        lines['c'] += [np.linalg.norm(E)]

        r += (E/ 10)


        #if np.abs(x) < 1 and np.abs(y) < 1 and np.abs(z) < 1:
        #    reached = True
        #    num = j
        #if j > num * 2 and reached:
        #    done = True
        if j > 10:
            done = True

        j += 1

    i += 1


fig = go.Figure()


i = 0
while i < len(lines['x']):
    j = i

    lineDict = {'x':[],'y':[],'z':[]}
    colorOld = 'blue'
    first = True
    while j < len(lines['x']) and lines['lineNum'][j] == lines['lineNum'][i]:
        if lines['c'][j] > 0:
            color = 'blue'
        else:
            color = 'red'

        if not first and colorOld != color:
            #print('oh snap')
            pass

        lineDict['x'] += [lines['x'][j]]
        lineDict['y'] += [lines['y'][j]]
        lineDict['z'] += [lines['z'][j]]

        j += 1
        #print(j)
        colorOld = color
        first = False

    fig.add_trace(go.Scatter3d(x=lineDict['x'],y=lineDict['y'], z=lineDict['z'], mode='lines', line={'color': color}, showlegend=False))

    i = j+1


fig.add_trace(go.Scatter3d(x = seedPoints['x'], y = seedPoints['y'], z = seedPoints['z']))
fig.write_html("plot.html")
fig.show()
