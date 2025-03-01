import plotly.express as px
import pandas as pd
import math

seeds = 50

def fibonacci_sphere(samples):

    points = [(0,0,0)]
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

def planeSeeds(seeds, seed_distance):
    seeds = int(seeds**.5)
    xT = seeds * seed_distance
    yT = seeds * seed_distance

    

    l = []

    for x in range(xT):
        for y in range(yT):
            if x * seed_distance != 0 or y * seed_distance != 0:
                l += [[(x * seed_distance) - (xT / 2), (y * seed_distance) - (yT / 2), 0]]

    return l

seedPoints = pd.DataFrame(planeSeeds(seeds, 1), columns=['x','y','z'])

print(len(seedPoints))
fig = px.scatter_3d(seedPoints, x = 'x', y = 'y', z = 'z')
fig.show()
