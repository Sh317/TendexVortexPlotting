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

def fibonacci_sphere(samples,a):

    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x* seed_distance + a, y * seed_distance, z * seed_distance))

    return points

#Point seperation
a = 100

#Seed params
seeds = 10
seed_distance = 1
seedPoints = fibonacci_sphere(seeds,a)
seed_point_dict = pd.DataFrame(seedPoints, columns=['x','y','z'])

#Spin
s = .5

#Variable setup
xA = np.array([a,0,0])
xB = np.array([-a,0,0])
SA = [0,0,s]
SB = [0,0,-s]
lines = {'i':[],'x':[],'y':[],'z':[]}

#Dimenisons
jTot = 3
kTot = 3

i = 0
for row in seedPoints:
    x = np.array(row)
    
    prev_vec = []
    first = True
    n = 0
    while n < 1000:

        rA = np.linalg.norm(x - xA)
        rB = np.linalg.norm(x - xB)
        
        nA = (x - xA) / rA
        nB = (x - xB) / rB

        SA_dot_nA = np.dot(SA,nA)
        SB_dot_nB = np.dot(SB,nB)
        BMatrix = []

        for j in range(jTot):
            temp_row = []
            for k in range(kTot):
                
                delt = 0
                if j == k:
                    delt = 1

                temp_rowA = (3 / (rA**4)) * (2 * SA[j] + SA_dot_nA * (delt - 5 * nA[j] * nA[k]))
                temp_rowB = (3 / (rB**4)) * (2 * SB[j] + SB_dot_nB * (delt - 5 * nB[j] * nB[k]))
                temp_row += [temp_rowA + temp_rowB]

            BMatrix += [temp_row]

        B_Eig = np.linalg.eig(BMatrix)
        B_Eig_Val = B_Eig[0]
        B_Eig_Vec = B_Eig[1]
        
        min_dot_ind = 0
        vecInd = 0
        if first:
            vecInd = np.argmax(B_Eig_Val)
            first = False

            prev_vec = B_Eig_Vec[vecInd]
            val = B_Eig_Val[vecInd]
        else: # TODO max with positive vorticity
            dot0 = np.dot(B_Eig_Vec[0], prev_vec)
            dot1 = np.dot(B_Eig_Vec[1], prev_vec)
            dot2 = np.dot(B_Eig_Vec[2], prev_vec)
            min_dot_ind = np.argmin([dot0, dot1, dot2])

            
            val = B_Eig_Val[min_dot_ind]
        
        

        if val < 0:
            print(i)
            prev_vec = B_Eig_Vec[min_dot_ind]
            prev_vec = np.real(prev_vec)

            lines['i'] += [i]
            lines['x'] += [x[0]]
            lines['y'] += [x[1]]
            lines['z'] += [x[2]]
            x += prev_vec# * np.real(val)

        
        #print('________________________')
        #print(B_Eig_Val[0], B_Eig_Vec[0])
        #print(B_Eig_Val[1], B_Eig_Vec[1])
        #print(B_Eig_Val[2], B_Eig_Vec[2])
        #print('________________________')
        n += 1
    i += 1

    #print('-----')
                

        



print('Processing Time: %s' %(dt.now() - start))

fig = px.line_3d(lines, x = 'x', y = 'y', z = 'z', line_group='i')

fig.update_layout(scene = dict(xaxis = dict(range=[-200,200],),yaxis = dict(range=[-200,200]),zaxis = dict(range=[-200,200])))
fig.update_layout(scene_aspectmode='cube')

fig.show()