from ctypes import *
import numpy as np
from datetime import datetime as dt

a_I = 1

def rka(r_vec, icity = 1):
    delta_0 = .0004
    safety = .95
    h = .0001
    a = {'2':1/5,         '3':3/10,     '4':3/5,         '5':1,            '6': 7/8}
    b = {'21':1/5, 
        '31':3/40,       '32':9/40, 
        '41':3/10,       '42':-9/10,   '43':6/5, 
        '51':-11/54,     '52':5/2,     '53':-70/27,    '54': 35/27, 
        '61':1631/55296, '62':175/512, '63':575/13824, '64':44275/110592, '65': 253/4096}
    c_dif = {'1':37/378 - 2825/27648, '3':250/621 - 18575/48384, '4':125/621 - 13525/55296, '5':-277/14336, '6':512/1771 - 1/4}
    c = {'1':37/378, '3':250/621, '4':125/621, '6':512/1771}
    step_sizing = True
    while step_sizing:
        #k1
        E = f(r_vec)
        k1 = h * e_solve(E, icity)

        #k2
        E_temp = f(r_vec + k1 * b['21'])
        k2 = h * e_solve(E, icity)

        #k3
        E_temp = f(r_vec + k1 * b['31'] + k2 * b['32'])
        k3 = h * e_solve(E, icity)

        #k4
        E_temp = f(r_vec + k1 * b['41'] + k2 * b['42'] + k3 * b['43'])
        k4 = h * e_solve(E, icity)

        #k5
        E_temp = f(r_vec + k1 * b['51'] + k2 * b['52'] + k3 * b['53'] + k4 * b['54'])
        k5 = h * e_solve(E, icity)

        #k6
        E_temp = f(r_vec + k1 * b['61'] + k2 * b['62'] + k3 * b['63'] + k4 * b['64'] + k5 * b['65'])
        k6 = h * e_solve(E, icity)

        delta = c_dif['1'] * k1 + c_dif['3'] * k3 + c_dif['4'] * k4 + c_dif['5'] * k5 + c_dif['6'] * k6 #c2 =0
        
        if (np.abs(delta) > np.abs(delta_0)).any():
            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .25
        else:
            h *= safety * np.abs(np.linalg.norm(delta_0) / np.linalg.norm(delta)) ** .2
            r_change = c['1'] * k1 + c['3'] * k3 + c['4'] * k4 + c['6'] * k6
            step_sizing = False

    return r_change

def f(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
    I = np.array([[a_I, 0, 0], 
              [0, -a_I, 0],
              [0, 0, 0]])
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1**2 - x_2**2)))

    E_1 = I

    E_2 = np.array([[2 * x_1**2 * a_I, 0, 0], 
                    [0, -2 * a_I * x_1 * x_2, 0], 
                    [0, 0, 0]])
    
    
    E_3 = np.identity(3) * IpqXpXq
    E_4 = np.zeros((3,3))
    for i in range(len(I)):
        for j in range(len(I)):
            E_4[j,i] = (r_V[i] * r_V[j])

    E_4 *= IpqXpXq

    E = (-6 * E_1) + (30 * E_2 / (r**2)) + (15 * E_3 / (r**2)) + (-105 * E_4 / (r**4))
    
    if give:
        return E, E_1, E_2, E_3, E_4
    return E

def e_solve(E, icity, val_return = False, s = None, r_vec = None): # r_past dotted with e vects
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    maxe_val = -10e10
    maxindex = -1

    for i in range(3):
        val = e_vals[i]

        if icity * val > 0 and icity * val > maxe_val:
            maxindex = i
            maxe_val = val

    r_change = e_vecs[maxindex]

    if maxindex == -1:
        print(e_vals)

    if val_return: # Needed for saving eigenvalue
        mag = e_vals[index]
        #with open("output.txt", "a") as file:
        #    file.write(str(e_vals) + ' ' + str(index) + ' ' + str(r_vec) + ' ' + str(e_valcheck) + ' ' + str(e_vals[index]) + ' ' + str(s) + '\n')
        return r_change, mag
    else:
        return r_change

class vect(Structure):
    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double)]

rka_c = CDLL("main").rka  # Use "example.dll" on Windows
rka_c.argtypes = [c_float, c_float, c_float, c_int]
rka_c.restype = vect


l = 10000
rand_x = np.random.rand(l) * 100
rand_y = np.random.rand(l) * 100
rand_z = np.zeros(l)

start = dt.now()

for x,y,z in zip (rand_x, rand_y, rand_z):
    vect_c = rka_c(x, y, z, 1)
    vect = rka(np.array([x, y, z]), 1)  
    print(vect- np.array([vect_c.x, vect_c.y, vect_c.z]))


c_time = dt.now() - start
start = dt.now()

for x,y,z in zip (rand_x, rand_y, rand_z):
    vect = rka(np.array([x, y, z]), 1)

py_time = dt.now() - start

print(c_time, py_time)