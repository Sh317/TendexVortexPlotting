import matplotlib.pyplot as plt
import numpy as np
a_I = 1
I = np.array([[a_I, 0, 0], 
              [0, -a_I, 0],
              [0, 0, 0]])

def f(r_V, give = False): #Assuming off diagonal = 0, assuming I33 = 0
    r = np.linalg.norm(r_V)
    x_1 = r_V[0]
    x_2 = r_V[1]
    x_3 = r_V[2]
    IpqXpXq = ((a_I * (x_1**2 - x_2**2)))

    E_1 = I
    E_2 = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            E_2[i,j] = I[0,i] * x_1 * r_V[j] + I[1,i] * x_2 * r_V[j] + I[0,j] * x_1 * r_V[i] + I[1,j] * x_2 * r_V[i]

    #E_2 = np.array([[x_1**2 * a_I, 1/2 * a_I * x_1 * x_2, 1/2 * a_I * x_1 * x_3], 
    #                    [1/2 * a_I * x_1 * x_2, -1 * x_2**2 * a_I, -1/2 * a_I * x_2 * x_3], 
    #                    [1/2 * a_I * x_1 * x_3, -1/2 * a_I * x_2 * x_3, 0]])
    

    E_3 = np.identity(3) * IpqXpXq
    E_4 = np.zeros((3,3))
    for i in range(len(I)):
        for j in range(len(I)):
            E_4[j,i] = (r_V[i] * r_V[j])

    E_4 *= IpqXpXq


    E = (-6 * E_1) + (30 * E_2 / (r**2)) + (15 * E_3 / (r**2)) + (-105 * E_4 / (r**4))
    
    if give:
        return E
    return E

def e_solve(E, r_vec, r_past, icity, mag_return = False): # r_past dotted with e vects
    # Find e vec/vals
    e_vals, e_vecs = np.linalg.eig(E)

    #dots = np.zeros(3)
    #dots_check = np.zeros(3)
    

    #for i in range(3):
    #    e_vec = e_vecs[:, i]
    #    dots_check[i] = True

    #    dot = np.dot(r_past, e_vec)
    #    dot = int(dot * 10e9) / 10e9

    #    if e_vals[i] * icity > 0:
    #        dots[i] = np.abs(dot)
    #        dots_check[i] = False

    #dots = ma(dots, dots_check)
    
    #index = np.argmax(dots)
    e_valcheck = np.ma.greater(e_vals * -1 * icity, 0)

    index = np.argmin(e_valcheck)
    if e_vals[index] < 0:
        print('sfsfgasfga')

    r_change = e_vecs[:, index]

    if mag_return:
        #mag = e_vals[index]
        return r_change, e_vals, index, e_vecs
    else:
        return r_change

def get_perpendicular_points(x1, y1, dx, dy, n, distance):
    # Compute the perpendicular vector
    perp_dx, perp_dy = -dy, dx
    
    x = []
    y = []
    
    # Generate points on either side
    for i in range(1, n + 1):
        offset = i * distance
        x += [x1 + perp_dx * offset, x1 - perp_dx * offset]
        y += [y1 + perp_dy * offset, y1 - perp_dy * offset]

    z = np.zeros(len(x))
    return x, y, z

def generate_grid_from_function(f, x_start, x_end, n, distance):
    x_list = []
    y_list = []
    z_list = []

    x_values = np.linspace(x_start, x_end, n)
    
    for x in x_values:
        y = f(x)
        
        # Compute numerical derivative (dx, dy)
        epsilon = 1e-10
        dx = epsilon
        dy = f(x + epsilon) - y
        length = np.hypot(dx, dy)
        dx, dy = dx / length, dy / length
        
        x0, y0, z0 = get_perpendicular_points(x, y, dx, dy, n, distance)

        x_list.extend(x0)
        y_list.extend(y0)
        z_list.extend(z0)

    
    return x_list, y_list, z_list

def grid(xSeeds, fStart, fStop, fSeeds, seed_distance): #0

    x0 = []
    y0 = []
    z0 = []
    x = np.linspace(fStart, fStop, fSeeds)
    y = g1P(x)

    for i,j in zip(x,y):
        x0 += [i]
        y0 += [j]
        z0 += [0]
        x_new1 = np.linspace(i, i+xSeeds*seed_distance, xSeeds)
        x_new2 = np.linspace(i-xSeeds*seed_distance, i, xSeeds)
        y_new = np.linspace(j, j, xSeeds)
        for n in x_new1:
            x0 += [n]
        for n in x_new2:
            x0 += [n]
        for n in y_new:
            y0 += [n]
            y0 += [n]
            z0 += [0]
            z0 += [0]

    return x0, y0, z0

def g1P(x):
    return (2.2749) * x


if __name__ == '__main__':

    #fig, ax = plt.subplots()
    #ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(-30,30,1000)
    y = g1P(x)
    #z_list = np.zeros(1000)

    plt.plot(x, y, color = 'black', label = 'line where x and y components of eigenvectors equal')
    x_list,y_list,z = grid(5,1,3,10,.1)#generate_grid_from_function(g1P, 1, 2, 10, .1)

    U = []
    V = []
    for x, y in zip(x_list, y_list):
        r = np.array([x,y,0])
        E = f(r)
        e_vec = e_solve(E, r, 0, 1)
        print(e_vec, x, y, g1P(e_vec[0]))
        U += [e_vec[0]]
        V += [e_vec[1]]

    plt.quiver(x_list, y_list, U, V, color = 'blue')
    #plt.scatter(x_list, y_list)

    #ax.scatter(x,y,z)




    plt.legend()
    plt.show()



