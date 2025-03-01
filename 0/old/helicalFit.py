import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Given Cartesian points (x, y, z)
points_cartesian = [
    (1.009509, -2.147538, 0),
    (-2.760612, 2.575502, 0),
    (2.764213, 3.647096, 0),
    (-1.982757, -3.996672, 0)
]

# Convert Cartesian points to cylindrical coordinates (r, theta, z)
r_vals = []
theta_vals = []
for x, y, z in points_cartesian:
    r_vals += [np.sqrt(x**2 + y**2)]
    theta_vals += [np.arctan(y/x)]

r_vals = np.array(r_vals)
theta_vals = np.array(theta_vals)

# Define linear helical equation
def linear_r(theta, k,a,b):
    return k/(theta)

# Fit r = k * theta for radial component
params_k, _ = curve_fit(linear_r, theta_vals, r_vals, p0=(1,1,1))

# Display the results
k = params_k[0]
a = params_k[1]
b = params_k[2]
print(params_k)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta_vals, r_vals, label ='real')
ax.plot(theta_vals, linear_r(theta_vals, k, a, b))
ax.legend()
plt.show()

print(f"Helical equation parameters:")
print(f"r = {k:.4f} * theta")
