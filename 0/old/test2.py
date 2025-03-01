import numpy as np
import plotly.graph_objects as go

# Constants
omega = 2 # Angular frequency
t = 0 # Start time
c = 3e8 #light
eps = 8.85e-12
k = omega / c #wave #

# Time variables for animation
num_frames = 30
time_steps = np.linspace(0, 2 * np.pi / omega, num_frames)


def EField():
    # Initialize field component arrays
    Ex = np.zeros_like(X)
    Ey = np.zeros_like(Y)
    Ez = np.zeros_like(Z)

    # Calculate electric field at each point
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                r_vec = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]])
                r_mag = np.linalg.norm(r_vec)
                
                # Avoid division by zero at the origin
                #if r_mag != 0:
                n_vec = r_vec / r_mag
                E_vec = (3 * np.dot(p, n_vec) * n_vec - p) / r_mag**3
                Ex[i, j, k] = E_vec[0]
                Ey[i, j, k] = E_vec[1]
                Ez[i, j, k] = E_vec[2]
    return Ex, Ey, Ez

# Funcs to generate seeds
def planeSeeds(seeds, seed_distance):
    seeds = int(seeds**.5)
    width = seed_distance * seeds

    x0 = []
    y0 = []
    z0 = []

    for x in range(seeds):
        for y in range(seeds):
            if x * seed_distance != 0 or y * seed_distance != 0:
                x0 += [(x * seed_distance) - (width / 2)]
                y0 += [0]
                z0 += [(y * seed_distance) - (width / 2)]

    return x0, y0, z0

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

#Model setup
spacing = .05 # Spacing between points
p0 = 10
seeds = 50
x0, y0, z0 = planeSeeds(seeds, .5)
seeds = len(x0)

#Frame list for animation
frames = []
for t in time_steps:
    # Initialize lists for field line points
    lines_segments = {'x':[],'y':[],'z':[],'c':[],'w':[]}

    # Calculate the dipole moment at time t
    p = np.array([0, p0 * np.cos(omega * t), 0])
    segments_x = []
    segments_y = []
    segments_z = []
    segment_colors = []
    segment_widths = []
    for i in range(seeds):
        # Starting point of each field line
        x, y, z = x0[i], y0[i], z0[i]
        x00, y00, z00= x, y, z

        for n in range(50):  # Adjust for longer/shorter field lines
            # Interpolate field at this point
            r_vec = np.array([x, y, z])
            r_mag = np.linalg.norm(r_vec)
            
            if r_mag != 0:
                n_vec = r_vec / r_mag
                # Radiation component
                E_1 = np.cross(n_vec, p)
                E_1 = np.cross(E_1, n_vec)
                E_1 *= k**2 * np.cos(k * r_mag) / r_mag
                
                # Static dipole field component
                E_2 = (3 * np.dot(n_vec, p) * n_vec - p) / r_mag**3
                
                # Total field
                E_vec = E_2 + E_1
                
                # Compute magnitude and direction for coloring and thickness
                E_mag = np.linalg.norm(E_vec)
                E_vec /= E_mag
                dir_fac = np.sign(E_vec[0]) + np.sign(E_vec[1]) + np.sign(E_vec[2])  # Use x-component direction for color


                if y > 0 and z > 0 and np.sign(E_vec[2]) > 0:
                    color = "blue"
                elif y < 0 and z > 0 and np.sign(E_vec[2]) < 0:
                    color = "blue"
                elif y > 0 and z < 0 and np.sign(E_vec[2]) < 0:
                    color = "blue"
                elif y < 0 and z < 0 and np.sign(E_vec[2]) > 0:
                    color = "blue"
                else:
                    color = "red"
                #if z > 0 and np.sign(E_vec[2]) < 0:
                #    if dir_fac > 1:
                #        color = "blue"
                #    elif dir_fac < 1:
                #        color = "red"
                #    else:
                #        color = "grey"

                width = E_mag # Scale line width with magnitude
                
                # Next point in the field direction
                next_x = x + E_vec[0] * spacing
                next_y = y + E_vec[1] * spacing
                next_z = z + E_vec[2] * spacing
                
                # Append this segment
                segments_x.append([next_x])
                segments_y.append([next_y])
                segments_z.append([next_z])
                segment_colors.append(color)
                segment_widths.append(width)
                
                # Update current point
                x, y, z = next_x, next_y, next_z

                if np.abs(x) < x00 + spacing and np.abs(y) < y00 + spacing and np.abs(z) < z00 + spacing:
                    break

    # Add all segments for this time step as a frame
    l = 0
    for i in segment_colors:
        if i == "blue":
            l += 1
        elif i == "red":
            l -= 1

    if l > 0:
        color = "blue"
    if l < 0:
        color = "red"
    else:
        color = "grey"

    frame_data = [
        go.Scatter3d(
            x= segments_x,
            y= segments_y,
            z= segments_z,
            mode="lines",
            line=dict(width=np.average(segment_widths)+5, color=color),
            showlegend=False
        )
    ]
    frames.append(go.Frame(data=frame_data, name=f"frame{t}"))
print('s')
# Create figure and add the initial frame for display


print(len(frames))
fig = go.Figure(
    layout=go.Layout(
        title="Animated Electric Field Lines of an Oscillating Dipole",
        scene=dict(
            xaxis_title="x (m)",
            yaxis_title="y (m)",
            zaxis_title="z (m)",
            aspectratio=dict(x=1, y=1, z=1),
        ),
        updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, 30],
                        "label": "&#9654;", # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], 0],
                        "label": "&#9724;", # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
         ],
    ),
    frames=frames

)
fig.add_trace(frames[0].traces())
# Show the animation
fig.show()

if False:
    # Add each segment to the plot
    for i in range(len(lines_segments)):
        fig.add_trace(
            go.Scatter3d(
                x=lines_segments[i]['x'],
                y=lines_segments[i]['y'],
                z=lines_segments[i]['z'],
                mode="lines",
                line=dict(width=lines_segments[i]['w'], color=lines_segments[i]['c']),
                showlegend=False
            )
        )

    #fig.add_trace(go.Scatter3d(x = x0, y = y0, z = z0))

    # Set layout for the 3D plot
    fig.update_layout(
        title="Electric Field Lines of an Oscillating Dipole",
        scene=dict(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title="z",
            aspectratio=dict(x=1, y=1, z=1)
        )
    )

    # Show the plot
    fig.write_html("file.html")
    fig.show()