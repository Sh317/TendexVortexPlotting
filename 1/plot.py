import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from datetime import datetime as dt


def plot(data):
    """
    Create a 3D line chart with varying widths and colors at each point.

    Parameters:
        data (dict): A dictionary with keys 'x', 'y', 'z', 'm', and 'c'.
                     - 'x', 'y', 'z' are lists of coordinates.
                     - 'm' is a list of line widths at each point.
                     - 'c' is a list of colors at each point (as RGB or hex values).
    """
    interp_amount = 2
    # Extract data
    x_list = data['x']
    y_list = data['y']
    z_list = data['z']
    m_list = data['m']
    c_list = data['c']

    # Create the figure and axis
    fig = plt.figure()
    ax = fig.add_subplot(111)#, projection='3d')

    for x, y, z, m, c in zip(x_list, y_list, z_list, m_list, c_list):
        # Draw segments with varying widths and colors
        for i in range(len(x) - 1):
            # Interpolate points for smooth rendering
            xs = np.linspace(x[i], x[i+1], interp_amount)
            ys = np.linspace(y[i], y[i+1], interp_amount)
            #zs = np.linspace(z[i], z[i+1], interp_amount)
            
            # Interpolate width and color
            interp_widths = np.linspace(m[i], m[i+1], interp_amount)
            interp_colors = np.linspace(
                np.array(plt.cm.colors.to_rgba(c[i])),
                np.array(plt.cm.colors.to_rgba(c[i+1])),
                interp_amount
            )

            for j in range(len(xs) - 1):
                ax.plot(xs[j:j+2], ys[j:j+2],
                        linewidth=interp_widths[j],
                        color=interp_colors[j])

    # Show the plot
    plt.show()