
# Funcs to generate seeds
def plane_seeds(xSeeds, ySeeds, zSeeds, seed_distance): #0
    totX_dist = (xSeeds - 1) / 2
    totY_dist = (ySeeds - 1) / 2
    totZ_dist = (zSeeds - 1) / 2

    x0 = []
    y0 = []
    z0 = []

    for x in range(xSeeds):
        for y in range(ySeeds):
            for z in range(zSeeds):
                if (x - totX_dist) * seed_distance != 0 or (y - totY_dist) * seed_distance != 0 or (z - totZ_dist) * seed_distance != 0:
                    x0 += [(x - totX_dist) * seed_distance]
                    y0 += [(y - totY_dist) * seed_distance]
                    z0 += [(z - totZ_dist) * seed_distance]

    return x0, y0, z0

def circular_seeds(seeds, starting_radius, radius_steps, radius_increment): #1
    x_list = []
    y_list = []
    z_list = []
    phi = np.pi * (np.sqrt(5.) - 1.)  # golden angle in radians
    for n in range(radius_steps):
        radius = starting_radius + (radius_increment * n)
        for i in range(seeds):
            theta = phi * i  # golden angle increment

            x = np.cos(theta) * radius
            z = np.sin(theta) * radius
            if x != 0 or z != 0:
                x_list += [x]
                y_list += [0]
                z_list += [z]

    return x_list, y_list, z_list

def spherical_seeds(seeds, spacing): #2
    x_list = []
    y_list = []
    z_list = []
    phi = np.pi * (np.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(seeds):
        y = 1 - (i / float(seeds - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        if x != 0 or y != 0 or z != 0:
            x_list += [x * spacing]
            y_list += [y * spacing]
            z_list += [z * spacing]

    return x_list, y_list, z_list

def helical_seeds(start_rad, spacing, seeds): #3
    x_list = []
    y_list = []
    z_list = []
 
    for i in range(1,seeds+1):
        a = (0.78643418 * i) + 3.76697152
        x_list += [a * np.cos(i)]
        y_list += [a * np.sin(i)]
        z_list += [0]

    return x_list, y_list, z_list

def seed(seed_type):
    rotate = False
    if seed_type == 0: # Plane
            
        xSeeds = 7
        ySeeds = 5
        zSeeds = 1
        seed_distance = 30 #Spacing between seeds

        if rotate:
            title = "Seeding: (%s, %s, %s). Spacing: %s. RK Adaptive" % (xSeeds, ySeeds, zSeeds, seed_distance)
        else:
            title = "Seeding: (%s, %s, %s). Spacing: %s. RK Adaptive" % (xSeeds, ySeeds, zSeeds, seed_distance)

        x0, y0, z0 = plane_seeds(xSeeds, ySeeds, zSeeds, seed_distance)
    elif seed_type == 1: # Circular
        starting_radius = 2
        radius_steps = 2
        radius_increment = 5
        seeds = 5
        if rotate:
            title = "Circular seeding. Spacing: %s. RK Adaptive" % (radius_increment)
        else:
            title = "Circular seeding. Spacing: %s. RK Adaptive" % (radius_increment)
        
        x0, z0, y0 = circular_seeds(seeds, starting_radius, radius_steps, radius_increment)
    elif seed_type == 2: # Spherical
        radius = 10
        seeds = 15
        if rotate:
            title = "Spherical seeding. Radius: %s. RK Adaptive" % (radius)
        else:
            title = "Spherical seeding. Radius: %s. RK Adaptive" % (radius)

        x0, y0, z0 = spherical_seeds(seeds, radius)
    elif seed_type == 3: # Helical
        start_rad = 0
        spacing = 2
        seeds = 5

        if rotate:
            title = "Rotating dipole. Helical seeding. RK Adaptive"
        else:
            title = "Two dipole. Helical seeding. RK Adaptive"

        x0, y0, z0 = helical_seeds(start_rad, spacing, seeds)

    return title, x0, y0, z0