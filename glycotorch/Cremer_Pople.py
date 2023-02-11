import math
import numpy as np


# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# import numpy as np
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')


# list_of_ring_coords = o5, c1, c2, c3, c4, c5
def cremer_pople(list_of_ring_coords):
    """
    A function to calculate the Cremer-Pople ring puckering co-ordinants of a six-membered ring
    :param list_of_ring_coords: list of x,y,z co-ordinants for o5, c1, c2, c3, c4, c5
                                Cremer-Pople co-ordinants require these atoms to be in this order

    :return: a dictionary in the form {"phi": Phi, "theta": Theta, "Q": Q, "q_2": q_2, "q_3": q_3}
    """

    #  Make two empty arrays for the vectors R' and R''
    R1 = np.array([0.0, 0.0, 0.0])
    R2 = np.array([0.0, 0.0, 0.0])
    N = len(list_of_ring_coords)

    #  Make an empty list to store the atom positions
    Rjs = []
    for index in range(N):
        Rj = np.array(list_of_ring_coords[index])
        Rjs.append(Rj)

    #  Make a new vector to translate the centre of the atoms to the point (0, 0, 0)
    center = np.array([0.0, 0.0, 0.0])
    for x in Rjs:
        center += x
    center = center / 6

    #  Make a new list to place in the centered atom positions
    centered = []
    for index in range(N):
        Rj = np.array(list_of_ring_coords[index])
        centered.append(Rj - center)

    for index in range(N):
        R1 += centered[index] * np.sin(2 * math.pi * index / 6)
        R2 += centered[index] * np.cos(2 * math.pi * index / 6)

    R1xR2 = np.cross(R1, R2)
    unit_vector_n = R1xR2 * (1 / np.linalg.norm(R1xR2))

    z = []
    for index in range(N):
        dot = 0
        for i, j in zip(unit_vector_n, centered[index]):
            dot += i * j
        z.append(dot)

    sum_for_eq12 = 0
    sum_for_eq13 = 0
    sum_for_eq14 = 0
    Q = 0
    for index in range(N):
        #  equations 12, 14 refer to the original Cremer Pople paper
        sum_for_eq12 += z[index] * np.cos(2 * math.pi * 2 * index / 6)
        sum_for_eq13 += z[index] * np.sin(2 * math.pi * 2 * index / 6)
        sum_for_eq14 += np.cos(index * math.pi) * z[index]
        Q += z[index] ** 2

    sqr2_div_6 = math.sqrt(2 / 6)
    q_2_cos0_m = sqr2_div_6 * sum_for_eq12
    q_2_sin0_m = -1 * sqr2_div_6 * sum_for_eq13
    q_3 = 1/math.sqrt(6) * sum_for_eq14
    q_2 = math.atan(q_2_sin0_m/q_2_cos0_m)

    Q = math.sqrt(Q)
    theta_radians = math.acos(q_3/Q)
    theta = np.rad2deg(math.acos(q_3/Q))

    # if q_2 < 0:
    #     phi += math.pi
    #     q_2 = q_2 * -1
    #
    # phi = np.rad2deg(phi)

    phi_radians = math.atan(q_2_sin0_m / q_2_cos0_m)

    if q_2_cos0_m > 0:
        if q_2_sin0_m > 0:

            phi = np.rad2deg(math.atan(q_2_sin0_m / q_2_cos0_m))
        else:
            phi = 360 - abs(np.rad2deg(math.atan(q_2_sin0_m / q_2_cos0_m)))
    else:
        if q_2_sin0_m > 0:
            phi = 180 - abs(np.rad2deg(math.atan(q_2_sin0_m / q_2_cos0_m)))
        else:
            phi = 180 + abs(np.rad2deg(math.atan(q_2_sin0_m / q_2_cos0_m)))

    return {"phi": phi, "theta": theta, "Q": Q, "q_2": q_2, "q_3": q_3,
            "phi_radians": phi_radians, "theta_radians": theta_radians }


# cremer_pople([[0, 1.3839, 0.1976],
#               [1.997,0.7624, -0.2106],
#               [1.2356, -0.704, 0.2393],
#               [0.0110, -1.4564, -0.2550],
#               [-1.23, -0.7208, 0.2420],
#               [-1.2164, 0.7350, -0.2133]])
