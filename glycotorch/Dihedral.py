import numpy as np


def dihedral(p0, p1, p2, p3, negative=False):
    """


    :param p0:
    :param p1:
    :param p2:
    :param p3:
    :param negative:
    :return:
    """
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

    y = np.dot(b0xb1_x_b1xb2, b1) * (1.0 / np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)

    ans = np.degrees(np.arctan2(y, x))

    if negative == True:
        return ans

    if ans < 0:
        return ans + 360
    return ans
