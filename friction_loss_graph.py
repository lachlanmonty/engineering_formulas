import nump as np
from scipy.optimize import fslove


def darcy_weisbach(f_d, q, d, g=9.81):
    """Darcy Weisbach equation
    https://en.wikipedia.org/wiki/Darcy-Weisbach_equation

    Args:
        f_d (float): _description_
        q (float): Volumetric Flow [m3/s]
        d (float): Hydraulic Diameter [m]
        g (float, optional): _description_. Defaults to 9.81.
    """
    return f_d * (8 / np.pi) * ()


def colebrook_white(f):

    return 1 / np.sqrt(f)

