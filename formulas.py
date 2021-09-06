import numpy as np


def bingham_darcy_friction_factor(d_m, v, rho, t_0, mu_b):
    
    r'''Calculates the Darcy friction factor for a Bingham plastic.
    Calculated for any Reynolds number, from laminar through turbulent. [1]
       
    Parameters
    ----------
    d_m : float
        Inner Diameter of Pipe, [m]
    v : float
        Velocity, [m/s]
    rho : float
        Density of Slurry, [kg/m^3]
    t_0 : float
        Yield stress, [Pa]
    mu_b : float
        Plastic Viscocity [Pa.s]
        
    Returns
    -------
    f_d : float
        Darcy Friction Factor, [-]
    
    Examples
    --------
    >>> bingham_darcy_friction_factor(0.254, 2.3, 1300, 6, 0.02)
    0.01905007708620241
    
    References
    ----------
    .. [1] Ron Darby, Chemical Engineering Fluid Mechanics. 2nd edition. 2001.

    '''

    N_re = d_m*v*rho/mu_b

    N_he = d_m**2*rho*t_0/mu_b**2

    # Laminar
    f_L = 16/N_re*(1+1/6*(N_he/N_re))
    
    for i in range(10):
        f_L = 16/N_re*(1+1/6*(N_he/N_re)-((N_he**4)/(3*f_L**3*N_re**7)))  
        
    f_L = 4 * f_L       # Converts to Darcy Friction Factor
    
    # Turbulent
    a = -1.47*(1+0.146*np.exp(-2.9*10**-5*N_he))
    f_T = 4 * (10**a/N_re**0.193)
    
    # Combine
    m = 1.7+40000/N_re
        
    f_d = (f_L**m+f_T**m)**(1/m)

    return f_d
