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


def oroskar_settling_velocity(D, d, S, C_c, p_f, u_f, X = 1, g = 9.81):
    
    r'''Calculates the settling or critical velocity for a slurry with a 
    significant portion of fines and coarse particles. [1]
       
    Parameters
    ----------
    D : float
        Pipe inner diameter [m]
    d : float
        Particle diameter [m]
    S : float
        Ratio of the coarse solid density to carrier fluid density, ps/pf [-]
    C_c : float
        Coarse particle volume fraction (i.e. particles exceeding 74 um) [-]
    p_f : float
        Carrier fluid, including fines density [kg/m]
    u_f : float
        Carrier fluid dynamic viscosity [Pa.s]
    X : float
        Hindered settling factor (1) [-]
    g : float
        Gravitational acceleration (9.81) [m/s2]

    Returns
    -------
    V_ot : float
        Oroskar and Turian (1980) critical velocity [m/s]
    
    Examples
    --------
    >>> oroskar_settling_velocity(TBC)
    TBC
    
    References
    ----------
    .. [1] AP Poloski et al., Deposition Velocities of Newtonian and
    Non-Newtonian Slurries in Pipelines. 2009.
    https://www.pnnl.gov/rpp-wtp/documents/WTP-RPT-175.pdf

    '''
        
    V_ot = ((g*d*(S-1))**0.5)*\
    (\
        (1.85*C_c**0.1536)*\
        ((1-C_c)**0.3567)*\
        ((D/d)**0.378)*\
        (((p_f*D*((g*d*(S-1))**0.5))/u_f)**0.09)*\
        (X**0.3)
    )
    return V_ot


def pump_power(q_m3h, rho, head, efficiency, g = 9.81):

    P_kw = round((((q_m3h/60/60)* rho * g * head)/ efficiency)/1000,2)

    return P_kw
