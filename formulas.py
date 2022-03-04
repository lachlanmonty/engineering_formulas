import math


def bingham_darcy_friction_factor(d_m, v, rho, t_0, mu_b):
    
    r"""Calculates the Darcy friction factor for a Bingham plastic.
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

    """

    N_re = d_m*v*rho/mu_b

    N_he = d_m**2*rho*t_0/mu_b**2

    # Laminar
    f_L = 16/N_re*(1+1/6*(N_he/N_re))
    
    for i in range(10):
        f_L = 16/N_re*(1+1/6*(N_he/N_re)-((N_he**4)/(3*f_L**3*N_re**7)))  
        
    f_L = 4 * f_L       # Converts to Darcy Friction Factor
    
    # Turbulent
    a = -1.47*(1+0.146*math.exp(-2.9*10**-5*N_he))
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



def api_standard_520_part1_sizing(d_m, v, rho, t_0, mu_b):
    
    r'''Calculates required size for rupture disc with liquid medium. [1]
       
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

    return


def cema_area(BW, beta, phi):
    r'''### CEMA Standard Capacity Cross Sectional Area, $A_s$

    Equation 4.15 is used to calculate, As, for standard CEMA three equal roll 
    roughing idlers based on the average CEMA center roll length circular 
    surcharge surface and the CEMA standard belt edge. [1]

    Parameters
    ----------
    - BW : float
        - Belt Width, [mm]
    - $\beta$ : float
        - Idler trough angle, [deg]
    - $\phi_s$ : float
        - Material surcharge angle, [deg]
        
    Returns
    -------
    - $A_s$ : float
        - CEMA Standard Cross Sectional Area, area based on three equal roll 
        CEMA troughing idler, the surcharge angle with circular top surface, 
        and standard edge distance, [m²]

    '''
    beta =  math.radians(beta)
    phi  =  math.radians(phi)

    bc = (0.371 * BW + 6.35)/BW         # Equation 4.11
    bw = (1-bc)/2                       # Equation 4.13
    bwe = (0.055 * BW + 22.9)/BW        # Equation 4.12
    bwmc = bw - bwe
    rsch = ((bc/2)/(math.sin(phi)))+((math.cos(beta)*bwmc)/(math.sin(phi)))

    A_s = 2*BW**2 * \
        (rsch**2*((phi/2)-((math.sin(phi)*math.cos(phi))/2)))+\
        ((bc/2) * bwmc * math.sin(beta)) + \
        (bwmc**2 * ((math.sin(beta) * math.cos(beta))/2))

    A_s = A_s / 1000**2                 # Convert to m²

    return A_s

    
def req_conv_area(bulk_density, tonnage, belt_velocity):
    r'''### Required Area for Mass Flow

    The conveyed bulk material cross sectional area, A, can be calculated from
    the design inputs for tonnage, Q, belt speed, V, bulk density, $\gamma_m$,
    and $\Theta $ = 0 degrees.  [1]

    Formula
    ----------
    $$ A = \frac{Q}{V \gamma_m} $$
    
    Parameters
    ----------
    - $\gamma_m$ : float
        - Conveyed bulk density, [kg/m³]
    - $Q$ : float
        - Bulk material tonage, [t/h]
    - $V$ : float
        - Belt speed, [m/s]
    
    Returns
    -------
    - $A$ : float
        - Required conveyed cross sectional area, [m²]
    '''
    return (tonnage * 1000/60/60) / (belt_velocity * bulk_density)