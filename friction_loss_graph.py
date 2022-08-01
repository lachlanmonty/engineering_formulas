import numpy as np
from scipy.optimize import fsolve
import plotly.express as px


# def colebrook(f):
#     fun = 1 / np.sqrt(f) + 2 * np.log((e / (3.7 * d)) + (2.51 / (re * np.sqrt(f))))
#     return fun


# res = fsolve(colebrook, 0.001)


# calculate and plot frisction loss as a function of pipe size @ 3 m/s
def swamee_jain(e, D, re):
    f = 0.25 / (np.log((e / (3.7 * D)) + (5.74 / (re ** 0.9)))) ** 2
    return f


def reynolds(rho, mu, D, V):
    re = (rho * V * D) / mu
    return re


def friction_loss(f, L, D, V, g=9.81):
    h_l = f * (L / D) * ((V ** 2) / (2 * g))
    return h_l


dn = np.linspace(0.02, 0.2, 1000)

rho = 1000  # kg/m3
mu = 0.001  # Pa.s
V = 3  # m/s
e = 0.001  # mm?
l = 1000  # 100 m

re = reynolds(rho=rho, mu=mu, D=dn, V=V)

f = swamee_jain(e=e, D=dn, re=re)

h = friction_loss(f=f, L=l, D=dn, V=V)


# plot

fig = px.scatter(x=dn * 1000, y=h)
fig.show()

