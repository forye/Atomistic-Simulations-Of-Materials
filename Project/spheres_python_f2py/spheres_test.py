__author__ = 'Idan'
import spheres as sp
import numpy as np

Iterations = 1
Ncolissions = 30001
Ns = 32
rv = 1.8
n = int((Ns/4.) ** (1./3))

n = int(np.round(125**(1/3.)))

kin, mom, delta_vel, tx, y1, kinetic = sp.spheres(rv, Ncolissions=Ncolissions, n=n, Iterations=Iterations)
