__author__ = 'Idan'
import spheres as sp

Ncolissions = 3001
Ns = 4
rv = 1.2
n = int((Ns/4.) ** (1./3))


kin, mom, delta_vel, tx, y1, kinetic = sp.spheres(rv, Ncolissions=Ncolissions, n=n)
