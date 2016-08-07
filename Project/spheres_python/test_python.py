__author__ = 'Idan'
import spheres as sp


Ncolissions = 30000
print "colissions", Ncolissions
for rv in [1.2, 1.8, 3.6]:
    print "rv:", rv
    for Ns in [4, 32, 108, 256, 500]:
        n = int((Ns / 4.) ** (1. / 3))
        print "Ns:", Ns
        kin, mom, delta_vel, tx, y1, kinetic = sp.spheres(rv, Ncolissions=Ncolissions, n=n)
