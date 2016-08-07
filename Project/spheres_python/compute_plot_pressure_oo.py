__author__ = 'Idan'
import time as tm

import matplotlib.pylab as plt

import spheres as sp
from Project.spheres_python_f2py.backup.course_resources import compute_pressure as lagacy

kin_fix = 0
rv0 = 1.1
drv = 0.15
jumps = 10

Ncolissions = 30001

Ns = 500

n = int((Ns/4.) ** (1./3))
speedScale = 1.0


x = []
y = []
rvs = [rv0 + drv * i for i in range(jumps)]
for rv in rvs:
    start = tm.time()
    kin, mom, delta_vel, tx, y1, kinetic = sp.spheres(rv, Ncolissions=Ncolissions, n=n, speedScale=1.0).run()
    result = lagacy.calc_pressure()
    x.append(result[0])
    y.append(result[1])

plt.plot(x, y)
plt.xlabel("V/Vo")
plt.ylabel("Pressure")
plt.show()
