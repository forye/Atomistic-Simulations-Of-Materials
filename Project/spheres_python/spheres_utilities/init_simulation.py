__author__ = 'Idan'

import math

import numpy as np

import testing_utils as test
import configuration as conf


def init_parameters(n=0, Ncollisions = 1000 , Rvol = 0, L=1.0):
    L = float(L)
    test.icheck(type(n) is not int, "ERROR n must be int!.../n the number of spheres in a dimention (4n^3)=N!!")
    test.icheck(type(Ncollisions) is not int, "ERROR Ncolissions must be int!.../n the number of spheres in a dimention (4n^3)=N!!")
    Ns = (n ** 3) * 4
    Rv = Rvol
    if not Rvol:
        Rv = conf.VOLUME_RATIO
    sigma = compute_diameter(Ns,Rv) * L
    n = int(round((Ns / 4.0) ** (1./3)))
    return Ns, Rv, Ncollisions, sigma, n, L


def validate_input(nspheres,rvolume):
    '''
    ! This subroutine makes sure that a number of spheres is suitable for
    ! packing them in a fcc arrangement inside a cube, and that a reduced volume
    ! is larger than 1.
    True is Ok
    '''
    if conf.debug_mode: print'ENTERING SUBROUTINE validate_input'
    n = int(round( (nspheres/4.0) ** (1./3.) ))
    if 4*n**3 is not nspheres:
        print'ERROR: Wrong nspheres!'
        return False
    else:
        print' Number of spheres: OK'
    if rvolume < 1.0:
        print'ERROR: Wrong rvolume!'
        return False
    else:
        print' Reduced volume: OK'
    return True



def compute_diameter(nspheres,rvolume):
    '''
    ! Given a number of spheres (nspheres) occupying a reduced volume (rvolume)
    ! in a cube of unit volume, this subroutine returns the value of the diameter
    ! of the spheres (sigma). It uses the formula of homework 02 (exercise 1).
    '''
    if conf.debug_mode: print'ENTERING SUBROUTINE compute_diameter'
    sigma = ( math.sqrt(2.) / (nspheres*float(rvolume)) ) ** (1./3.)
    if conf.debug_mode: print' Diameter of spheres:', sigma
    return sigma


def init_positions(L, n, Ns):
    pos = np.zeros([Ns, 3])
    a =float(L)/n
    i=0
    for z in xrange(n):
        for y in xrange(n):
            for x in xrange(n):
                ref =np.array([x,y,z])*a
                pos[i:i+4] = ref + fcc_mask(a)
                i+=4
    return pos

def fcc_mask(a):
    #create FCC
    pos = np.zeros([4,3])
    i = 0
    pos[i] = np.array([0,  0     ,0   ])
    i += 1
    pos[i] = np.array([a/2.,a/2. ,0   ])
    i += 1
    pos[i] = np.array([a/2.,0    ,a/2.])
    i += 1
    pos[i] = np.array([0   ,a/2. ,a/2.])
    return pos


def init_velocities(Ns,kT,speedScale):
    mass = conf.MASS
    speed = speedScale * math.sqrt(3.)
    u = np.random.random(Ns)
    v = np.random.random(Ns)
    theta = 2.*np.pi * u
    phi =  np.arccos(2. * v - 1)
    vel = speed *np.array([np.sin(phi) * np.cos(theta),np.sin(phi) * np.sin(theta),np.cos(phi)])
    if conf.ZERO_AVG_MOM:
        avg_mom = np.mean(vel)
        vel = vel-avg_mom
    if (conf.debug_mode):
        print' Average momentum before setting it to 0:', avg_mom
    # avg_kin = np.sqrt(np.sum(vel**2))/ (2.* Ns)# was a mistake!!
    avg_kin =mass* np.sum(vel**2)/ (2.* Ns)
    alpha = math.sqrt(kT * 3./(2. * avg_kin) )# scale to temperture and other constants
    vel = alpha * np.array(vel)
    vel = np.transpose(vel)
    return vel,avg_kin,avg_mom


def init_collision_time(pos,vel,pos_trans,sigma,L=1.0):
    """
    deprecated due to numpy implementation, kept for testing
    """
    Ns = np.shape(pos)[0]
    ctime = conf.INF* np.ones((Ns,Ns))
    for mol in range(Ns):
        for j in range(mol+1,Ns):
            tc_minimum=conf.INF
            for trans in pos_trans:
                rij = pos[mol] + trans - pos[j] #27x3 for posj in flatten(pos_trans)
                uij = vel[mol] - vel[j]
                am = sum(uij**2)
                bm = np.sum(rij*uij)
                cm = sum(rij**2) - sigma**2
                if bm < 0:
                    dm = bm**2 - am*cm
                    if dm > 0:
                        tc_new = (-bm - np.sqrt(dm))/am
                        if tc_new < tc_minimum:
                            tc_minimum = tc_new
                            test.icheck(np.any(tc_minimum < 0), "error! negative time!", False)
            ctime[mol, j] = tc_minimum
            ctime[j, mol] = tc_minimum
    test.icheck(np.any(ctime) < 0, "initial collision time with negtive t_cols !!!", True)
    return ctime