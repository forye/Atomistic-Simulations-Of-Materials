__author__ = 'Idan'

import numpy as np
import configuration as conf
import testing_utils as test

def get_next_collision_time(ctime , Ns ,nc=0):
    min_t = np.argmin(ctime)
    jcol , icol = min_t % Ns , min_t / Ns

    test.icheck(conf.debug_mode,' Minimum collision time:'+str(ctime[icol,jcol] ),False)
    test.icheck( conf.debug_mode ,' Minimum indices:'+ str( icol )+str( jcol) , False)
    test.icheck( ctime[icol,jcol] != ctime[jcol,icol] , ", at col = "+str(nc) +" time matrix is wrong\n not simetric!" + str([icol,jcol]) ,False)
    test.icheck( ctime[icol,jcol]<= 0 ,  "ERROR! bad collision time, at col = "+str(nc) +" t_col=" +  str(ctime[icol,jcol]) )
    test.icheck( ctime[icol,jcol]> 1000000 , "ERROR! bad collision time, at col = "+str(nc) +" t_col=" +  str(ctime[icol,jcol]))

    tcol = ctime[icol,jcol]
    return tcol,icol,jcol


pos_trans = []
def get_pos_period_bc_trans(L=1.0):
    '''
    boundry conditions
    :param L: L is the scale
    :return: 3x3 cube of 9 possible positions
    '''
    global pos_trans
    if len(pos_trans)==0:
        # pos_trans = np.arange(1)
        trans=[]
        for x in xrange(-1,2):
            for y in xrange(-1,2):
                for z in xrange(-1,2):
                    trans.append(float(L) * np.array([x,y,z]) )
        pos_trans = np.array(np.array(trans))
    return pos_trans


def get_boundry_condition(L=1.0):
    # if not trans:
    arng = np.arange(-1, 2, 1)
    trans_ = np.array(
        [[ix, iy, iz] for ix in arng for iy in arng for iz in arng])
    trans_ = trans_.reshape(1, 27, 3)
    return float(L) * trans_


def update_collision_row(i,Ns,icol,jcol,pos,sigma_sq,L,sigma,nc,vel,c_time):
    v_rij =pos[i] - pos[:,np.newaxis] + pos_trans  # pos[i] - pos + pos_trans # check , index of trans is 2- have 3 dimentions- sphere,xyz,trans?
    v_vij = vel[i] - vel
    v_aij = np.sum(v_vij**2, 1)
    v_inv__aij = np.reshape( np.tile( 1 / np.sum(v_vij**2, 1) , 27),(108,27) )
    # v_vij[:,np.newaxis]
    v_vij_t = np.reshape( np.tile( v_vij , 27),(108,27,3) )
    v_bij =np.sum(v_rij*v_vij_t,2 )
    # v_bij = v_rij.dot(v_vij, 1)  # check!
    v_cij = np.sum(v_rij**2, 2) - sigma_sq

    sq_b =  v_bij**2   # check
    ac = v_aij[:,np.newaxis] * v_cij  #  np.sum( v_aij * v_cij , 2)  # check
    v_cond1 = np.where(v_bij < 0 , sq_b > ac , False )  # np.where(v_bij < 0 , sq_b > ac , False )  #  & v_cij < 0)
    v_cond1[i,:] = np.zeros((1,27),bool)
# ERROR! bad collision time, at col = 2 t_col=0.0
    index_of_trans = 2
    tcols_new =np.where(v_cond1, (-v_bij - np.sqrt(sq_b - ac ))/ ( v_inv__aij ) , np.inf)  # check- that matrixes are aligned, that the inverse is done
    res = np.min(tcols_new,1)
    c_time[i,:] = res
    c_time[:,i] = res
    test.icheck(c_time[i,i]< 10^20 ,"computation error ! in mol ="+str(i)+ ", error- in nc = "+str(nc)+"!\n sphere collides with itself :  c_time[i,i] <inf "+ str(i)  )
    return c_time
    # tcols_new <- minimal in trans
    # c_time[i,:] = tcols_new
    # c_time[:,i] = tcols_new
    # for j in range(Ns):
    #     if j in [icol,jcol]: continue
    #     tc_minimum = conf.INF
    #     rij_base = pos[i] - pos[j]
    #     for trans in pos_trans:
    #         rij = rij_base + trans  #27x3 for posj in flatten(pos_trans)
    #         if conf.debug_mode:
    #             if np.abs( np.sum( rij**2 ) -sigma_sq  ) <  10.**(-15)*L :
    #                 print " after collidion of "+str([icol,jcol])+". distance between spheres on nc="\
    #                       +str(nc) + " i=" + str(i) + " and j=" +str(j) + ", is d = " \
    #                       + str(np.abs( np.sum( rij**2 ) -sigma**2  ))
    #         uij = vel[i] - vel[j]
    #         am = sum( uij**2)
    #         bm = np.sum( rij *uij )
    #         cm = sum( rij**2) - sigma_sq
    #         if cm< 0 :
    #             print "cm is negative: spheres are colliding each other in n=" + str(nc)
    #         if bm<0:#approaching -> vel makes distance gets smaller
    #             dm = bm**2 -am *cm
    #             if dm > 0 :# there are solutions to where the spheres meet
    #                 test.icheck(-bm - np.sqrt(float(dm)) < 0 , " somthing is wrong, got a negative collision time. got -bm - np.sqrt(float(dm)) =" + str(-bm - np.sqrt(float(dm))) ,False)#extra condition: make sure the collision is in the future
    #                 tc_new = (-bm - np.sqrt(float(dm)))/am
    #                 test.icheck( tc_new <0, "negative predicted collision time!! t=" +str( tc_new) + " nc = " +str(nc) + " col index = " +str([icol,jcol])+ "  [i,j] = " +str([i,j]) , False)
    #                 test.icheck( cm < -1*10.**(-2), "error in col "+str(nc)+ ", spheres are crushing in to each other")
    #                 if tc_new < tc_minimum :#sooner collision i
    #                     tc_minimum = tc_new
    #                     test.icheck(tc_minimum <0 ,"error! negative time!",False)
    #                     test.icheck(tc_new > L*3600 ,"too big prdeiced collision time, t=" +str( tc_minimum),False)
    #
    #     c_time[i,j] = tc_minimum
    #     c_time[j,i] = tc_minimum
    # test.icheck(c_time[i,i]< 10^20 ,"computation error ! in mol ="+str(i)+ ", error- in nc = "+str(nc)+"!\n sphere collides with itself :  c_time[i,i] <inf "+ str(i)  )
    # return c_time




def update_collision_row_old(i,Ns,icol,jcol,pos,sigma_sq,L,sigma,nc,vel,c_time):
    '''

    :param i: colliding sphere
    :param Ns:
    :param icol:
    :param jcol:
    :param pos:
    :param sigma_sq:
    :param L:
    :param sigma:
    :param nc:
    :param vel:
    :param c_time:
    :return:
    the updated collsion time matrix by row
    '''
    '''
    :param i: row to be updated
    :param Ns: total number of rows
    :param icol: last collided sphere #1
    :param jcol: last collided sphere #2
    :param pos: all spheres positions
    :param sigma_sq: square of diameter
    :param L: total box axis side length
    :param sigma:diameter
    :param nc: the undex of the collision
    :param vel: all spheres velocities
    :param c_time: a matrix that indicate the time of next collision for each sphere
    :return: c_time
    '''
    for j in range(Ns):
        if j in [icol,jcol]: continue
        tc_minimum = conf.INF
        rij_base = pos[i] - pos[j]
        for trans in pos_trans:
            rij = rij_base + trans  #27x3 for posj in flatten(pos_trans)
            if conf.debug_mode:
                if np.abs( np.sum( rij**2 ) -sigma_sq  ) <  10.**(-15)*L :
                    print " after collidion of "+str([icol,jcol])+". distance between spheres on nc="\
                          +str(nc) + " i=" + str(i) + " and j=" +str(j) + ", is d = " \
                          + str(np.abs( np.sum( rij**2 ) -sigma**2  ))
            uij = vel[i] - vel[j]
            am = sum( uij**2)
            bm = np.sum( rij *uij )
            cm = sum( rij**2) - sigma_sq
            if cm< 0 :
                print "cm is negative: spheres are colliding each other in n=" + str(nc)
            if bm<0:#approaching -> vel makes distance gets smaller
                dm = bm**2 -am *cm
                if dm > 0 :# there are solutions to where the spheres meet
                    test.icheck(-bm - np.sqrt(float(dm)) < 0 , " somthing is wrong, got a negative collision time. got -bm - np.sqrt(float(dm)) =" + str(-bm - np.sqrt(float(dm))) ,False)#extra condition: make sure the collision is in the future
                    tc_new = (-bm - np.sqrt(float(dm)))/am
                    test.icheck( tc_new <0, "negative predicted collision time!! t=" +str( tc_new) + " nc = " +str(nc) + " col index = " +str([icol,jcol])+ "  [i,j] = " +str([i,j]) , False)
                    test.icheck( cm < -1*10.**(-2), "error in col "+str(nc)+ ", spheres are crushing in to each other")
                    if tc_new < tc_minimum :#sooner collision i
                        tc_minimum = tc_new
                        test.icheck(tc_minimum <0 ,"error! negative time!",False)
                        test.icheck(tc_new > L*3600 ,"too big prdeiced collision time, t=" +str( tc_minimum),False)

        c_time[i,j] = tc_minimum
        c_time[j,i] = tc_minimum
    test.icheck(c_time[i,i]< 10^20 ,"computation error ! in mol ="+str(i)+ ", error- in nc = "+str(nc)+"!\n sphere collides with itself :  c_time[i,i] <inf "+ str(i)  )
    return c_time