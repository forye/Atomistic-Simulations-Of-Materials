__author__ = 'Idan'

import numpy as np
import visual as v

import configuration as conf
import init_simulation as init
import testing_utils as test
import os

# sudo ~/anaconda/bin/pythonw personal/hagasha/Project/spheres_python/visualize.py


def visualize_last_result(dt=0.0000005, rate=1000 ,L=0, secconds=0, result_path=None):
    print "os.getwd()", os.getcwd(), "conf.RESULTS", conf.RESULTS
    if 'spheres_python' in os.getcwd():
        result = open('results.txt').readlines()
    else:
        if not result_path : result = open('personal/hagasha/Project/spheres_python/results.txt').readlines()
        else: result =open(result_path).readlines()

    if L == 0:
        L = conf.L

    pos=[]
    vel=[]
    # tcol = []


    cur = 0
    Ns = int(result[cur])
    cur+=1
    rv = float(result[cur])
    cur+=1
    Ncollisions =int(result[cur])
    for ns in range(Ns):
        cur+=1
        row = result[cur].split()
        pos.append( np.array( [ float(row[1])/L,float(row[2])/L,float(row[3])/L ]))
        vel.append( np.array(  [ float(row[4])/L,float(row[5])/L,float(row[6])/L]))

    Collisions =[]
    time=0.0

    while cur < len(result)-1:
        d = dict()
        cur+=1
        d['n_col'] = int(result[cur])
        cur+=1
        d['row'] = result[cur].split()
        mom = np.array ( [ float(row[0]),float(row[1]),float(row[2])] )
        cur+=1
        d['kin'] = float(result[cur])
        cur+=1
        d['tcol'] = float(result[cur])
        time += d['tcol']
        d['time'] = time
        cur+=1
        row= result[cur].split()
        d['icol'],d['jcol'] = int(row[0]) ,int(row[1])
        cur+=1
        row = result[cur].split()
        d['delta_vel'] = np.array ( [ float(row[0])/L,float(row[1])/L,float(row[2])/L] )
        Collisions.append(d)
    test.icheck(Ncollisions == len(Collisions) , " not enough collision in file while reading to visualize.\n Ncollisions = "
                + str(Ncollisions) + "len(Collisionsol) = " + str(len(Collisions)),False)


    if secconds != 0 :
        total_time = Collisions[-1]['time']
        dt = total_time / (secconds * 24)
        rate = int(24)
        # dt = step size
        # dt = total_time / (secconds * 1000)
    visualize_pheres_col(L, Ns, init.compute_diameter(Ns, rv), Collisions, pos, vel, dt, rate)
    return



def visualize_pheres_col(L,Ns,sigma ,collisions,pos,vel,dt,_rate):
    v.scene.width = 400
    v.scene.height = 400
    screen_2_lattice = 3.0
    X=screen_2_lattice
    Y=screen_2_lattice
    Z=screen_2_lattice
    cent= screen_2_lattice* 0.25*0.25
    v.scene.range = (X,Y,Z)
    v.scene.center = (cent,cent,cent)

    Lo=L
    if Lo is not 1.0: print "L is : " + str(Lo)
    # min_t = min([col["tcol"] for col in collisions])

    RATE =_rate#* Ns
    balls =[]*Ns
    for ns in range(Ns):#create list
        balls.append(v.sphere( pos = (pos[ns][0],pos[ns][1],pos[ns][2]),radius = sigma/2.0 ))

    tend = float( collisions[-1]['time'])
    # dt = dt * tend / len(collisions)
    timeline = np.arange(dt,tend,dt)

    col_idx=0

    ic=0
    for t in timeline:
        # '''
        # t_0: the time that is left untill the next continues time sample
        # RATE: number of itteration per sec (uppper bound)
        # '''
        v.rate( RATE )
        t_0 = dt
        ex_t = 0
        #handle collisions
        # while col_idx != len( collisions)-1 -1 and collisions[ col_idx][ 'time'] <= t + t_0:#dt:
        while col_idx < len( collisions)-1  and collisions[ col_idx ][ 'time'] <= t + dt - ex_t:#dt:t?
            ex_t += collisions[ col_idx][ 'time'] -t
            tm = collisions[col_idx]['time'] -t
        # while col_idx != len(collisions)-1 -1 and collisions[col_idx ]['tcol'] <=  t_0:# + t_0:
            #bad collisions-> in case of 2 collision, what is the updatoing while rule
            #the velocity is not correct!!! first collision yeilds 0,0,0 detla vel- added 1 to col_idx and reduced 1 frp, col idx. i dont need last collision without delta vel
            for ns in range(Ns):
                # balls[ns].pos = v.vector(pos[ns])+ v.vector(vel[ns])*collisions[col_idx]['tcol']
                # balls[ns].pos = (balls[ns].pos + v.vector(vel[ns])*collisions[col_idx]['tcol'])%L
                # temp =  balls[ns].pos+ v.vector(vel[ns])* ( collisions[col_idx]['time'] - t)#t_0
                pos_new = balls[ns].pos + v.vector( vel[ns] * tm) #) collisions[col_idx]['tcol'] )#
                balls[ns].pos = v.vector([ pos_new[0] % Lo , pos_new[1] % Lo , pos_new[2] % Lo ] )
            i_col = collisions[col_idx]['icol']
            j_col = collisions[col_idx]['jcol']
            vel[ i_col ] = np.array( vel[ i_col ] + collisions[ col_idx+1 ]['delta_vel'])
            vel[ j_col ] = np.array( vel[ j_col ] - collisions[ col_idx+1 ]['delta_vel'])
            # t_0=t_0 -collisions[col_idx ]['tcol']
            t_0=t_0 -tm

            # t_0=t_0 -collisions[col_idx ]['tcol'] #(collisions[col_idx]['time']-t)
            # vel_new = np.array(vel[ i_col ] + collisions[ col_idx + 1]['delta_vel'])
            # vel[ i_col ] = vel_new
            # vel_new = np.array( vel[ j_col ] - collisions[ col_idx + 1]['delta_vel'])
            # vel[j_col] = vel_new
            # t - collisions[col_idx]['time']

            # t_0=t_0 - ( t - collisions[col_idx]['time'])
            if t_0 < 0 :
                print"error t_0 < 0"
                test.icheck(t_0 < 0, "error in visualization. after calculating  collision, t_0 < 0\n in col=" +str(col_idx),False)
                col_idx+=1
                continue
            col_idx+=1
            #paint colliding next
            if ns in [collisions[col_idx]['icol'],collisions[col_idx]['jcol']]:balls[ns].color = v.color.red
            else: balls[ns].color = v.color.white

        #t_o = dt-ex_t
        #advance regular
        for ns in range(Ns):#finish time from the last collision to th end
            # balls[ns].pos = v.vector(pos[ns])+ v.vector(vel[ns])*t_0
            pos_new =balls[ns].pos + v.vector(vel[ns]) * t_0

            for q in range(Ns):
                if q is ns : continue
                cond21 = np.sum( (np.array(pos_new) - np.array(balls[q].pos)) **2)  < sigma**2
                if cond21:
                    test.icheck( cond21 , "error in simuilation, spheres are collapsing to each other, on " + str([q,ns]),False)
                    break
                if cond21 : break

            balls[ns].pos =v.vector([pos_new[0] % Lo, pos_new[1] % Lo, pos_new[2] % Lo])
            # if ns in [collisions[col_idx]['icol'], collisions[col_idx]['jcol']]:
            #     balls[ns].color = v.color.red
            # else: balls[ns].color = v.color.white
        #
        # ic +=1
        # if not ic %10:
        #     for p in range(Ns):
        #         for q in range(Ns):
        #             if q is p : continue
        #             cond21 = np.sum( (np.array( balls[q].pos) - np.array(balls[p].pos)) **2)  < sigma**2
        #             if cond21:
        #                 test.icheck( cond21 , "error in simuilation, spheres are collapsing to each other, on " + str([q,p]),False)
        #                 break
        #         if cond21 : break



    print "done visualization!"





# class sp:
#     n_col
#
#     def __init__(self, cur , result):
#         if cur < len(result):
#             cur+=1
#             n_col = int(result[cur])
#             cur+=1
#             row = result[cur].split()
#             mom = np.array ( [ float(row[0]),float(row[1]),float(row[2])] )
#             cur+=1
#             kin = float(result[cur])
#             cur+=1
#             tcol = float(result[cur])
#             cur+=1
#             row = result[cur].split()
#             icol,jcol = int(row[0]) ,int(row[1])
#             cur+=1
#             row = result[cur].split()
#             delta_vel = np.array ( [ float(row[0]),float(row[1]),float(row[2])] )
   # # visualize_spheres_vel(rep.pos_t,rep.vel_t,rep.pos_t[0],t_range,rep.t_totals ,_sigma = rep.sigma,
    # #                                                    L=rep.L,_rate = rate,icol=rep.icols,jcol=rep.jcols)
    # print "visualizing ended at: ", str(t_range[-1])



    # t_range =[]
    # ti=0
    # t_range.append(ti)

    # while ti <= rep.t_totals[-1]:
    #     ti+=dt
    #     t_range.append(ti)

    # #make sure that the time vector is
    # if len(t_range) is not np.shape( rep.pos_t )[0]:
    #     while len(t_range) < np.shape( rep.pos_t )[0]:
    #         t_range.append(t_range[-1] + dt)