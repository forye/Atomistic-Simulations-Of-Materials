__author__ = 'Idan'
import numpy as np
import Project.spheres_python.spheres_utilities.testing_utils as test
import Project.spheres_python.spheres_utilities.init_simulation as it
import Project.spheres_python.spheres_utilities.simulation_utils as sim
import Project.spheres_python.spheres_utilities.io_utils as io
import Project.spheres_python.spheres_utilities.configuration as conf

# class spheres(object):

def spheres(Rvt=None, n =3 , Ncolissions = 100,r= None,L=1.0,kT=1,speedScale = 1.0 , original_write = True):#0.0001

    vel_t = []
    kin_evol = []
    time_total, t_now, delta_vel, mom = [], 0, np.array([0, 0, 0]), np.array([0, 0, 0])

    notification_interval = 4
    collisions_to_report = range(Ncolissions / notification_interval, Ncolissions, Ncolissions / notification_interval)
    collisions_to_report.append(Ncolissions - 1)
    pos_trans = sim.get_boundry_condition(L)
    boundary_condition = sim.get_pos_period_bc_trans(L)
    # boundary_condition = pos_trans.reshape(1, 27, 3)/


    Ns, Rv, Ncollisions, sigma, n, L = it.init_parameters(n, Ncolissions, r, Rvt, L)
    pos = it.init_positions(L, n, Ns)
    vel, avg_kin, avg_mom = it.init_velocities(Ns, kT, speedScale)
    c_time = it.init_collision_time(pos, vel, boundary_condition, sigma)
    io.write_initial(Ns, Rv, Ncolissions, pos, vel)

    kin = (0.5/Ns) * np.sum(vel ** 2)
    kin_old = kin
    sigma_sq = sigma ** 2

    for nc in range(Ncollisions):
        if nc in collisions_to_report:
            print "v/v0 ={0}: Collisions simulation is at {1}% of {2} collisions" \
                .format(str(Rv), (nc + 1) * 100.0 / Ncolissions, str(Ncolissions))

        # tcol, jcol, icol = retrieve_collision_info(c_time)

        min_t = np.argmin(c_time)
        jcol, icol = min_t % Ns, min_t / Ns
        tcol = c_time[icol, jcol]

        t_now += tcol
        c_time = c_time - tcol

        if conf.debug_mode:
            test.icheck(conf.debug_mode, ' Minimum collision time:' + str(c_time[icol, jcol] ), False)
            test.icheck(conf.debug_mode, ' Minimum indices:' + str(icol) + str(jcol), False)
            test.icheck(c_time[icol, jcol] != c_time[jcol, icol], ", at col = "+str(nc) +
                                " time matrix is wrong\n not simetric!" + str([icol, jcol]), False)
            test.icheck(c_time[icol, jcol] <= 0, "ERROR! bad collision time, at col = " +
                                str(nc) + " t_col=" + str(c_time[icol, jcol]))
            test.icheck(c_time[icol, jcol] > 3600 * L * speedScale,
                                "ERROR! bad collision time, at col = "+str(nc) + " t_col=" + str(c_time[icol, jcol]))
            test.icheck(np.any(c_time < 0), "negative collision time, a collision was passed on nc=" + str(nc), True)

        #save itteration results
        if original_write: io.write_properties(nc, mom, kin, tcol, icol, jcol, delta_vel)
        # else: sr_old.write_status(pos,vel,tcol,icol,jcol,t_now)

        # advance simulation :
        #   move the spheres


        # mom, kin, pos, vel, delta_vel = advance_simulation(tcol, icol, jcol, nc, pos,
        #                                     vel, L, boundary_condition, sigma_sq, Ns, sigma)
        pos, vel, delta_vel, mom, kin = advance_simulation(pos, vel, tcol, icol, jcol, sigma_sq, sigma, Ns, nc,
                                                           boundary_condition, conf, L)



        # update collision table
        c_time[icol, jcol] = conf.INF
        c_time[jcol, icol] = conf.INF
        for i in [icol, jcol]:
            time = update_ctime_row(i, pos, vel, sigma, pos_trans)
            c_time[i, :] = time
            c_time[:, i] = time

        if conf.debug_mode:
            test.icheck(np.any(c_time < 0),
                    "c_time contaiins negative collision time at tcol=" + str(np.min(c_time)) + " nc=" + str(nc), True)
        vel_t.append(vel)
        time_total.append(t_now)
        if conf.debug_mode:
            test.icheck(abs(1 - kin / kin_old) > 0.05, " somthing is worong. diffrence in kinetic energy")
        kin_old = kin
        kin_evol.append(kin)

        if original_write: io.write_properties(nc, mom, kin, tcol, icol, jcol, delta_vel)


    #final test
    if conf.debug_mode:
        val1 = np.sum((pos - pos[:, np.newaxis]) ** 2, 2) - sigma_sq < -sigma_sq * 10 ** (-2)
        val2 = np.sum((pos - pos[:, np.newaxis]) ** 2, 2) > 0
        for i, u in enumerate(val1):
            for j, v1 in enumerate(u):
                test.icheck(v1 and val2[i, j]
                            , "final test:\nproblem with the simulation: in indexes [i,j]=" + str([i, j])
                            + "\n cought some spheres that have collapesd to each other. "
                            , False)
    return kin, mom, delta_vel, time_total, vel_t, kin_evol
#
#
# def advance_simulation(tcol, icol, jcol, nc, pos, vel, L, boundary_condition, sigma_sq, Ns, sigma):
#
#     pos = (pos + tcol * vel) % int(L)
#     # boundary_condition = boundary_condition  # sim.get_pos_period_bc_trans(L)
#     uij = vel[icol] - vel[jcol]
#     rij = pos[icol] - pos[jcol]  # the current colliding spheres distance
#
#     rij_image = (boundary_condition + rij)[np.argmin(np.sum((boundary_condition + rij) ** 2, 1))]
#     bij = rij_image.dot(uij)  # np.linalg.norm(rij_image) #bij = radial velocity component
#     delta_vel = (-rij_image / sigma_sq) * bij  # sigma**2 ==rij_image**2 @ collision. safer then sigma
#     vel[icol] += delta_vel
#     vel[jcol] -= delta_vel
#
#     mom = (1. / Ns) * np.sum(vel, 0)
#     kin = (0.5 / Ns) * np.sum(vel ** 2)
#     if conf.debug_mode:
#         print' Average linear mom:', mom
#         print' Average kinetic e:', kin
#
#     if conf.debug_mode:
#         test.icheck(sum(rij_image ** 2) - sigma ** 2 > L * 10 ** (-6),
#                     " ERROR distance in collision is not the diameter! " + str(np.sqrt(sum(rij_image ** 2))) + "  i,j="
#                     + str([icol, jcol]) + " \n diameter is :" + str(sigma) + "\n iter: "
#                     + str(nc) + "\nat t_col=" + str(tcol), True)
#
#
#     return mom, kin, pos, vel, delta_vel

def retrieve_collision_info(c_time):
    Ns = c_time.shape[0]
    min_t = np.argmin(c_time)
    jcol, icol = min_t % Ns, min_t / Ns
    tcol = c_time[icol, jcol]
    return tcol, jcol, icol


def update_ctime_row(i, pos, vel, sigma, trans_):

    Ns = np.shape(pos)[0]
    sigma_sq = sigma ** 2

    pos_ = pos.reshape(Ns, 1, 3)
    rij_ = pos_[i] + trans_ - pos_
    rij_[i] = np.zeros([27, 3])
    uij = np.tile(vel[i], [Ns, 1]) - vel
    uij_ = np.tile(uij, [1, 27, ]).reshape(Ns, 27, 3)

    aij = np.zeros([Ns, 27])
    bij = np.zeros([Ns, 27])
    cij = np.zeros([Ns, 27])
    disc = np.zeros([Ns, 27])
    time_ = np.ones([Ns, 27]) * np.inf

    bij = np.sum(uij_ * rij_, axis=2)
    col_mask = bij < 0  # mask only aproaching speres

    col_mask[i] = np.zeros((27), dtype=bool)
    am_ = np.sum(uij ** 2, axis=1)

    # aij[col_mask] = np.tile(am_, [1, 27]).reshape(Ns, 27)[col_mask]
    # aij = np.tile(am_, [1, 27]).reshape(Ns, 27)
    aij = np.sum(uij_**2, axis =2)

    cij[col_mask] = np.sum(rij_ ** 2, axis=2)[col_mask] - sigma_sq

    disc[col_mask] = bij[col_mask] ** 2 - aij[col_mask] * cij[col_mask]

    col_mask2 = np.array(col_mask)
    col_mask2[col_mask] = disc[col_mask] > 0  # mask only apraching spheres, that are colliding spheres

    col_mask3 = np.array(col_mask2)
    col_mask3[col_mask2] = aij[col_mask2] != 0
    time_[col_mask3] = -(bij[col_mask3] + np.sqrt(disc[col_mask3])) / aij[
        col_mask3]  # stops at col=30 neg t- numeric

    time = np.min(time_, axis=1)
    return time


def advance_simulation(pos, vel, tcol, icol, jcol, sigma_sq, sigma, Ns, nc, boundary_condition, conf, L):
    pos = (pos + tcol * vel) % int(L)
    # boundary_condition = sim.get_pos_period_bc_trans(L)
    uij = vel[icol] - vel[jcol]
    rij = pos[icol] - pos[jcol]  # the current colliding spheres distance

    rij_image = (boundary_condition + rij)[np.argmin(np.sum((boundary_condition + rij) ** 2, 1))]
    bij = rij_image.dot(uij)  # np.linalg.norm(rij_image) #bij = radial velocity component
    delta_vel = (-rij_image / sigma_sq) * bij  # sigma**2 ==rij_image**2 @ collision. safer then sigma
    vel[icol] += delta_vel
    vel[jcol] -= delta_vel

    if conf.debug_mode:
        test.icheck(sum(rij_image ** 2) - sigma ** 2 > L * 10 ** (-6),
                    " ERROR distance in collision is not the diameter! " + str(np.sqrt(sum(rij_image ** 2))) + "  i,j="
                    + str([icol, jcol]) + " \n diameter is :" + str(sigma) + "\n iter: "
                    + str(nc) + "\nat t_col=" + str(tcol), True)

    mom = (1. / Ns) * np.sum(vel, 0)
    kin = (0.5 / Ns) * np.sum(vel ** 2)
    if conf.debug_mode:
        print' Average linear mom:', mom
        print' Average kinetic e:', kin
    return pos, vel, delta_vel, mom, kin
