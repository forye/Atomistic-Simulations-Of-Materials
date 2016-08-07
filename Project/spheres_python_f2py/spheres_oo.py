__author__ = 'Idan'
import numpy as np
import spheres_utilities.testing_utils as test
import spheres_utilities.init_simulation as it
import spheres_utilities.simulation_utils as sim
import spheres_utilities.io_utils as io
import spheres_utilities.configuration as conf

class spheres(object):
    def __init__(self, Rvt=None, n =3 , Ncolissions = 100,r= None,L=1.0,kT=1,speedScale = 1.0 , original_write = True):#0.0001):
        self.Ncolissions = Ncolissions
        self.n = n
        notification_interval = 4
        self.collisions_to_report = range(Ncolissions / notification_interval, Ncolissions,
                                     Ncolissions / notification_interval)
        self.collisions_to_report.append(Ncolissions - 1)
        self.pos_trans = sim.get_boundry_condition(L)
        self.boundary_condition = sim.get_pos_period_bc_trans(L)
        self.Ns, self.Rv, self.Ncollisions, self.sigma, self.n, self.L = it.init_parameters(n, Ncolissions, r, Rvt, L)
        self.sigma_sq = self.sigma ** 2
        self.L = L
        self.kT = kT
        self.speedScale = speedScale
        self.original_write = original_write


    def run(self):#0.0001

        vel_t = []
        kin_evol = []
        time_total, t_now, delta_vel, mom = [], 0, np.array([0, 0, 0]), np.array([0, 0, 0])
        self.pos = it.init_positions(self.L, self.n, self.Ns)
        self.vel, avg_kin, avg_mom = it.init_velocities(self.Ns, self.kT, self.speedScale)
        self.c_time = it.init_collision_time(self.pos, self.vel, self.boundary_condition, self.sigma)
        io.write_initial(self.Ns, self.Rv, self.Ncolissions, self.pos, self.vel)
        kin = (0.5/self.Ns) * np.sum(self.vel ** 2)
        kin_old = kin

        for nc in range(self.Ncollisions):
            if nc in self.collisions_to_report:
                print "v/v0 ={0}: Collisions simulation is at {1}% of {2} collisions" \
                    .format(str(self.Rv), (nc + 1) * 100.0 / self.Ncolissions, str(self.Ncolissions))

            min_t = np.argmin(self.c_time)
            jcol, icol = min_t % self.Ns, min_t / self.Ns
            tcol = self.c_time[icol, jcol]
            t_now += tcol
            self.c_time = self.c_time - tcol

            if conf.debug_mode:
                test.icheck(conf.debug_mode, ' Minimum collision time:' + str(self.c_time[icol, jcol] ), False)
                test.icheck(conf.debug_mode, ' Minimum indices:' + str(icol) + str(jcol), False)
                test.icheck(self.c_time[icol, jcol] != self.c_time[jcol, icol], ", at col = "+str(nc) +
                                    " time matrix is wrong\n not simetric!" + str([icol, jcol]), False)
                test.icheck(self.c_time[icol, jcol] < 0, "ERROR! bad collision time, at col = " +
                                    str(nc) + " t_col=" + str(self.c_time[icol, jcol]))
                test.icheck(self.c_time[icol, jcol] > 3600 * self.L * self.speedScale,
                                    "ERROR! bad collision time, at col = "+str(nc) + " t_col=" + str(self.c_time[icol, jcol]))
                test.icheck(np.any(self.c_time < 0), "negative collision time, a collision was passed on nc=" + str(nc), True)

            #save itteration results
            if self.original_write: io.write_properties(nc, mom, kin, tcol, icol, jcol, delta_vel)
            # else: sr_old.write_status(pos,vel,tcol,icol,jcol,t_now)

            delta_vel, mom, kin = self.advance_simulation(tcol, icol, jcol, nc)
            # update collision table
            self.c_time[icol, jcol] = conf.INF
            self.c_time[jcol, icol] = conf.INF
            for i in [icol, jcol]:
                self.update_ctime_row(i)

            if conf.debug_mode:
                test.icheck(np.any(self.c_time < 0),
                        "c_time contaiins negative collision time at tcol=" + str(np.min(self.c_time)) + " nc=" + str(nc), True)
            vel_t.append(self.vel)
            time_total.append(t_now)
            if conf.debug_mode:
                test.icheck(abs(1 - kin / kin_old) > 0.05, " somthing is worong. diffrence in kinetic energy")
            kin_old = kin
            kin_evol.append(kin)

            if self.original_write: io.write_properties(nc, mom, kin, tcol, icol, jcol, delta_vel)


        #final test
        if conf.debug_mode:
            val1 = np.sum((self.pos - self.pos[:, np.newaxis]) ** 2, 2) - self.sigma_sq < -self.sigma_sq * 10 ** (-2)
            val2 = np.sum((self.pos - self.pos[:, np.newaxis]) ** 2, 2) > 0
            for i, u in enumerate(val1):
                for j, v1 in enumerate(u):
                    test.icheck(v1 and val2[i, j]
                                , "final test:\nproblem with the simulation: in indexes [i,j]=" + str([i, j])
                                + "\n cought some spheres that have collapesd to each other. "
                                , False)
        return kin, mom, delta_vel, time_total, vel_t, kin_evol


    def update_ctime_row(self,i):

        trans_ = self.pos_trans
        Ns = np.shape(self.pos)[0]
        sigma_sq = self.sigma_sq

        pos_ = self.pos.reshape(Ns, 1, 3)
        rij_ = pos_[i] + trans_ - pos_
        rij_[i] = np.zeros([27, 3])
        uij = np.tile(self.vel[i], [Ns, 1]) - self.vel
        uij_ = np.tile(uij, [1, 27, ]).reshape(Ns, 27, 3)
        cij = np.zeros([Ns, 27])
        disc = np.zeros([Ns, 27])
        time_ = np.ones([Ns, 27]) * np.inf
        bij = np.sum(uij_ * rij_, axis=2)
        col_mask = bij < 0  # mask only aproaching speres
        col_mask[i] = np.zeros((27), dtype=bool)
        aij = np.sum(uij_**2, axis =2)
        cij[col_mask] = np.sum(rij_ ** 2, axis=2)[col_mask] - sigma_sq
        disc[col_mask] = bij[col_mask] ** 2 - aij[col_mask] * cij[col_mask]
        col_mask2 = np.array(col_mask)
        col_mask2[col_mask] = disc[col_mask] > 0  # mask only apraching spheres, that are colliding spheres
        col_mask3 = np.array(col_mask2)
        col_mask3[col_mask2] = aij[col_mask2] != 0
        time_[col_mask3] = -(bij[col_mask3] + np.sqrt(disc[col_mask3])) / aij[col_mask3]
        time = np.min(time_, axis=1)
        self.c_time[i, :] = time
        self.c_time[:, i] = time

    def advance_simulation(self, tcol, icol, jcol, nc):
        self.pos = (self.pos + tcol * self.vel) % int(self.L)
        uij = self.vel[icol] - self.vel[jcol]
        rij = self.pos[icol] - self.pos[jcol]  # the current colliding spheres distance
        rij_image = (self.boundary_condition + rij)[np.argmin(np.sum((self.boundary_condition + rij) ** 2, 1))]
        bij = rij_image.dot(uij)
        delta_vel = (-rij_image / self.sigma_sq) * bij
        self.vel[icol] += delta_vel
        self.vel[jcol] -= delta_vel

        if conf.debug_mode:
            test.icheck(sum(rij_image ** 2) - self.sigma ** 2 > self.L * 10 ** (-6),
                        " ERROR distance in collision is not the diameter! " + str(np.sqrt(sum(rij_image ** 2))) +
                        "  i,j=" + str([icol, jcol]) + " \n diameter is :" + str(self.sigma) + "\n iter: " +
                        str(nc) + "\nat t_col=" + str(tcol), True)

        mom = (1. / self.Ns) * np.sum(self.vel, 0)
        kin = (0.5 / self.Ns) * np.sum(self.vel ** 2)
        if conf.debug_mode:
            print' Average linear mom:', mom
            print' Average kinetic e:', kin
        return delta_vel, mom, kin
