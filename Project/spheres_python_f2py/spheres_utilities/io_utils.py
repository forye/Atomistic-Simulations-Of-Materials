__author__ = 'Idan'

import configuration as conf

def write_initial(nspheres,rvolume,ncollisions,pos,vel):
    '''
    ! This subroutine writes the initial positions and velocities of a
    ! set of spheres to an output file.
    '''

    filename = conf.RESULTS# "results.txt"
    if (conf.debug_mode): print'ENTERING SUBROUTINE write_initial'
    res = open(filename,'w')
    res.write( str(nspheres)+'\n')
    res.write( str(rvolume) +'\n')
    res.write(str(ncollisions) +'\n')
    for i in range(nspheres):
        res.write( str(i) +' '+ str( pos[i][0])+' '+ str( pos[i][1])+' '+ str( pos[i][2])+' ' + str(vel[i][0])+' ' + str(vel[i][1])+' ' + str(vel[i][2]) + '\n')
    return filename

def write_properties(nc,mom,kin,tcol,icol,jcol,last_delta_vel):
    '''
    ! This subroutine writes to a file the relevant information it receives:
    ! simulation step (c), average linear momentum (mom), average kinetic
    ! energy (kin), time to collision (tcol), indices of balls that collided
    ! (icol, jcol), and change in velocities at the collision (delta_vel).
    '''
    filename = conf.RESULTS
    if (conf.debug_mode): print'ENTERING SUBROUTINE write_properties 2'
    res = open(filename,'a')
    res.write(str(nc))
    res.write('\n')
    res.write(str( mom[0]) +' ' +str( mom[1])+' ' +str( mom[2]) )
    res.write('\n')
    res.write(str( kin))
    res.write('\n')
    res.write(str( tcol))
    res.write('\n')
    res.write( str(icol)+' '+str( jcol))
    res.write('\n')
    res.write(str(last_delta_vel[0])+' ' + str(last_delta_vel[1])+' ' + str(last_delta_vel[2]) )
    res.write('\n')
    res.close()
