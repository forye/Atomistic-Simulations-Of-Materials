
subroutine assign_velocities(vel)
!  This subroutine fills up array vel with three-dimensional vectors that
!  represent the velocities of spheres (as many spheres as the size of the
!  second dimension of vel). The module of the vectors is set to sqrt(3)
!  so that the energy per sphere is 3/2 in reduced units (kT = 1). The
!  direction of the vector is random, with equal probability of lying
!  in any area element of a sphere.

    implicit none
    double precision, intent(inout) :: vel(:,:)
    ! start here
    double precision :: alpha
    double precision :: avg_kin
    double precision :: avg_vel(1:3)
    integer :: i
    integer :: nspheres
    double precision :: phi
    double precision, parameter :: pi = 3.1415926535897932
    double precision :: speed
    double precision :: theta

    double precision :: u
    double precision :: v
    
    nspheres = size(vel,2)
    speed = sqrt(3.0)
    avg_vel(1:3) = 0.0d0
    call random_seed
        
    ! init velocoties and calculate velocity
    do i = 1, nspheres
        call random_number(u)
        call random_number(v)
        theta = 2.0d0 * pi * u
        phi = acos(2.0d0 * v - 1.0d0)
        vel(1,i) = speed * sin(phi) * cos(theta)
        vel(2,i) = speed * sin(phi) * sin(theta)
        vel(3,i) = speed * cos(phi)
        avg_vel(1:3) = avg_vel(1:3) + vel(1:3,i)
    end do    
    avg_vel(1:3) = avg_vel(1:3) / dble(nspheres)
    
    ! standatizize volocities
    do i = 1, nspheres
        vel(1:3,i) = vel(1:3,i) - avg_vel(1:3)
    end do     
    
end subroutine assign_velocities