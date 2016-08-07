! Exercise 2: Implementing Subroutine advance simulation

subroutine advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
!  Given the time to the next collision (tcol) of a system of spheres of given
!  diameter (sigma), and the indices of the two spheres involved in that
!  collision (icol,jcol), this subroutine modifies the positions (pos)
!  and velocities (vel) of all the spheres to reflect the situation
!  immediately after that collision. It also returns the change in velocities
!  of the spheres that collide (delta_vel). The periodic boundary conditions of
!  the simulation cell containing the spheres are taken into account.
    implicit none
    double precision, intent(in) :: tcol
    integer, intent(in) :: icol
    integer, intent(in) :: jcol
    double precision, intent(in) :: sigma
    double precision, intent(inout) :: pos(:,:)
    double precision, intent(inout) :: vel(:,:)
    double precision, intent(out) :: delta_vel(1:3)   

    double precision :: bij
    double precision :: distance2
    integer :: i
    double precision, parameter :: infinity = huge(0.0d0)
    integer :: jx
    integer :: jy
    integer :: jz
    double precision :: min_distance2
    integer :: nspheres
    double precision :: pos_image(1:3)
    double precision :: rij(1:3)
    double precision :: rij_image(1:3)
    double precision :: translat(1:3)
    double precision :: uij(1:3)
    

    nspheres = size(pos,2)
    ! Update positions of all spheres:
    do i = 1, nspheres
        ! Compute new position vectors:
        pos(1:3,i) = pos(1:3,i) + vel(1:3,i) * tcol
        ! If new position is outside simulation cell, bring it inside:
        pos(1,i) = pos(1,i) - floor(pos(1,i))
        pos(2,i) = pos(2,i) - floor(pos(2,i))
        pos(3,i) = pos(3,i) - floor(pos(3,i))
    end do

    ! Update the new velocities of the colliding spheres spheres:
    uij(1:3) = vel(1:3,icol) - vel(1:3,jcol)        
    min_distance2 = infinity
    do jx = -1, 1
        do jy = -1, 1
            do jz = -1, 1
                translat = (/ jx, jy, jz /)
                pos_image(1:3) = pos(1:3,jcol) + translat(1:3)
                rij_image(1:3) = pos(1:3,icol) - pos_image(1:3)
                distance2 = dot_product(rij_image(1:3),rij_image(1:3))
                if (distance2 < min_distance2) then
                    rij(1:3) = rij_image(1:3)
                    bij = dot_product(rij(1:3),uij(1:3))
                    min_distance2 = distance2
                end if
            end do
        end do
    end do
    delta_vel(1:3) = - (bij/(sigma*sigma)) * rij(1:3)
    vel(1:3,icol) = vel(1:3,icol) + delta_vel(1:3)
    vel(1:3,jcol) = vel(1:3,jcol) - delta_vel(1:3)
end subroutine advance_simulation


