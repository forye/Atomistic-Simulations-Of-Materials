! Exercise 4: Implementing Subroutine write properties

subroutine write_properties(c,mom,kin,tcol,icol,jcol,delta_vel)
!  This subroutine writes to a file the relevant information it receives:
!  simulation step (c), average linear momentum (mom), average kinetic
!  energy (kin), time to collision (tcol), indices of balls that collided
!  (icol, jcol), and change in velocities at the collision (delta_vel).
    implicit none
    integer, intent(in) :: c
    double precision, intent(in) :: mom(1:3)
    double precision, intent(in) :: kin
    double precision, intent(in) :: tcol
    integer, intent(in) :: icol
    integer, intent(in) :: jcol
    double precision, intent(in) :: delta_vel(1:3)
    character(11), parameter :: filename="results.txt"
    open (unit=1, file=filename, access="append", status="old")
    write (1,*) c
    write (1,*) mom
    write (1,*) kin
    write (1,*) tcol
    write (1,*) icol, jcol
    write (1,*) delta_vel
    close (1)
end subroutine write_properties

