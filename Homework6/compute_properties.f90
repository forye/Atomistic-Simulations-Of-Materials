! Exercise 3: Implementing Subroutine compute properties

subroutine compute_properties(vel,mom,kin)
!  Given the velocity vectors (vel) of a set of spheres of unit mass,     
!  this subroutine returns their average linear momentum and their average
!  kinetic energy.
    implicit none
    double precision, intent(in) :: vel(:,:)
    double precision, intent(out) :: mom(1:3)
    double precision, intent(out) :: kin
    integer :: i
    integer :: nspheres

    nspheres = size(vel,2)
    ! Compute average linear momentum:
    mom = (/ 0.0d0 , 0.0d0 , 0.0d0 /)
    do i = 1, nspheres
        mom(1:3) = mom(1:3) + vel(1:3,i)
    end do
    mom(1:3) = (1.0d0/nspheres) * mom(1:3)

    ! Compute average kinetic energy:
    kin = 0.0d0
    do i = 1, nspheres
        kin = kin + dot_product(vel(1:3,i), vel(1:3,i))
    end do

    
    kin = (0.5d0/nspheres) * kin
    !write(*,*) ' Average kinetic energy (where KT=1):', kin
    
end subroutine compute_properties


