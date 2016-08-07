! Exercise 5: Implementing Subroutine update collisions table
subroutine update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)
!  Given the positions (pos) and velocities (vel) of a set of spheres of given
!  diameter (sigma), and given the latest time in between collisions
!  (tcol) and the indices of the spheres involved in it (icol,jcol),
!  this subroutine updates the value of the table that stores the times
!  till next collision for each pair of spheres.
    implicit none
    double precision, intent(in) :: pos(:,:)
    double precision, intent(in) :: vel(:,:)
    double precision, intent(in) :: sigma
    double precision, intent(in) :: tcol
    integer, intent(in) :: icol
    integer, intent(in) :: jcol
    double precision, intent(inout) :: ctime(:,:)

    double precision :: bij
    double precision :: cij
    double precision :: discriminant
    integer :: i
    double precision, parameter :: infinity = huge(0.0d0)
    integer :: j
    integer :: jx
    integer :: jy
    integer :: jz
    integer :: nspheres
    double precision :: pos_image(1:3)
    double precision :: time
    double precision :: translat(1:3)
    double precision :: rij(1:3)
    double precision :: uij(1:3)
    double precision :: uij2

    nspheres = size(pos,2)
    do i = 1, nspheres-1
        do j = i+1, nspheres
            ctime(i,j) = ctime(i,j) - tcol
            if ( (i .eq. icol) .or. (i .eq. jcol) .or. (j .eq. icol) .or. (j .eq. jcol) ) then 
                ctime(i,j) = infinity
                uij(1:3) = vel(1:3,i) - vel(1:3,j)
                do jx = -1, 1
                    do jy = -1, 1
                        do jz = -1, 1
                            translat = (/ jx, jy, jz /)
                            pos_image(1:3) = pos(1:3,j) + translat(1:3)
                            rij(1:3) = pos(1:3,i) - pos_image(1:3)
                            bij = dot_product(rij,uij)
                            cij = dot_product(rij,rij) - sigma * sigma
                            if (bij < 0.0) then
                                uij2 = dot_product(uij,uij)
                                discriminant = bij * bij - uij2 * cij
                                if (discriminant > 0.0) then
                                    time = - ( bij + sqrt(discriminant)) / uij2
                                    if (time < ctime(i,j)) then
                                        ctime(i,j) = time
                                    end if
                                end if
                            end if
                        end do
                    end do                
                end do
            end if
        end do    
    end do
end subroutine update_collisions_table