! Exercise 6: Implementing Subroutine initialize_collisions_table()

subroutine initialize_collisions_table(pos,vel,sigma,ctime)
!  Given the positions (pos) and velocities (vel) of a number of spheres
!  (deduced from the size of these arrays) of given diameter (sigma)
!  enclosed in a cube of unit volume, this subroutine computes
!  how long it takes for every pair of spheres to collide, and returns these
!  times in a two-dimensional array (ctime). The periodic boundary conditions
!  of the box are taken into account.
    implicit none
    double precision, intent(in) :: pos(:,:)
    double precision, intent(in) :: vel(:,:)
    double precision, intent(in) :: sigma
    double precision, intent(out) :: ctime(:,:)

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
    double precision :: rij(1:3)
    double precision :: time
    double precision :: periodic_mirror(1:3)
    double precision :: uij(1:3)
    double precision :: uij2
    
    ! claculate the minimal colusion time, next collision, for each reflection of all spheres
    nspheres = size(pos,2)
    do i = 1, nspheres-1
        do j = i+1, nspheres
            ctime(i,j) = infinity
            uij(1:3) = vel(1:3,i) - vel(1:3,j)
            do jx = -1, 1
                do jy = -1, 1
                    do jz = -1, 1
                        periodic_mirror = (/ jx, jy, jz /)
                        pos_image(1:3) = pos(1:3,j) + periodic_mirror(1:3)
                        rij(1:3) = pos(1:3,i) - pos_image(1:3)
                        bij = dot_product(rij,uij)
                        cij = dot_product(rij,rij)-sigma*sigma
                        if (bij < 0.0) then
                            uij2 = dot_product(uij, uij)
                            discriminant = bij * bij - uij2 * cij
                            if (discriminant > 0.0) then
                                !init with next collision
                                time = (-bij - sqrt(discriminant)) / uij2
                                if (time < ctime(i,j)) then
                                    ctime(i,j) = time
                                end if                            
                            end if
                        end if                
                    end do
                end do
            end do
        end do
    end do
    write(*,*) 'ctime dimensions:', size(ctime,1), size(ctime,2)
end subroutine initialize_collisions_table