! Exercise 4: Implementing Subroutine assign_positions()

subroutine assign_positions(pos)
!  Given an array pos, this subroutine returns in it the positions of the
!  sites of a fcc lattice uniformly distributed in a cube of unit volume
!  It computes how many sites there are from the size of array pos.

    implicit none
    double precision, intent(inout) :: pos(:,:)

    double precision :: a
    double precision :: img_dist
    integer :: ix, iy, iz
    integer :: n
    integer :: nspheres
    integer :: sphere_idx
    double precision :: periodic_mirror(1:3)
    
    nspheres = size(pos,2)
    n = nint((nspheres/4)**(1.0/3.0))
    a = 1.0 / dble(n)
    sphere_idx = 0
    img_dist = a / 2.0
    do ix = 1, n
        do iy = 1, n
            do iz = 1, n
                periodic_mirror(1:3) = (/ dble(ix-1) * a, dble(iy - 1) * a, dble(iz - 1) * a /)
                sphere_idx = sphere_idx + 1
                
                pos(1:3, sphere_idx) = periodic_mirror(:)
                sphere_idx = sphere_idx + 1                
                pos(1:3, sphere_idx) = (/ img_dist, img_dist, 0.0d0 /) + periodic_mirror(:)
                sphere_idx = sphere_idx + 1
                pos(1:3, sphere_idx) = (/ img_dist, 0.0d0, img_dist /) + periodic_mirror(:)
                sphere_idx = sphere_idx + 1
                pos(1:3, sphere_idx) = (/ 0.0d0, img_dist, img_dist /) + periodic_mirror(:)
            end do
        end do
    end do
end subroutine assign_positions

