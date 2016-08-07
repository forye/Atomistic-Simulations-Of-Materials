! Exercise 1: Implementing Subroutine retrieve_collision_info
! test with test_retrieve_collision_info.f90

subroutine retrieve_collision_info(array,minimum,imin,jmin)
!  This subroutine returns the minimum value (minimum) of an squared array
!  (array) and the elements (imin,jmin) for which it occurs. It only searches
!  in the part of the array where the first index is lower than the second
!  index.
    implicit none
    double precision, intent(in) :: array(:,:)
    double precision, intent(out) :: minimum
    integer, intent(out) :: imin
    integer, intent(out) :: jmin

    integer :: i
    double precision :: infinity = 9999999.0! huge(0.0d0)
    integer :: j
    integer :: n

    n = size(array,1)
    minimum = infinity
    do i = 1, n-1
        do j = i+1, n
            if (array(i,j) < minimum) then
                minimum = array(i,j)
                imin = i
                jmin = j
            end if
        end do
    end do
end subroutine retrieve_collision_info