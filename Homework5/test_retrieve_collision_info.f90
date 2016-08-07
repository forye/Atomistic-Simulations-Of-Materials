program test_retrieve_collision_info

   implicit none

   interface
      subroutine retrieve_collision_info(array,minimum,imin,jmin)
         implicit none
         double precision, intent(in) :: array(:,:)
         double precision, intent(out) :: minimum
         integer, intent(out) :: imin
         integer, intent(out) :: jmin
      end subroutine retrieve_collision_info
   end interface

   double precision, allocatable :: array(:,:)
   integer :: i
   integer :: imin
   integer :: imin_correct
   integer :: j
   integer :: jmin
   integer :: jmin_correct
   double precision :: minimum
   double precision :: minimum_correct
   integer :: narray

   ! Decide on a particular case:
   narray = 1000
   allocate(array(1:narray,1:narray))
   array = huge(1.0d0)
   do i = 1 , narray-1
      do j = i+1 , narray
         array(i,j) = dble(i*j)
      end do
   end do
   minimum_correct = 2.0d0
   imin_correct = 1
   jmin_correct = 2

   ! Call the subroutine being tested:
   call retrieve_collision_info(array,minimum,imin,jmin)

   ! Check if the results look correct:
   write(*,*) 'Expected mimimum element value:', minimum_correct
   write(*,*) 'Computed mimimum element value:', minimum
   write(*,*) 'Expected mimimum element indices:', imin_correct, jmin_correct
   write(*,*) 'Computed mimimum element indices:', imin, jmin

   deallocate(array)

contains

end program test_retrieve_collision_info
