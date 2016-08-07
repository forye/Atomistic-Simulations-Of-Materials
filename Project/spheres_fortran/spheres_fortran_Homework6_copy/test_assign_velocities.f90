program test_assign_velocities

   implicit none

   interface
      subroutine assign_velocities(vel)
         implicit none
         double precision, intent(inout) :: vel(:,:)
      end subroutine assign_velocities
   end interface

   double precision :: height
   integer :: npoints_1
   integer :: npoints_2
   integer :: npoints_3
   integer :: npoints_4
   integer :: npoints_5
   integer :: npoints_6
   integer :: npoints_correct
   integer :: nspheres
   double precision, allocatable :: vel(:,:)

   ! Decide on a particular case:
   nspheres = 4000000
   height = 0.1

   ! Call the subroutine being tested:
   allocate(vel(1:3,1:nspheres))
   call assign_velocities(vel)

   ! To check if the results look correct, we count points generated in a
   ! spherical cap located around (100), (010), (001) and their opposites:
   npoints_correct = int(dble(nspheres)*height/2.0d0)
   npoints_1 = counted_points(vel, 1, 1, height)
   npoints_2 = counted_points(vel, 2, 1, height)
   npoints_3 = counted_points(vel, 3, 1, height)
   npoints_4 = counted_points(vel, 1, -1, height)
   npoints_5 = counted_points(vel, 2, -1, height)
   npoints_6 = counted_points(vel, 3, -1, height)
   write(*,*) ' Expected number of points per cap...', npoints_correct
   write(*,*) '    Number of points per cap (x) ...', npoints_1
   write(*,*) '    Number of points per cap (y) ...', npoints_2
   write(*,*) '    Number of points per cap (z) ...', npoints_3
   write(*,*) '    Number of points per cap (-x)...', npoints_4
   write(*,*) '    Number of points per cap (-y)...', npoints_5
   write(*,*) '    Number of points per cap (-z)...', npoints_6

   deallocate(vel)

contains

   integer function counted_points(points, idirection, idirsign, height)
   
      double precision, intent(in) :: points(:,:)
      integer, intent(in) :: idirection
      integer, intent(in) :: idirsign
      double precision, intent(in) :: height

      integer :: i
      double precision :: magnitude
      integer :: npoints

      npoints = size(points,2)
      counted_points = 0
      do i=1, npoints
         magnitude = points(idirection,i)/norm2(points(:,i))
         if ( idirsign*magnitude > (1.0d0 - height) ) then
            counted_points = counted_points + 1
         end if
      end do

   end function counted_points

end program test_assign_velocities
