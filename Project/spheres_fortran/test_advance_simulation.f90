program test_advance_simulation

   implicit none

   interface
      subroutine initialize_collisions_table(pos,vel,sigma,ctime)
         implicit none
         double precision, intent(in) :: pos(:,:)
         double precision, intent(in) :: vel(:,:)
         double precision, intent(in) :: sigma
         double precision, intent(out) :: ctime(:,:)
      end subroutine initialize_collisions_table
      subroutine retrieve_collision_info(array,minimum,imin,jmin)
         implicit none
         double precision, intent(in) :: array(:,:)
         double precision, intent(out) :: minimum
         integer, intent(out) :: imin
         integer, intent(out) :: jmin
      end subroutine retrieve_collision_info
      subroutine advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
         double precision, intent(in) :: tcol
         integer, intent(in) :: icol
         integer, intent(in) :: jcol
         double precision, intent(in) :: sigma
         double precision, intent(inout) :: pos(:,:)
         double precision, intent(inout) :: vel(:,:)
         double precision, intent(out) :: delta_vel(1:3)
      end subroutine advance_simulation
   end interface

   double precision, allocatable :: ctime(:,:)
   double precision :: delta_vel(1:3)
   double precision :: delta_vel_correct(1:3)
   integer :: icol
   integer :: jcol
   integer :: nspheres
   double precision, allocatable :: pos(:,:)
   double precision, allocatable :: pos_correct(:,:)
   double precision, allocatable :: vel(:,:)
   double precision, allocatable :: vel_correct(:,:)
   double precision :: sigma
   double precision :: tcol
   double precision, parameter :: tol = 1.0d-10

   ! Decide on a particular case with simple known result:
   nspheres = 4
   sigma = sqrt(2.0d0)/10.0d0
   allocate(pos(1:3,1:nspheres))
   pos(1:3,1) = (/ 0.0d0 , 0.0d0 , 0.0d0 /)
   pos(1:3,2) = (/ 0.0d0 , 0.5d0 , 0.5d0 /)
   pos(1:3,3) = (/ 0.5d0 , 0.0d0 , 0.5d0 /)
   pos(1:3,4) = (/ 0.5d0 , 0.5d0 , 0.0d0 /)
   allocate(vel(1:3,1:nspheres))
   vel(1:3,1) = (/ 0.0d0 , 0.0d0 , 0.0d0 /)
   vel(1:3,2) = (/ 0.0d0 , 1.0d0 , 1.0d0 /)
   vel(1:3,3) = (/ 0.0d0 ,-1.0d0 , 1.0d0 /)
   vel(1:3,4) = (/ 0.0d0 , 1.0d0 ,-1.0d0 /)
   allocate(pos_correct(1:3,1:nspheres))
   pos_correct(1:3,1) = (/ 0.0d0 , 0.0d0 , 0.0d0 /)
   pos_correct(1:3,2) = (/ 0.0d0 , 0.7d0 , 0.7d0 /)
   pos_correct(1:3,3) = (/ 0.5d0 , 0.8d0 , 0.7d0 /)
   pos_correct(1:3,4) = (/ 0.5d0 , 0.7d0 , 0.8d0 /)
   allocate(vel_correct(1:3,1:nspheres))
   vel_correct(1:3,1) = (/ 0.0d0 , 0.0d0 , 0.0d0 /)
   vel_correct(1:3,2) = (/ 0.0d0 , 1.0d0 , 1.0d0 /)
   vel_correct(1:3,3) = (/ 0.0d0 , 1.0d0 ,-1.0d0 /)
   vel_correct(1:3,4) = (/ 0.0d0 ,-1.0d0 , 1.0d0 /)
   delta_vel_correct = (/ 0.0d0 , 2.0d0 ,-2.0d0 /) 

   ! Call the subroutine being tested:
   allocate(ctime(1:nspheres,1:nspheres))
   call initialize_collisions_table(pos,vel,sigma,ctime)
   call retrieve_collision_info(ctime,tcol,icol,jcol)
   call advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)

   ! Check if the results are as expected:
   write(*,*) '  Position of sphere 1 (correct):', pos_correct(1:3,1)
   write(*,*) ' Position of sphere 1 (computed):', pos(1:3,1)
   write(*,*) '  Position of sphere 2 (correct):', pos_correct(1:3,2)
   write(*,*) ' Position of sphere 2 (computed):', pos(1:3,2)
   write(*,*) '  Position of sphere 3 (correct):', pos_correct(1:3,3)
   write(*,*) ' Position of sphere 3 (computed):', pos(1:3,3)
   write(*,*) '  Position of sphere 4 (correct):', pos_correct(1:3,4)
   write(*,*) ' Position of sphere 4 (computed):', pos(1:3,4)
   write(*,*) '  Position of sphere 1 (correct):', vel_correct(1:3,1)
   write(*,*) ' Position of sphere 1 (computed):', vel(1:3,1)
   write(*,*) '  Position of sphere 2 (correct):', vel_correct(1:3,2)
   write(*,*) ' Position of sphere 2 (computed):', vel(1:3,2)
   write(*,*) '  Position of sphere 3 (correct):', vel_correct(1:3,3)
   write(*,*) ' Position of sphere 3 (computed):', vel(1:3,3)
   write(*,*) '  Position of sphere 4 (correct):', vel_correct(1:3,4)
   write(*,*) ' Position of sphere 4 (computed):', vel(1:3,4)
   write(*,*) '    Change in velocity (correct):', delta_vel_correct(1:3)
   write(*,*) '   Change in velocity (computed):', delta_vel(1:3)

   deallocate(pos)
   deallocate(vel)
   deallocate(ctime)

end program test_advance_simulation
