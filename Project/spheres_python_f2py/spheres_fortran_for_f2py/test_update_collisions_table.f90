program test_initialize_collisions_table

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
      subroutine update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)
         implicit none
         double precision, intent(in) :: pos(:,:)
         double precision, intent(in) :: vel(:,:)
         double precision, intent(in) :: sigma
         double precision, intent(in) :: tcol
         integer, intent(in) :: icol
         integer, intent(in) :: jcol
         double precision, intent(inout) :: ctime(:,:)
      end subroutine update_collisions_table
   end interface

   double precision, allocatable :: ctime(:,:)
   double precision, allocatable :: ctime_correct(:,:)
   double precision :: delta_vel(1:3)
   integer :: icol
   integer :: jcol
   integer :: nspheres
   double precision, allocatable :: pos(:,:)
   double precision, allocatable :: vel(:,:)
   double precision :: sigma
   double precision :: tcol

   ! Decide on a particular case:
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
   allocate(ctime_correct(1:nspheres,1:nspheres))
   ctime_correct(:,:) = huge(1.0d0)
   ctime_correct(1,2) = 0.2d0
   ctime_correct(3,4) = 0.4d0

   ! Call the subroutine being tested:
   allocate(ctime(1:nspheres,1:nspheres))
   call initialize_collisions_table(pos,vel,sigma,ctime)
   call retrieve_collision_info(ctime,tcol,icol,jcol)
   call advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
   call update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)

   ! Check if the results look correct:
   write(*,*) '  Correct time of collision (1,2)...', ctime_correct(1,2)
   write(*,*) ' Computed time of collision (1,2)...', ctime(1,2)
   write(*,*) '  Correct time of collision (3,4)...', ctime_correct(3,4)
   write(*,*) ' Computed time of collision (3,4)...', ctime(3,4)

   deallocate(pos)
   deallocate(vel)
   deallocate(ctime)
   deallocate(ctime_correct)

end program test_initialize_collisions_table
