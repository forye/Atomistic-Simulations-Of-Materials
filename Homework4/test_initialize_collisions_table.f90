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
   end interface

   integer :: nspheres
   double precision, allocatable :: ctime(:,:)
   double precision, allocatable :: ctime_correct(:,:)
   double precision, allocatable :: pos(:,:)
   double precision, allocatable :: sigma
   double precision :: tol = 1.0d-10
   double precision, allocatable :: vel(:,:)

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
   ctime_correct(1,2) = 0.4d0
   ctime_correct(3,4) = 0.2d0

   ! Call the subroutine being tested:
   allocate(ctime(1:nspheres,1:nspheres))
   call initialize_collisions_table(pos,vel,sigma,ctime)

   ! Check if the results look correct:
   if (abs(ctime(1,2) - ctime_correct(1,2)) < tol) then
      write(*,*) 'Particular test case passed'
      write(*,*) ' Expected time of collision (1,2)...', ctime_correct(1,2)
      write(*,*) ' Computed time of collision (1,2)...', ctime(1,2)
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) ' Expected time of collision (1,2)...', ctime_correct(1,2)
      write(*,*) ' Computed time of collision (1,2)...', ctime(1,2)
   end if
   if (abs(ctime(3,4) - ctime_correct(3,4)) < tol) then
      write(*,*) 'Particular test case passed'
      write(*,*) ' Expected time of collision (3,4)...', ctime_correct(3,4)
      write(*,*) ' Computed time of collision (3,4)...', ctime(3,4)
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) ' Expected time of collision (3,4)...', ctime_correct(3,4)
      write(*,*) ' Computed time of collision (3,4)...', ctime(3,4)
   end if

   deallocate(pos)
   deallocate(vel)
   deallocate(ctime)
   deallocate(ctime_correct)

contains

end program test_initialize_collisions_table
