program test_assign_positions

   implicit none

   interface
      subroutine assign_positions(pos)
         implicit none
         double precision, intent(inout) :: pos(:,:)
      end subroutine assign_positions
   end interface

   double precision :: angle
   double precision :: angle_correct
   integer :: i
   integer :: j
   integer :: k
   integer :: l
   double precision :: length1
   double precision :: length1_correct
   double precision :: length2
   double precision :: length2_correct
   integer :: nspheres
   double precision, allocatable :: pos(:,:)
   double precision, parameter :: tol = 1.0d-10

   ! Decide on a particular case with simple known result:
   nspheres = 32
   length1_correct = 0.353553390593274d0
   length2_correct = 0.5d0
   angle_correct = 60.0d0

   ! Call the subroutine being tested:
   allocate(pos(1:3,1:nspheres))
   call assign_positions(pos)

   ! Check if the results are as expected:
   i = 17
   j = 18
   k = 21
   l = 19
   length1 = computed_distance(pos(1:3,i), pos(1:3,j))
   length2 = computed_distance(pos(1:3,i), pos(1:3,k))
   angle = computed_angle(pos(1:3,j), pos(1:3,i), pos(1:3,l))
   if (abs(length1 - length1_correct) < tol) then
      write(*,*) 'Particular test case passed'
      write(*,*) ' Expected length...', length1_correct
      write(*,*) ' Computed length...', length1
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) ' Expected length...', length1_correct
      write(*,*) ' Computed length...', length1
      stop('Exiting...')
   end if
   if (abs(length2 - length2_correct) < tol) then
      write(*,*) 'Particular test case passed'
      write(*,*) ' Expected length...', length2_correct
      write(*,*) ' Computed length...', length2
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) ' Expected length...', length2_correct
      write(*,*) ' Computed length...', length2
      stop('Exiting...')
   end if
   if (abs(angle - angle_correct) < tol) then
      write(*,*) 'Particular test case passed'
      write(*,*) ' Expected angle...', angle_correct
      write(*,*) ' Computed angle...', angle
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) ' Expected angle...', angle_correct
      write(*,*) ' Computed angle...', angle
      stop('Exiting...')
   end if

   deallocate(pos)

contains

   double precision function computed_distance(p1, p2)

      implicit none

      double precision, intent(in) :: p1(1:3)
      double precision, intent(in) :: p2(1:3)

      computed_distance = norm2(p1-p2)

   end function computed_distance

   double precision function computed_angle(p1, p2, p3)

      implicit none

      double precision, intent(in) :: p1(1:3)
      double precision, intent(in) :: p2(1:3)
      double precision, intent(in) :: p3(1:3)

      double precision :: v1(1:3)
      double precision :: v2(1:3)

      v1 = p1 - p2
      v2 = p3 - p2
      computed_angle = acos(dot_product(v1,v2) / (norm2(v1)*norm2(v2)))
      computed_angle = computed_angle/(4.0d0*atan(1.0d0)/180.0d0)

   end function computed_angle

end program test_assign_positions
