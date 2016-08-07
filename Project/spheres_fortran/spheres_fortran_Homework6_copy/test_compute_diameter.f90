program test_compute_diameter

   implicit none

   interface
      subroutine compute_diameter(nspheres,rvolume,sigma)
         implicit none
         integer, intent(in) :: nspheres
         double precision, intent(in) :: rvolume
         double precision, intent(out) :: sigma
      end subroutine compute_diameter
   end interface

   integer :: nspheres
   integer :: nspheres_correct
   double precision :: rvolume
   double precision :: sigma
   double precision :: sigma_correct
   double precision, parameter :: tol = 1.0d-10

   ! Test a particular case with simple known result:
   nspheres = 4
   rvolume = 1
   sigma_correct = 1.0d0/sqrt(2.0d0)
   call compute_diameter(nspheres,rvolume,sigma)
   if (abs(sigma - sigma_correct) < tol) then
      write(*,*) 'Particular test case passed'
   else
      write(*,*) 'ERROR: Particular test case NOT passed'
      write(*,*) '         sigma =', sigma
      write(*,*) ' sigma_correct =', sigma_correct
      stop('Exiting...')
   end if

end program test_compute_diameter
