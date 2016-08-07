! Exercise 2: Implementing Subroutine validate_input()

subroutine validate_input(nspheres,rvolume)
!  This subroutine makes sure that a number of spheres is suitable for
!  packing them in a fcc arrangement inside a cube, and that a reduced volume
!  is larger than 1.
   implicit none
   integer, intent(in) :: nspheres
   double precision, intent(in) :: rvolume
   
   integer :: n
   n = NINT((nspheres/4)**(1.0/3.0))
   if ((4*n*n*n /= nspheres) .or. (n <= 0)) then
        stop ('Bad Number of spheres (Ns)')
   end if
   
   if (rvolume < 1.0) then
        stop ('Bad wrong volumes ratio (rv=v/v_0)')
   end if
   
end subroutine validate_input
