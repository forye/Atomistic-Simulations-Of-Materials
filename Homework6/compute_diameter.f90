! Exercise 3: Implementing Subroutine compute_diameter()


subroutine compute_diameter(nspheres,rvolume,sigma)
!  Given a number of spheres (nspheres) occupying a reduced volume (rvolume)
!  in a cube of unit volume, this subroutine returns the value of the diameter
!  of the spheres (sigma).
   implicit none
   integer, intent(in) :: nspheres
   double precision, intent(in) :: rvolume
   double precision, intent(out) :: sigma   
   sigma = (sqrt(2.0d0)/(nspheres*rvolume))**(1.0d0/3.0d0)

end subroutine compute_diameter