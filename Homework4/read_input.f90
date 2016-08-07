! Exercise 1: Implementing Subroutine read input()

subroutine read_input(nspheres,rvolume,ncollisions)
!  This subroutine reads three values from standard input, prints them to
!  standard output, and passes them as arguments to the main program.
   implicit none
   integer, intent(out) :: nspheres
   double precision, intent(out) :: rvolume
   integer, intent(out) :: ncollisions
   write(*,*) 'Choose number of spheres:  Ns = 4*(n**3), n=1, 2, 3...'
   read(*,*) nspheres
   write(*,*) 'Enter volumes ratio: v / v0 >= 1'
   read(*,*) rvolume
   write(*,*) 'Enter the number of collisions to simulate'
   read(*,*) ncollisions
   write(*,*) ' Number of spheres:', nspheres
   write(*,*) ' Reduced volume:', rvolume
   write(*,*) ' Number of collisions:', ncollisions
      
end subroutine read_input
