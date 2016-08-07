program spheres

implicit none

    interface
        subroutine spheres(nspheres,rvolume,ncollisions)
            implicit none
            integer, intent(in) :: nspheres
            double precision, intent(in) :: rvolume
            integer, intent(in) :: ncollisions
    end subroutine spheres
    
    integer :: Iterations = 10
    integer :: i    
    integer :: nspheres,ncollisions
    
    double precision :: rvolume
    double precision :: time_diff = 0.0
    double precision :: t_stop, t_start
    
    nspheres = 500
    rvolume = 1.2
    ncollisions = 30000
    
    call cpu_time(t_start)
    do i = 1, Iterations
        call spheres(nspheres,rvolume,ncollisions)
    end do
    call cpu_time(t_stop)
    
    print *, Iterations, "Iterations running time:", t_stop- t_start, " [seconds]"

end program spheres