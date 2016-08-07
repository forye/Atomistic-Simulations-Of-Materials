! Exercise 5: Implementing Subroutine write_initial()

subroutine write_initial(nspheres,rvolume,ncollisions,pos,vel)
!  This subroutine writes the input given and the initial positions and
!  velocities of a set of spheres to an output file.    
    implicit none
    integer, intent(in) :: nspheres
    double precision, intent(in) :: rvolume
    integer, intent(in) :: ncollisions
    double precision, intent(in) :: pos(:,:)
    double precision, intent(in) :: vel(:,:)
    character(11), parameter :: filename = "results.txt"
    integer :: i
    
    open(unit=1, file=filename, action="write", status="replace")    
!    write(*,*) 'number of spheres: Ns=', nspheres
    write(1,*) nspheres
    write(1,*) rvolume
    write(1,*) ncollisions
    do i = 1, nspheres        
!        write(*,*) i, pos(1,i), pos(2,i), pos(3,i), vel(1,i), vel(2,i), vel(3,i)
        write(1,*) i, pos(:,i), vel(:,i)
        ! write(1,*) pos(2,i)
        ! write(1,*) pos(3,i)
        !write(1,*) vel(1,i)
        !write(1,*) vel(2,i)
        !write(1,*) vel(3,i)
    end do
    close(1)
    
!    write(*,*) 'Finished writing lines: #', nspheres
!    write(*,*) 'Results written to: ', filename
    
end subroutine write_initial