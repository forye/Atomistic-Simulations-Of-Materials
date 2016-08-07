! Exercise 2: Complete Simulations

!gfortran -o spheres_test spheres_test.f90 advance_simulation.f90 assign_positions.f90 assign_velocities.f90 compute_diameter.f90 compute_properties.f90 initialize_collisions_table.f90 read_input.f90 retrieve_collision_info.f90 update_collisions_table.f90 validate_input.f90 write_initial.f90 write_properties.f90
!./spheres_test

program spheres

    implicit none

    interface
        subroutine read_input(nspheres,rvolume,ncollisions)
            implicit none
            integer, intent(out) :: nspheres
            double precision, intent(out) :: rvolume
            integer, intent(out) :: ncollisions
        end subroutine read_input

        subroutine validate_input(nspheres,rvolume)
            implicit none
            integer, intent(in) :: nspheres
            double precision, intent(in) :: rvolume
        end subroutine validate_input

        subroutine compute_diameter(nspheres,rvolume,sigma)
            implicit none
            integer, intent(in) :: nspheres
            double precision, intent(in) :: rvolume
            double precision, intent(out) :: sigma
        end subroutine compute_diameter

        subroutine assign_positions(pos)
            implicit none
            double precision, intent(inout) :: pos(:,:)
        end subroutine assign_positions

        subroutine assign_velocities(vel)
            implicit none
            double precision, intent(inout) :: vel(:,:)
        end subroutine assign_velocities

        subroutine write_initial(nspheres,rvolume,ncollisions,pos,vel)
            implicit none
            integer, intent(in) :: nspheres
            double precision, intent(in) :: rvolume
            integer, intent(in) :: ncollisions
            double precision, intent(in) :: pos(:,:)
            double precision, intent(in) :: vel(:,:)
        end subroutine write_initial

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

        subroutine compute_properties(vel,mom,kin)
            implicit none
            double precision, intent(in) :: vel(:,:)
            double precision, intent(out) :: mom(1:3)
            double precision, intent(out) :: kin
        end subroutine compute_properties

        subroutine write_properties(c,mom,kin,tcol,icol,jcol,delta_vel)
            implicit none
            integer, intent(in) :: c
            double precision, intent(in) :: mom(1:3)
            double precision, intent(in) :: kin
            double precision, intent(in) :: tcol
            integer, intent(in) :: icol
            integer, intent(in) :: jcol
            double precision, intent(in) :: delta_vel(1:3)
        end subroutine write_properties

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


    integer :: c
    double precision, allocatable :: ctime(:,:)
    double precision :: delta_vel(1:3)
    integer :: icol
    integer :: jcol
    double precision :: kin
    double precision :: mom(1:3)
    integer :: nspheres
    integer :: ncollisions
    double precision, allocatable :: pos(:,:)
    double precision :: rvolume
    double precision :: sigma
    double precision :: tcol
    double precision, allocatable :: vel(:,:)
    integer :: rv
    integer :: ns
    double precision :: t_stop
    double precision :: t_start
    
    double precision :: rvolumes(1:3)
    double precision :: nsphereses(1:5)
    
    
    rvolumes(1) = 1.2
    rvolumes(2) = 1.8
    rvolumes(3) = 3.6
    nsphereses(1) = 4
    nsphereses(2) = 32
    nsphereses(3) = 108
    nsphereses(4) = 256
    nsphereses(5) = 500
    ncollisions = 30001    
    write(*,*) 'Running spheres test on ', ncollisions, 'collisions'
    write(*,*) 'volume ratios: ', rvolumes, '.  numbers of spheres: ', nsphereses
    
    do rv = 1,3
        rvolume = rvolumes(rv)
        do ns = 1,5 
            nspheres = nsphereses(ns)     
!    call read_input(nspheres,rvolume,ncollisions)
!    write(*,*) 'ENTERING valiudate input'
!    call validate_input(nspheres,rvolume)
!    write(*,*) 'ENTERING compute diameter'
            call compute_diameter(nspheres,rvolume,sigma)
            allocate(pos(1:3,1:nspheres))
!            write(*,*) 'alocated pos'
            call assign_positions(pos)
!            write(*,*) 'assigned pos'
            allocate(vel(1:3,1:nspheres))
!            write(*,*) 'alocated vel'
            call assign_velocities(vel)
!            write(*,*) 'assigned vel, writing initial'
            call write_initial(nspheres,rvolume,ncollisions,pos,vel)
!            write(*,*) 'initial written'
            allocate(ctime(1:nspheres,1:nspheres))
            call initialize_collisions_table(pos,vel,sigma,ctime)

            call cpu_time(t_start)
            
            do c = 1, ncollisions
                call retrieve_collision_info(ctime,tcol,icol,jcol)
                call advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
                call compute_properties(vel,mom,kin)
                call write_properties(c,mom,kin,tcol,icol,jcol,delta_vel)
                call update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)
            end do

            call cpu_time(t_stop)    

            deallocate(pos)
            deallocate(vel)
            deallocate(ctime)
            write(*,*) 'done with nspheres:',nspheres,' spheres. volumes ratio: ',rvolume
            write(*,*) 'ncollisions:', ncollisions, ' took: ', t_stop-t_start, 'seconds'
        end do
    end do

end program spheres



