# run on cli: sh  ./compile_and_test_all.sh

# run on cli: sh  ./compile_and_test_all.sh
#
#gfortran -o test_assign_positions assign_positions.f90 test_assign_positions.f90
#
#rm test_assign_positions
#
#./test_assign_positions
#
#gfortran -o test_compute_diameter compute_diameter.f90 test_compute_diameter.f90
#
#./test_compute_diameter
#
#rm test_compute_diameter
#
#gfortran -o test_initialize_collisions_table initialize_collisions_table.f90 test_initialize_collisions_table.f90
#
#./test_initialize_collisions_table
#
#rm test_initialize_collisions_table
#
#
#
## test assign_velocities
#gfortran -o test_assign_velocities assign_velocities.f90 test_assign_velocities.f90
#
#./assign_velocities
#
#rm assign_velocities
#
#
#
## run on cli: sh  ./compile_and_test_all.sh
#
#gfortran -o test_advance_simulation advance_simulation.f90 test_advance_simulation.f90 initialize_collisions_table.f90 retrieve_collision_info.f90
#
#./test_advance_simulation
#
#rm test_advance_simulation
#
#gfortran -o test_retrieve_collision_info retrieve_collision_info.f90 test_retrieve_collision_info.f90
#
#./test_retrieve_collision_info
#
#rm test_retrieve_collision_info
#
#gfortran -o test_update_collisions_table update_collisions_table.f90 test_update_collisions_table.f90 initialize_collisions_table.f90 retrieve_collision_info.f90 advance_simulation.f90
#
#./test_update_collisions_table
#
#rm test_update_collisions_table
#
#
#
#
#gfortran -o test_assign_velocities assign_velocities.f90 test_assign_velocities.f90 
#
#./test_assign_velocities
#
#rm test_assign_velocities
#
#
#


gfortran -o spheres_fortran spheres.f90 advance_simulation.f90 assign_positions.f90 assign_velocities.f90 compute_diameter.f90 compute_properties.f90 initialize_collisions_table.f90 read_input.f90 retrieve_collision_info.f90 update_collisions_table.f90 validate_input.f90 write_initial.f90 write_properties.f90

./spheres_fortran

