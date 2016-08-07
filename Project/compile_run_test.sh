
# procedure for testing execution time of spheres_fortran

# for each scenario, edit the file spheres_test.f90 and spheres_test.py  updating the changes in the variables accroding to this scenario (Itterations, rvolume, Ncollisions, number of spheres (4n**3)) and run

gfortran -o test_sphere spheres_test.f90 advance_simulation.f90 assign_positions.f90 assign_velocities.f90 compute_diameter.f90 compute_properties.f90 initialize_collisions_table.f90 read_input.f90 retrieve_collision_info.f90 update_collisions_table.f90 validate_input.f90 write_initial.f90 write_properties.f90

./spheres_fortran/test_sphere

./spheres_python/python spheres_test.py

