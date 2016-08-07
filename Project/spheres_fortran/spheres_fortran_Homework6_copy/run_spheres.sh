gfortran -o spheres_fortran spheres.f90 advance_simulation.f90 assign_positions.f90 assign_velocities.f90 compute_diameter.f90 compute_properties.f90 initialize_collisions_table.f90 read_input.f90 retrieve_collision_info.f90 update_collisions_table.f90 validate_input.f90 write_initial.f90 write_properties.f90

./spheres_fortran


echo "write to screen the reduced volume used and the value of pv/kT computed."
python compute_pressure.py results.txt

echo "produce a graph similar to that of slide ->24<- in Lecture <08> (if everything is correct with your code)."
python graph_order_parameter.py results.txt

echo "produce a graph similar to that of slide ->27<- in Lecture <08> (if everything is correct with your code)."
python graph_velocity_distribution.py results.txt
