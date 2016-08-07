# run on cli: sh  ./compile_and_test_all.sh

gfortran -o test_advance_simulation advance_simulation.f90 test_advance_simulation.f90 initialize_collisions_table.f90 retrieve_collision_info.f90

./test_advance_simulation

rm test_advance_simulation

gfortran -o test_retrieve_collision_info retrieve_collision_info.f90 test_retrieve_collision_info.f90

./test_retrieve_collision_info

rm test_retrieve_collision_info

gfortran -o test_update_collisions_table update_collisions_table.f90 test_update_collisions_table.f90 initialize_collisions_table.f90 retrieve_collision_info.f90 advance_simulation.f90

./test_update_collisions_table

rm test_update_collisions_table
