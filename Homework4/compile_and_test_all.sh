# run on cli: sh  ./compile_and_test_all.sh

gfortran -o test_assign_positions assign_positions.f90 test_assign_positions.f90

rm test_assign_positions

./test_assign_positions

gfortran -o test_compute_diameter compute_diameter.f90 test_compute_diameter.f90

./test_compute_diameter

rm test_compute_diameter

gfortran -o test_initialize_collisions_table initialize_collisions_table.f90 test_initialize_collisions_table.f90

./test_initialize_collisions_table

rm test_initialize_collisions_table