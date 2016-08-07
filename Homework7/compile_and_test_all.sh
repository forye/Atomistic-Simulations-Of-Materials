gfortran -o  eigen_vals.o eigenvalues.f90 -llapack

./eigen_vals.o

gfortran -o tbxu.o tbxu.f90 -lblas -llapack

./tbxu.o