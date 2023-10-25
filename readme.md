# Fortran vs Numpy for classical MD
### Code in Different Language Project comparing Fortran implementation and NumPy implementation of classical Molecular dynamics code. 

This project provides a Python alternative to the Fortran code, simulating the colliding spheres model with periodic boundary conditions. The model is described in [link_to_article](#link_to_article).

The NumPy implementation achieves the same order of magnitude in speed as Fortran even though its interpreted and not compiled

## Introduction

Python is an object-oriented interpreted language. It serves as a generic wrapper for a plethora of modules and packages (some written in C and Fortran) for various purposes and disciplines. While Python is highly versatile, making it apt for web, cloud computing, string manipulations, and machine learning, it can be slow and inefficient if misused. However, when used correctly, Python can be an incredible tool.

The crux of numerical computations in Python lies in the **Numpy** package, which offers `numpy.array` (and `f2py`). With its roots tied to BLAS and LAPCK and being implemented in C (and occasionally Fortran), Numpy is optimized to minimize iteration overhead, resulting in rapid calculations.

To exploit Python's numerical computation capabilities, iterative loop operations should be rewritten using numpy arrays. This will involve data structures, array indexing, and masked array operations to concentrate vector operations where required.

The intuitive expectation is that Numpy in Python might run slower than pure Fortran due to inherent overheads. However, for larger data structures, the performance difference should shrink, making Python a more favorable choice.

## Step 1: Preparing the Code

To compare Fortran and Python, five subroutines were selected. These subroutines are repeated with each collision and have been implemented in Numpy in a Python script named `spheres`.

Additional functionalities include:
- Testing if spheres collide
- Ensuring positive collision time
- Checking for a symmetric collision time matrix
- A configuration file

There are different versions of the code:

1. **Object-oriented version**
2. **Aligned with the Fortran implementation version**


Python functions signitures:
```python
    def spheres.retrieve_collision_info(c_time):        
        returns tcol, jcol, icol
    def spheres.advance_simulation(pos, vel, tcol, icol, jcol, sigma_sq, sigma, Ns, nc,pos_trans[0], conf, L)        
        returns pos, vel, delta_vel, mom, kin 
    def spheres.compute_properties(vel)
        returns mom, kin
    def utilities.io_utils.write_properties(nc, mom, kin, tcol, icol, jcol, delta_vel)
        returns None
    def spheres.update_collisions_table(pos, vel, sigma, tcol, icol, jcol, c_time, pos_trans)
        returns c_time
```      
* original code at spheres_python folder

in fortran, I've canceled the UI and the testing the main function
```fortran

Fortran subroutiens to be tested:
        call retrieve_collision_info(ctime,tcol,icol,jcol)
        call advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
        call compute_properties(vel,mom,kin)
        call write_properties(c,mom,kin,tcol,icol,jcol,delta_vel)
        call update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)
```
* original code at spheres_fortran/spheres_fortran_Homework6_copy folder 
        
## step 2 - time measurement

The Python implementation was first tested using the provided script from the class and
then modified for benchmarks:

1. execute
expirement parameteres: 300,000 times to increse meausurement in the tendency for equalibirium
    4 cases:

        1. rv = 1.2, n=2 ( number of spheres = 32)
        python: 439.53 secs fortran: 50.58 secs
        
        2. rv = 1.8, n=2 
        python: 422.47 secs fortran: 58.44 secs
        
        3. rv = 1.2, n=5 ( number of spheres = 500)
        python: 1393.77 secs fortran: 1455.4049540000001 secs
        
        4. rv = 1.8, n=5
        python: 2265.97 secs fortran: 1436.77 secs
        
2. execute
test for 30,000 collisions
rv = 1.2, 1.8, 3.6
Ns = 4, 32, 108, 256, 500
fortran:{
    1.8: {32: 3.913207, 108: 9.610123000000002, 4: 2.249561, 500: 71.317622, 256: 27.273652}, 
    1.2: {32: 4.012453, 108: 9.551079999999999, 4: 2.378793, 500: 72.441602, 256: 28.153621}, 
    3.6: {32: 3.945790, 108: 9.659133999999995, 4: 2.369048, 500: 71.455186, 256: 27.323331}}
python:{    
    1.2: {32: 32.6938040257, 108: 55.231487035, 4: 23.444005012, 500: 103.05965685, 256: 54.6719648838},
    1.8: {32: 32.3645670414, 108: 54.926554203, 4: 23.145853042, 500: 101.40656399, 256: 54.8189430237}, 
    3.6: {32: 32.4476668835, 108: 55.170861959, 4: 22.943297147, 500: 104.00603389, 256: 54.5032930374}} 
    
made by:
    read_results.py

watch :
    spheres_compare.png

from the graph, it is seen that both implemtations have an order of magnitude difference in a small number of sphres, but converges in large scale simulations (high number of spheres), with a pseudo constant overhead and regardles of the volumes ratio


possible future work:
    add 2nd diameter (d2, sigma2) of the colliding spheres (realizing the potential function)
    fix the visualization tool
    realize the entire program in f2py and benchmark each function
    repeat process for bcc
    different radiuses
    






out of scope work:    

-------------------------------------------------------------------------------------------------

step 3 3D visualization


the code is in visual_spheres module


note: for running with this tool, you require Anaconda distrubution

sudo ~/anaconda/bin/conda clean --lock
 6124  sudo ~/anaconda/bin/conda install -c https://conda.binstar.org/mwcraig vpython


~/anaconda/bin/pythonw personal/hagasha/Project/spheres_python/visualize.py


...

---------------------------------------------------------------------------------------------------------------------
step 4 f2py realization and testing 5x4 tests in hybrid f2py vs numpy
this section is not finished, but its out of scope

first test simple function, compute diamter:
    lets run
>>> cd spheres_fortran_for_f2py    
>>>f2py -c -m compute_diameter_f2py compute_diameter.f90
>>>python
>>>import test_f2py
>>> dir(test_f2py)
['__doc__', '__file__', '__name__', '__package__', '__version__', 'compute_diameter']
>>> print test_f2py.compute_diameter(32, 1.2)
0.332706478681
>>> print test_f2py.compute_diameter(500, 1.2)
0.133082591473
    
    looks good
    
lets compile the 5 functions to a spheres_fortran module:
>>>f2py -c -m spheres_fortran retrieve_collision_info.f90 advance_simulation.f90 compute_properties.f90 write_properties.f90 update_collisions_table.f90

we've got a warning:
    In: :spheres_fortran:advance_simulation.f90:advance_simulation
    get_parameters: got "invalid syntax (<string>, line 1)" on 'huge(0.0d0)'

so I changed huge(0.0d0) to be 9999999.9

lets check the interface:
>>>f2py -h stdout -m spheres_fortran retrieve_collision_info.f90 advance_simulation.f90 compute_properties.f90 write_properties.f90 update_collisions_table.f90

lets check the module:
@python>>> dir(spheres_fortran)
['__doc__', '__file__', '__name__', '__package__', '__version__', 'advance_simulation', 'compute_properties', 'retrieve_collision_info', 'update_collisions_table', 'write_properties']

use of the f2py code is in sphres_python_f2py
