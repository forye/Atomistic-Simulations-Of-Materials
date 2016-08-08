Code in diffrent language project:

This paper is to provied a python alternative to the Fortran code, simulating the coliding spheres model with periodinc boundry conditions, described in #link_to_article

Python is an object oriented interpreted language. It can be seen as a generic wrapper for a vast selection of modules and packages (some written in C and Fortran) for various porpeouses and dicipilnes. It's good for Web and clowd computing, great with strings and files, most common language for Machine learning,  This can make it very slow and ineffician if miss-used or a great tool if used well. 

The main package in the python is Numpy, a computational package that implements numpy.array ( and f2py...) leaned to BLAS and LAPCK with C (and somtimes fortran) implementation. By the use of vectors opperation, Numpy is optimized to skip the overhead necsary for itterations, and can result really fast calculations. 

In order to utilize python best numeric calculation efficiency, the loops itterative operations need to be rewriten as vector opperation using numpy arrays data structures, array indexing and masked arrays operations will be used to focus vector operations to only where needed

I inutivly expect Numpy Python should run slower due to the overhead it contains in compairisment to pure Fortran, but in a magnitude order, and in large data structures operations, performance should converge. This should make Python faverable due to additional options it provides, the ease of writing and the fact that its an object oriented language.


Step 1 - preparing the code

In order to comapre fortran and Python implementation, I've selected the 5 subroutines that are repeated each collision and implemented them in Numpy in a python script named spheres. 

The Python scrip imitates the execution of the fortran code with use of numpy functions. In addition additional tests has been added, for example: test if the spheres are coliding with each other, if the collision time is positive, if the collision time matrix is simetric (instead of triangular matrix, due to realization threw vector operations)
also a configuration file has been added.

I've also realized of the python script as a python obejct oriented module.
Then I have modified the implementations to 2 equivavlent version of the same algorithm, one in Python-Numpy and the other in Fortran and time each. The program allows to itterate several time on the same exectution.

Ive created several versions of the code
    
    1st :Is the Object oriented version
    
    2nd: Is an aligned to the fortran implemetation version


Python functions signitures:

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
        
* original code at spheres_python folder

in fortran, I've canceled the UI and the testing the main function

Fortran subroutiens to be tested:
        call retrieve_collision_info(ctime,tcol,icol,jcol)
        call advance_simulation(tcol,icol,jcol,sigma,pos,vel,delta_vel)
        call compute_properties(vel,mom,kin)
        call write_properties(c,mom,kin,tcol,icol,jcol,delta_vel)
        call update_collisions_table(pos,vel,sigma,tcol,icol,jcol,ctime)

* original code at spheres_fortran/spheres_fortran_Homework6_copy folder 
        
step 2 - time measurement

I've tested the execution of the python implementation with the python script given in class

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
