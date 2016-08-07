
supposing elastic collision, between 2 sigma diametered spheres.
cosidering the velocity change, resulting from a colludion

from :
	(dv[i] = -dv[j], where dv is on the vector of velocity canges for all spheres)

show:	
	dv[i] = -dv[j] = - dr[i,j] * b[i,j] / sigma ** 2

where:
	r[i,j] = r[i] - r[j]
	b[i,j] = r[i,j] * (u[i] - u[j])


dataas structures:
r - location of each spher [N_spehres x 3]
V - velocity of each sphere [N_spehres x 3]
dv - velocity change of spehere in collission between t0 and t0+1 
	(in [t0, t0+1] duration there was only 1 collision)


solution:

on collusion where, 2 spehres with same mass m, and temporal speed V [ N x T_collision_time_samples], 
considering linear momentum conservation :

on the direction of the dv verctors which is determained by a line pulled from each center of sphere, 
	in the exact moment of collision, is sparse except the indexes of collision

for example: coll at i=3, j = 5, dv would look like: [0, 0, dv_i, 0, dv_j...]

t0 < t_collision and to + 1 > t_collision

	m * V[i,t0] + m * V[j,t0]  = m * ( V[i,t1] )  + m * V[j,t1]  =>

	 							 m * ( V[j,t0] + dv[i])  + m * (V[j,t0] + dv[j])  =

	 							 m * ( V[j,t0] + dv[i])  + m * (V[j,t0] + dv[j])  => 
	 					 dv[i] = -dv[j]

conservetion of energy in ellastic collision:

	0.5* m * v[i, t0] **2 +
	0.5* m * v[j, t0] ** 2 = 0.5 * m * v[i, t0 + 1] **2 
								+0.5 * m * v[j, t0 + 1] **2 

						= 0.5 * m * (v[i, t0] + dv[i])**2 
						 +0.5 * m * (v[j, t0] + dv[j])**2 =

								 = 0.5 * m * (v[i, t0] + dv[i])**2 
								 +0.5 * m * (v[j, t0] - dv[i])**2

								 => ( v[i, t0] - v[j, t0] ) * dv[i] = - dv[i]**2 =>
								 abs(dv[i]) = - ( v[i, t0] - v[j, t0]) [1]



dv[i] is a 3d vector - it can be represented polarly: 
	
	dv[i] = abs(dv[i]) * r[i,j] / sigma  ( r_hat = r[i,j] / sigma)
	
	
	abs(dv[i]) = -(V[i, t0] - V[j, t0]) * r[i,j] /sigma 

	dv[i] = abs(dv[i]) * r[i,j] / sigma 
		  = ( -(V[i, t0] - V[j, t0]) * r[i,j] /sigma ) * r[i,j] / sigma

				[  b[i,j] = r[i,j] * (u[i] - u[j]) =r[i,j] *  V[i, t0] - V[j, t0]  ]
		  
		   =  -b[i,j] * r[i,j] /sigma**2 *


were done!

Exercise2: Reduced Volume:

	show thatL

		v/v_0 - sqrt(2) * L**3 / (N * sigma**3)

	where
		N= number of spheres
		sigma= diameter of each sphere
		L= edge length of cubic box

		v/v_0= reduced volume - 
		V=L**3,
		v = V / N  

	V_0 is the volume available where the spheres are in dense fcc structure
		V_0=N * v_0



we wil inspect a cube unit cell:
	4 spheres in a unic cell, intercecting inth ecubic diagonals. 
	dense fcc structure results that each unit size a, of the cube is euqal to: 
		a = sqrt(2) * sigma
volume per spehre in dense fcc: 
	v_0 = sigma**3 / (sqrt(2))

	v / v_0 = (L**3 / N) /  (sigma**3 / (sqrt(2)))  

		= sqrt(2) * L**3 / ( N * sigma**3)
were done!



mkdir yad/yad/yada/gedit hello.f90


gfortran -o hello hello.f90
./hello



To create the executable:
$ gfortran -o spheres spheres.f90
To run the executable, reading from file input.txt:
$ ./spheres < input.txt
which writes to the screen; if we want to write to a file output.txt:
●
●
$ ./spheres < input.txt > output.txt


sudo apt-get install -y build-essential git python-pip libfreetype6-dev libxft-dev libncurses-dev libopenblas-dev gfortran python-matplotlib libblas-dev liblapack-dev libatlas-base-dev python-dev python-pydot linux-headers-generic linux-image-extra-virtual unzip python-numpy swig python-pandas python-sklearn unzip wget pkg-config zip g++ zlib1g-dev
