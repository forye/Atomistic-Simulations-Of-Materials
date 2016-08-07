import math
import matplotlib.pyplot as plt
import numpy as np
import sys


def read_file_lines(filename):
   lines = []
   try:
      file = open(filename)
      while 1:
         line = file.readline()
         if not line:
            break
         lines.append(line)
      print 'Read lines from file...', filename
      return lines
   except:
      sys.exit( 'Problem reading file. EXITING.' )

def extract_ncollisions(lines):
   ncollisions = int(lines[2].split()[0])
   return ncollisions

def extract_initial_positions(lines):
   npos = int(lines[0].split()[0])
   pos = []
   for i in range(npos):
      line_index = i + 3
      pos.append(np.array(map(float,lines[line_index].split()[1:4])))
   return pos

def extract_initial_velocities(lines):
   nvel = int(lines[0].split()[0])
   vel = []
   for i in range(nvel):
      line_index = i + 3
      vel.append(np.array(map(float,lines[line_index].split()[4:7])))
   return vel

def compute_fcc_latpar(npoints):
   box_side = 1.0
   latpar = box_side / float((npoints/4)**(1.0/3.0))
   return latpar

def extract_time(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 3
   ctime = float(lines[line_index].split()[0])
   return ctime

def update_positions(pos, time, vel):
   for i in range(len(pos)):
      pos[i] = pos[i] + time*vel[i]

def compute_order_parameter(pos,latpar):
   npos = len(pos)
   gamma_x = 0.0
   gamma_y = 0.0
   gamma_z = 0.0
   for i in range(npos):
      (x, y, z) = pos[i]
      gamma_x = gamma_x + math.cos(4.0*math.pi*x/latpar)
      gamma_y = gamma_y + math.cos(4.0*math.pi*y/latpar)
      gamma_z = gamma_z + math.cos(4.0*math.pi*z/latpar)
   gamma = (1.0/float(npos))*(gamma_x+gamma_y+gamma_z)/3.0
   return gamma

def extract_indices(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 4
   i, j = map(int, lines[line_index].split())
   return i, j

def extract_delta_vel(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 5
   delta_vel = np.array(map(float, lines[line_index].split()))
   return delta_vel

def update_velocities(vel, i, j, delta_vel):
   vel[i-1] = vel[i-1] + delta_vel
   vel[j-1] = vel[j-1] - delta_vel

def plot_gamma_vs_iterations(lcollisions,lgamma):
   plt.plot(lcollisions,lgamma)
   plt.xlabel('number of collisions')
   plt.ylabel('gamma order parameter')
   plt.show()

nstep = 10
filelines = read_file_lines('../results.txt')
ncollisions = extract_ncollisions(filelines)
pos = extract_initial_positions(filelines)
vel = extract_initial_velocities(filelines)
a = compute_fcc_latpar(len(pos))
collisions_list = [0]
gamma_list = [compute_order_parameter(pos, a)]
for c in range(ncollisions):
   time = extract_time(c, filelines)
   update_positions(pos, time, vel)
   if c % nstep == 0:
      collisions_list.append(c+1)
      gamma_list.append(compute_order_parameter(pos, a))
   (i,j) = extract_indices(c, filelines)
   delta_vel = extract_delta_vel(c, filelines)
   update_velocities(vel, i, j, delta_vel)
plot_gamma_vs_iterations(collisions_list,gamma_list)
