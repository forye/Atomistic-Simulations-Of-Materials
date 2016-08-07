import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import spheres_utilities.configuration as conf

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

def extract_nspheres(lines):
   nspheres = int(lines[0].split()[0])
   return nspheres

def extract_rvolume(lines):
   rvolume = float(lines[1].split()[0])
   return rvolume

def extract_ncollisions(lines):
   ncollisions = int(lines[2].split()[0])
   return ncollisions

def extract_ncollisions(lines):
   ncollisions = int(lines[2].split()[0])
   return ncollisions

def extract_initial_velocities(lines):
   nvel = int(lines[0].split()[0])
   vel = []
   for i in range(nvel):
      line_index = i + 3
      vel.append(np.array(map(float,lines[line_index].split()[4:7])))
   return vel

def extract_indices(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 4
   i, j = map(int, lines[line_index].split())
   return i, j

def extract_delta_vel(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 5
   # delta_vel = np.array(map(float, lines[line_index].split()))
   try: delta_vel = np.array(map(float, re.sub( "[\[\]]","",lines[line_index] ) .split()))
   except:
       print "error! col = " , c,line_index
       return np.array([0,0,0])
   return delta_vel

def extract_time(c, lines):
   nspheres= int(lines[0].split()[0])
   line_index = 3 + nspheres + 6*c + 3
   ctime = float(lines[line_index].split()[0])
   return ctime

# cmin = 9999
# cmax = 14999-5
# filelines = read_file_lines('results.txt')
# nspheres = extract_nspheres(filelines)
# rvolume = extract_rvolume(filelines)
# sigma = (math.sqrt(2.0) / (rvolume*float(nspheres)))**(1.0/3.0)
# ncollisions = extract_ncollisions(filelines)
# t = 0.0
# sum_abs_delta_vel = 0.0
# for c in range(cmin, cmax):
#    delta_vel = extract_delta_vel(c, filelines)
#    ctime = extract_time(c, filelines)
#    t = t + ctime
#    sum_abs_delta_vel = sum_abs_delta_vel + np.linalg.norm(delta_vel)
# z = 1 + sigma*sum_abs_delta_vel / (nspheres*3.0*t)
# print rvolume, z


def calc_pressure():

    cmin = 1#9999
    cmax = 145#14999-5
    # filelines = read_file_lines('results.txt')
    filelines = read_file_lines(conf.RESULTS)
    nspheres = extract_nspheres(filelines)
    rvolume = extract_rvolume(filelines)
    sigma = (math.sqrt(2.0) / (rvolume*float(nspheres)))**(1.0/3.0)
    ncollisions = extract_ncollisions(filelines)
    t = 0.0
    sum_abs_delta_vel = 0.0
    for c in range(cmin, cmax):
       delta_vel = extract_delta_vel(c, filelines)
       ctime = extract_time(c, filelines)
       t = t + ctime
       sum_abs_delta_vel = sum_abs_delta_vel + np.linalg.norm(delta_vel)
    z = 1 + sigma*sum_abs_delta_vel / (nspheres*3.0*t)
    return rvolume, z
