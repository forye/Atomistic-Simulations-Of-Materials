import math
import matplotlib.pyplot as plt
import numpy as np

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
   delta_vel = np.array(map(float, lines[line_index].split()))
   return delta_vel

def update_velocities(vel, i, j, delta_vel):
   vel[i-1] = vel[i-1] + delta_vel
   vel[j-1] = vel[j-1] - delta_vel

def initialize_lists(list_values, list_count, bin_size, max_value):
   nbins = int(max_value/bin_size)
   for ibin in range(nbins):
      list_values.append(bin_size*(0.5+float(ibin)))
      list_count.append(0)

def count_vmag_in_bins(vel, count_list, bin_size, max_value):
   nvel = len(vel)
   for i in range(nvel):
      vmag = np.linalg.norm(vel[i])
      if vmag < max_value:
         count_list[int(vmag/bin_size)] = count_list[int(vmag/bin_size)] + 1

def plot_count_vs_vmag(xlist, ylist, bin_size, max_value):
   ylist_scaled = [x / (float(sum(ylist))*bin_size) for x in ylist]
   plt.plot(xlist, ylist_scaled,'ro')
   plt.xlabel('velocity size')
   plt.ylabel('distribution function')
   x = (np.arange(0.0, max_value, 0.1)).tolist()
   y = []
   for element in x:
      y.append(fmb(element))
   plt.plot(x,y)
   plt.show()

def fmb(x):
   return math.sqrt(2.0/math.pi)*x*x*math.exp(-x*x/2.0)

bin_size = 0.1
max_value = 5.0
cmin = 9999
cmax = 14999
filelines = read_file_lines('results.txt')
ncollisions = extract_ncollisions(filelines)
vel = extract_initial_velocities(filelines)
counting_list = []
vmag_list = []
initialize_lists(vmag_list, counting_list, bin_size, max_value)
for c in range(ncollisions):
   (i,j) = extract_indices(c, filelines)
   delta_vel = extract_delta_vel(c, filelines)
   update_velocities(vel, i, j, delta_vel)
   if c >= cmin and c <= cmax:
      count_vmag_in_bins(vel, counting_list, bin_size, max_value)
plot_count_vs_vmag(vmag_list, counting_list, bin_size, max_value)
