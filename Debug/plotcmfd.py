import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from itertools import cycle
import numpy as np
import os

geometries = ['geometry_c5g7_refl.xml']
geometry = geometries[0]

# plot color 
plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])
#num_plots = 4
#colormap = plt.cm.gist_ncar
#plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
#colors = ('b','r','g','c', 'm', 'y', 'k')
#colors=('red', 'green', 'blue', 'yellow')
#colorcycle = cycle(colors)

# line style
lines = ['-o']
linecycle = cycle(lines)

# font sizes 
fontP = FontProperties()
fontP.set_size('small')

max_num_lines = 0

# list of l2_norm files
l2_norm_files = []

# get all l2_norm file names in directory
for file in os.listdir("."):
    if file.startswith(geometry[0:-4] + "_l2_norm_64_0.50",) and file.endswith(".txt"):
        l2_norm_files.append(file)
		
# parse output files
for file in l2_norm_files:
    logfile = open(file, "r").readlines()
        
    for c_num, c in enumerate(file):
        if c == '_':
            method = file[-8:-4]
            break

    for c_num, c in enumerate(file):
        if c == '_':
            update = file[-15:-9]
            break

    for c_num, c in enumerate(file):
        if c == '_':
            damp = file[-19:-16]
            break

    for c_num, c in enumerate(file):
        if c == '_':
            ts = file[-24:-20]
            break

    for c_num, c in enumerate(file):
        if c == '_':
            na = file[-27:-25]
            break

    # find number of lines in file
    for num_lines, l in enumerate(logfile):
        pass

    max_num_lines = max(num_lines, max_num_lines)

    # create numpy arras
    iteration = np.zeros(num_lines)
    cell_l2   = np.zeros(num_lines)
    fsr_linf  = np.zeros(num_lines)
    fsr_l2    = np.zeros(num_lines)
    num       = np.zeros(num_lines)

    for i, line in enumerate(logfile):
        if i is not 0:
            iteration[i-1] = line.split()[0]
            cell_l2[i-1]   = line.split()[1]
            fsr_linf[i-1]  = line.split()[2]
            fsr_l2[i-1]    = line.split()[3]
            num[i-1]       = line.split()[5]

	# plot l2 norm
    var = []
    var.append(cell_l2);
    var.append(fsr_linf);
    var.append(fsr_l2);
    var.append(num);

    for i in range(4):  
        plt.figure(i)
        plt.semilogy(iteration, var[i], next(linecycle), 
        label = ("%s boundary %s damp %s" %(method, update, damp)), 
                     markersize=5)
        plt.xlim(0, max_num_lines + 1)
        plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                   bbox_to_anchor=(0.5,-0.1),fancybox=True)
    # end of one file

# save figure including different configuration of the same geometries.
plt.figure(0)
plt.xlabel('# MOC iteration')
plt.ylabel('Mesh Cell L2 Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
		  + ' #angles: %s'%str(na[0]))
plt.savefig(geometry[9:-4] + '_cell_l2.png', bbox_inches='tight')
plt.clf()

plt.figure(1)
plt.xlabel('# MOC iteration')
plt.ylabel('FSR L-infinity Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
	  + ' #angles: %s'%str(na[0]))
plt.savefig(geometry[9:-4] + '_fsr_linf.png', bbox_inches='tight')
plt.clf()

plt.figure(2)
plt.xlabel('# MOC iteration')
plt.ylabel('FSR L2 Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
			  + ' #angles: %s'%str(na[0]))
plt.savefig(geometry[9:-4] + '_fsr_l2.png', bbox_inches='tight')
plt.clf()

plt.figure(3)
plt.xlabel('# MOC iteration')
plt.ylabel('# Acceleration Iterations Taken at Each MOC Iteration')
plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
  + ' #angles: %s'%str(na[0]))
plt.savefig(geometry[9:-4] + '_num_acc.png', bbox_inches='tight')
plt.clf()
