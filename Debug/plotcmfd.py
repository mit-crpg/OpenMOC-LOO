import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from itertools import cycle
import numpy as np
import os
import os.path 
import sys


flag = 'loo1'
geometries = ['geometry_c5g7.xml']
geometry = geometries[0]

# line style
lines = ['-o']
linecycle = cycle(lines)

# font sizes 
fontP = FontProperties()
fontP.set_size('small')

max_num_lines = 0
l2_norm_files = []
num = 1;
# get all l2_norm file names in directory
for file in os.listdir("."):
    if file.startswith(geometry[0:-4] + "_l2_norm_64_0.50",) and file.endswith(flag + ".txt"):
        l2_norm_files.append(file)
        num = num+1

    basepath = os.path.dirname('./')

# parse output files
counter = 0;
for file in l2_norm_files:
    counter = counter + 1
    filepath = os.path.abspath(os.path.join(basepath, file))
    logfile = open(filepath, "r").readlines()

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
            bi = file[-21:-20]
            break

    

    na = file[-31:-29]
    # find number of lines in file
    for num_lines, l in enumerate(logfile):
        pass

    max_num_lines = max(num_lines, max_num_lines)
    
    # create numpy arras
    iteration = np.zeros(num_lines)
    cell_l2   = np.zeros(num_lines)
    fsr_l2    = np.zeros(num_lines)
    rho       = np.zeros(num_lines)

    for i, line in enumerate(logfile):
        if i is not 0:
            iteration[i-1] = line.split()[0]
            cell_l2[i-1]   = line.split()[1]
            fsr_l2[i-1]    = line.split()[3]
            rho[i-1]       = line.split()[7]

    # plot l2 norm
    var = []
    var.append(cell_l2);
    var.append(fsr_l2);
    var.append(rho);

    for i in range(3):
        plt.figure(i)
        plt.semilogy(iteration, var[i], 
                     next(linecycle),
                     color=cm.jet(1.*counter/num), 
                     label = ("%s bi %s %s damp %s" 
                              %(method, bi, update, damp)), 
                     markersize=5)
        plt.xlim(0, max_num_lines + 1)
        plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                   bbox_to_anchor=(0.5,-0.1),fancybox=True)


# save figure including different configuration of the same geometries.
plt.figure(0)
plt.xlabel('# MOC iteration')
plt.ylabel('Mesh Cell L2 Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(geometry[9:-4]))
plt.savefig(geometry[9:-4] + '_cell_l2_' + flag + '.png', bbox_inches='tight')
plt.clf()

plt.figure(1)
plt.xlabel('# MOC iteration')
plt.ylabel('FSR L2 Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(geometry[9:-4]) )
plt.savefig(geometry[9:-4] + '_fsr_l2_' + flag +  '.png', bbox_inches='tight')
plt.clf()

plt.figure(2)
plt.xlabel('# MOC iteration')
plt.ylabel('Spectral Radius (numerical approximation)')
plt.title('Geometry: %s,'%(geometry[9:-4]))
plt.savefig(geometry[9:-4] + '_rho_' + flag + '.png', bbox_inches='tight')
plt.clf()
