import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from itertools import cycle
import numpy as np
import os
import os.path 
import sys

geometries = ['geometry_c5g7.xml']
geometry = geometries[0]

# line style
lines = ['-o']
linecycle = cycle(lines)

# font sizes 
fontP = FontProperties()
fontP.set_size('small')


figures = ['fig11', 'fig12'];

for fig in figures:
    max_num_lines = 0
    l2_norm_files = []

    # get all l2_norm file names in directory
    for file in os.listdir("./"+fig):
        if file.startswith("l2_norm_64_0.50",) and file.endswith(".txt"):
            l2_norm_files.append(file)

    basepath = os.path.dirname('./')

    # parse output files
    for file in l2_norm_files:
        filepath = os.path.abspath(os.path.join(basepath, fig, file))
        logfile = open(filepath, "r").readlines()
        method = file[-8:-4]
        damp = file[-19:-16]
        print(damp)
        ts = file[-29:-25]
        na = file[-31:-29]
        # find number of lines in file
        for num_lines, l in enumerate(logfile):
            pass

        max_num_lines = max(num_lines, max_num_lines)

        # create numpy arras
        iteration = np.zeros(num_lines)
        fsr_l2    = np.zeros(num_lines)

        for i, line in enumerate(logfile):
            if i is not 0:
                iteration[i-1] = line.split()[0]
                fsr_l2[i-1]    = line.split()[3]

        # plot l2 norm
        var = []
        var.append(fsr_l2);

        plt.figure(0)
        plt.semilogy(iteration, var[0], next(linecycle), 
                     label = ("%s damping factor %s" %(method, damp)), 
                     markersize=5)
        plt.xlim(0, max_num_lines + 1)
        plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                   bbox_to_anchor=(0.5,-0.1),fancybox=True)


    # save figure including different configuration of the same geometries.
    plt.figure(0)
    plt.xlabel('# MOC iteration')
    plt.ylabel('FSR L2 Norm on Fission Source Relative Change')
    plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
              + ' #angles: %s'%str(na[0]))
    plt.savefig(geometry[9:-4] + '_fsr_l2_' + fig + '.png', bbox_inches='tight')
    plt.clf()
