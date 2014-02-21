import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import numpy as np
import os


# Control variables: where C4 default is 0.5cm and 64 azimuthal angle
geometries = ['xml-sample/geometry_11x11_vacuum.xml']
materials = ['xml-sample/material_c5g7.xml']
ts = [0.1, 0.05] 
na = [32, 64]
fc = [1e-5]

# Parameters for plotting
ls = ['--o', '-s', '-.v', '-<', '-^', '-.>', '--s', '-v']
fontP = FontProperties()
fontP.set_size('small')
max_num_lines = 0
num = 1;
counter = 0;

for i, geometry in enumerate(geometries):

    # Generates the geometry names
    geometry_no_slash = geometry.replace("/", "_")
    geometry_name = geometry_no_slash[25:-4]

    # Runs OpenMOC
    for spacing in ts:
        for angle in na:
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry 
                      + ' -na ' + str(angle) 
                      + ' -ts ' + str(spacing) 
                      + ' -fc ' + str(fc[0]) 
                      + ' -wc -df 1.0')

            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry 
                      + ' -na ' + str(angle) 
                      + ' -ts ' + str(spacing) 
                      + ' -fc ' + str(fc[0]) 
                      + ' -wc -df 0.7')

    l2_norm_files = []

    # Obtain and sorts l2_norm file names in directory
    for file in os.listdir("../"):
        if file.startswith(geometry_no_slash) and file.endswith(".txt"):
            l2_norm_files.append(file)
            num = num + 1
    l2_norm_files.sort()

    # parse output files
    for file in l2_norm_files:
        counter = counter + 1
        logfile = open('../'+file, "r").readlines()
        
        method = file[-8:-4]
        update = file[-15:-9]
        damp = file[-19:-16]
        bi = file[-21:-20]
        print("ts = %s, na = %s"%(ts, na))

        # find number of lines in file
        for num_lines, l in enumerate(logfile):
            pass

        max_num_lines = max(num_lines, max_num_lines)

        # create numpy arrays
        iteration = np.zeros(num_lines)
        fsr_l2    = np.zeros(num_lines)
        rho       = np.zeros(num_lines)

        # collect data together
        for k, line in enumerate(logfile):
            if k is not 0:
                iteration[k-1] = line.split()[0]
                fsr_l2[k-1]    = line.split()[3]
                rho[k-1]       = line.split()[7]

        var = []
        var.append(fsr_l2);
        var.append(rho);

        # plotting :)
        for j in range(2):  
            plt.figure(j)
            plt.semilogy(iteration, var[j], ls[counter - 1],
                         color=cm.jet(1.*counter / num),
                         label = ("%s %s"%(method, update)), markersize=5)
            plt.xlim(0, max_num_lines + 1)
            plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                       bbox_to_anchor=(0.5,-0.1),fancybox=True)
        # end of one file

    # save figure including different configuration of the same geometries.
    plt.figure(0)
    plt.xlabel('# MOC iteration')
    plt.ylabel('FSR L2 Norm on Fission Source Relative Change')
    plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
              + ' #angles: %s'%str(na[0]))
    plt.savefig(geometry_name + '_fsr_l2.png', bbox_inches='tight')
    plt.clf()

    plt.figure(1)
    plt.xlabel('# MOC iteration')
    plt.ylabel('Spectral Radius (numerical approximation)')
    plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
              + ' #angles: %s'%str(na[0]))
    plt.savefig(geometry_name + '_rho.png', bbox_inches='tight')
    plt.clf()
