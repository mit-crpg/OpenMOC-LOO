import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import numpy as np
import os


# Control variables: where C4 default is 0.5cm and 64 azimuthal angle
geometries = ['xml-sample/geometry_c5g7_fc.xml']
materials = ['xml-sample/material_c5g7.xml']

ts = [0.05] 
na = [32]
fc = [1e-8]

# Parameters for plotting
ls = ['--o', '-s', '-.v', '-<', '-^', '-.>', '--s', '-v', '-o', '-.<', '--<', '-->', '-.s']
fontP = FontProperties()
fontP.set_size('small')
max_num_lines = 0
num = 1;
counter = 0;

for i, geometry in enumerate(geometries):

    # Generates the geometry names
    geometry_no_slash = geometry.replace("/", "_")
    geometry_name = geometry_no_slash[20:-4]

    # Runs OpenMOC
    for spacing in ts:
        for angle in na:
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry 
                      + ' -na ' + str(angle) 
                      + ' -ts ' + str(spacing) 
                      + ' -fc ' + str(fc[0])
                      + ' -wc -df 0.7 -ro -ub')
    l2_norm_files = []

    # Obtain and sorts l2_norm file names in directory
    for file in os.listdir("../"):
        if file.startswith("%s_%s_%s"%(geometry_no_slash, na[0], ts[0])) and file.endswith(".txt"):
            l2_norm_files.append(file)
            num = num + 1
    l2_norm_files.sort()
    
    counter = 0;
    # parse output files
    for file in l2_norm_files:
        counter = counter + 1
        logfile = open('../'+file, "r").readlines()
        
        method = file[-8:-4]
        upscat = file[-15:-9]
        update = file[-22:-16]
        damp = file[-26:-23]
        bi = file[-28:-27]

        print("geometry = %s, ts = %s, na = %s, upscat = %s"%(geometry_name, ts, na, upscat))

        # find number of lines in file
        for num_lines, l in enumerate(logfile):
            pass

        max_num_lines = max(num_lines, max_num_lines)

        # create numpy arrays
        iteration = np.zeros(num_lines)
        fsr_l2    = np.zeros(num_lines)
        k_l2 = np.zeros(num_lines)

        # collect data together
        for k, line in enumerate(logfile):
            if k is not 0:
                iteration[k-1] = line.split()[0]
                fsr_l2[k-1]    = line.split()[1]
                k_l2[k-1] = line.split()[3]
                k_l2[k-1] = abs(k_l2[k-1])
        var = []
        var.append(fsr_l2);
        var.append(k_l2);

        # plotting :)
        for j in range(2):  
            plt.figure(j)
            plt.semilogy(iteration, 
                         var[j], 
                         ls[counter - 1],
                         color=cm.jet(1.*counter / num),
                         label = ("%s DF = %s"%(method, damp)), markersize=5)
            plt.xlim(0, max_num_lines + 1)
            plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                       bbox_to_anchor=(0.5,-0.1),fancybox=True)
        # end of one file

    # save figure including different configuration of the same geometries.
    plt.figure(0)
    plt.xlabel('# MOC iteration')
    plt.ylabel('L2 Norm on Cell Fission Source Relative Change')
    plt.title('Geometry: %s,'%(geometry_name) + ' spacing: %s,'%str(ts[0]) 
              + ' #angles: %s'%str(na[0]))
    plt.savefig(geometry_name + '_l2.png', bbox_inches='tight')
    plt.clf()

    plt.figure(1)
    plt.xlabel('# MOC iteration')
    plt.ylabel('k Relative Change')
    plt.title('Geometry: %s,'%(geometry_name) + ' spacing: %s,'%str(ts[0])
              + ' #angles: %s'%str(na[0]))
    plt.savefig(geometry_name + '_kl2.png', bbox_inches='tight')
    plt.clf()

