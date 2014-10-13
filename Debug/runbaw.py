import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import numpy as np
import os

# Control variables: where C4 default is 0.5cm and 64 azimuthal angle
geometries = ['xml-sample/geometry_1484_1_8.xml',
              'xml-sample/geometry_1484_2_8.xml',
              'xml-sample/geometry_1484_3_8.xml',
              'xml-sample/geometry_1484_9_8.xml',
              'xml-sample/geometry_1810_1_8.xml',
              'xml-sample/geometry_1810_2_8.xml',
              'xml-sample/geometry_1810_3_8.xml', 
              'xml-sample/geometry_1810_4_8.xml',
              'xml-sample/geometry_1810_5_8.xml',
              'xml-sample/geometry_1810_10_8.xml',
              'xml-sample/geometry_1810_12_8.xml',
              'xml-sample/geometry_1810_13_8.xml',
              'xml-sample/geometry_1810_14_8.xml',
              'xml-sample/geometry_1810_15_8.xml', 
              'xml-sample/geometry_1810_16_8.xml',
              'xml-sample/geometry_1810_17_8.xml',
              'xml-sample/geometry_1810_18_8.xml',
              'xml-sample/geometry_1810_19_8.xml',
              'xml-sample/geometry_1810_20_8.xml']

materials = ['xml-sample/material_1484_1_8g.xml', 
             'xml-sample/material_1484_2_8g.xml',
             'xml-sample/material_1484_3_8g.xml',
             'xml-sample/material_1484_9_8g.xml',
             'xml-sample/material_1810_1_8g.xml',
             'xml-sample/material_1810_2_8g.xml',
             'xml-sample/material_1810_3_8g.xml',
             'xml-sample/material_1810_4_8g.xml', 
             'xml-sample/material_1810_5_8g.xml',
             'xml-sample/material_1810_10_8g.xml',
             'xml-sample/material_1810_12_8g.xml',
             'xml-sample/material_1810_13_8g.xml',
             'xml-sample/material_1810_14_8g.xml',
             'xml-sample/material_1810_15_8g.xml',
             'xml-sample/material_1810_16_8g.xml', 
             'xml-sample/material_1810_17_8g.xml',
             'xml-sample/material_1810_18_8g.xml',
             'xml-sample/material_1810_19_8g.xml',
             'xml-sample/material_1810_20_8g.xml']

ts = [0.4] 
na = [16]
fc = [1e-5]

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
                      + ' -wc -df 0.5')
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry
                      + ' -na ' + str(angle)
                      + ' -ts ' + str(spacing)
                      + ' -fc ' + str(fc[0])
                      + ' -wc -lp')
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry
                      + ' -na ' + str(angle)
                      + ' -ts ' + str(spacing)
                      + ' -fc ' + str(fc[0])
                      + ' -wl1')
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry
                      + ' -na ' + str(angle)
                      + ' -ts ' + str(spacing)
                      + ' -fc ' + str(fc[0])
                      + ' -wl1 -lp')
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry
                      + ' -na ' + str(angle)
                      + ' -ts ' + str(spacing)
                      + ' -fc ' + str(fc[0])
                      + ' -wl2')
            os.system('cd .. && ./bin/openmoc'
                      + ' -m ' + materials[i]
                      + ' -g ' + geometry
                      + ' -na ' + str(angle)
                      + ' -ts ' + str(spacing)
                      + ' -fc ' + str(fc[0])
                      + ' -wl2 -lp')

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

        # collect data together
        for k, line in enumerate(logfile):
            if k is not 0:
                iteration[k-1] = line.split()[0]
                fsr_l2[k-1]    = line.split()[1]

        var = []
        var.append(fsr_l2);

        # plotting :)
        for j in range(1):  
            plt.figure(j)
            plt.semilogy(iteration, 
                         var[j], 
                         ls[counter - 1],
                         color=cm.jet(1.*counter / num),
                         label = ("%s %s"%(method, upscat)), markersize=5)
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
