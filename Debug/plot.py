import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import numpy as np
import os
import os.path 

# Control variables
flag = ''
#printgeometry='2D C5G7 (c5g7_cc_full)'
geometries = ['geometry_1810_2_8.xml']
"""geometries = ['geometry_1484_1_8.xml',
              'geometry_1484_2_8.xml',
              'geometry_1484_3_8.xml',
              'geometry_1484_9_8.xml',
              'geometry_1810_1_8.xml',
              'geometry_1810_2_8.xml',
              'geometry_1810_3_8.xml', 
              'geometry_1810_4_8.xml',
              'geometry_1810_5_8.xml',
              'geometry_1810_10_8.xml',
              'geometry_1810_12_8.xml',
              'geometry_1810_13_8.xml',
              'geometry_1810_14_8.xml',
              'geometry_1810_15_8.xml', 
              'geometry_1810_16_8.xml',
              'geometry_1810_17_8.xml',
              'geometry_1810_18_8.xml',
              'geometry_1810_19_8.xml',
              'geometry_1810_20_8.xml']
"""

# Define styles: four types of linestyles, 7 types of markers. 
ls = ['--o', '-s', '-.v', ':o', '-<', '-p', ':s','-^', '--s','-.p',
      '-.>', '--p', '-v', '-o', '-.<', ':p','--<', '-.o', ':v', '-->',
      '-.s','--v','--^', '-.^', '->', ':>', ':<',':^']
fontP = FontProperties()
fontP.set_size('small')

# Find the string of geometry name after "geometry_" before ".xml"
for geometry in geometries:
    index = geometry.find("geometry_") + 9
    geometry_name = geometry[index:-4]
    max_num_lines = 0
    l2_norm_files = []
    num = 1;

    for file in os.listdir("../"): # ../
        if file.startswith(geometry) and file.endswith(flag+".txt"):
            print("parsed file %s" % file)
            l2_norm_files.append(file)
            num = num+1
            basepath = os.path.dirname('../') # ../

    l2_norm_files.sort()

    # parse output files
    counter = 0;
    for file in l2_norm_files:
        counter = counter + 1
        filepath = os.path.abspath(os.path.join(basepath, file))
        logfile = open(filepath, "r").readlines()

        # sample: geometry_c5g7_cc.xml_1_64_0.05_1.0_np_loo1.txt
        method = file[-8:-4]
        update = file[-11:-9]
        damp = file[-15:-12]
        ts = file[-20:-16]
        na = file[-23:-21]
        index = file.find(".xml_") + 5
        closure = file[index:-24]
        print("plotted geometry %s, method = %s, closure = %s, update = %s, damp = %s, ts = %s, na = %s"%(geometry_name, method, closure, update, damp, ts, na))

        # find number of lines in file
        for num_lines, l in enumerate(logfile):
            pass

        max_num_lines = max(num_lines, max_num_lines)

        if (max_num_lines > 30):
            max_num_lines = 30

        # create numpy arras
        iteration = np.zeros(num_lines)
        fsr_l2    = np.zeros(num_lines)
        fsr_ratio = np.zeros(num_lines)

        for i, line in enumerate(logfile):
            if i is not 0:
                iteration[i-1] = line.split()[0]
                fsr_l2[i-1]    = line.split()[1]
                fsr_ratio[i-1] = line.split()[2]

        # plot l2 norm
        var = [];
        var.append(fsr_l2);
        var.append(fsr_ratio);
        
        for j in range(2):
            plt.figure(j)
            plt.semilogy(iteration, 
                         var[j], 
                         ls[counter], 
                         color=cm.jet(1.*counter/num), 
                         label = ("%s #%s %s"%(method, closure, update)), markersize=5)
            plt.xlim(0, max_num_lines + 1)
            plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                       bbox_to_anchor=(0.5,-0.1),fancybox=True)

    # save figure including different configuration of the same geometries.
    plt.figure(0)
    plt.xlabel('# MOC iteration')
    plt.ylabel('RMS of Relative Change in Cell Fission Source')
    plt.title('Geometry: %s,'%(geometry_name) + ' spacing: %s cm,'%str(ts) 
              + ' #angles: %s'%str(na))
    plt.savefig(geometry_name +  '_eps.png', bbox_inches='tight')
    plt.clf()

    plt.figure(1)
    plt.yscale('linear', nonposy='clip')
    plt.xlabel('# MOC iteration')
    plt.ylabel('Ratio of Successive Iteration RMS of Relative Change in Cell Fission Source')
    plt.title('Geometry: %s,'%(geometry_name) + ' spacing: %s cm,'%str(ts) 
              + ' #angles: %s'%str(na))
    plt.savefig(geometry_name +  '_eps_ratio.png', bbox_inches='tight')
    plt.clf()