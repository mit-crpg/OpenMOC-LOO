import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import numpy as np
import os
import os.path 

# Control variables
flag = ''
printgeometry='2D C5G7 Benchmark'
geometries = ['xml-sample_geometry_c5g7_cc.xml']

# Find the string of geometry name after "geometry_" before ".xml"
geometry = geometries[0]
index = geometry.find("geometry_") + 9
geometry_name = geometry[index:-4]

ls = ['--o', '-s', '-.v', '-<', '-^', '-.>', '--s', '-v']
fontP = FontProperties()
fontP.set_size('small')

max_num_lines = 0
l2_norm_files = []
num = 1;

for file in os.listdir("../"):
    if file.startswith(geometry[0:-4]) and file.endswith(flag+".txt"):
        print("parsed file %s" % file)
        l2_norm_files.append(file)
        num = num+1

    basepath = os.path.dirname('../')

l2_norm_files.sort()

# parse output files
counter = 0;
for file in l2_norm_files:
    counter = counter + 1
    filepath = os.path.abspath(os.path.join(basepath, file))
    logfile = open(filepath, "r").readlines()

    # sample: xml-sample_geometry_c5g7_cc.xml_64_0.05_1.0_update_loo1.txt
    method = file[-8:-4]
    update = file[-15:-9]
    damp = file[-19:-16]
    ts = file[-24:-20]
    na = file[-27:-25]
    print("ts = %s, na = %s"%(ts, na))

    # find number of lines in file
    for num_lines, l in enumerate(logfile):
        pass

    max_num_lines = max(num_lines, max_num_lines)
    
    # create numpy arras
    iteration = np.zeros(num_lines)
    fsr_l2    = np.zeros(num_lines)
    rho       = np.zeros(num_lines)

    for i, line in enumerate(logfile):
        if i is not 0:
            iteration[i-1] = line.split()[0]
            fsr_l2[i-1]    = line.split()[3]
            rho[i-1]       = line.split()[7]

    # plot l2 norm
    var = [];
    var.append(fsr_l2);
    var.append(rho);

    for j in range(2):
        plt.figure(j)
        plt.semilogy(iteration, var[j], ls[counter-1], 
                     color=cm.jet(1.*counter/num), 
                     label = ("%s DF: %s"%(method, damp)), markersize=5)
        plt.xlim(0, max_num_lines + 1)
        plt.legend(loc='upper center', ncol=3, prop = fontP, shadow=True, 
                   bbox_to_anchor=(0.5,-0.1),fancybox=True)


# save figure including different configuration of the same geometries.
plt.figure(0)
plt.xlabel('# MOC iteration')
plt.ylabel('FSR L2 Norm on Fission Source Relative Change')
plt.title('Geometry: %s,'%(printgeometry) + ' spacing: %s cm,'%str(ts) 
          + ' #angles: %s'%str(na))

if flag == '':
    plt.savefig(geometry_name + '_fsr_l2' +  '.png', bbox_inches='tight')
else:
    plt.savefig(geometry_name + '_fsr_l2_' + flag + '.png', bbox_inches='tight')
plt.clf()

plt.figure(1)
plt.xlabel('# MOC iteration')
plt.ylabel('Apparent Spectral Radius')
plt.title('Geometry: %s,'%(printgeometry) + ' spacing: %s cm,'%str(ts) 
          + ' #angles: %s'%str(na))
if flag == '':
    plt.savefig(geometry_name + '_rho' + '.png', bbox_inches='tight')
else:
    plt.savefig(geometry_name + '_rho_' + flag + '.png', bbox_inches='tight')
plt.clf()
