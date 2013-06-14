import matplotlib.pyplot as plt
import numpy as np
import os

# run OpenMOC
geometries = ['geometry_pin.xml']
materials = ['material_simple.xml']

#geometries = ['geometry_c5g7_wo_refl.xml', 'geometry_c5g7.xml', 'geometry_UO2.xml', 'geometry_MOX.xml', 'geometry_pin.xml']
#materials = ['material_c5g7.xml', 'material_c5g7.xml', 'material_c5g7.xml', 'material_c5g7.xml', 'material_simple.xml']
#geometries = ['geometry_LRA.xml']
#materials = ['material_LRA.xml']
#ts = [.2, .1, .05, .025, .01]

ts = [.2, .05]
na = [32, 128]

for i, geometry in enumerate(geometries):
    for spacing in ts:
        for angle in na:
            if geometry != 'geometry_pin.xml':
                os.system('../bin/openmoc'
						  + ' -m ../xml-sample/Cmfd/' + materials[i]
						  + ' -g ../xml-sample/Cmfd/' + geometry 
						  + ' -na ' + str(angle) 
						  + ' -ts ' + str(spacing) 
						  + ' -fc 1e-5 -uk -wc -mg -cl 2')
            else:
                os.system('../bin/openmoc'
						  + ' -m ../xml-sample/Cmfd/' + materials[i]
						  + ' -g ../xml-sample/Cmfd/' + geometry 
						  + ' -na ' + str(angle) 
						  + ' -ts ' + str(spacing) 
						  + ' -fc 1e-5 -uk -wl -mg')

    # list of l2_norm files
    l2_norm_files = []

    # make figure objects
    fig_l2_norm = plt.figure(1)
    fig_keff = plt.figure(2)

    # get all l2_norm files in directory
    for file in os.listdir("."):
        if file.startswith("l2_norm") and file.endswith(".txt"):
            l2_norm_files.append(file)

    # parse output files
    for file in l2_norm_files:
        logfile = open(file, "r").readlines()
        os.rename(file, geometry[:-4] + file)
        
        for c_num, c in enumerate(file):
            if c == '_':
                angle_space = file[8:-4]
                break

        # find number of lines in file
        for num_lines, l in enumerate(logfile):
            pass

        # create numpy arras
        iteration = np.zeros(num_lines)
        l2_norm   = np.zeros(num_lines)
        keff      = np.zeros(num_lines)


        for i, line in enumerate(logfile):

            if i is not 0:
                iteration[i-1] = line.split()[0]
                l2_norm[i-1]   = line.split()[1]
                keff[i-1]      = line.split()[2]

        # plot l2 norm
        plt.figure(1)
        fit_coeff = np.polyfit(iteration, np.log(l2_norm), 1)
        fit_vals  = np.polyval(fit_coeff, iteration)
        plt.semilogy(iteration, l2_norm, '+', label=str(angle_space))
        plt.semilogy(iteration, np.exp(fit_vals))
        plt.xlim(0, iteration[num_lines-1] + 1)
        plt.xlabel('iteration')
        plt.ylabel('l2 norm')
        plt.title('l2 norm vs iteration')

        plt.figure(2)
        plt.plot(iteration, keff, '+', label=str(angle_space))
        plt.xlim(0, iteration[num_lines-1] + 1)
        plt.xlabel('iteration')
        plt.ylabel('keff')
        plt.title('keff vs iteration')


    plt.figure(1)
    plt.savefig(geometry[9:-4] + '_l2_norm.png')
    plt.figure(2)
    plt.savefig(geometry[9:-4] + '_keff.png')
