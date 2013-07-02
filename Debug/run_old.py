import matplotlib.pyplot as plt
import numpy as np
import os

# run OpenMOC
geometries = ['geometry_LRA.xml', 
			  'geometry_UO2.xml', 
			  'geometry_c5g7_wo_refl.xml']
materials = ['material_LRA.xml', 
			 'material_c5g7.xml',
			 'material_c5g7.xml' ]
ts = [0.01]
na = [128]

for i, geometry in enumerate(geometries):
    for spacing in ts:
        for angle in na:
			os.system('../bin/openmoc'
					  + ' -m ../xml-sample/Cmfd/' + materials[i]
					  + ' -g ../xml-sample/Cmfd/' + geometry 
					  + ' -na ' + str(angle) 
					  + ' -ts ' + str(spacing) 
					  + ' -fc 1e-12')
			os.system('../bin/openmoc'
					  + ' -m ../xml-sample/Cmfd/' + materials[i]
					  + ' -g ../xml-sample/Cmfd/' + geometry 
					  + ' -na ' + str(angle) 
					  + ' -ts ' + str(spacing) 
					  + ' -fc 1e-12 -wc -cl 2')
			os.system('../bin/openmoc'
					  + ' -m ../xml-sample/Cmfd/' + materials[i]
					  + ' -g ../xml-sample/Cmfd/' + geometry 
					  + ' -na ' + str(angle) 
					  + ' -ts ' + str(spacing) 
					  + ' -fc 1e-12 -wl -cl 2')
			
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
        os.rename(file, geometry[:-4] + '_' + file)
        
        for c_num, c in enumerate(file):
            if c == '_':
                angle_space = file[:-4]
                break

        # find number of lines in file
		for num_lines, l in enumerate(logfile):
			pass

        # create numpy arras
		iteration = np.zeros(num_lines)
		l2_norm   = np.zeros(num_lines)
		linf      = np.zeros(num_lines)
		keff      = np.zeros(num_lines)
		num       = np.zeros(num_lines)

		for i, line in enumerate(logfile):
			if i is not 0:
				iteration[i-1] = line.split()[0]
				l2_norm[i-1]   = line.split()[1]
				linf[i-1]      = line.split()[2]
				keff[i-1]      = line.split()[4]
				num[i-1]       = line.split()[5]
				
				
		# plot l2 norm
		plt.figure(1)
		fit_coeff = np.polyfit(iteration, np.log(l2_norm), 1)
		fit_vals  = np.polyval(fit_coeff, iteration)
		plt.semilogy(iteration, l2_norm, '+', label = "damp = %s" %(angle_space))
		plt.semilogy(iteration, np.exp(fit_vals))
		plt.xlim(0, iteration[num_lines-1] + 1)
		plt.legend(loc=1, ncol=3, shadow=False)
		
		plt.figure(2)
		plt.semilogy(iteration, linf, '+', label= "damp = %s" %(angle_space))
		plt.xlim(0, iteration[num_lines-1] + 1)
		plt.legend(loc=1, ncol=3, shadow=False)
		
		plt.figure(3)
		plt.plot(iteration, keff, '+', label= "damp = %s" %(angle_space))
		plt.xlim(0, iteration[num_lines-1] + 1)
		plt.legend(loc=1, ncol=3, shadow=False)
		
		plt.figure(4)
		plt.plot(iteration, num, '+', label = "damp = %s" %(angle_space))
		plt.xlim(0, iteration[num_lines-1] + 1)
		plt.legend(loc=1, ncol=3, shadow=False)
		
	plt.figure(1)
	plt.xlabel('# MOC iteration')
	plt.ylabel('Mesh Cell L2 Norm on Fission Source')
	plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
			  + ' #angles: %s'%str(na[0]))
	plt.savefig(geometry[9:-4] + '_l2_norm.png')

	plt.figure(2)
	plt.xlabel('# MOC iteration')
	plt.ylabel('FSR Max Relative Change on Fission Source')
	plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
			  + ' #angles: %s'%str(na[0]))
	plt.savefig(geometry[9:-4] + '_linf_norm.png')

	plt.figure(3)
	plt.xlabel('# MOC iteration')
	plt.ylabel('keff relative change')
	plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
			  + ' #angles: %s'%str(na[0]))
	plt.savefig(geometry[9:-4] + '_keff.png')

	plt.figure(4)
	plt.xlabel('# MOC iteration')
	plt.ylabel('# Acceleration Iterations Taken at Each MOC Iteration')
	plt.title('Geometry: %s,'%(geometry[9:-4]) + ' spacing: %s,'%str(ts[0]) 
			  + ' #angles: %s'%str(na[0]))
	plt.savefig(geometry[9:-4] + '_num_acc.png')
