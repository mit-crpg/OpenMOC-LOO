/*
 * Options.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#include "Options.h"

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)


/**
 * Options constructor
 * @param argc the number of command line arguments from console
 * @param argv a char array of command line arguments from console
 */
Options::Options(int argc, char **argv) {
	_cmfd = false; 					
	_loo = true;
	_plot_quad_flux = true; 
	/* Checks the working directory to set the relative path for input files
	 * This is important so that default input files work when program is run
	 * from both eclipse and the console */
	if (std::string(getenv("PWD")).find("Release") != std::string::npos)
		_relative_path = "../";
	else
		_relative_path = "";

	_geometry_file = _relative_path + "xml-sample/SimpleLattice/geometry.xml"; 	 /* Default geometry input file */
	_material_file = _relative_path + "xml-sample/SimpleLattice/material.xml";    /* Default material input file */
	_track_spacing = 0.1;				/* Default track spacing */
	_num_azim = 16;						/* Default number of azimuthal angles */
	_bit_dimension = 1000;				/* y dimension of tracks and segments plots */
	_verbosity = "NORMAL";				/* Default logging level */
	_dump_geometry = false;				/* Default will not dump geometry */
	_extension = "png";					/* Default will plot png */
	_plot_specs = false;            	/* Default will not plot materials, cells, FSRs, tracks, or segments */
	_plot_fluxes = false;				/* Default will not plot fluxes */
	_compute_pin_powers = false;		/* Default will not compute pin powers */
	_compress_cross_sections = false;	/* Default will not compress cross-sections */
	_update_flux = false;  				/* Default will not use CMFD to update flux */
	_plot_current = false;				/* Default will not plot cmfd currents */
	_plot_diffusion = false;			/* Default will not plot diffusion flux */
	_keff_conv_thresh = 1e-6;			/* Default will set keff conv thresh to 1e-6 */
	_multigroup = false;				/* Default sets CMFD to one group structure */
	_print_matrices = false;			/* Default will not print matrices */
	_cmfd_level = 1;					/* Default cmfd level is 1 (hightest level) */
	_plot_keff = false;					/* Default will not plot keff */
	_transient = false;					/* Default will not solve transient problem */
	_time_end = 1.0;					/* Default end time is 1 second */
	_dt_outer = 1.0e-4;					/* Default outer time step is 1.0e-4 seconds */
	_dt_inner = 5.0e-7;					/* Default inner time step is 5.0e-7 seconds */
	_diffusion = false;					/* Default will not solve diffusion problem */

	for (int i = 0; i < argc; i++) {
		if (i > 0) {
			if (LAST("--geometryfile") || LAST("-g")) {
				_geometry_file = argv[i];
			}
			else if (LAST("--materialfile") || LAST("-m")) {
				_material_file = argv[i];
			}
			else if (LAST("--trackspacing") || LAST("-ts"))
				_track_spacing = atof(argv[i]);
			else if (LAST("--numazimuthal") || LAST("-na"))
				_num_azim = atoi(argv[i]);
			else if (LAST("--bitdimension") || LAST("-bd"))
							_bit_dimension = atoi(argv[i]);
			else if (LAST("--verbosity") || LAST("-v"))
				_verbosity = strdup(argv[i]);
			else if (strcmp(argv[i], "-dg") == 0 ||
					strcmp(argv[i], "--dumpgeometry") == 0)
				_dump_geometry = true;
			else if (LAST("--extension") || LAST("-ex"))
							_extension = argv[i];
			else if (strcmp(argv[i], "-ps") == 0 ||
					strcmp(argv[i], "--plotspecs") == 0)
				_plot_specs = true;
			else if (strcmp(argv[i], "-pf") == 0 ||
					strcmp(argv[i], "--plotfluxes") == 0)
				_plot_fluxes = true;
			else if (strcmp(argv[i], "-cp") == 0 ||
					strcmp(argv[i], "--computepowers") == 0)
				_compute_pin_powers = true;
			else if (strcmp(argv[i], "-cxs") == 0 ||
					strcmp(argv[i], "--compressxs") == 0)
				_compress_cross_sections = true;
			else if (strcmp(argv[i], "-uf") == 0 ||
					strcmp(argv[i], "--updateflux") == 0)
				_update_flux = true;
			else if (strcmp(argv[i], "-nc") == 0 ||
					strcmp(argv[i], "--nocmfd") == 0)
				_cmfd = false;
			else if (strcmp(argv[i], "-wc") == 0 ||
					strcmp(argv[i], "--withcmfd") == 0)
				_cmfd = true;
			else if (strcmp(argv[i], "-nl") == 0 ||
					strcmp(argv[i], "--noloo") == 0)
				_loo = false;
			else if (strcmp(argv[i], "-wl") == 0 ||
					strcmp(argv[i], "--withloo") == 0)
				_loo = true;
			else if (strcmp(argv[i], "-pc") == 0 ||
					strcmp(argv[i], "--plotcurrent") == 0)
				_plot_current = true;
			else if (strcmp(argv[i], "-pk") == 0 ||
					strcmp(argv[i], "--plotkeff") == 0)
				_plot_keff = true;
			else if (strcmp(argv[i], "-tr") == 0 ||
					strcmp(argv[i], "--transient") == 0)
				_transient = true;
			else if (strcmp(argv[i], "-df") == 0 ||
					strcmp(argv[i], "--diffusion") == 0)
				_diffusion = true;
			else if (strcmp(argv[i], "-pd") == 0 ||
					strcmp(argv[i], "--plotdiffusion") == 0)
				_plot_diffusion = true;
			else if (LAST("--keffconv") || LAST("-kc"))
				_keff_conv_thresh = atof(argv[i]);
			else if (LAST("--timeend") || LAST("-te"))
				_time_end = atof(argv[i]);
			else if (LAST("--timestepouter") || LAST("-dto"))
				_dt_outer = atof(argv[i]);
			else if (LAST("--timestepinner") || LAST("-dti"))
				_dt_inner = atof(argv[i]);
			else if (strcmp(argv[i], "-mg") == 0 ||
					strcmp(argv[i], "--multigroup") == 0)
				_multigroup = true;
			else if (strcmp(argv[i], "-pm") == 0 ||
					strcmp(argv[i], "--printmatrices") == 0)
				_print_matrices = true;
			else if (LAST("--cmfdlevel") || LAST("-cl"))
				_cmfd_level = atoi(argv[i]);
		}
	}
}

Options::~Options(void) { }

/**
 * Returns a character array with the path to the geometry input file. By default this
 * will return the path to /xml-sample/1/geometry.xml if not set at runtime from the
 * console
 * @return path to the geometry input file
 */
const char *Options::getGeometryFile() const {
    return (const char*)_geometry_file.c_str();
}

/**
 * Returns a character array with the path to the material input file. By default this
 * will return the path to /xml-sample/1/material.xml if not set at runtime from the
 * console
 * @return path to the geometry input file
 */
const char *Options::getMaterialFile() const {
    return (const char*)_material_file.c_str();
}


/**
 * Returns a boolean representing whether or not to dump the geometry to the
 * console. If true, the geometry will be printed out after parsing is complete
 * @return whether or not to dump the geometry to the console
 */
bool Options::dumpGeometry(){
	return _dump_geometry;
}

/**
 * Returns the number of azimuthal angles. By default this will return 128 angles if
 * not set at runtime from the console
 * @return the number of azimuthal angles
 */
double Options::getNumAzim(){
    return _num_azim;
}

/**
 * Returns the y dimension of plots. By default this will return 1000 bits
 * (or pixels) if not set at runtime from the console
 * @return the y dimension of plots.
 */
int Options::getBitDimension(){
	return _bit_dimension;
}



/**
 * Returns the track spacing. By default this will return 0.05 if not set at runtime
 * from the console
 * @return the track spacing
 */
double Options::getTrackSpacing(){
    return _track_spacing;
}


/**
 * Returns the verbosity logging level. By default this will return NORMAL if not set
 * at runtime from the console
 * @return the verbosity
 */
char* Options::getVerbosity(){
    return (char*)_verbosity.c_str();
}

/**
 * Returns the image files extension. By default this will return .png if not set
 * at runtime from the console
 * @return the image files extension
 */
std::string Options::getExtension(){
    return _extension;
}

/**
 * Returns a boolean representing whether or not to plot the specs.
 *  If true, the specs will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotSpecs(){
	return _plot_specs;
}

/**
 * Returns a boolean representing whether or not to plot the cells.
 *  If true, the cells will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotFluxes(){
	return _plot_fluxes;
}


/**
 * Returns a boolean representing whether or not to compute the powers
 * in each pin. If true, txt files with the pin powers will be created
 * in a new directory called "PinPowers"
 * @return whether or not to compute the pin powers
 */
bool Options::computePinPowers(){
	return _compute_pin_powers;
}


/**
 * Returns a boolean representing whether or not to compress the
 * cross-sections for each material. If true, the starting and ending
 * index for the non-zero cross-section values will be computed to
 * help speed up fixed source iteraiton. Note: this will only speed
 * up fixed source iteration for materials with many zeroes in their
 * cross-section values
 * @return whether or not to compute the pin powers
 */
bool Options::compressCrossSections(){
	return _compress_cross_sections;
}

/**
 * Returns a boolean representing whether or not to perform CMFD acceleration
 * @return whether or not to perform CMFD acceleration
 */
bool Options::cmfd(){
	return _cmfd;
}

/**
 * Returns a boolean representing whether or not to plot the cmfd fluxes
 * at each step. If true, the net current will be plotted in a file of 
 * _extension type
 * @return whether or not to plot net current
 */
bool Options::plotCurrent(){
	return _plot_current;
}

/**
 * Returns a boolean representing whether or not to plot the loo surface
 * averaged quadrature fluxes at each step. If true, the quadrature fluxes 
 * will be plotted in a file of _extension type
 * @return whether or not to plot quadrature flux
 */
bool Options::plotQuadFlux(){
	return _plot_quad_flux;
}

/**
 * Returns a boolean representing whether or not to plot the diffusion flux.
 *  If true, the net current will be plotted in a file of _extension type
 * @return whether or not to plot net current
 */
bool Options::plotDiffusion(){
	return _plot_diffusion;
}

/**
 * Returns a boolean representing whether or not to plot keff.
 *  If true, the net current will be plotted in a file of _extension type
 * @return whether or not to plot net current
 */
bool Options::plotKeff(){
	return _plot_keff;
}

/**
 * Returns a boolean representing whether or not to use CMFD to update flux
 * @return whether or not to use CMFD to update flux
 */
bool Options::updateFlux(){
	return _update_flux;
}

/**
 * Returns the keff convergence threshold
 * @return keff convergence threshold
 */
double Options::getKeffConvThresh() {
	return _keff_conv_thresh;
}

/**
 * Returns bool telling us cmfd group structure
 * @return bool telling us cmfd group structure
 */
bool Options::getGroupStructure() {
	return _multigroup;
}


/**
 * Returns bool telling us cmfd group structure
 * @return bool telling us cmfd group structure
 */
bool Options::getPrintMatrices() {
	return _print_matrices;
}

int Options::getCmfdLevel(){
	return _cmfd_level;
}

bool Options::getCmfd(){
	return _cmfd;
}

bool Options::getLoo(){
	return _loo;
}

bool Options::getTransient(){
	return _transient;
}

bool Options::getDiffusion(){
	return _diffusion;
}

double Options::getTimeEnd(){
	return _time_end;
}

double Options::getTimeStepOuter(){
	return _dt_outer;
}

double Options::getTimeStepInner(){
	return _dt_inner;
}




