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
Options::Options(int argc, const char **argv) {

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
	_cmfd = true; 						/* Default will not perform CMFD acceleration */
	_update_flux = false;  				/* Default will not use CMFD to update flux */
	_plot_current = false;				/* Default will not plot net current */
	_keff_conv_thresh = 1e-6;			/* Default will set keff conv thresh to 1e-6 */
	_multigroup = false;				/* Default sets CMFD to one group structure */
	_print_matrices = false;			/* Default will not print matrices */
	_cmfd_level = 1;					/* Default cmfd level is 1 (hightest level) */


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
			else if (strcmp(argv[i], "-pc") == 0 ||
					strcmp(argv[i], "--plotcurrent") == 0)
				_plot_current = true;
			else if (LAST("--keffconv") || LAST("-kc"))
				_keff_conv_thresh = atof(argv[i]);
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
    return _geometry_file.c_str();
}

/**
 * Returns a character array with the path to the material input file. By default this
 * will return the path to /xml-sample/1/material.xml if not set at runtime from the
 * console
 * @return path to the geometry input file
 */
const char *Options::getMaterialFile() const {
    return _material_file.c_str();
}


/**
 * Returns a boolean representing whether or not to dump the geometry to the
 * console. If true, the geometry will be printed out after parsing is complete
 * @return whether or not to dump the geometry to the console
 */
bool Options::dumpGeometry() const {
	return _dump_geometry;
}

/**
 * Returns the number of azimuthal angles. By default this will return 128 angles if
 * not set at runtime from the console
 * @return the number of azimuthal angles
 */
double Options::getNumAzim() const {
    return _num_azim;
}

/**
 * Returns the y dimension of plots. By default this will return 1000 bits
 * (or pixels) if not set at runtime from the console
 * @return the y dimension of plots.
 */
int Options::getBitDimension() const{
	return _bit_dimension;
}



/**
 * Returns the track spacing. By default this will return 0.05 if not set at runtime
 * from the console
 * @return the track spacing
 */
double Options::getTrackSpacing() const {
    return _track_spacing;
}


/**
 * Returns the verbosity logging level. By default this will return NORMAL if not set
 * at runtime from the console
 * @return the verbosity
 */
const char* Options::getVerbosity() const {
    return _verbosity.c_str();
}

/**
 * Returns the image files extension. By default this will return .png if not set
 * at runtime from the console
 * @return the image files extension
 */
std::string Options::getExtension() const {
    return _extension.c_str();
}

/**
 * Returns a boolean representing whether or not to plot the specs.
 *  If true, the specs will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotSpecs() const {
	return _plot_specs;
}

/**
 * Returns a boolean representing whether or not to plot the cells.
 *  If true, the cells will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotFluxes() const {
	return _plot_fluxes;
}


/**
 * Returns a boolean representing whether or not to compute the powers
 * in each pin. If true, txt files with the pin powers will be created
 * in a new directory called "PinPowers"
 * @return whether or not to compute the pin powers
 */
bool Options::computePinPowers() const {
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
bool Options::compressCrossSections() const {
	return _compress_cross_sections;
}

/**
 * Returns a boolean representing whether or not to perform CMFD acceleration
 * @return whether or not to perform CMFD acceleration
 */
bool Options::cmfd() const {
	return _cmfd;
}

/**
 * Returns a boolean representing whether or not to plot the net current.
 *  If true, the net current will be plotted in a file of _extension type
 * @return whether or not to plot net current
 */
bool Options::plotCurrent() const {
	return _plot_current;
}

/**
 * Returns a boolean representing whether or not to use CMFD to update flux
 * @return whether or not to use CMFD to update flux
 */
bool Options::updateFlux() const {
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

