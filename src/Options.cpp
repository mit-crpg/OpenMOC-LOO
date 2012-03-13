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
	_track_spacing = 0.1;					 /* Default track spacing */
	_num_azim = 16;					 /* Default number of azimuthal angles */
	_bit_dimension = 1000;					 /* y dimension of tracks and segments plots */
	_verbosity = "NORMAL";				 /* Default logging level */
	_dump_geometry = false;				/* Default will not dump geometry */
	_extension = "png";				/* Default will plot png */
	_plot_materials = false;			/* Default will not plot materials */
	_plot_cells = false;				/* Default will not plot cells */
	_plot_fluxes = false;				/* Default will not plot fluxes */
	_compute_pin_powers = false;		/* Default will not compute pin powers */


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
			else if (strcmp(argv[i], "-pm") == 0 ||
					strcmp(argv[i], "--plotmaterials") == 0)
				_plot_materials = true;
			else if (strcmp(argv[i], "-pc") == 0 ||
					strcmp(argv[i], "--plotcells") == 0)
				_plot_cells = true;
			else if (strcmp(argv[i], "-pf") == 0 ||
					strcmp(argv[i], "--plotfluxes") == 0)
				_plot_fluxes = true;
			else if (strcmp(argv[i], "-cp") == 0 ||
					strcmp(argv[i], "--computepowers") == 0)
				_compute_pin_powers = true;
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
 * Returns a boolean representing whether or not to plot the materials.
 *  If true, the materials will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotMaterials() const {
	return _plot_materials;
}

/**
 * Returns a boolean representing whether or not to plot the cells.
 *  If true, the cells will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotCells() const {
	return _plot_cells;
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
