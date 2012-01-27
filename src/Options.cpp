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
#include <string.h>
#include <stdlib.h>

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)


/**
 * Options constructor
 * @param argc the number of command line arguments from console
 * @param argv a char array of command line arguments from console
 */
Options::Options(int argc, const char **argv) {

	_geometry_file = strdup("xml-sample/1/geometry.xml"); 	 /* Default geometry input file */
	_track_spacing = 0.05;									 /* Default track spacing */
	_num_azim = 128;										 /* Default number of azimuthal angles */
	_num_sectors = 0;										 /* Default number of sectors */
	_num_rings = 0;											 /* Default number of rings */
	_sector_offset = 0;										 /* Default sector offset */
	_verbosity = strdup("NORMAL");						 	 /* Default logging level */

	for (int i = 0; i < argc; i++) {
		if (i > 0) {
			if (LAST("--geometryfile") || LAST("-g")) {
				if (_geometry_file != NULL)
					free(_geometry_file);
				_geometry_file = strdup(argv[i]);
			}
			else if (LAST("--trackspacing") || LAST("-ts"))
				_track_spacing = atof(argv[i]);
			else if (LAST("--numazimuthal") || LAST("-na"))
				_num_azim = atoi(argv[i]);
			else if (LAST("--numrings") || LAST("-nr"))
				_num_rings = atoi(argv[i]);
			else if (LAST("--numsectors") || LAST("-ns"))
				_num_sectors = atoi(argv[i]);
			else if (LAST("--sectoroffset") || LAST("--so"))
				_sector_offset = atof(argv[i]);
			else if (LAST("--verbosity") || LAST("-v"))
				_verbosity = strdup(argv[i]);
		}
	}
}

Options::~Options(void) {
	if (this->_geometry_file != NULL)
		free(this->_geometry_file);
	if (this->_verbosity != NULL)
		free(this->_verbosity);
}

/**
 * Returns a character array with the path to the geometry input file. By default this
 * will return the path to /xml-sample/1/geometry.xml if not set at runtime from the
 * console
 * @return path to the geometry input file
 */
char *Options::getGeometryFile() const {
    return _geometry_file;
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
 * Returns the track spacing. By default this will return 0.05 if not set at runtime
 * from the console
 * @return the track spacing
 */
double Options::getTrackSpacing() const {
    return _track_spacing;
}


/**
 * Returns the numbe of rings used in subdividing flat source regions
 * @return the number of rings
 */
int Options::getNumRings() const {
    return _num_rings;
}


/**
 * Returns the number of sectors used in subdividing flat source regions
 * @return the number of sectors
 */
int Options::getNumSectors() const {
    return _num_sectors;
}


/**
 * Returns the angular offset of the sectors used in subdividing the flat
 * source regions
 * @return the angular offset in degrees
 */
double Options::getSectorOffset() const {
    return _sector_offset;
}

/**
 * Returns the verbosity logging level. By default this will return NORMAL if not set
 * at runtime from the console
 * @return the verbosity
 */
char* Options::getVerbosity() const {
    return _verbosity;
}
