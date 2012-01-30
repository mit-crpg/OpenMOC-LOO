/*
 * Options.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include "log.h"

class Options {
private:
	char* _geometry_file;
	char* _material_file;
	double _track_spacing;
	int _num_azim;
	int _num_sectors;
	int _num_rings;
	double _sector_offset;
	char* _verbosity;
	bool _dump_geometry;
public:
    Options(int argc, const char **argv);
    ~Options(void);
    char *getGeometryFile() const;
    char *getMaterialFile() const;
    bool dumpGeometry() const;
    double getNumAzim() const;
    double getTrackSpacing() const;
    char* getVerbosity() const;
    int getNumRings() const;
    int getNumSectors() const;
    double getSectorOffset() const;
};

#endif
