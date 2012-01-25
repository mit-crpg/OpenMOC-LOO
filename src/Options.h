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
	double _track_spacing;
	int _num_azim;
	char* _verbosity;
public:
    Options(int argc, const char **argv);
    char *getGeometryFile() const;
    double getNumAzim() const;
    double getTrackSpacing() const;
    char* getVerbosity() const;
};

#endif
