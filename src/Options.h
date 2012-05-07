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

#include <string.h>
#include <stdlib.h>
#include "log.h"

class Options {
private:
	std::string _relative_path;
	std::string _geometry_file;
	std::string _material_file;

	std::string _track_input_file;

	std::string _extension;
	double _track_spacing;
	int _num_azim;
	int _bit_dimension;
	std::string _verbosity;
	bool _dump_geometry;
	bool _plot_specs;
	bool _plot_fluxes;
	bool _compute_pin_powers;
	bool _compress_cross_sections;
	bool _cmfd;
public:
    Options(int argc, const char **argv);
    ~Options(void);
    const char *getGeometryFile() const;
    const char *getMaterialFile() const;
    bool dumpGeometry() const;
    double getNumAzim() const;
    int getBitDimension() const;
    double getTrackSpacing() const;
    const char* getVerbosity() const;
    std::string getExtension() const;
    bool plotSpecs() const;
    bool plotFluxes() const;
    bool computePinPowers() const;
    bool compressCrossSections() const;
	bool cmfd() const;
};

#endif
