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
	std::string _extension;
	double _track_spacing;
	int _num_azim;
	int _bit_dimension;
	int _num_sectors;
	int _num_rings;
	double _sector_offset;
	std::string _verbosity;
	bool _dump_geometry;
	bool _plot_materials;
	bool _plot_cells;
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
    int getNumRings() const;
    int getNumSectors() const;
    double getSectorOffset() const;
    bool plotMaterials() const;
    bool plotCells() const;
};

#endif
