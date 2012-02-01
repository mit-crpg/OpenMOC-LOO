/*
 * openmoc.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  This file defines the executable for OpenMOC.
 *
 */

#include "Geometry.h"
#include "TrackGenerator.h"
#include "Parser.h"
#include "log.h"
#include "Options.h"
#include "LocalCoords.h"
#include "Plotting.h"
#include "Geometry.h"

// FIXME: These should be removed when main() is properly implemented
#pragma GCC diagnostic ignored "-Wunused"
#pragma GCC diagnostic ignored "-Wunused-variable"

// TODO: This is just stubbed out for now
int main(int argc, const char **argv) {
	log_printf(NORMAL, "Starting OpenMOC...");

	/* Create an options class to parse command line options */
	Options opts(argc, argv);

	/* Set the verbosity */
	log_setlevel(opts.getVerbosity());

	/* Initialize the parser */
	Parser parser(&opts);

	/* Create an empty geometry */
	Geometry geometry(opts.getNumSectors(), opts.getNumRings(),
			  opts.getSectorOffset(), &parser);
	Plotting plotter(&geometry);

//	/* Print out geometry to console if requested at runtime*/
	if (opts.dumpGeometry())
		geometry.printString();

	/* Adjust the indices for each geometry class to use uids */
//	geometry.adjustKeys();

	/* Generate the neighbor cells for each surface in geometry */
//	geometry.buildNeighborsLists();

	/* Initialize the trackgenerator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
									opts.getTrackSpacing());

	track_generator.generateTracks();
	plotter.plotTracksTiff(&track_generator);
	track_generator.makeReflective();
	track_generator.segmentize();
	plotter.plotSegments();

	/* Testing findCell method */
//	LocalCoords* test_coords = new LocalCoords(0, 0);
//	test_coords->setUniverse(0);
//	Cell* result = geometry.findCell(test_coords);
//	log_printf(DEBUG, "Found cell for point %s: %s", test_coords->toString().c_str(), result->toString().c_str());
//	result = geometry.findNextCell(test_coords, M_PI/2);
//	log_printf(DEBUG, "Found next cell for point %s: %s", test_coords->toString().c_str(), result->toString().c_str());

	log_printf(INFO, "Program complete");
}
