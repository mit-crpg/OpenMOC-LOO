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
#include "Geometry.h"
#include "Surface.h"

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

	/* Initialize the geometry with surfaces, cells & materials */
	Geometry geometry(opts.getNumSectors(), opts.getNumRings(),
			  opts.getSectorOffset(), &parser);

	/* Prints out all universes */
//	geometry.checkUniverse();

	/* Plot geometry */
//	Plotting plotter(&geometry, opts.getBitDimension());

//	/* Print out geometry to console if requested at runtime*/
	if (opts.dumpGeometry())
		geometry.printString();

	/* Adjust the indices for each geometry class to use uids */
//	geometry.adjustKeys();

	/* Generate the neighbor cells for each surface in geometry */
//	geometry.buildNeighborsLists();

	/* Initialize the trackgenerator */
	TrackGenerator track_generator(&geometry, opts.getNumAzim(),
				       opts.getTrackSpacing(), opts.getBitDimension());

	track_generator.generateTracks();
	track_generator.makeReflective();
	track_generator.segmentize();

	log_printf(NORMAL, "Program complete");
}
