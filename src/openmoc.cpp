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
#include "Options.h"
#include "Solver.h"
#include "Timer.h"
#include "log.h"
#include "Plotter.h"

// FIXME: These should be removed when main() is properly implemented
#pragma GCC diagnostic ignored "-Wunused"
#pragma GCC diagnostic ignored "-Wunused-variable"

// TODO: This is just stubbed out for now
int main(int argc, const char **argv) {
	log_printf(NORMAL, "Starting OpenMOC...");

	Timer timer;

	/* Create an options class to parse command line options */
	Options opts(argc, argv);

	/* Set the verbosity */
	log_setlevel(opts.getVerbosity());

	/* Initialize the parser and time the parser */
	timer.start();
	Parser parser(&opts);
	timer.stop();
	timer.recordSplit("Parsing input files");

	/* Initialize the geometry with surfaces, cells & materials */
	timer.reset();
	timer.start();
	Geometry geometry(opts.getNumSectors(), opts.getNumRings(),
			  opts.getSectorOffset(), &parser);
	timer.stop();
	timer.recordSplit("Geomery initialization");

//	/* Print out geometry to console if requested at runtime*/
	if (opts.dumpGeometry())
		geometry.printString();

	Plotter plotter(&geometry, opts.getBitDimension(), opts.getExtension());

	/* Adjust the indices for each geometry class to use uids */
//	geometry.adjustKeys();

	/* Generate the neighbor cells for each surface in geometry */
//	geometry.buildNeighborsLists();

	/* Initialize the trackgenerator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
				       opts.getTrackSpacing());

	/* Generate tracks */
	timer.reset();
	timer.start();
	track_generator.generateTracks();
	track_generator.makeReflective();
	timer.stop();
	timer.recordSplit("Generating tracks");

	/* Segment tracks */
	timer.reset();
	timer.start();
	track_generator.segmentize();
	timer.stop();
	timer.recordSplit("Segmenting tracks");

	Solver solver(&geometry, &track_generator);
	solver.zeroTrackFluxes();
	solver.oneFSRFluxes();

	/* Print timer splits to console */
	log_printf(NORMAL, "Program complete");
	timer.printSplits();
}
