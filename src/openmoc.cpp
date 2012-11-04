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
#include "configurations.h"
#include "Plotter.h"
#include "Cmfd.h"

// FIXME: These should be removed when main() is properly implemented
#pragma GCC diagnostic ignored "-Wunused"
#pragma GCC diagnostic ignored "-Wunused-variable"

int main(int argc, const char **argv) {
	log_printf(NORMAL, "Starting OpenMOC...");

	double k_eff;
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
	Geometry geometry(&parser);
	timer.stop();
	timer.recordSplit("Geometry initialization");

	/* Print out geometry to console if requested at runtime*/
	if (opts.dumpGeometry())
		geometry.printString();

	/* Compress cross-sections if requested at runtime */
	if (opts.compressCrossSections())
		geometry.compressCrossSections();

	Plotter plotter(&geometry, opts.getBitDimension(), opts.getExtension(),
			opts.plotSpecs(), opts.plotFluxes(), opts.plotCurrent());

#if !CMFD_ONLY
	/* Initialize the trackgenerator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
				       opts.getTrackSpacing());
#endif

	/* create CMFD Mesh */
#if CMFD_ACCEL
		geometry.makeCMFDMesh(opts.getNumAzim());
		if (opts.plotSpecs()){
			plotter.plotCMFDMesh(geometry.getMesh());
		}
#endif

#if !CMFD_ONLY
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
#endif

	/* create CMFD class */
	Cmfd cmfd(&geometry, &plotter);

	/* Fixed source iteration to solve for k_eff */
#if !CMFD_ONLY
	Solver solver(&geometry, &track_generator, &plotter, &cmfd, opts.updateFlux());
	timer.reset();
	timer.start();
	k_eff = solver.computeKeff(MAX_ITERATIONS);
	timer.stop();
	timer.recordSplit("Fixed source iteration");
#else
	timer.reset();
	timer.start();
	k_eff = cmfd.computeKeff(MAX_ITERATIONS);
	timer.stop();
	timer.recordSplit("cmfd solve");
#endif

	log_printf(RESULT, "k_eff = %f", k_eff);

	/* Print timer splits to console */
	log_printf(NORMAL, "Program complete");
	timer.printSplits();
}
