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

// FIXME: These should be removed when main() is properly implemented
#pragma GCC diagnostic ignored "-Wunused"
#pragma GCC diagnostic ignored "-Wunused-variable"

// TODO: This is just stubbed out for now
int main(int argc, const char **argv) {

	/* Create an options class to parse command line options */
	Options opts = Options(argc, argv);

	/* Set the verbosity */
	log_setlevel(opts.getVerbosity());
	log_printf(NORMAL, "Starting OpenMOC...\n");

	/* Create an empty geometry */
	Geometry* geometry = new Geometry();

	/* Initialize the parser */
	Parser* parser = new Parser(&opts);

	/* Initialize the trackgenerator */
	TrackGenerator* trackGenerator = new TrackGenerator(geometry,
			opts.getTrackSpacing(), opts.getTrackSpacing());

	/* Parse the input geometry and materials files and register input
	 * with the geometry object
	 */
	parser->parseMaterials();
	parser->parseGeometry();
}
