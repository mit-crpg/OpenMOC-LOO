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

	// FIXME: Put this inside the parser to hide it from the main program
	/* Initialize the parser */
	Parser* parser = new Parser(&opts);
	parser->each_surface([geometry](Surface *s) -> void
			     {
				     geometry->addSurface(s);
				     return;
			     });
	parser->each_cell([geometry](Cell *c) -> void
			  {
				  geometry->addCell(c);
				  return;
			  });
	

	/* Print out geometry to console */
	log_printf(INFO, geometry->toString());

	/* Initialize the trackgenerator */
	TrackGenerator* trackGenerator = new TrackGenerator(geometry,
			opts.getTrackSpacing(), opts.getTrackSpacing());
}
