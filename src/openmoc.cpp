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

// TODO: This is just stubbed out for now
int main(int argc, const char **argv) {
	Options opts = Options(argc, argv);
	log_printf(NORMAL, "Starting OpenMOC...\n");

// Variable is unused, causes warning with -Wall
#if 0
	Geometry* geometry = new Geometry();
	TrackGenerator* trackGenerator = new TrackGenerator();
#endif

	Parser* parser = new Parser(&opts);
	parser->parseMaterials();
	parser->parseGeometry();
}
