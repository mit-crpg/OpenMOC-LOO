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
#include "Parser.h"
#include "log.h"

// TODO: This is just stubbed out for now
int main(int argc, char *argv[]) {

	log_printf(NORMAL, "Starting OpenMOC...\n");


#if 0
	// The two command line arguments: materials.xml and geometry.xml
	if (argc != 3)
		log_printf(WARNING, "OpenMOC takes 2 input arguments, but only %d were given.\n", argc-1);
#endif

// Variable is unused, causes warning with -Wall
#if 0
	Geometry* geometry = new Geometry();
#endif

	Parser* parser = new Parser("../xml-sample/1/geometry.xml");
	parser->parseMaterials();
	parser->parseGeometry();
}
