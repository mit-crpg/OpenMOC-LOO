/*
 * openmoc.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 *
 *  This file defines the executable for OpenMOC.
 *
 */

#include <stdlib.h>

// TODO: This is just stubbed out for now
int main(int argc, char *argv[]) {

	// The two command line arguments: materials.xml and geometry.xml
	if (argc != 3)
		std::cout << "OpenMOC takes 2 input arguments, but only " << argc - 1
				<< " were given. " << std::endl;

	Geometry* geometry = new Geometry();
	Parser* parser = new Parser(argv[0], argv[1]);

	parser->parseMaterials();
	parser->parseGeometry();

}

