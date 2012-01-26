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

	log_printf(NORMAL, "Starting OpenMOC...\n");

	/* Create an options class to parse command line options */
	Options opts = Options(argc, argv);

	/* Set the verbosity */
	log_setlevel(opts.getVerbosity());

	/* Create an empty geometry */
	Geometry* geometry = new Geometry();

	/* Set the geomery's number of rings, sectors, and sector offset */
	geometry->setNumRings(opts.getNumRings());
	geometry->setNumSectors(opts.getNumSectors());
	geometry->setSectorOffset(opts.getSectorOffset());


	// FIXME: PUT THIS INSIDE THE PARSER TO HIDE IT FROM THE MAIN PROGRAM
	/* We only want general functional calls inside the main program
	 * The logic going on here is too specific and should be hidden from main
	 */
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
	

	parser->each_lattice([geometry](Lattice *l) -> void
			  {
				  geometry->addLattice(l);
				  return;
			  });
	
	/* Adjust the indices for each geometry class to use uids */
	geometry->adjustKeys();

	/* Generate the neighbor cells for each surface in geometry */
	geometry->buildNeighborsLists();

	/* Print out geometry to console */
	log_printf(INFO, geometry->toString());

	/* Initialize the trackgenerator */
	TrackGenerator* trackGenerator = new TrackGenerator(geometry,
			opts.getTrackSpacing(), opts.getTrackSpacing());
}
