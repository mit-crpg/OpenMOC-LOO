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
#include "Transient.h"
#include "configurations.h"
#include "Plotter.h"
#include "Cmfd.h"
#include "petsc.h"
#include "mpi.h"
#include <petscmat.h>

// FIXME: These should be removed when main() is properly implemented
#pragma GCC diagnostic ignored "-Wunused"
#pragma GCC diagnostic ignored "-Wunused-variable"

int main(int argc, char **argv) {
	log_printf(NORMAL, "Starting OpenMOC...");

	double k_eff;
	Timer timer;

	/* initialize Petsc */
	int petsc_err = 0;
	CHKERRQ(petsc_err);
	PetscInitialize(&argc, &argv, 0, 0);
	CHKERRQ(petsc_err);

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

	/* Initialize plotter */
	Plotter plotter(&geometry, opts.getBitDimension(), opts.getExtension(),
			opts.plotSpecs(), opts.plotFluxes(), opts.plotCurrent(), opts.plotDiffusion(), opts.plotKeff());

	/* Initialize track generator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
				       opts.getTrackSpacing());

	/* Make CMFD mesh */
	geometry.makeCMFDMesh(geometry.getMesh(), opts.getNumAzim(), opts.getGroupStructure(), opts.getPrintMatrices(), opts.getCmfdLevel());

	/* make FSR map for plotting */
	if (opts.plotCurrent() || opts.plotDiffusion() || opts.plotFluxes() || opts.plotSpecs())
		plotter.makeFSRMap();

	/* plot CMFD mesh */
	if (opts.plotSpecs()){
		plotter.plotCMFDMesh(geometry.getMesh());
	}

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

	/* create CMFD class */
	Cmfd cmfd(&geometry, &plotter, geometry.getMesh(), opts.updateFlux());

	/* Fixed source iteration to solve for k_eff */
	Solver solver(&geometry, &track_generator, &plotter, &cmfd, opts.updateFlux(), opts.getKeffConvThresh(), opts.computePinPowers(), opts.getCmfd(), opts.getDiffusion());

	if (opts.getTransient()){
		/* Solve the transient problem */
		Transient transient(&geometry, &cmfd, geometry.getMesh(), &solver, &plotter, opts.getTimeEnd(), opts.getTimeStepOuter(), opts.getTimeStepInner());
		timer.reset();
		timer.start();
		transient.solve();
		timer.stop();
		timer.recordSplit("Time dependent problem");
	}
	else{
		/* solve steady state problem */
		timer.reset();
		timer.start();
		k_eff = solver.computeKeff(MAX_ITERATIONS);
		timer.stop();
		timer.recordSplit("Fixed source iteration");
		log_printf(RESULT, "k_eff = %f", k_eff);
	}

	/* Finalize petsc */
	PetscFinalize();
	CHKERRQ(petsc_err);

	/* Print timer splits to console */
	log_printf(NORMAL, "Program complete");
	timer.printSplits();
}
