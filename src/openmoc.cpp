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

	/* Create an options class to parse command line options */
	Options opts(argc, argv);

	/* initialize Petsc */
	int petsc_err = 0;
	CHKERRQ(petsc_err);
	PetscInitialize(&(opts.extra_argc), &(opts.extra_argv), (char*)0, NULL);
	CHKERRQ(petsc_err);

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
					opts.plotSpecs(), opts.plotFluxes(), opts.plotCurrent(), 
					opts.plotDiffusion(), opts.plotKeff(), opts.plotQuadFlux());

	/* Initialize track generator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
				       opts.getTrackSpacing());

	/* Tell geometry whether CMFD is on/off */
	geometry.setCmfd(opts.getCmfd());
	geometry.setLoo(opts.getLoo());

	/* Make CMFD mesh */
	if (opts.getCmfd() || opts.getLoo())
		geometry.makeCMFDMesh(geometry.getMesh(), opts.getNumAzim(), 
							  opts.getGroupStructure(), opts.getPrintMatrices(),
							  opts.getCmfdLevel());

	/* make FSR map for plotting */
	if (opts.plotCurrent() || opts.plotDiffusion() || opts.plotFluxes() || 
		opts.plotSpecs())
		plotter.makeFSRMap();

	/* plot CMFD mesh */
	if (opts.plotSpecs() && opts.getCmfd())
	{
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

	/* Create CMFD class */
	Cmfd cmfd(&geometry, &plotter, geometry.getMesh(), 
			  opts.getCmfd(), opts.getLoo(), opts.getLoo1(), opts.getLoo2(),
			  opts.getDiffusionCorrection(), opts.plotProlongation(), 
			  opts.getL2NormConvThresh(), opts.getDampFactor(),
			  &track_generator);

	/* Creat Solver class */
	Solver solver(&geometry, &track_generator, &plotter, &cmfd, &opts);

	cmfd.setFSRs(solver.getFSRs());

	/* Solve steady state problem */
	timer.reset();
	timer.start();
	k_eff = solver.kernel(MAX_ITERATIONS);
	timer.stop();
	timer.recordSplit("Fixed source iteration");
	log_printf(RESULT, "k_eff = %.10f", k_eff);

	/* Finalize petsc */
	PetscFinalize();
	CHKERRQ(petsc_err);

	/* Print timer splits to console */
	log_printf(NORMAL, "Program complete");
	timer.printSplits();
}
