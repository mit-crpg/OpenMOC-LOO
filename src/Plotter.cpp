/*
 * Plotter.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: samuelshaner
 */

#include "Plotter.h"


/**
 * Plotting constructor
 * @param geom a pointer to a geometry object
 */
Plotter::Plotter(Geometry* geom, const int bitDim, std::string extension, bool specs, bool fluxes, bool netCurrent) {

	/* extension for plotting files */
	_extension = extension;
	_geom = geom;
	_specs = specs;
	_fluxes = fluxes;
	_net_current = netCurrent;

	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;

	/* set pixel dimensions of plots */
	_bit_length_x = int (bitDim*ratio);
	_bit_length_y = bitDim;

	/* make multiplier for translating geometry coordinates
	 * to bitmap coordinates.
	 */
	_x_pixel = double(_bit_length_x)/_width;
	_y_pixel = double(_bit_length_y)/_height;
}

/**
 * Plotter destructor frees memory for all tracks
 */
Plotter::~Plotter() {
}

/**
 * Return the x bit length of plots
 * @return the x bit length of plots
 */
int Plotter::getBitLengthX(){
	return _bit_length_x;
}

/**
 * Return the y bit length of plots
 * @return the y bit length of plots
 */
int Plotter::getBitLengthY(){
	return _bit_length_y;
}

/**
 * Return the x dimension pixel multiplier for plots
 * @return the x dimension pixel multiplier for plots
 */
double Plotter::getXPixel(){
	return _x_pixel;
}

/**
 * Return the y dimension pixel multiplier for plots
 * @return the y dimension pixel multiplier for plots
 */
double Plotter::getYPixel(){
	return _y_pixel;
}

/**
 * Return the plotting extension
 * @return the plotting extension
 */
std::string Plotter::getExtension(){
	return _extension;
}

/**
 * Return boolean to decide whether to plot specs
 * @return boolean to decide whether to plot specs
 */
bool Plotter::plotSpecs(){
	return _specs;
}

/**
 * Return boolean to decide whether to plot fluxes
 * @return boolean to decide whether to plot fluxes
 */
bool Plotter::plotFlux(){
	return _fluxes;
}

/**
 * Return boolean to decide whether to plot net current
 * @return boolean to decide whether to plot net current
 */
bool Plotter::plotCurrent(){
	return _net_current;
}


/* PLOTTING HELPER FUNCTIONS */
/* These functions manipulate data and call the generic plotting functions
 */

/**
 * Convert an x value our from geometry coordinates to Bitmap coordinates.
 */
double Plotter::convertToGeometryX(int x){
	return double((x - _bit_length_x/2.0) / _x_pixel);
}

/**
 * Convert an y value our from geometry coordinates to Bitmap coordinates.
 */
double Plotter::convertToGeometryY(int y){
	return double(-(y - _bit_length_y/2.0) / _y_pixel);
}

/**
 * Convert an x value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToPixelX(double x){
	return int(((x + _width / 2.0) / _width) * _bit_length_x);
}

/**
 * Convert an y value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToPixelY(double y){
	return int(((y + _height / 2.0) / _height) * _bit_length_y);
}



/**
 * plot given track and numReflect reflected tracks
 */
void Plotter::plotTracksReflective(Track* track, int numReflect){
	log_printf(NORMAL, "Writing tracks reflect bitmap...");

	/* initialize variables */
	double sin_phi, cos_phi, phi;
	Track *track2;
	bool get_out = TRUE;
	double x0, y0, x1, y1;
	int num_segments;

	/* create BitMap for plotting */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _width;
	bitMap->geom_y = _height;
	bitMap->color_type = RANDOM;


	/* loop through tracks and write to pixMap array */
	for (int i = 0; i < (numReflect + 1); i++){
		log_printf(DEBUG, "plotting reflective track: %d", i);
		log_printf(DEBUG, "x_start, y_start, x_end, y_end: %f, %f, %f, %f",
				track->getStart()->getX(),track->getStart()->getY(),
				track->getEnd()->getX(),track->getEnd()->getY());
		phi = track->getPhi();
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		x0 = track->getStart()->getX();
		y0 = track->getStart()->getY();
		num_segments = track->getNumSegments();
		for (int k=0; k < num_segments; k++){
			x1 = x0 + cos_phi*track->getSegment(k)->_length;
			y1 = y0 + sin_phi*track->getSegment(k)->_length;
			drawLine(bitMap, x0, y0, x1, y1, track->getSegment(k)->_region_id);
			x0 = x1;
			y0 = y1;
		}

		/* Get next track */
		track2 = track;
		if (get_out){
			track = track->getTrackOut();
		}
		else {
			track = track->getTrackIn();
		}

		/*determine whether we want TrackIn or TrackOut of next track */
		if (track->getTrackOut() == track2){
			get_out = FALSE;
		}
		else{
			get_out = TRUE;
		}
	}

	/* plot pixMap array */
	plot(bitMap, "reflect", _extension);

	/* release memory */
	deleteBitMap(bitMap);
}

/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion) in either a
 * scaled Magick Plot or a silo plot.
 */
void Plotter::makeRegionMap(int* pixMapFSR, int* pixMap, int* regionMap){

	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = (int)regionMap[pixMapFSR[y * _bit_length_x + x]];
		}
	}
}


/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion) in either a
 * scaled Magick Plot or a silo plot.
 */
void Plotter::makeRegionMap(int* pixMapFSR, float* pixMap, double* regionMap){

	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = (float)regionMap[pixMapFSR[y * _bit_length_x + x]];
		}
	}
}


/**
 * Loops over pixels in pixMap array and finds the corresponding
 * FSR in geometry.
 */
void Plotter::makeFSRMap(int* pixMap){
	log_printf(NORMAL, "Generating FSR maps...");

	/* initialize variables */
	double x_global;
	double y_global;

	/* loop over pixels */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);

			log_printf(DEBUG, "finding cell for bit x: %i, bit y: %i, "
					"global x: %f, global %f", x, y, x_global, y_global);

			/* create point located in universe 0 */
			LocalCoords point(x_global,y_global);
			point.setUniverse(0);

			/* find which cell the point is in */
			_geom->findCell(&point);

			/* Store FSR id in pixMap */
			pixMap[y * _bit_length_x + x] = _geom->findFSRId(&point);

			/* Remove all allocated localcoords */
			point.prune();
		}
	}
}

/* plot CMFD mesh */
void Plotter::plotCMFDMesh(Mesh* mesh){
	log_printf(NORMAL, "plotting CMFD mesh...");

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _width;
	bitMap->geom_y = _height;
	bitMap->color_type = SCALED;

	double x_global;
	double y_global;

	/* find meshCell for each pixel */
	for (int y=0;y < _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = mesh->findMeshCell(x_global, y_global);
		}
	}

	plot(bitMap, "cmfd", _extension);
	deleteBitMap(bitMap);
}

void Plotter::plotNetCurrents(Mesh* mesh){
	log_printf(NORMAL, "plotting net currents...");

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _width;
	bitMap->geom_y = _height;
	bitMap->color_type = SCALED;

	double x_global;
	double y_global;

	/* find meshCell for each pixel */
	for (int y=0;y < _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = mesh->findMeshCell(x_global, y_global);
		}
	}

	double x_mid, y_mid;
	MeshCell* meshCell;
	std::stringstream text_stream;
	std::string text;

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++){
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++){
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){

				/* SIDE 0 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[0]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(0)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid + 20, y_mid + 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 1 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[1]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(1)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid + 10 * (NUM_ENERGY_GROUPS - group + 1));
				text.clear();

				/* SIDE 2 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[2]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(2)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 80, y_mid + 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 3 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[3]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(3)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid - 10 * (group + 1));
				text.clear();
			}

		}
	}


	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, "cmfd_current", _extension);
	}
	else{
		log_printf(WARNING, "Currents can only be plotted in tiff, jpg, and png. Plotting CMFD currents as png...");
		plot(bitMap, "cmfd_current", "png");
	}

	deleteBitMap(bitMap);
}


void Plotter::plotSurfaceFlux(Mesh* mesh){
	log_printf(NORMAL, "plotting surface flux...");

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _width;
	bitMap->geom_y = _height;
	bitMap->color_type = SCALED;

	double x_global;
	double y_global;

	/* find meshCell for each pixel */
	for (int y=0;y < _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = mesh->findMeshCell(x_global, y_global);
		}
	}

	double x_mid, y_mid;
	MeshCell* meshCell;
	std::stringstream text_stream;
	std::string text;

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++){
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++){
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){

				/* SIDE 0 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[0]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(0)->getFlux(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid + 20, y_mid + 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 1 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[1]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(1)->getFlux(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid + 10 * (NUM_ENERGY_GROUPS - group + 1));
				text.clear();

				/* SIDE 2 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[2]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(2)->getFlux(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 80, y_mid + 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 3 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[3]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(3)->getFlux(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid - 10 * (group + 1));
				text.clear();
			}

		}
	}


	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, "cmfd_flux", _extension);
	}
	else{
		log_printf(WARNING, "Suface fluxes can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
		plot(bitMap, "cmfd_flux", "png");
	}

	deleteBitMap(bitMap);
}




