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
Plotter::Plotter(Geometry* geom, const int bitDim, std::string extension, 
				 bool specs, bool fluxes, bool netCurrent, bool plotDiffusion, 
				 bool plotKeff, bool plotQuadFluxFlag) {

	/* extension for plotting files */
	_extension = extension;
	_geom = geom;
	_specs = specs;
	_fluxes = fluxes;
	_net_current = netCurrent;
	_plot_diffusion = plotDiffusion;
	_plot_keff = plotKeff;
	_plot_quad_flux_flag = plotQuadFluxFlag;
	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;


	/* set pixel dimensions of plots */
	_bit_length_x = int (bitDim*ratio);
	_bit_length_y = bitDim;

	_FSR_map = new int[_bit_length_x*_bit_length_y];

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
 * Return boolean to decide whether to plot cmfd flux at each step
 * @return boolean to decide whether to plot cmfd flux at each step
 */
bool Plotter::plotDiffusion(){
	return _plot_diffusion;
}

/**
 * Return boolean to decide whether to plot keff
 * @return boolean to decide whether to plot keff
 */
bool Plotter::plotKeff(){
	return _plot_keff;
}

/**
 * Return boolean to decide whether to plot net current
 * @return boolean to decide whether to plot net current
 */
bool Plotter::plotCurrent(){
	return _net_current;
}

/**
 * Return boolean to decide whether to plot quadrature flux
 * @return boolean to decide whether to plot quadrature flux
 */
bool Plotter::plotQuadFluxFlag(){
	return _plot_quad_flux_flag;
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
	return int((1.0 - (y + _height / 2.0) / _height) * _bit_length_y);
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
void Plotter::makeFSRMap(){
	log_printf(NORMAL, "Generating FSR maps...");

	/* initialize variables */
	double x_global;
	double y_global;

	/* loop over pixels */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);

			/* create point located in universe 0 */
			LocalCoords point(x_global,y_global);
			point.setUniverse(0);

			/* find which cell the point is in */
			_geom->findCell(&point);

			/* Store FSR id in pixMap */
			_FSR_map[y * _bit_length_x + x] = _geom->findFSRId(&point);

			/* Remove all allocated localcoords */
			point.prune();
		}
	}
}

int *Plotter::getFSRMap(){
	return _FSR_map;
}

void Plotter::copyFSRMap(int *pixels){

	for (int x = 0; x < _bit_length_x; x++){
		for (int y = 0; y < _bit_length_y; y++){
			pixels[y * _bit_length_x + x] = _FSR_map[y * _bit_length_x + x];
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
	bitMap->color_type = RANDOM;

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

	BitMap<int>* bitMap2 = new BitMap<int>;
	bitMap2->pixel_x = _bit_length_x;
	bitMap2->pixel_y = _bit_length_y;
	initialize(bitMap2);
	bitMap2->geom_x = _width;
	bitMap2->geom_y = _height;
	bitMap2->color_type = SCALED;

	double x_global;
	double y_global;

	/* find meshCell for each pixel */
	for (int y=0;y < _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = mesh->findMeshCell(x_global, y_global);
			bitMap2->pixels[y * _bit_length_x + x] = mesh->findMeshCell(x_global, y_global);
		}
	}

	double x_mid, y_mid;
	MeshCell* meshCell;
	std::stringstream text_stream;
	std::string text;
	double current0, current1, current2, current3;

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++){
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++){
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);
			current0 = 0;
			current1 = 0;
			current2 = 0;
			current3 = 0;
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
				current0 += meshCell->getMeshSurfaces(0)->getCurrent(group);
				current1 += meshCell->getMeshSurfaces(1)->getCurrent(group);
				current2 += meshCell->getMeshSurfaces(2)->getCurrent(group);
				current3 += meshCell->getMeshSurfaces(3)->getCurrent(group);

				/* SIDE 0 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[0]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(0)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid + 20, y_mid - 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 1 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[1]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(1)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid - 10 * (NUM_ENERGY_GROUPS - group));
				text.clear();

				/* SIDE 2 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX(meshCell->getBounds()[2]);
				y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(2)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 80, y_mid - 10 * (NUM_ENERGY_GROUPS / 2.0 - group));
				text.clear();

				/* SIDE 3 */
				/* get midpoint of mesh surface */
				x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
				y_mid = convertToPixelY(meshCell->getBounds()[3]);

				/* create string and draw on bitMap */
				text_stream << meshCell->getMeshSurfaces(3)->getCurrent(group);
				text = text_stream.str();
				text_stream.str("");
				drawText(bitMap, text, x_mid - 20, y_mid + 10 * (group + 1.5));
				text.clear();
			}

			/* SIDE 0 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[0]);
			y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

			/* create string and draw on bitMap */
			text_stream << "tally: " << current0;
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap2, text, x_mid + 20, y_mid);
			text.clear();

			/* SIDE 1 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[1]);

			/* create string and draw on bitMap */
			text_stream << "tally: " << current1;
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap2, text, x_mid - 20, y_mid - 10);
			text.clear();

			/* SIDE 2 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[2]);
			y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

			/* create string and draw on bitMap */
			text_stream << "tally: " << current2;
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap2, text, x_mid - 80, y_mid);
			text.clear();

			/* SIDE 3 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[3]);

			/* create string and draw on bitMap */
			text_stream << "tally: " << current3;
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap2, text, x_mid - 20, y_mid + 15);
			text.clear();


		}
	}

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, "current", _extension);
		plot(bitMap2, "current_1G", _extension);
	}
	else{
		log_printf(WARNING, "Currents can only be plotted in tiff, jpg, and png. Plotting CMFD currents as png...");
		plot(bitMap, "current", "png");
		plot(bitMap2, "current_1G", "png");
	}

	deleteBitMap(bitMap);
	deleteBitMap(bitMap2);
}

void Plotter::plotDHats(Mesh* mesh, int iter_num){
	log_printf(NORMAL, "plotting D Hats...");

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

			/* SIDE 0 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[0]);
			y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

			/* create string and draw on bitMap */
			text_stream << "DHat0: " << meshCell->getMeshSurfaces(0)->getDHat()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid + 20, y_mid);
			text.clear();
			text_stream << "DTilde: " << meshCell->getMeshSurfaces(0)->getDTilde()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid + 20, y_mid - 15.0);
			text.clear();


			/* SIDE 1 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[1]);

			/* create string and draw on bitMap */
			text_stream << "DHat1: " << meshCell->getMeshSurfaces(1)->getDHat()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 20, y_mid - 20.0);
			text.clear();
			text_stream << "DTilde: " << meshCell->getMeshSurfaces(1)->getDTilde()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 20, y_mid - 35.0);
			text.clear();

			/* SIDE 2 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[2]);
			y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

			/* create string and draw on bitMap */
			text_stream << "DHat2: " << meshCell->getMeshSurfaces(2)->getDHat()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 80, y_mid);
			text.clear();
			text_stream << "DTilde: " << meshCell->getMeshSurfaces(2)->getDTilde()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 80, y_mid - 15.0);
			text.clear();

			/* SIDE 3 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[3]);

			/* create string and draw on bitMap */
			text_stream << "DHat3: " << meshCell->getMeshSurfaces(3)->getDHat()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 20, y_mid + 35.0);
			text.clear();
			text_stream << "DTilde: " << meshCell->getMeshSurfaces(3)->getDTilde()[0];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 20, y_mid + 20.0);
			text.clear();
		}
	}

	std::stringstream string;
	string << "cmfd_dhats_i_" << iter_num;
	std::string title_str = string.str();

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, title_str, _extension);
	}
	else{
		log_printf(WARNING, "D Hats can only be plotted in tiff, jpg, and png. Plotting CMFD D Hats as png...");
		plot(bitMap, title_str, "png");
	}

	deleteBitMap(bitMap);
}

void Plotter::plotGeometry(Mesh* mesh)
{
	log_printf(INFO, "plotting geometry for debugging ...");

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
	for (int y = 0; y < _bit_length_y; y++){
		for (int x =  0 ; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = 
				mesh->findMeshCell(x_global, y_global);
		}
	}

	double x_mid, y_mid;
	MeshCell* meshCell;
	std::stringstream text_stream;
	std::string text;

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++)
{
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++)
		{
			/* FIXME: should be cellX * mesh->getCellHeight() + cellY? */
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);

			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] + 
									meshCell->getBounds()[2]) / 2.0);

			y_mid = convertToPixelY((meshCell->getBounds()[1] + 
									 meshCell->getBounds()[3]) / 2.0);

			/* create string and draw on bitMap */
			text_stream << "Mesh Cell Index: (" << 
				x_mid << y_mid << ")";
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid, y_mid);
			text.clear();
		}
	}

	std::stringstream string;
	string << "mesh-cell-index";
	std::string title_str = string.str();

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, title_str, _extension);
	}
	else{
		log_printf(WARNING, "Mesh cell index can only be plotted in tiff, jpg,"
				   " and png. Plotting LOO Quad Flux as png...");
		plot(bitMap, title_str, "png");
	}

	deleteBitMap(bitMap);
}


void Plotter::plotQuadFlux(Mesh* mesh, int iter_num){
	log_printf(NORMAL, "plotting quadrature fluxes for %d-th iteration ...", 
			   iter_num);

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
	for (int y = 0; y < _bit_length_y; y++){
		for (int x =  0 ; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = 
				mesh->findMeshCell(x_global, y_global);
		}
	}

	double x_mid, y_mid;
	double x_min = 0.0, y_min = 0.0;
	double x_max, y_max;
	MeshCell* meshCell;
	std::stringstream text_stream;
	std::string text;

	text_stream.precision(10);
	
	double scale = 1;
	/* y_max correponds to (0, max_y)'s bottom surface (surface 3) */
	y_max = convertToPixelX(mesh->getCells(mesh->getCellWidth() - 1)
							->getBounds()[3]);
	/* x_max corresponds to (max_x, 0)'s right surface (surface 2) */
	x_max = y_max;
		/*convertToPixelX(mesh->getCells((mesh->getCellHeight() - 1)
										   * mesh->getCellWidth())
										   ->getBounds()[2]); */
	log_printf(NORMAL, "x_max = %f, y_max = %f", x_max, y_max);

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++){
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++){
			/* FIXME: should be cellX * mesh->getCellHeight() + cellY? */
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);

			/* SIDE 0 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[0]) - 170.0;
			y_mid = convertToPixelY((meshCell->getBounds()[1] + 
									 meshCell->getBounds()[3]) / 2.0);

			
			if (x_mid < x_min)
				x_mid = convertToPixelX(meshCell->getBounds()[0]) + 20.0;

			/* create string and draw on bitMap */
			/* getQuadCurrent(group, index) */
			text_stream << "surf[0].flux[1]: " << 
				scale * meshCell->getMeshSurfaces(0)->getQuadCurrent(0,1);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid, y_mid - 20.0);
			text.clear();

			text_stream << "surf[0].flux[0]: " << 
				scale *meshCell->getMeshSurfaces(0)->getQuadCurrent(0,0);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid, y_mid + 20.0);
			text.clear();

			/* SIDE 2 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX(meshCell->getBounds()[2]) + 40.0;
			y_mid = convertToPixelY((meshCell->getBounds()[1] 
									 + meshCell->getBounds()[3]) / 2.0);

			if (x_mid > x_max)
				x_mid = convertToPixelX(meshCell->getBounds()[2]) - 180.0;

			/* create string and draw on bitMap */
			text_stream << "surf[2].flux[0]: " << 
				scale *meshCell->getMeshSurfaces(2)->getQuadCurrent(0,0);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid, y_mid - 20);
			text.clear();

			text_stream << "surf[2].flux[1]: " << 
				scale *meshCell->getMeshSurfaces(2)->getQuadCurrent(0,1);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid, y_mid + 20);
			text.clear();


			/* SIDE 1 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] 
									 + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[1]) + 20.0;

			if (y_mid > y_max)
				y_mid -= 40.0;

			/* create string and draw on bitMap */
			text_stream << "surf[1].flux[0]: " << 
				scale *meshCell->getMeshSurfaces(1)->getQuadCurrent(0,0);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 200, y_mid);
			text.clear();

			text_stream << "surf[1].flux[1]: " << 
				scale *meshCell->getMeshSurfaces(1)->getQuadCurrent(0,1);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid + 20, y_mid);
			text.clear();

			/* SIDE 3 */
			/* get midpoint of mesh surface */
			x_mid = convertToPixelX((meshCell->getBounds()[0] 
									 + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY(meshCell->getBounds()[3]) - 20.0;

			if (y_mid < y_min)
				y_mid += 40.0;

			/* create string and draw on bitMap */
			text_stream << "surf[3].flux[1]: " << 
				scale *meshCell->getMeshSurfaces(3)->getQuadCurrent(0,1);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 200, y_mid);
			text.clear();

			text_stream << "surf[3].flux[0]: " << 
				scale *meshCell->getMeshSurfaces(3)->getQuadCurrent(0,0);
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid + 20, y_mid);
			text.clear();
		}
	}

	std::stringstream string;
	string << "loo_quad_flux_iter_" << iter_num;
	std::string title_str = string.str();

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, title_str, _extension);
	}
	else{
		log_printf(WARNING, "Quad Flux can only be plotted in tiff, jpg, and"
				   " png. Plotting LOO Quad Flux as png...");
		plot(bitMap, title_str, "png");
	}

	deleteBitMap(bitMap);
}


void Plotter::plotXS(Mesh* mesh, int iter_num){
	log_printf(NORMAL, "plotting cross sections...");

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

	int e = 6;

	/* plot mesh currents next to surface */
	for (int cellY = 0; cellY < mesh->getCellHeight(); cellY++){
		for (int cellX = 0; cellX < mesh->getCellWidth(); cellX++){
			meshCell = mesh->getCells(cellY * mesh->getCellWidth() + cellX);

			/* find middle of cell */
			x_mid = convertToPixelY((meshCell->getBounds()[0] + meshCell->getBounds()[2]) / 2.0);
			y_mid = convertToPixelY((meshCell->getBounds()[1] + meshCell->getBounds()[3]) / 2.0);

			/* Cell Index */
			text_stream << "Cell Index: ("<< cellX << " , " << cellY << ")";
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 50, y_mid - 30);
			text.clear();			

			/* Sigma A */
			text_stream << "SigmaA: " << meshCell->getSigmaA()[e];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 50, y_mid + 30);
			text.clear();

			/* Sigma S */
			text_stream << "SigmaS: " << meshCell->getSigmaS()[e];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 50, y_mid + 15);
			text.clear();

			/* Nu Sigma F */
			text_stream << "NuSigmaF: " << meshCell->getNuSigmaF()[e];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 50, y_mid);
			text.clear();

			/* Diffusivity */
			text_stream << "Diffusivity: " << meshCell->getDiffusivity()[e];
			text = text_stream.str();
			text_stream.str("");
			drawText(bitMap, text, x_mid - 50, y_mid - 15);
			text.clear();
		}
	}

	std::stringstream string;
	string << "xs_iter_" << iter_num;
	std::string title_str = string.str();

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, title_str, _extension);
	}
	else{
		log_printf(WARNING, "Cross sections can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
		plot(bitMap, title_str, "png");
	}

	deleteBitMap(bitMap);
}


void Plotter::plotCMFDflux(Mesh* mesh, std::string title1, int iter_num){
	log_printf(NORMAL, "plotting flux...");

	/* set up bitMap */
	BitMap<double>* bitMap = new BitMap<double>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _width;
	bitMap->geom_y = _height;
	bitMap->color_type = SCALED;

	double x_global;
	double y_global;
	std::stringstream string;
	std::stringstream num;
	std::string title_str;

	int ng = 1;
	if (mesh->getMultigroup() == true){
		ng = NUM_ENERGY_GROUPS;
	}

	/* PLOT OLD FLUX */

	for (int e = 0; e < ng; e++){

		/* find meshCell for each pixel */
		for (int y=0;y < _bit_length_y; y++){
			for (int x = 0; x < _bit_length_x; x++){
				x_global = convertToGeometryX(x);
				y_global = convertToGeometryY(y);
				bitMap->pixels[y * _bit_length_x + x] = mesh->getCells(mesh->findMeshCell(x_global, y_global))->getOldFlux()[e];
			}
		}


		num.str("");
		num << std::setw(5) << std::setfill('0') << iter_num;
		string.str("");
		string << title1 << "_old_i_" << num.str() << "_g_" << e+1;
		title_str = string.str();


		/* create filename with correct extension */
		if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
			plot(bitMap, title_str, _extension);
		}
		else{
			log_printf(WARNING, "CMFD flux can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
			plot(bitMap, title_str, "png");
		}
	}

	/* PLOT NEW FLUX */

	for (int e = 0; e < ng; e++){

		/* find meshCell for each pixel */
		for (int y=0;y < _bit_length_y; y++){
			for (int x = 0; x < _bit_length_x; x++){
				x_global = convertToGeometryX(x);
				y_global = convertToGeometryY(y);
				bitMap->pixels[y * _bit_length_x + x] = mesh->getCells(mesh->findMeshCell(x_global, y_global))->getNewFlux()[e];
			}
		}

		string.str("");
		string << title1 << "_new_i_" << num.str() << "_g_" << e+1;
		title_str = string.str();

		/* create filename with correct extension */
		if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
			plot(bitMap, title_str, _extension);
		}
		else{
			log_printf(WARNING, "CMFD flux can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
			plot(bitMap, title_str, "png");
		}

	}

	/* PLOT FLUX DIFFERENCE */

	for (int e = 0; e < ng; e++){

		/* find meshCell for each pixel */
		for (int y=0;y < _bit_length_y; y++){
			for (int x = 0; x < _bit_length_x; x++){
				x_global = convertToGeometryX(x);
				y_global = convertToGeometryY(y);
				bitMap->pixels[y * _bit_length_x + x] = mesh->getCells(mesh->findMeshCell(x_global, y_global))->getNewFlux()[e]
				   - mesh->getCells(mesh->findMeshCell(x_global, y_global))->getOldFlux()[e];
			}
		}

		string.str("");
		string << "res" << "_new_i_" << num.str() << "_g_" << e+1;
		title_str = string.str();

		/* create filename with correct extension */
		if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
			plot(bitMap, title_str, _extension);
		}
		else{
			log_printf(WARNING, "CMFD flux can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
			plot(bitMap, title_str, "png");
		}

	}


	deleteBitMap(bitMap);
}

void Plotter::plotCMFDKeff(Mesh* mesh, int num_iter){
	log_printf(NORMAL, "plotting CMFD keff...");

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _bit_length_x;
	bitMap->geom_y = _bit_length_y;
	bitMap->color_type = BLACKWHITE;

	std::stringstream text_stream;
	std::string text;
	std::string text2;

	/* draw and label axes */
	drawLine(bitMap, _bit_length_x / 10, 10, _bit_length_x / 10, _bit_length_y - 10);
	drawLine(bitMap, 10, 9 * _bit_length_y / 10, _bit_length_x - 10, 9 * _bit_length_y / 10);
	text = "keff";
	drawText(bitMap, text, 30, _bit_length_y / 2 - 30);
	text = "iteration";
	drawText(bitMap, text, _bit_length_x / 2, _bit_length_y - 30);
	text_stream << num_iter;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, 9.5 * _bit_length_x / 10, 9 * _bit_length_y / 10 + 20);


	/* create x axis scale */
	double keff_max = 0, keff_min = 2;
	for (int i = 1; i < num_iter; i++){
		if (mesh->getKeffCMFD(i) > keff_max){
			keff_max = mesh->getKeffCMFD(i);
		}
		if (mesh->getKeffMOC(i) > keff_max){
			keff_max = mesh->getKeffMOC(i);
		}
		if (mesh->getKeffCMFD(i) < keff_min){
					keff_min = mesh->getKeffCMFD(i);
		}
		if (mesh->getKeffMOC(i) < keff_min){
			keff_min = mesh->getKeffMOC(i);
		}
	}

	text_stream << keff_min;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, 9 * _bit_length_y / 10 - 10);
	text_stream << keff_max;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 10);
	double keff_avg = (keff_max - keff_min)/2 + keff_min;
	double keff_range = keff_max - keff_min;
	text_stream << keff_avg;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 2);

	text = "blue";
	text2 = "white";

	/* draw CMFD keff */
	for (int i = 1; i < num_iter; i++){

		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 2 - (mesh->getKeffCMFD(i) - keff_avg) / (keff_range / 2) * (4 * _bit_length_y / 10), 3);

	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30, 3);
	text = "CMFD";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 5);


	text = "black";

	/* draw MOC keff */
	for (int i = 1; i < num_iter; i++){

		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 2 - (mesh->getKeffMOC(i) - keff_avg) / (keff_range / 2) * (4 * _bit_length_y / 10), 3);

	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30 + 15, 3);
	text = "MOC";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 20);


	text = "cmfd_keff";

	/* create filename with correct extension */
	if (_extension == "tiff" || _extension == "jpg" || _extension == "png"){
		plot(bitMap, text, _extension);
	}
	else{
		log_printf(WARNING, "Cross sections can only be plotted in tiff, jpg, and png. Plotting CMFD flux as png...");
		plot(bitMap, text, "png");
	}

	deleteBitMap(bitMap);

}










