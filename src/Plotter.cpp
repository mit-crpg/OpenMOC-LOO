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
Plotter::Plotter(Geometry* geom, const int bitDim, std::string extension) {

	/* extension for plotting files */
	_extension = extension;
	_geom = geom;

	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;

	/* set pixel dimensions of plots */
	_bit_length_x = int (bitDim*ratio) + 1;
	_bit_length_y = bitDim + 1;

	/* make multiplier for translating geometry coordinates
	 * to bitmap coordinates.
	 */
	_x_pixel = double(_bit_length_x)/_width;
	_y_pixel = double(_bit_length_y)/_height;

	/* make _color_map for plotting */
	_color_map.insert(std::pair<int, std::string>(0,"indigo"));
	_color_map.insert(std::pair<int, std::string>(1,"black"));
	_color_map.insert(std::pair<int, std::string>(2,"blue"));
	_color_map.insert(std::pair<int, std::string>(3,"green"));
	_color_map.insert(std::pair<int, std::string>(4,"magenta"));
	_color_map.insert(std::pair<int, std::string>(5,"orange"));
	_color_map.insert(std::pair<int, std::string>(6,"maroon"));
	_color_map.insert(std::pair<int, std::string>(7,"orchid"));
	_color_map.insert(std::pair<int, std::string>(8,"blue violet"));
	_color_map.insert(std::pair<int, std::string>(9,"crimson"));
	_color_map.insert(std::pair<int, std::string>(10,"salmon"));
	_color_map.insert(std::pair<int, std::string>(11,"gold"));
	_color_map.insert(std::pair<int, std::string>(12, "DarkSlateGray"));
	_color_map.insert(std::pair<int, std::string>(13,"orange red"));
	_color_map.insert(std::pair<int, std::string>(14,"red"));
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


/* GENERIC PLOTTING FUNCTIONS */
/* These functions write pixMap array data to image files */


/**
 * Generic function for plotting int pixMap in png, tiff, or jpg file
 * using Magick++
 */
void Plotter::plotMagick(int* pixMap, std::string type){
	log_printf(NORMAL, "Writing Magick bitmap...");

	/* pixel color */
	int color_int;

	/* create image and open for modification */
	Magick::Image image(Magick::Geometry(_bit_length_x,_bit_length_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,_bit_length_x,_bit_length_y);

	/* Write pixMap array to Magick cache */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			if (pixMap[y * _bit_length_x + x] != -1){
				color_int = pixMap[y * _bit_length_x + x] % 15;
				*(pixels+(y * _bit_length_x + x)) =
						Magick::Color(_color_map.at(color_int));
			}
		}
	}

	/* Sync pixel cache with Magick image */
	pixel_cache.sync();

	/* create filename with correct extension */
	std::stringstream string;
	string << type << "." << _extension;
	std::string title = string.str();

	/* write Magick image to file */
	image.write(title);
}

/**
 * Generic function for plotting float pixMap in png, tiff, or jpg file
 * using Magick++
 */
void Plotter::plotMagick(float* pixMap, std::string type){
	log_printf(NORMAL, "Writing Magick bitmap...");

	/* pixel color */
	int color_int;

	/* create image and open for modification */
	Magick::Image image(Magick::Geometry(_bit_length_x,_bit_length_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,_bit_length_x,_bit_length_y);

	/* Write pixMap array to Magick cache */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			if (pixMap[y * _bit_length_x + x] != -1){
				color_int = int(floor(pixMap[y * _bit_length_x + x]*1E8)) % 15;
				*(pixels+(y * _bit_length_x + x)) = Magick::Color(_color_map.at(color_int));
			}
		}
	}

	/* Sync pixel cache with Magick image */
	pixel_cache.sync();

	/* create filename with correct extension */
	std::stringstream string;
	string << type << "." << _extension;
	std::string title = string.str();

	/* write Magick image to file */
	image.write(title);
}

/**
 * Generic function for plotting double inline RGB pixMap in png, tiff, or jpg file
 * using Magick++
 */
void Plotter::plotMagickScaled(double* pixMapRGB, std::string type){
	log_printf(NORMAL, "Writing Magick bitmap...");

	/* declare variables */
	double red, green, blue;

	/* create image and open for modification */
	Magick::Image image(Magick::Geometry(_bit_length_x,_bit_length_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,_bit_length_x,_bit_length_y);

	/* Write pixMapRGB array to Magick pixel_cache */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			red = pixMapRGB[3 * y * _bit_length_x + 3 * x];
			green = pixMapRGB[3 * y * _bit_length_x + 3 * x + 1];
			blue = pixMapRGB[3 * y * _bit_length_x + 3 * x + 2];
			*(pixels+(y * _bit_length_x + x)) = Magick::ColorRGB(red, green, blue);
		}
	}

	/* make color bar */
	double* colors = new double[3];
	double value;
	/* pixel coordinates for color bar */
	int x_start = _bit_length_x - 40;
	int y_start = 10;
	int x_end = _bit_length_x - 20;
	int y_end = 90;

	/* draw color bar on pixel map */
	for (int y = y_start;y < y_end; y++){
		value = (double(y_end) - y)/(y_end - y_start);
		colors = getScaledColors(value, 0.0, 1.0, colors);
		for (int x = x_start; x < x_end; x++){
			*(pixels+(y * _bit_length_x + x)) = Magick::ColorRGB(colors[0], colors[1], colors[2]);
		}
	}

	delete [] colors;

	/* Sync pixel cache with Magick image */
	pixel_cache.sync();

	/* Draw black box around color bar */
	std::list<Magick::Drawable> drawList;
	drawList.push_back(Magick::DrawableStrokeColor("black"));    // Outline color
	drawList.push_back(Magick::DrawableStrokeWidth(2)); 		 // Stroke width
	drawList.push_back(Magick::DrawableStrokeAntialias(false));  // Don't antialias

	drawList.push_back(Magick::DrawableLine(x_start,y_start,x_end,y_start));  // top
	drawList.push_back(Magick::DrawableLine(x_start,y_end,x_end,y_end));      // bottom
	drawList.push_back(Magick::DrawableLine(x_start,y_start,x_start,y_end));  // left
	drawList.push_back(Magick::DrawableLine(x_end,y_start,x_end,y_end));      // right

	/* draw box on image */
	image.draw(drawList);

	/* create filename with correct extension */
	std::stringstream string;
	string << type << "." << _extension;
	std::string title = string.str();

	/* write Magick image to file */
	image.write(title);
}



/**
 * Generic function for plotting pixMap array on a structured mesh array
 * in a pdb file using silo I/O library
 */
void Plotter::plotSilo(int* pixMap, std::string type){
	log_printf(NORMAL, "plotting silo mesh...");

	/* Create file pointer */
    DBfile *file;

    /* create filename with correct extension */
	std::stringstream string;
	string << type << "." << _extension;
	std::string title_str = string.str();
	const char* title = title_str.c_str();

	/* Create pdb file */
	if (_extension == "h5"){
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_HDF5);
	}
	else{
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_PDB);
	}

    /* create mesh point arrays */
	double mesh_x[_bit_length_x + 1];
	double mesh_y[_bit_length_y + 1];

	/* generate structured mesh */
	for (int i = 0; i < (_bit_length_x + 1); i++){
		mesh_x[i] = (double(i) - double(_bit_length_x)/2.0 + 1.0) * (_width/double(_bit_length_x));
	}
	for (int i = 0; i < (_bit_length_y + 1); i++){
		mesh_y[i] = (double(i) - double(_bit_length_y)/2.0) * (_height/double(_bit_length_y));
	}

	/* descriptions of mesh */
	double *coords[] = {mesh_x, mesh_y};
	int dims[] = {_bit_length_x + 1, _bit_length_y + 1};
	int ndims = 2;

	/* Write structured mesh to pdb file */
	DBPutQuadmesh(file, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);

	/* dimensions of mesh */
	int dimsvar[] = {_bit_length_x, _bit_length_y};

	/* description of what is being plotted */
	const char* type_char = type.c_str();

	/* flip pixMap from Bitmap coordinates to cartesian coordinates */
	FlipBitmap(pixMap);

	/* write pixMap data to pdb file */
	DBPutQuadvar1(file, type_char, "quadmesh", pixMap, dimsvar, ndims, NULL, 0, DB_INT, DB_ZONECENT, NULL);

	/* flip pixMap from cartesian coordinates back to Bitmap coordinates */
	FlipBitmap(pixMap);

	/* close pdb file */
    DBClose(file);
	log_printf(NORMAL, "done plotting silo mesh...");
}

void Plotter::plotSilo(float* pixMap, std::string type){
	log_printf(NORMAL, "plotting silo mesh...");

	/* Create file pointer */
    DBfile *file;

    /* create filename with correct extension */
	std::stringstream string;
	string << type << "." << _extension;
	std::string title_str = string.str();
	const char* title = title_str.c_str();

	/* Create file */
	if (_extension == "h5"){
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_HDF5);
	}
	else{
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_PDB);
	}

    /* create mesh point arrays */
	double mesh_x[_bit_length_x + 1];
	double mesh_y[_bit_length_y + 1];

	/* create pixmap mesh */
	for (int i = 0; i < (_bit_length_x + 1); i++){
		mesh_x[i] = (double(i) - double(_bit_length_x)/2.0 + 1.0) * (_width/double(_bit_length_x));
	}
	for (int i = 0; i < (_bit_length_y + 1); i++){
		mesh_y[i] = (double(i) - double(_bit_length_y)/2.0) * (_height/double(_bit_length_y));
	}

	/* generate structured mesh */
	double *coords[] = {mesh_x, mesh_y};
	int dims[] = {_bit_length_x + 1, _bit_length_y + 1};
	int ndims = 2;

	/* Write structured mesh to pdb file */
	DBPutQuadmesh(file, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);

	/* dimensions of mesh */
	int dimsvar[] = {_bit_length_x, _bit_length_y};

	/* description of what is being plotted */
	const char* type_char = type.c_str();

	/* flip pixMap from Bitmap coordinates to cartesian coordinates */
	FlipBitmap(pixMap);

	/* write pixMap data to pdb file */
	DBPutQuadvar1(file, type_char, "quadmesh", pixMap, dimsvar, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);

	/* flip pixMap from cartesian coordinates back to Bitmap coordinates */
	FlipBitmap(pixMap);

	/* close pdb file */
    DBClose(file);
	log_printf(NORMAL, "done plotting silo mesh...");
}


/* PLOTTING HELPER FUNCTIONS */
/* These functions manipulate data and call the generic plotting functions
 */


/**
 * Based on user specified file extension, choose function to plot int pixMap
 */
void Plotter::plot(int* pixMap, std::string type){
	if (_extension == "png" || _extension == "tiff" || _extension == "jpg"){
		plotMagick(pixMap, type);
	}
	else if (_extension == "pdb" || _extension == "h5"){
		plotSilo(pixMap, type);
	}
}

/**
 * Based on user specified file extension, choose function to plot float pixMap
 */
void Plotter::plot(float* pixMap, std::string type){
	if (_extension == "png" || _extension == "tiff" || _extension == "jpg"){
		plotMagick(pixMap, type);
	}
	else if (_extension == "pdb" || _extension == "h5"){
		plotSilo(pixMap, type);
	}
}

/**
 * Vertically flip in int pixMap array
 */
void Plotter::FlipBitmap(int* pixMap){

	/* allocate temporary array */
	int* pixMapTemp = new int[_bit_length_x * _bit_length_y];

	/* Invert pixMap and store in pixMapTemp */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMapTemp[(_bit_length_y - 1 - y) * _bit_length_x + x] = pixMap[y * _bit_length_x + x];
		}
	}

	/* Write pixMapTemp to pixMap */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = pixMapTemp[y * _bit_length_x + x];
		}
	}

	/* release memory */
	delete [] pixMapTemp;
}

/**
 * Convert an x value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToBitmapX(double x){
	return int(x*_x_pixel + _bit_length_x/2);
}

/**
 * Convert an y value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToBitmapY(double y){
	return int(-y*_y_pixel + _bit_length_y/2);
}

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
 * Vertically flip in float pixMap array
 */
void Plotter::FlipBitmap(float* pixMap){

	/* allocate temporary array */
	float* pixMapTemp = new float[_bit_length_x * _bit_length_y];

	/* Invert pixMap and store in pixMapTemp */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMapTemp[(_bit_length_y - 1 - y) * _bit_length_x + x] = pixMap[y * _bit_length_x + x];
		}
	}

	/* Write pixMapTemp to pixMap */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = pixMapTemp[y * _bit_length_x + x];
		}
	}

	/* release memory */
	delete [] pixMapTemp;
}

/**
 * Initialize pixMap to -1. This allows us to see which array values
 * have not been edited.
 */
void Plotter::initializePixMap(int* pixMap){
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = -1;
		}
	}
}

/**
 * Initialize pixMap to -1. This allows us to see which array values
 * have not been edited.
 */
void Plotter::initializePixMap(float* pixMap){
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = -1;
		}
	}
}

/**
 * Plots track segments in a pixMap array
 */
void Plotter::plotSegments(Track* track, double sin_phi,
		double cos_phi, int* pixMap){

	/* Initialize variables */
	double x0, y0, x1, y1;
	int num_segments;

	/* Set first segment start point and get the number of tracks*/
	x0 = track->getStart()->getX();
	y0 = track->getStart()->getY();
	num_segments = track->getNumSegments();

	log_printf(DEBUG, "looping over segments in plotter endy: %f...", track->getEnd()->getY());
	/* loop over segments and write to pixMap array */
	for (int k=0; k < num_segments; k++){
		x1 = x0 + cos_phi*track->getSegment(k)->_length;
		y1 = y0 + sin_phi*track->getSegment(k)->_length;

		/* "draw" segment on pixMap array */
		log_printf(DEBUG, "calling linefct x0: %f, y0: %f, x1: %f, y1: %f,length: %f, sin_phi: %f ...", x0, y0, x1, y1, track->getSegment(k)->_region_id, track->getSegment(k)->_length, sin_phi);
		LineFct(x0, y0, x1, y1, pixMap, track->getSegment(k)->_region_id);
		log_printf(DEBUG, "done calling linefct...");

		x0 = x1;
		y0 = y1;
	}
}

/**
 * plot given track and numReflect reflected tracks
 */
void Plotter::plotTracksReflective(Track* track, int numReflect){
	log_printf(NORMAL, "Writing tracks reflect bitmap...");

	/* allocate pixMap array */
	int* pixMap = new int[_bit_length_x*_bit_length_y];

	/* initialize variables */
	double sin_phi, cos_phi, phi;
	Track *track2;
	bool get_out = TRUE;

	/* loop through tracks and write to pixMap array */
	for (int i = 0; i < (numReflect + 1); i++){
		log_printf(DEBUG, "plotting reflective track: %d", i);
		log_printf(DEBUG, "x_start, y_start, x_end, y_end: %f, %f, %f, %f",
				track->getStart()->getX(),track->getStart()->getY(),
				track->getEnd()->getX(),track->getEnd()->getY());
		phi = track->getPhi();
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		plotSegments(track, sin_phi, cos_phi, pixMap);

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
	plot(pixMap, "reflect");

	/* release memory */
	delete [] pixMap;
}

/**
 * Bresenham's line drawing algorithm. Takes in the start and end coordinates
 * of line (in geometry coordinates), pointer to pixMap array, and line color.
 * "Draws" the line on pixMap array.
 * Taken from "Simplificaiton" code at link below
 * http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
*/
void Plotter::LineFct(double xIn, double yIn, double xOut, double yOut, int* pixMap, int color){

	/* initialize variables */
	int x0, y0, x1,y1;

	/* convert geometry coordinates to bitmap coordinates */
	x0 = convertToBitmapX(xIn);
	y0 = convertToBitmapY(yIn);
	x1 = convertToBitmapX(xOut);
	y1 = convertToBitmapY(yOut);

	log_printf(DEBUG, "linefct bitmap x0: %i, y0: %i, x1: %i, y1: %i", x0, y0, x1, y1);

	/* "draw" line on pixMap array */
	int dx = abs(x1-x0);
	int dy = abs(y1-y0);
	int sx, sy;
	if (x0 < x1){
		sx = 1;
	}
	else{
		sx = -1;
	}
	if (y0 < y1){
		sy = 1;
	}
	else{
		sy = -1;
	}
	int error = dx - dy;
	pixMap[y0 * _bit_length_x + x0] = color;
	pixMap[y1 * _bit_length_x + x1] = color;
	while (x0 != x1 && y0 != y1){
		pixMap[y0 * _bit_length_x + x0] = color;
		int e2 = 2 * error;
		if (e2 > -dy){
			error = error - dy;
			x0 = x0 + sx;
		}
		if (e2 < dx){
			error = error + dx;
			y0 = y0 + sy;
		}
	}
}

/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion).
 */
void Plotter::plotRegion(int* pixMap, int* regionMap, std::string regionName){

	int regionValue;

	/* allocate memory for pixMapRegion array */
	int* pixMapRegion = new int[_bit_length_x * _bit_length_y];

	/* translate FSR id's stored in pixMap to pixMapRegion using regionMap */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			regionValue = regionMap[pixMap[y * _bit_length_x + x]];
			pixMapRegion[y * _bit_length_x + x] = regionValue;
		}
	}

	/* plot pixMapRegion array */
	plot(pixMapRegion, regionName);

	/* release memory */
	delete [] pixMapRegion;
}

/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion) in either a
 * scaled Magick Plot or a silo plot.
 */
void Plotter::plotRegion(int* pixMap, double* regionMap, std::string regionName){

	/* determine whethere color map will be scaled */
	if (_extension == "png" || _extension == "jpg" || _extension == "tiff"){

		/* allocate memory for pixMapRegion array */
		double* pixMapRegionRGB = new double[9 * _bit_length_x * _bit_length_y];

		pixMapRegionRGB = makeScaledMap(regionMap, pixMap, pixMapRegionRGB);

		plotMagickScaled(pixMapRegionRGB, regionName);

		delete [] pixMapRegionRGB;
	}
	else{

		float regionValue;

		/* allocate memory for pixMapRegion array */
		float* pixMapRegion = new float[_bit_length_x * _bit_length_y];

		/* translate FSR id's stored in pixMap to pixMapRegion using regionMap */
		for (int y=0;y< _bit_length_y; y++){
			for (int x = 0; x < _bit_length_x; x++){
				regionValue = regionMap[pixMap[y * _bit_length_x + x]];
				pixMapRegion[y * _bit_length_x + x] = regionValue;
			}
		}

		/* plot pixMapRegion array */
		plot(pixMapRegion, regionName);

		/* release memory */
		delete [] pixMapRegion;

	}
}

/**
 * Takes in a FSR pixMap array (pixMap), a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a new blank pixMap (pixMapRegion).
 * Writes RGB color data to pixMapRegion based on the intensity of the variable
 * being plotted.
 */
double* Plotter::makeScaledMap(double* regionMap, int* pixMap, double* pixMapRegionRGB){

	/* find max, min and range of regionMap */
	double max = regionMap[pixMap[0]];
	double min = regionMap[pixMap[0]];

	/* find min and max values in regionMap */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			max = std::max(max, regionMap[pixMap[y * _bit_length_x + x]]);
			min = std::min(min, regionMap[pixMap[y * _bit_length_x + x]]);
		}
	}

	log_printf(DEBUG, "pixel map min: %f, max: %f", min, max);

	/* allocate memory for RGB array */
	double* colors = new double[3];

	/* translate FSR id's stored in pixMap to pixMapRegion (RGB using regionMap */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			colors = getScaledColors(regionMap[pixMap[y * _bit_length_x + x]], min, max, colors);
			pixMapRegionRGB[3 * y * _bit_length_x + 3 * x] = colors[0];
			pixMapRegionRGB[3 * y * _bit_length_x + 3 * x + 1] = colors[1];
			pixMapRegionRGB[3 * y * _bit_length_x + 3 * x + 2] = colors[2];
		}
	}

	delete [] colors;

	return pixMapRegionRGB;
}


/**
 * Assigns a RGB color array (colors) to a pixel based on the value to be plotted
 * and the min and max of the variable being plotted (e.g. flux).
 */
double* Plotter::getScaledColors(double value, double min, double max, double* colors){

	/* initialize variables */
	double red;
	double green;
	double blue;
	double range = max - min;
	double mid = range / 2.0;

	/* assign a value to red */
	red = (value - mid) / range;

	/* if value is above mid => red-green */
	if (red > 0.0){
		blue = 0.0;
		green = 1.0 - 2 * red;
		red = 2 * red;
	    }
	/* if value is below mid => blue-green */
	else {
		blue = -2 * red;
		green = 1.0 + 2 * red;
		red = 0.0;
	}

	/* write color data to colors array */
	colors[0] = red;
	colors[1] = green;
	colors[2] = blue;

	return colors;
}


/**
 * Loops over pixels in pixMap array and finds the corresponding
 * FSR in geometry.
 */
void Plotter::plotFSRs(int* pixMap){
	log_printf(NORMAL, "Generating FSR maps...");

	/* initialize variables */
	double x_global;
	double y_global;

	/* loop over pixels */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){

			/* If pixel is blank, find what FSR it is in */
			if (pixMap[y * _bit_length_x + x] == -1){

				x_global = convertToGeometryX(x);
				y_global = convertToGeometryY(y);

				/* create point located in universe 0 */
				LocalCoords point(x_global,y_global);
				point.setUniverse(0);

				/* find which cell the point is in */
				_geom->findCell(&point);
				_geom->findCell(&point);

				/* Store FSR id in pixMap */
				pixMap[y * _bit_length_x + x] = _geom->findFSRId(&point);

			}
		}
	}

	/* plot pixMap array */
	plot(pixMap, "FSRs");
}
