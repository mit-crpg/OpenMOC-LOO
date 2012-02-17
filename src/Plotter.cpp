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
 * @param plotter a pointer to a plotting object
 * @param num_azim number of azimuthal angles
 * @param spacing track spacing
 */
Plotter::Plotter(Geometry* geom, const int bitDim, std::string extension,
		bool plotMaterials, bool plotCells) {

	_extension = extension;
	_geom = geom;
	_plot_materials = plotMaterials;
	_plot_cells = plotCells;

	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;

	_bit_length_x = int (bitDim*ratio) + 1;
	_bit_length_y = bitDim + 1;

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

	try{
		_pix_map_segments = new int[_bit_length_x*_bit_length_y];
		_pix_map_tracks = new int[_bit_length_x*_bit_length_y];
		_pix_map_FSR = new int[_bit_length_x*_bit_length_y];

		if (_plot_materials == true){
			_pix_map_materials = new int[_bit_length_x*_bit_length_y];
		}
		if (_plot_cells == true){
			_pix_map_cells = new int[_bit_length_x*_bit_length_y];
		}


	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to allocate memory needed to generate pixel maps"
				".Backtrace:\n%s", e.what());
	}
	for (int i=0;i<_bit_length_x; i++){
		for (int j = 0; j < _bit_length_y; j++){
			_pix_map_segments[i * _bit_length_x + j] = -1;
			_pix_map_tracks[i * _bit_length_x + j] = -1;
			_pix_map_FSR[i * _bit_length_x + j] = -1;

			if (_plot_materials == true){
				_pix_map_materials[i * _bit_length_x + j] = -1;;
			}
			if (_plot_cells == true){
				_pix_map_cells[i * _bit_length_x + j] = -1;;
			}
		}
	}

	_x_pixel = double(_bit_length_x)/_width;
	_y_pixel = double(_bit_length_y)/_height;
}

/**
 * Plotter destructor frees memory for all tracks
 */
Plotter::~Plotter() {
	delete [] _pix_map_segments;
	delete [] _pix_map_tracks;
}


/**
 * Plot _pix_map (cartesian coordinates)
 * in png, tiff, or jpg file using Magick++
 */
void Plotter::plotMagick(int* pixMap, std::string type){
	log_printf(NORMAL, "Writing Magick bitmap...");

	int color_int;

	Magick::Image image(Magick::Geometry(_bit_length_x,_bit_length_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,_bit_length_x,_bit_length_y);

	/* Convert _pix_map_tracks bitmap array to Magick bitmap pixel array. */
	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			color_int = pixMap[y * _bit_length_x + x] % 15;

			if (color_int != -1){
				*(pixels+(y * _bit_length_x + x)) =
						Magick::Color(_color_map.at(color_int));
			}
		}
	}

	/* Close pixel viewing/changing */
	pixel_cache.sync();

	/* Write pixel bitmap to png file */

	std::stringstream string;
	string << type << "." << _extension;
	std::string title = string.str();

	image.write(title);
}


/**
 * Plots track segments in a _pix_map_segments bitmap array on the fly
 */
void Plotter::plotSegmentsBitMap(Track* track, double sin_phi,
		double cos_phi, int* pixMap){

	/* Initialize variables */
	double start_x, start_y, end_x, end_y;
	int num_segments;

	/* Set first segment start point and get the number of tracks*/
	start_x = track->getStart()->getX();
	start_y = track->getStart()->getY();
	num_segments = track->getNumSegments();

	/* loop over segments and write to _pix_map_segments bitmap array */
	for (int k=0; k < num_segments; k++){
		end_x = start_x + cos_phi*track->getSegment(k)->_length;
		end_y = start_y + sin_phi*track->getSegment(k)->_length;

		/*
		 * Add line to _pix_map segments bitmap array.
		 * Note conversion from geometric coordinate system to
		 * bitmap coordinates.
		 */
		LineFct(start_x*_x_pixel + _bit_length_x/2,
				-start_y*_y_pixel + _bit_length_y/2,
				end_x*_x_pixel + _bit_length_x/2,
				-end_y*_y_pixel + _bit_length_y/2,
				pixMap, track->getSegment(k)->_region_id);


		start_x = end_x;
		start_y = end_y;
	}
}


/**
 * plot flat source regions in pdb file using segments bitmap
 */
void Plotter::plotSilo(int* pixMap, std::string type){
	log_printf(NORMAL, "plotting silo mesh...");

	/* Create pdb file */
    DBfile *pdb_file;

	std::stringstream string;
	string << type << "." << _extension;
	std::string title_str = string.str();
	const char* title = title_str.c_str();


    pdb_file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_PDB);

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

	/* descriptions of mesh */
	double *coords[] = {mesh_x, mesh_y};
	int dims[] = {_bit_length_x + 1, _bit_length_y + 1};
	int ndims = 2;

	/* create structured mesh bit map in pdb_file */
	DBPutQuadmesh(pdb_file, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);

	/* dimensions of mesh */
	int dimsvar[] = {_bit_length_x, _bit_length_y};

	/* Save FSR bit map (_bit_map_visit) to the pdb file */

	const char* type_char = type.c_str();

	DBPutQuadvar1(pdb_file, type_char, "quadmesh", pixMap, dimsvar, ndims, NULL, 0, DB_INT, DB_ZONECENT, NULL);

	/* close pdb file */
    DBClose(pdb_file);
}

/**
 * plot given track and numReflect reflected tracks
 */
void Plotter::plotTracksReflective(Track* track, int numReflect){
	log_printf(NORMAL, "Writing tracks reflect bitmap...");

	int* pix_map_reflect = new int[_bit_length_x*_bit_length_y];

	/* initialize variables */
	double sin_phi, cos_phi, phi;
	Track *track2;
	bool get_out = TRUE;

	/* loop through tracks and write to _pix_map_reflect bitmap */
	for (int i = 0; i < (numReflect + 1); i++){
		log_printf(DEBUG, "plotting reflective track: %d", i);
		log_printf(DEBUG, "x_start, y_start, x_end, y_end: %f, %f, %f, %f",
				track->getStart()->getX(),track->getStart()->getY(),
				track->getEnd()->getX(),track->getEnd()->getY());
		phi = track->getPhi();
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		plotSegmentsBitMap(track, sin_phi, cos_phi, pix_map_reflect);

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

	if (_extension == "png" || _extension == "tiff" || _extension == "jpg"){
		plotMagick(pix_map_reflect, "reflect");
	}
	else if (_extension == "pdb"){
		plotSilo(pix_map_reflect, "reflect");
	}

	delete [] pix_map_reflect;

}

/**
 * Bresenham's line drawing algorithm. Takes in the start and end coordinates
 * of line (gometric coordinates), pointer to _pix_map bitmap array (pixMap), and line color.
 * "Draws" the line on _pix_map bitmap array.
 * Taken from "Simplificaiton" code at link below
 * http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
*/
void Plotter::LineFct(int x0, int y0, int x1, int y1, int* pixMap, int color){
	log_printf(DEBUG, "Drawing segment...");

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


int* Plotter::getPixMap(std::string type){
	if (type == "tracks"){
		return _pix_map_tracks;
	}
	else if (type == "segments"){
		return _pix_map_segments;
	}
	else if (type == "FSR"){
		return _pix_map_FSR;
	}
	else {
		log_printf(ERROR, "Invalid pixMap name give to Plotter::getPixMap");
		return NULL;
	}
}

/*
 * Generate bitmap of geometry with pixels mapped to FSR ids
 */
void Plotter::generateFsrMap(){
	log_printf(NORMAL, "Generating FSR maps...");

	double x_global;
	double y_global;

	for (int y=0;y<_bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){

			x_global = (double(x) - _bit_length_x/2)/_x_pixel;
			y_global = -(double(y) - _bit_length_y/2)/_y_pixel;

			LocalCoords point(x_global,y_global);
			point.setUniverse(0);

			Cell* curr = _geom->findCell(&point);

			_pix_map_FSR[y * _bit_length_x + x] = _geom->findFSRId(&point);

			if (_plot_materials == true){
				_pix_map_materials[y * _bit_length_x + x] = static_cast<CellBasic*>(curr)->getMaterial();
			}
			if (_plot_cells == true){
				_pix_map_cells[y * _bit_length_x + x] = curr->getUid();
			}
		}
	}

	if (_extension == "png" || _extension == "tiff" || _extension == "jpg"){
		plotMagick(_pix_map_FSR, "fsr");
		if (_plot_materials == true){
			plotMagick(_pix_map_materials, "materials");
		}
		if (_plot_cells == true){
			plotMagick(_pix_map_cells, "cells");
		}
	}
	else if (_extension == "pdb"){
		plotSilo(_pix_map_FSR, "fsr");
		if (_plot_materials == true){
			plotSilo(_pix_map_materials, "materials");

		}
		if (_plot_cells == true){
			plotSilo(_pix_map_cells, "cells");
		}
	}
}

std::string Plotter::getExtension(){
	return _extension;
}

int Plotter::getBitLengthX(){
	return _bit_length_x;
}

int Plotter::getBitLengthY(){
	return _bit_length_y;
}

double Plotter::getXPixel(){
	return _x_pixel;
}

double Plotter::getYPixel(){
	return _y_pixel;
}

