/*
 * plotterNew.cpp
 *
 *  Created on: Apr 5, 2012
 *      Author: samuelshaner
 */

#include "plotterNew.h"

void plotSilo(BitMap* bitMap, std::string name, std::string extension){
	printf("plotting silo mesh...");

	/* create array to store color values */
	float* bitMapRGB = new float[bitMap->pixel_x * bitMap->pixel_y];

	/* Create file pointer */
    DBfile *file;

    /* create filename with correct extension */
	std::stringstream string;
	string << name << "." << extension;
	std::string title_str = string.str();
	const char* title = title_str.c_str();

	/* Create file */
	if (extension == "h5"){
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_HDF5);
	}
	else{
		file = DBCreate(title, DB_CLOBBER, DB_LOCAL, "structured mesh bitmap", DB_PDB);
	}

	/* copy bitMap to bitMapRGB */
	copyBitMap(bitMap, bitMapRGB);

	/* if color_type is RANDOM, randomize bitMapRGB */
	if (bitMap->color_type == RANDOM){
		normalize(bitMapRGB);
		randomize(bitMapRGB);
	}

    /* create mesh point arrays */
	double mesh_x[bitMap->pixel_x + 1];
	double mesh_y[bitMap->pixel_y + 1];

	/* generate structured mesh */
	for (int i = 0; i < (bitMap->pixel_x + 1); i++){
		mesh_x[i] = (double(i) - double(bitMap->pixel_x)/2.0 + 1.0) * (bitMap->geom_x/double(bitMap->pixel_x));
	}
	for (int i = 0; i < (bitMap->pixel_y + 1); i++){
		mesh_y[i] = (double(bitMap->pixel_y)/2.0 - double(i)) * (bitMap->geom_y/double(bitMap->pixel_y));
	}

	/* descriptions of mesh */
	double *coords[] = {mesh_x, mesh_y};
	int dims[] = {bitMap->pixel_x + 1, bitMap->pixel_y + 1};
	int ndims = 2;

	/* Write structured mesh to file */
	DBPutQuadmesh(file, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, NULL);

	/* dimensions of mesh */
	int dimsvar[] = {bitMap->pixel_x, bitMap->pixel_y};

	/* description of what is being plotted */
	const char* type_char = name.c_str();

	/* write pixMap data to file */
	DBPutQuadvar1(file, type_char, "quadmesh", bitMap->pixels, dimsvar, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);

	/* close file */
    DBClose(file);
	printf("done plotting silo mesh...");

	delete [] bitMapRGB;
}

/**
 * Generic function for plotting pixMap in png, tiff, or jpg file
 * using Magick++
 */
void plotMagick(BitMap* bitMap, std::string name, std::string extension){
	printf("Writing Magick bitmap...");

	/* declare variables */
	double* color = new double[3];
	float* bitMapRGB = new float[bitMap->pixel_x * bitMap->pixel_y];

	/* copy bitMap to bitMapRGB and normalize */
	copyBitMap(bitMap, bitmapRGB);
	normalize(bitMapRGB);

	/* if color_type is RANDOM, randomize numbers */
	if (bitMap->color_type == RANDOM){
		randomize(bitMapRGB);
	}

	/* create image and open for modification */
	Magick::Image image(Magick::Geometry(bitMap->pixel_x,bitMap->pixel_y), "white");
	image.modifyImage();

	/* Make pixel cache */
	Magick::Pixels pixel_cache(image);
	Magick::PixelPacket* pixels;
	pixels = pixel_cache.get(0,0,bitMap->pixel_x,bitMap->pixel_y);

	/* Write pixMapRGB array to Magick pixel_cache */
	for (int y=0;y<bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			getColor(bitMapRGB[y * bitMap->pixel_x + x], color);
			*(pixels+(y * bitMap->pixel_x + x)) = Magick::ColorRGB(color[0], color[1], color[2]);
		}
	}

	/* Sync pixel cache with Magick image */
	pixel_cache.sync();

	/* create filename with correct extension */
	std::stringstream string;
	string << name << "." << extension;
	std::string title = string.str();

	/* write Magick image to file */
	image.write(title);

	delete [] bitMapRGB;
	delete [] color;
}


void copyBitMap(BitMap* bitMap, float* bitMapRGB){

	/* copy bitMap to bitMapRGB */
	for (int y=0;y<bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bitMapRGB[y * bitMap->pixel_x + x] = (float)bitMap->pixels[y * bitMap->pixel_x + x];
		}
	}

}


void normalize(float* bitMapRGB){

	double* bounds = new double[2];
	getBounds(bitMapRGB, bounds);

	/* copy bitMap to bitMapRGB and normalize */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x=0;x< bitMap->pixel_x; x++){
			bitMapRGB[y * bitMap->pixel_x + x] = bitMapRGB[y * bitMap->pixel_x + x] / bounds[1] + bounds[0] / bounds[1];
		}
	}

	delete [] bounds;
}


void getBounds(float* bitMapRGB, double* bounds){

	bounds[0] = bitMapRGB[0];
	bounds[1] = bitMapRGB[0];

	/* find regionMap */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bounds[0] = std::min(bounds[0], bitMap->pixels[y * bitMap->pixel_x + x]);
			bounds[1] = std::max(bounds[1], bitMap->pixels[y * bitMap->pixel_x + x]);
		}
	}
}


void getColor(float value, double* color){

	if (value < .333){
		color[0] = 0.0;
		color[1] = 3.0 * value;
		color[2] = 1.0;
	}
	else if (value < .666){
		color[0] = 3.0 * value - 1.0;
		color[1] = 1.0;
		color[2] = -3.0 * value + 1.0;
	}
	else {
		color[0] = 1.0;
		color[1] = -3.0 * value + 3.0;
		color[2] = 0.0;
	}
}


void randomize(float* bitMapRGB){

	/* make array to store random numbers */
	double* myRandoms = new double [bitMap->pixel_x * bitMap->pixel_y];

	/* make random numbers */
	srand(1);
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x=0;x< bitMap->pixel_x; x++){
			myRandoms[y * bitMap->pixel_x + x] = rand() / double(RAND_MAX);
		}
	}

	/* randomize bitMapRGB */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bitMapRGB[y * bitMap->pixel_x + x] = myRandoms[floor(bitMapRGB[y * bitMap->pixel_x + x]
			                                                                       * (bitMap->pixel_x * bitMap->pixel_y))];
		}
	}

	delete [] myRandoms;
}




