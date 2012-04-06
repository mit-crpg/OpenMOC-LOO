/*
 * plotterNew.h
 *
 *  Created on: Apr 5, 2012
 *      Author: samuelshaner
 */

#ifndef PLOTTERNEW_H_
#define PLOTTERNEW_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "Magick++.h"
#include "silo.h"

typedef enum datatypes {
	INTEGER,
	FLOAT,
	DOUBLE
}datatype;

typedef enum colortypes {
	SCALED,
	RANDOM,
	BLACKWHITE
}colortype;

/* define BitMap struct */
template <typename U, typename T>
struct BitMap {
	U* pixels;
	datatype pixel_type;
	colortype color_type;
	int pixel_x;
	int pixel_y;
	T geom_x;
	T geom_y;
};

/*
 *  function definitions
 */
template <typename U, typename T>
void plot(BitMap<U,T>* bitMap, std::string name, std::string extension);
template <typename U, typename T>
void plotSilo(BitMap<U,T>* bitMap, float* pixMap, std::string type, std::string extension);
template <typename U, typename T>
void plotMagick(BitMap<U,T>* bitMap, float* pixMap, std::string type, std::string extension);
template <typename U, typename T>
void copyBitMap(BitMap<U,T>* bitMap, float* pixMap);
template <typename U, typename T>
void normalize(BitMap<U,T>* bitMap, float* pixMap);
template <typename U, typename T>
void getBounds(BitMap<U,T>* bitMap, float* pixMap, float* bounds);
template <typename U, typename T>
void getColorRGB(BitMap<U,T>* bitMap, float value, float* color);
template <typename U, typename T>
void randomize(BitMap<U,T>* bitMap, float* pixMap);
template <typename U, typename T>
void initialize(BitMap<U,T>* bitMap);


/*
 * function declarations
 */

/* templated general plot function */
template <typename U, typename T>
void plot(BitMap<U,T>* bitMap, std::string name, std::string extension){

	/* create array to store color values */
	float* pixMap = new float[bitMap->pixel_x * bitMap->pixel_y];
	copyBitMap(bitMap, pixMap);

	/* decide which plot function to call */
	if (extension == "png" || extension == "tiff" || extension == "jpg"){
		plotMagick(bitMap, pixMap, name, extension);
	}
	else if (extension == "pdb" || extension == "h5"){
		plotSilo(bitMap, pixMap, name, extension);
	}

	delete [] pixMap;
}

template <typename U, typename T>
void plotSilo(BitMap<U,T>* bitMap, float* pixMap, std::string name, std::string extension){
	printf("plotting silo mesh...\n");

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

	/* if color_type is RANDOM, randomize bitMapRGB */
	if (bitMap->color_type == RANDOM){
		normalize(bitMap, pixMap);
		randomize(bitMap, pixMap);
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
	DBPutQuadvar1(file, type_char, "quadmesh", pixMap, dimsvar, ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);

	/* close file */
    DBClose(file);
	printf("done plotting silo mesh...\n");
}

/**
 * Generic function for plotting pixMap in png, tiff, or jpg file
 * using Magick++
 */
template <typename U, typename T>
void plotMagick(BitMap<U,T>* bitMap, float* pixMap, std::string name, std::string extension){
	printf("Writing Magick bitmap...\n");

	/* declare variables */
	float* color = new float[3];

	/* copy bitMap to bitMapRGB and normalize */
	normalize(bitMap, pixMap);

	/* if color_type is RANDOM, randomize numbers */
	if (bitMap->color_type == RANDOM){
		randomize(bitMap, pixMap);
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
			if (bitMap->pixels[y * bitMap->pixel_x + x] != -1){

				if (bitMap->color_type == BLACKWHITE){
					*(pixels+(y * bitMap->pixel_x + x)) = Magick::ColorRGB(0, 0, 0);
				}
				else{
					getColorRGB(bitMap, pixMap[y * bitMap->pixel_x + x], color);
					*(pixels+(y * bitMap->pixel_x + x)) = Magick::ColorRGB(color[0], color[1], color[2]);
				}
			}
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

	delete [] color;
}

/* copy elements in bitMap to bitMapRGB */
template <typename U, typename T>
void copyBitMap(BitMap<U,T>* bitMap, float* pixMap){

	/* copy bitMap to bitMapRGB */
	for (int y=0;y<bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			pixMap[y * bitMap->pixel_x + x] = (float)bitMap->pixels[y * bitMap->pixel_x + x];
		}
	}

}

/* normalize bitMapRGB to numbers between 0 and 1 */
template <typename U, typename T>
void normalize(BitMap<U,T>* bitMap, float* pixMap){

	float* bounds = new float[2];
	getBounds(bitMap, pixMap, bounds);

	/* copy bitMap to bitMapRGB and normalize */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x=0;x< bitMap->pixel_x; x++){
			pixMap[y * bitMap->pixel_x + x] = (pixMap[y * bitMap->pixel_x + x] - bounds[0]) /  (bounds[1] - bounds[0]);
		}
	}

	delete [] bounds;
}

/* get min and max bounds of bitMapRGB */
template <typename U, typename T>
void getBounds(BitMap<U,T>* bitMap, float* pixMap, float* bounds){

	bounds[0] = pixMap[0];
	bounds[1] = pixMap[0];

	/* find max */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bounds[1] = std::max(bounds[1], pixMap[y * bitMap->pixel_x + x]);
		}
	}

	bounds[0] = bounds[1];
	/* find min */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			if (pixMap[y * bitMap->pixel_x + x] != -1){
				bounds[0] = std::min(bounds[0], pixMap[y * bitMap->pixel_x + x]);
			}
		}
	}

}

/* write RGB triplet to color using HOT color scheme */
template <typename U, typename T>
void getColorRGB(BitMap<U,T>* bitMap, float value, float* color){

	if (value < 1.0/3.0){
		color[0] = 0.0;
		color[1] = 3.0 * value;
		color[2] = 1.0;
	}
	else if (value < 2.0/3.0){
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

/* pseudorandomize bitMapRGB with number between 0 and 1 */
template <typename U, typename T>
void randomize(BitMap<U,T>* bitMap, float* pixMap){

	/* make array to store random numbers */
	float* myRandoms = new float[97];

	/* make random numbers */
	srand(1);
	for (int i=0;i< 99; i++){
		myRandoms[i] = rand() / float(RAND_MAX);
	}

	/* randomize bitMapRGB */
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			pixMap[y * bitMap->pixel_x + x] = myRandoms[abs(int(pixMap[y * bitMap->pixel_x + x] / 1e-8)) % 97];
		}
	}

	delete [] myRandoms;
}

/* initialize values to -1 */
template <typename U, typename T>
void initialize(BitMap<U,T>* bitMap){
	for (int y=0;y< bitMap->pixel_y; y++){
		for (int x = 0; x < bitMap->pixel_x; x++){
			bitMap->pixels[y * bitMap->pixel_x + x] = -1;
		}
	}
}

#endif /* PLOTTERNEW_H_ */
