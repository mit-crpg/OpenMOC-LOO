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
#include "Magick++.h"
#include "silo.h"

typedef enum datatypes {
	INTEGER,
	FLOAT,
	DOUBLE
}datatype;

typedef enum colortypes {
	SCALED,
	RANDOM
}colortype;


struct BitMap {
	template <typename U>
	U* pixels;
	datatype pixel_type;
	colortype color_type;
	int pixel_x;
	int pixel_y;
	template <typename T>
	T geom_x;
	template <typename T>
	T geom_y;
};

/* function definitions */
void plotSilo(BitMap* bitMap, std::string type, std::string extension);
void plotMagick(BitMap* bitMap, std::string type, std::string extension);
void copyBitMap(BitMap* bitMap, float* bitMapRGB);
void normalize(float* bitMapRGB);
void getBounds(float* bitMapRGB, double* bounds);
void getColor(float value, double* color);
void randomize(float* bitMapRGB);

#endif /* PLOTTERNEW_H_ */
