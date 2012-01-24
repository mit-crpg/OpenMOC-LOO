/*
 * Options.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

class Options {
public:
	char *geometry_file;	/* TODO: make class attributes private with getter and setter methods */
    Options(int argc, const char **argv);
};

#endif
