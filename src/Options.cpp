/*
 * Options.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#include "Options.h"
#include <string.h>
#include <stdlib.h>

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)

Options::Options(int argc, const char **argv) {

	int i;
	this->geometry_file = NULL;

	for (i = 0; i < argc; i++) {
		if (i > 0) {
			if (LAST("--geometryfile") || LAST("-g")) {
				if (this->geometry_file != NULL)
				free(this->geometry_file);
				geometry_file = strdup(argv[i]);
			}
		}
	}
}
