/*
 * Parser.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <expat.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Surface.h"
#include "Options.h"

/**
 * Parses XML Files
 */
class Parser {
private:
//    std::string geoxml;
public:
	Parser(const Options *opts);
	virtual ~Parser();
	void parseMaterials(void);
	void parseGeometry(void);
};

#endif /* POINT_H_ */
