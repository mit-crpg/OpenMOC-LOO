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

/**
 * Parses XML Files
 */
class Parser {
private:
//    std::string geoxml;
public:
    Parser(const char* geoxml);
	virtual ~Parser();
	void parseMaterials(void);
	void parseGeometry(void);
};

#endif /* POINT_H_ */
