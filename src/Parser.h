/*
 * Parser.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <expat.h>

/**
 * Parses XML Files
 */
class Parser {
private:

public:
	Parser(const char *filename);
	virtual ~Parser();
	void parseMaterials(void);
	void parseGeometry(void);
};

#endif /* POINT_H_ */
