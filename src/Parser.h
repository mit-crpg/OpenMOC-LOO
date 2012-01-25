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
#include "Cell.h"
#include "Options.h"
#include "log.h"

#include <vector>
#include "Surface.h"
#include <functional>

/**
 * Parses XML Files
 */
class Parser {
private:
	std::vector<Surface *> surfaces;
	std::vector<Cell *> cells;

public:
	Parser(const Options *opts);
	virtual ~Parser();

	void each_surface(std::function<void(Surface *)> callback);
	void each_cell(std::function<void(Cell *)> callback);

	friend void XMLCALL Parser_XMLCallback_Start(void *context,
						     const XML_Char *name,
						     const XML_Char **atts);
};

#endif /* POINT_H_ */
