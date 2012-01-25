/*
 * Parser.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 */

#include "Parser.h"
#include <vector>
#include <cstring>
#include <string>
#include <stdexcept>


static void parse(XML_Parser parser, char c, int isFinal) {
	if (XML_STATUS_OK == XML_Parse(parser, &c, isFinal ^ 1, isFinal))
		return;
	
	fprintf( stderr, "ERROR: parsing XML failed at line %lu, pos %lu: %s\n",
		 (unsigned long)XML_GetCurrentLineNumber(parser),
		 (unsigned long)XML_GetCurrentColumnNumber(parser),
		 XML_ErrorString(XML_GetErrorCode(parser)) );
	exit(1);
}

typedef struct {
	FILE* out;
	int depth;
	int tagInd;
	int attrInd;
	std::vector<std::string *> *level;
} PContext;

/* callback for start element, e.g. <tag> */
static void XMLCALL startElementCallback( void *context,
					  const XML_Char *name,
					  const XML_Char **atts ) {
//    int is_key = 1;
	PContext *ctxt;
	/* surface related */
	Surface* surface;
	int surface_id;
	surfaceType surface_type;
	std::vector<double> coeffs;
	/* cell related */
	Cell* cell;
	int cell_id;
	cellType cell_type;
	int universe = 0;
	int num_surfaces;
	std::vector<int> surfaces;
	int material;
	int universe_fill = 0;

	ctxt = (PContext*)context;
	
#if 0
	fprintf(ctxt->out, "%*s%c%s\n", ctxt->depth, "", '>', name);
	ctxt->level->push_back(new std::string(name));	
	/* print level information to the screen  */
	for (std::vector<std::string *>::size_type i = 0;
	     i < ctxt->level->size(); i++) {
		fprintf(ctxt->out, " %s", ctxt->level->at(i)->c_str());
	}
	fprintf(ctxt->out, "\n");
#endif
	
	/* Parse surface types */
	if (strcmp(name, "surface") == 0) {
		while (*atts) {
			if (strcmp(*atts, "id") == 0) {
				++atts;
				surface_id = atoi(*atts);
			}
			
			if (strcmp(*atts, "type") == 0) {
				++atts;
				if (strcmp(*atts, "plane") == 0) {
					surface_type = PLANE;
				}
				if (strcmp(*atts, "x-plane") == 0) {
					surface_type = XPLANE;
				}
				if (strcmp(*atts, "y-plane") == 0) {
					surface_type = YPLANE;
				}
				if (strcmp(*atts, "circle") == 0) {
					surface_type = CIRCLE;
				}
			}
			
			if (strcmp(*atts, "coeffs") == 0) {
				++atts;
				char *long_str;
				char *tmp;
				char *result;
				
				long_str = strdup(*atts);
				result = strtok_r(long_str, " ,.", &tmp);
				coeffs.push_back(atof(result));
				while ((result = strtok_r(NULL, " ,.", &tmp))
				       != NULL) {
					coeffs.push_back(atof(result));
				}
			}		
			++atts;
		}
		
		/* TODO: check for number of coeffs */
		switch (surface_type) {
		case PLANE:
		        log_printf(NORMAL, "Parsed Surfs: plane, id = %d,"
				   " coeffs = %.1f, %.1f, %.1f.\n",
				   surface_id, coeffs.at(0), coeffs.at(1), 
				   coeffs.at(2));
			surface = new Plane(surface_id, coeffs.at(0),
					    coeffs.at(1), coeffs.at(2));
			break;
		case XPLANE:
		        log_printf(NORMAL, "Parsed Surfs: x-plane, id = %d,"
				   " coeffs = %.1f.\n",
				   surface_id, coeffs.at(0));
			surface = new XPlane(surface_id, coeffs.at(0));
			break;
		case YPLANE:
		        log_printf(NORMAL, "Parsed Surfs: y-plane, id = %d,"
				   " coeffs = %.1f.\n",
				   surface_id, coeffs.at(0));
			surface = new YPlane(surface_id, coeffs.at(0));
			break;
		case CIRCLE:
		        log_printf(NORMAL, "Parsed Surfs: circle, id = %d,"
				   " coeffs = %.1f, %.1f, %.1f\n",
				   surface_id, coeffs.at(0), coeffs.at(1), 
				   coeffs.at(2));
			surface = new Circle(surface_id, coeffs.at(0), 
					     coeffs.at(1), coeffs.at(2));
			break;
		case QUADRATIC:
			log_printf(ERROR, "quadratic is called; should not"
				   " be used in geometry.xml.\n");
			break;
		}
		
	}

	/* Parse Cell Types */
	if (strcmp(name, "cell") == 0) {
		while (*atts) {
			if (strcmp(*atts, "id") == 0) {
				++atts;
				cell_id = atoi(*atts);
			}
			
			if (strcmp(*atts, "universe") == 0) {
				++atts;
				universe = atoi(*atts);
			}

			if (strcmp(*atts, "material") == 0) {
				cell_type = MATERIAL;
				++atts;
				material = atoi(*atts);
			}

			if (strcmp(*atts, "fill") == 0) {
				cell_type = FILL;
				++atts;
				universe_fill = atoi(*atts);
			}

	 		if (strcmp(*atts, "surfaces") == 0) {
				++atts;
				char *long_str;
				char *tmp;
				char *result;
				
				long_str = strdup(*atts);
				result = strtok_r(long_str, " ,.", &tmp);
				surfaces.push_back(atof(result));
				while ((result = strtok_r(NULL, " ,.", &tmp))
				       != NULL) {
					surfaces.push_back(atof(result));
				}
				num_surfaces = surfaces.size();
			}
			++atts;
		}
		
		/* TODO: check for number of coeffs */
		switch (cell_type) {;
		case MATERIAL:
		        log_printf(NORMAL, "Parsed Cells: id = %d,"
				   " universe = %d, 1st surface = %d,"
				   " # of surfaces = %d, material = %d.\n",
				   cell_id, universe, surfaces.at(0),
				   num_surfaces, material);
			cell = new CellBasic(surface_id, universe, num_surfaces,
					     surfaces, material);
			break;
		case FILL:
		        log_printf(NORMAL, "Parsed Cells: id = %d,"
				   " universe = %d, 1st surface = %d,"
				   " # of surfaces = %d, u_fill = %d.\n",
				   cell_id, universe, surfaces.at(0),
				   num_surfaces, universe_fill);
			cell = new CellFill(surface_id, universe, num_surfaces,
					     surfaces, universe_fill);			
			break;
		}
	
	}

/* Enable the following for printing directly from the xml file */	
#if 0
	while (*atts) {
		if (is_key) {
			fprintf(ctxt->out, "%*c%s: ",
				ctxt->depth + ctxt->attrInd, ' ', *atts);
		}
		else {
			fprintf(ctxt->out, "%s\n", *atts);
		}
		is_key = !is_key;
		++atts;
	}
	ctxt->depth += ctxt->tagInd;
#endif
}


/* callback for end elements, e.g. </tag>,
 * it is called for empty elements, too
 */
static void XMLCALL
endElementCallback( void *context,
		    const XML_Char *name __attribute__((__unused__)) ) {
	PContext *ctxt = (PContext*)context;
	// ctxt->depth -= ctxt->tagInd;
	ctxt->level->pop_back();
}

/**
 * Default constructor
 */
Parser::Parser (const Options *opts) {
	const char *geoxml;
	XML_Parser parser;
	PContext ctxt;
	char c;
	FILE* geofile;
	
	if (opts->getGeometryFile() == NULL)
		throw std::runtime_error("No geometry file given");
	
	geoxml = opts->getGeometryFile();
	geofile = fopen(geoxml, "r");
	if (geofile == NULL)
		throw std::runtime_error("Given geometry file does not exist");
	
/* Create parser.
 * The only argument for XML_ParserCreate is encoding, and if it's NULL,
 * then encoding declared in the document is used.
 */
	parser = XML_ParserCreate(NULL);
	if (!parser) {
		fprintf(stderr, "Couldn't allocate memory for parser\n");
		exit(-1);
	}
	
/* Set context that will be passed by the parsers to all handlers */
	ctxt.out = stdout;
	ctxt.depth = 1;
	ctxt.tagInd = 4;
	ctxt.attrInd = 6;
	ctxt.level = new std::vector<std::string *>();
	
	XML_SetUserData(parser, &ctxt);
	
/* set callback for start element */
	XML_SetStartElementHandler(parser, &startElementCallback);
	
/* set callback for start element */
	XML_SetEndElementHandler(parser, &endElementCallback);
	
/* If you'd like to read input by large blocks, you can have a look at
 * XML_GetBuffer and XML_ParseBuffer functions.
 */
	while( EOF != (c = fgetc(geofile)) ) {
		parse(parser, c, 0);
	}
	
/* Finish parsing, note the last argument or XML_Parse */
	parse(parser, c, 1);
	
/* Free resource used by expat */
	XML_ParserFree(parser);
	
	fclose(geofile);
}


/**
 * Destructor
 */
Parser::~Parser() { }

void Parser::parseMaterials(void) {
}

void Parser::parseGeometry(void) {
}
