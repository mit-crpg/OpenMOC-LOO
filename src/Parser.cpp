/*
 * Parser.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
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
	int surface_id;
	surfaceType surface_type;
	Surface* surface;
	double coeff1, coeff2, coeff3;
	
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
				if (strcmp(*atts, "xplane") == 0) {
					surface_type = XPLANE;
				}
				if (strcmp(*atts, "yplane") == 0) {
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
				
				long_str = strdup(*atts);
				coeff1 = atof(strtok_r(long_str, " ,.", &tmp));
				coeff2 = atof(strtok_r(NULL, " ,.", &tmp));
				coeff3 = atof(strtok_r(NULL, " ,.", &tmp));
			}
			
			++atts;
		}
		
		switch (surface_type) {
		case PLANE:
			surface = new Plane(surface_id, coeff1, coeff2, coeff3);
			break;
		case XPLANE:
			surface = new XPlane(surface_id, coeff1);
			break;
		case YPLANE:
			surface = new YPlane(surface_id, coeff1);
			break;
		case CIRCLE:
			fprintf(ctxt->out, "%d %f %f %f\n",
				surface_id, coeff1, coeff2, coeff3);
			surface = new Circle(surface_id, coeff1, coeff2, coeff3);
			break;
		case QUADRATIC:
			fprintf(ctxt->out,
				"quadratic is called; should not"
				" be used in geometry.xml.\n");
			break;
		}
		
	}
	
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
	
	if (opts->geometry_file == NULL)
		throw std::runtime_error("No geometry file given");
	
	geoxml = opts->geometry_file;
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
