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

static void parse(XML_Parser parser, char c, int isFinal)
{
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
static void XMLCALL startElementCallback( void *context, const XML_Char *name, const XML_Char **atts )
{
    int is_key = 1;
    PContext *ctxt;
    
    ctxt = (PContext*)context;
    fprintf(ctxt->out, "%*s%c%s\n", ctxt->depth, "", '>', name);
    ctxt->level->push_back(new std::string(name)); /* copy name into level[key], need to free later */

    for (std::vector<std::string *>::size_type i = 0; i < ctxt->level->size(); i++)
    {
	fprintf(ctxt->out, " %s", ctxt->level->at(i)->c_str());
    }
    fprintf(ctxt->out, "\n");

    while (*atts)
    {
	if (is_key)
	{
	    fprintf(ctxt->out, "%*c%s: ",
		    ctxt->depth + ctxt->attrInd, ' ', *atts);
	}
	else
	{
	    fprintf(ctxt->out, "%s\n", *atts);
	}
	is_key = !is_key;
	++atts;
    }
    ctxt->depth += ctxt->tagInd;
}
 
/* callback for end elements, e.g. </tag>,
* it is called for empty elements, too
*/
static void XMLCALL
endElementCallback( void *context, const XML_Char *name __attribute__((__unused__)) )
{
    PContext *ctxt = (PContext*)context;
    ctxt->depth -= ctxt->tagInd;
    ctxt->level->pop_back();
}

/**
 * Default constructor
 */
Parser::Parser (const char* geoxml) {
    XML_Parser parser;
    PContext ctxt;
    char c;
    FILE* geofile;
    
    geofile = fopen(geoxml, "r");

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
    ctxt.tagInd = 1;
    ctxt.attrInd = 1;
    ctxt.level = new std::vector<std::string *>();

    XML_SetUserData(parser, &ctxt);
 
/* set callback for start element */
    XML_SetStartElementHandler(parser, &startElementCallback);
 
/* set callback for start element */
    XML_SetEndElementHandler(parser, &endElementCallback);
 
/* If you'd like to read input by large blocks, you can have a look at
 * XML_GetBuffer and XML_ParseBuffer functions.
 */
    while( EOF != (c = fgetc(geofile)) )
    {
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
