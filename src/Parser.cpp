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
#include <assert.h>

/* Verbose debugging of the parser, but doesn't follow log* format */
//#define DEBUG

/* These should really be static, but thanks to C++ that's impossible.  I've
 * still declared them up here though, as at least they can be made private!
 */
void XMLCALL Parser_XMLCallback_Start(void *context,
				      const XML_Char *name,
				      const XML_Char **atts);
void XMLCALL Parser_XMLCallback_End(void *context,
				    const XML_Char *name);
void XMLCALL Parser_XMLCallback_CData(void *context,
				      const XML_Char *s,
				      int len);

enum frame_type {
	NODE_TYPE_NONE,
	NODE_TYPE_GEOMETRY,
	NODE_TYPE_MIN = NODE_TYPE_GEOMETRY,
	NODE_TYPE_CELL,
	NODE_TYPE_LATTICE,
	NODE_TYPE_TYPE,
	NODE_TYPE_DIMENSION,
	NODE_TYPE_ORIGIN,
	NODE_TYPE_WIDTH,
	NODE_TYPE_UNIVERSES,
	NODE_TYPE_SURFACE,
	NODE_TYPE_MAX = NODE_TYPE_SURFACE /* Keep this in sync */
};

struct frame_geometry {
};

struct frame_cell {
	bool has_id;
	int id;

	bool has_fill;
	int fill;

	bool has_material;
	int material;

	bool has_universe;
	int universe;

	int surfaces_count;
	int *surfaces;
};

struct frame_ttype {
	char *data;
};

struct frame_dimension {
	char *data;
};

struct frame_origin {
	char *data;
};

struct frame_width {
	char *data;
};

struct frame_universes {
	char *data;
};

struct frame_lattice {
	bool has_id;
	int id;

	char *type;

	int *dimmensions;
	int dimmensions_count;

	int *universes;
	int universes_count;

	double *origin;
	int origin_count;

	double *width;
	int width_count;
};

struct frame_surface {
	bool has_id;
	int id;

	char *type;

	double *coeffs;
	int coeffs_count;
};

struct frame {
	struct frame *parent;
	enum frame_type type;
	unsigned int depth;

	union {
		struct frame_geometry geometry;
		struct frame_cell cell;
		struct frame_lattice lattice;
		struct frame_ttype ttype;
		struct frame_dimension dimmension;
		struct frame_origin origin;
		struct frame_width width;
		struct frame_universes universes;
		struct frame_surface surface;
	};
};

struct stack {
	struct frame *top;
	Parser *parser;
};

static inline const char *frame_type_string(enum frame_type type);

static inline struct frame *stack_push(struct stack *s, enum frame_type type);
static inline struct frame *stack_pop(struct stack *s);
static inline void stack_print(struct frame *f);
static inline void stack_print_help(struct frame *f);

static inline int *strtok_int(const char *str, int *count);
static inline double *strtok_double(const char *str, int *count);

static inline char *astrncat(char *orig, char *next, int len);

Parser::Parser (const Options *opts) {
	FILE* geofile;
	XML_Parser parser;
	struct stack stack;
	char c;

	/* Assures that the geometry file exists and is readable */
	geofile = fopen(opts->getGeometryFile(), "r");
	if (geofile == NULL) {
		log_printf(ERROR, "Given geometry file %s does not exist",
			   opts->getGeometryFile());
	}

	/* Sets up the parser */
	stack.top = NULL;
	stack.parser = this;
	parser = XML_ParserCreate(NULL); /* NULL -> system encoding */
	XML_SetUserData(parser, &stack);
	XML_SetStartElementHandler(parser, &Parser_XMLCallback_Start);
	XML_SetEndElementHandler(parser, &Parser_XMLCallback_End);
	XML_SetCharacterDataHandler(parser, &Parser_XMLCallback_CData);

	/* Passes single characters to the parser, which is quite slow but
	 * is the easiest for now. */
	while( EOF != (c = fgetc(geofile)) ) {
		if (XML_Parse(parser, &c, 1, false) != XML_STATUS_OK)
			log_printf(ERROR, "Expat error\n");
        }

	/* Tells the parse we've reached the end */
	XML_Parse(parser, NULL, 0, true);
	XML_ParserFree(parser);
	fclose(geofile);
}

Parser::~Parser() {
	unsigned int i;

	for (i = 0; i < this->surfaces.size(); i++)
		delete this->surfaces.at(i);
	for (i = 0; i < this->cells.size(); i++)
		delete this->cells.at(i);
	for (i = 0; i < this->lattices.size(); i++)
		delete this->lattices.at(i);
}

void Parser::each_surface(std::function<void(Surface *)> callback) {
	std::vector<Surface *>::size_type i;

	for (i = 0; i < this->surfaces.size(); i++)
		callback(this->surfaces.at(i));
}

void Parser::each_cell(std::function<void(Cell *)> callback) {
	std::vector<Cell *>::size_type i;

	for (i = 0; i < this->cells.size(); i++)
		callback(this->cells.at(i));
}

void Parser::each_lattice(std::function<void(Lattice *)> callback) {
	std::vector<Lattice *>::size_type i;

	for (i = 0; i < this->lattices.size(); i++)
		callback(this->lattices.at(i));
}

void XMLCALL Parser_XMLCallback_Start(void *context,
				      const XML_Char *name,
				      const XML_Char **attrs) {
	struct stack *s;
	struct frame *f;
	enum frame_type type;
	int i;

	s = (struct stack *)context;

	/* Checks what type of node this is */
	type = NODE_TYPE_NONE;
	for (i = (int)NODE_TYPE_MIN; i <= (int)NODE_TYPE_MAX; i++)
		if (strcmp(name, frame_type_string((enum frame_type)i)) == 0)
			type = (enum frame_type)i;

	/* Ensures that we know what type the node is */
	if (type == NODE_TYPE_NONE)
		log_printf(ERROR, "Unknown node type '%s'\n", name);

	/* Adds our item to the stack */
	f = stack_push(s, type);

	/* Parses every attribute */
	while (*attrs != NULL) {
		char *key, *value;
		
		/* Attributes are stored as a key-value pair, there are always
		 * two of them (one for the key, one for the value) so we can
		 * safely double-increment here.
		 */
		/* FIXME: Verify that a bad input file can't cause an odd
		 *        number of attributes.
		 */
		assert(sizeof(char) == sizeof(XML_Char));
		key = (char *)*attrs;
		attrs++;
		value = (char *)*attrs;
		attrs++;
		assert(key != NULL);
		assert(value != NULL);

		/* Does some type-specific parsing for some attributes */
		switch (f->type) {
		case NODE_TYPE_NONE:
			break;
		case NODE_TYPE_GEOMETRY:
			break;
		case NODE_TYPE_CELL:
			if (strcmp(key, "id") == 0) {
				if (f->cell.has_id == true)
					log_printf(ERROR, "Has 2 ids\n");

				f->cell.has_id = true;
				f->cell.id = atoi(value);
			} else if (strcmp(key, "fill") == 0) {
				if (f->cell.has_fill == true)
					log_printf(ERROR, "Has 2 fills\n");

				if (f->cell.has_material == true) {
					log_printf(ERROR,
						   "Has material and fill\n");
				}

				f->cell.has_fill = true;
				f->cell.fill = atoi(value);
			} else if (strcmp(key, "material") == 0) {
				if (f->cell.has_fill == true)
					log_printf(ERROR, "Has 2 material\n");

				if (f->cell.has_fill == true) {
					log_printf(ERROR,
						   "Has material and fill\n");
				}

				f->cell.has_material = true;
				f->cell.material = atoi(value);
			} else if (strcmp(key, "universe") == 0) {
				if (f->cell.has_universe == true)
					log_printf(ERROR, "Has 2 universes\n");

				f->cell.has_universe = true;
				f->cell.universe = atoi(value);
			} else if (strcmp(key, "surfaces") == 0) {
				if (f->cell.surfaces != NULL)
					log_printf(ERROR, "Has 2 surfaces\n");

				f->cell.surfaces =
					strtok_int(value,
						   &f->cell.surfaces_count);
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'\n",
					   key, value);
			}
			break;
		case NODE_TYPE_LATTICE:
			if (strcmp(key, "id") == 0) {
				if (f->lattice.has_id == true)
					log_printf(ERROR, "Has 2 ids\n");

				f->lattice.has_id = true;
				f->lattice.id = atoi(value);
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'\n",
					   key, value);
			}

			break;
		case NODE_TYPE_TYPE:
			break;
		case NODE_TYPE_DIMENSION:
		break;
		case NODE_TYPE_ORIGIN:
			break;
		case NODE_TYPE_WIDTH:
			break;
		case NODE_TYPE_UNIVERSES:
			break;
		case NODE_TYPE_SURFACE:
			if (strcmp(key, "id") == 0) {
				if (f->surface.has_id == true)
					log_printf(ERROR, "Has 2 ids\n");

				f->surface.has_id = true;
				f->surface.id = atoi(value);
			} else if (strcmp(key, "type") == 0) {
				if (f->surface.type != NULL)
					log_printf(ERROR, "Has 2 types\n");

				f->surface.type = strdup(value);
			} else if (strcmp(key, "coeffs") == 0) {
				if (f->surface.coeffs != NULL)
					log_printf(ERROR, "Has 2 coeffs\n");

				f->surface.coeffs =
					strtok_double(value,
						      &f->surface.coeffs_count);
			} else {
				log_printf(ERROR, "Unknown attribute '%s=%s'\n",
					   key, value);
			}
			break;
		}
	}
}

void XMLCALL Parser_XMLCallback_End(void *context,
				    const XML_Char *name) {
	struct stack *s;
	struct frame *f, *p;
	
	s = (struct stack *)context;
	f = stack_pop(s);
	p = s->top;

#ifdef DEBUG
	/* Prints out the stack */
	stack_print(f);
#endif

	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_CELL:
	{
		Cell *cell;

		cell = NULL;
		if (f->cell.has_fill) {
			cell = new CellFill(f->cell.id, 
					    f->cell.universe,
					    f->cell.surfaces_count,
					    f->cell.surfaces,
					    f->cell.fill);
		} else if (f->cell.has_material) {
			cell = new CellBasic(f->cell.id, 
					    f->cell.universe,
					    f->cell.surfaces_count,
					    f->cell.surfaces,
					    f->cell.material);
		} else {
			log_printf(ERROR, "Cell without material or fill\n");
		}

		if (cell != NULL)
			s->parser->cells.push_back(cell);
		else
			log_printf(ERROR, "Unknown cell type\n");

		if (f->cell.surfaces != NULL)
			free(f->cell.surfaces);
		break;
	}
	case NODE_TYPE_LATTICE:
		Lattice *lattice;

		lattice = NULL;
		if (f->lattice.has_id != true)
			log_printf(ERROR, "Lattice without id\n");
		if (f->lattice.dimmensions_count != 2)
			log_printf(ERROR, "Lattice without exactly 2 dimms\n");
		if (f->lattice.origin_count != 2)
			log_printf(ERROR, "Lattice without exactly 2 origin\n");
		if (f->lattice.width_count != 2)
			log_printf(ERROR, "Lattice without exactly 2 widths\n");
		if (f->lattice.universes == NULL)
			log_printf(ERROR, "Lattice without universes\n");

		lattice = new Lattice(f->lattice.id,
				      f->lattice.dimmensions[0],
				      f->lattice.dimmensions[1],
				      f->lattice.origin[0],
				      f->lattice.origin[1],
				      f->lattice.width[0],
				      f->lattice.width[1],
				      f->lattice.universes_count,
				      f->lattice.universes);

		s->parser->lattices.push_back(lattice);

		if (f->lattice.type != NULL)
			free(f->lattice.type);
		if (f->lattice.dimmensions != NULL)
			free(f->lattice.dimmensions);
		if (f->lattice.universes != NULL)
			free(f->lattice.universes);
		if (f->lattice.origin != NULL)
			free(f->lattice.origin);
		if (f->lattice.width != NULL)
			free(f->lattice.width);
		break;
	case NODE_TYPE_TYPE:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			p->lattice.type = f->ttype.data;
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_ORIGIN:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
			log_printf(ERROR, "Unexpected type subfield\n");
		}
		break;
	case NODE_TYPE_DIMENSION:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.dimmensions != NULL)
				log_printf(ERROR, "Has 2 dimmensions\n");

			p->lattice.dimmensions =
				strtok_int(f->dimmension.data,
					   &p->lattice.dimmensions_count);
			
			free(f->dimmension.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_ORIGIN:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
			log_printf(ERROR, "Unexpected dimmension subfield\n");
		}
		break;
	case NODE_TYPE_ORIGIN:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.origin != NULL)
				log_printf(ERROR, "Has 2 origins\n");

			p->lattice.origin =
				strtok_double(f->origin.data,
					      &p->lattice.origin_count);
			
			free(f->origin.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_ORIGIN:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
			log_printf(ERROR, "Unexpected dimmension subfield\n");
		}
		break;
	case NODE_TYPE_WIDTH:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.width != NULL)
				log_printf(ERROR, "Has 2 widths\n");

			p->lattice.width =
				strtok_double(f->width.data,
					      &p->lattice.width_count);
			
			free(f->width.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_ORIGIN:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
			log_printf(ERROR, "Unexpected dimmension subfield\n");
		}
		break;
	case NODE_TYPE_UNIVERSES:
		switch (p->type) {
		case NODE_TYPE_LATTICE:
			if (p->lattice.universes != NULL)
				log_printf(ERROR, "Has 2 universes\n");

			p->lattice.universes =
				strtok_int(f->universes.data,
					   &p->lattice.universes_count);
			
			free(f->universes.data);
			break;
		case NODE_TYPE_NONE:
		case NODE_TYPE_GEOMETRY:
		case NODE_TYPE_CELL:
		case NODE_TYPE_TYPE:	
		case NODE_TYPE_DIMENSION:
		case NODE_TYPE_ORIGIN:
		case NODE_TYPE_WIDTH:
		case NODE_TYPE_UNIVERSES:
		case NODE_TYPE_SURFACE:
			log_printf(ERROR, "Unexpected universes subfield\n");
		}
		break;
	case NODE_TYPE_SURFACE:
	{
		Surface *surface;

		surface = NULL;
		if (strcmp(f->surface.type, "plane") == 0) {
			if (f->surface.coeffs_count != 3)
				log_printf(ERROR, "Wrong number of coeffs\n");

			surface = new Plane(f->surface.id, 
					    f->surface.coeffs[0],
					    f->surface.coeffs[1],
					    f->surface.coeffs[2]);
		} else if (strcmp(f->surface.type, "x-plane") == 0) {
			if (f->surface.coeffs_count != 1)
				log_printf(ERROR, "Wrong number of coeffs\n");

			surface = new XPlane(f->surface.id, 
					     f->surface.coeffs[0]);
		} else if (strcmp(f->surface.type, "y-plane") == 0) {
			if (f->surface.coeffs_count != 1)
				log_printf(ERROR, "Wrong number of coeffs\n");

			surface = new YPlane(f->surface.id, 
					     f->surface.coeffs[0]);
		} else if (strcmp(f->surface.type, "circle") == 0) {
			if (f->surface.coeffs_count != 3)
				log_printf(ERROR, "Wrong number of coeffs\n");

			surface = new Circle(f->surface.id, 
					     f->surface.coeffs[0],
					     f->surface.coeffs[1],
					     f->surface.coeffs[2]);
		}
		
		if (surface != NULL)
			s->parser->surfaces.push_back(surface);
		else {
			log_printf(ERROR, "Unknown surface type '%s'\n",
				   f->surface.type);
		}

		if (f->surface.type != NULL)
			free(f->surface.type);
		if (f->surface.coeffs != NULL)
			free(f->surface.coeffs);
		break;
	}
	}

	free(f);
}

void XMLCALL Parser_XMLCallback_CData(void *context,
				      const XML_Char *str_uncast,
				      int len) {
	struct stack *s;
	struct frame *f;
	char *str;

	str = (char *)str_uncast;
	s = (struct stack *)context;
	f = s->top;

	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_CELL:
		break;
	case NODE_TYPE_LATTICE:
		break;
	case NODE_TYPE_TYPE:
		break;
	case NODE_TYPE_DIMENSION:
		f->dimmension.data = astrncat(f->dimmension.data, str, len);
		break;
	case NODE_TYPE_ORIGIN:
		f->origin.data = astrncat(f->origin.data, str, len);
		break;
	case NODE_TYPE_WIDTH:
		f->width.data = astrncat(f->width.data, str, len);
		break;
	case NODE_TYPE_UNIVERSES:
		f->universes.data = astrncat(f->universes.data, str, len);
		break;
	case NODE_TYPE_SURFACE:
		break;
	}
}

const char *frame_type_string(enum frame_type type) {
	switch (type) {
	case NODE_TYPE_NONE:
		return "none";
	case NODE_TYPE_GEOMETRY:
		return "geometry";
	case NODE_TYPE_CELL:
		return "cell";
	case NODE_TYPE_LATTICE:
		return "lattice";
	case NODE_TYPE_TYPE:
		return "type";
	case NODE_TYPE_DIMENSION:
		return "dimension";
	case NODE_TYPE_ORIGIN:
		return "origin";
	case NODE_TYPE_WIDTH:
		return "width";
	case NODE_TYPE_UNIVERSES:
		return "universes";
	case NODE_TYPE_SURFACE:
		return "surface";
	}
	
	abort();
	return NULL;
}

struct frame *stack_push(struct stack *s, enum frame_type type) {
	struct frame *f;

	/* Allocates a new stack frame (it gets added way down at the end) */
	f = (struct frame *)malloc(sizeof(*f));
	if (f == NULL)
		log_printf(ERROR, "malloc returned NULL!\n");

	/* Different node types get initialized differently */
	f->type = type;
	switch (type) {
	case NODE_TYPE_NONE:
		free(f);
		log_printf(ERROR, "Tried to push an unknown node type\n");
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_CELL:
		f->cell.has_id = false;
		f->cell.has_fill = false;
		f->cell.has_material = false;
		f->cell.has_universe = false;
		f->cell.surfaces = NULL;
		f->cell.surfaces_count = -1;
		break;
	case NODE_TYPE_LATTICE:
		f->lattice.has_id = false;
		f->lattice.type = NULL;
		f->lattice.dimmensions = NULL;
		f->lattice.universes = NULL;
		f->lattice.origin = NULL;
		f->lattice.width = NULL;
		break;
	case NODE_TYPE_TYPE:
		f->ttype.data = NULL;
		break;
	case NODE_TYPE_DIMENSION:
		f->dimmension.data = NULL;
		break;
	case NODE_TYPE_ORIGIN:
		f->origin.data = NULL;
		break;
	case NODE_TYPE_WIDTH:
		f->width.data = NULL;
		break;
	case NODE_TYPE_UNIVERSES:
		f->universes.data = NULL;
		break;
	case NODE_TYPE_SURFACE:
		f->surface.has_id = false;
		f->surface.type = NULL;
		f->surface.coeffs = NULL;
		break;
	}

	/* We always have one depth larger than our parent */
	if (s->top == NULL)
		f->depth = 0;
	else
		f->depth = s->top->depth + 1;

	/* Actually adds this to the stack */
	f->parent = s->top;
	s->top = f;
	return f;
}

struct frame *stack_pop(struct stack *s) {
	struct frame *f;

	f = s->top;
	if (s->top != NULL)
		s->top = s->top->parent;

	return f;
}

void stack_print(struct frame *f) {
	stack_print_help(f);
	fprintf(stderr, "\n");
}

void stack_print_help(struct frame *f) {
	unsigned int i;

	if (f == NULL)
		return;

	stack_print_help(f->parent);
	
	for (i = 0; i < f->depth; i++)
		fprintf(stderr, " ");
	fprintf(stderr, "%s", frame_type_string(f->type));
	
	switch (f->type) {
	case NODE_TYPE_NONE:
		break;
	case NODE_TYPE_GEOMETRY:
		break;
	case NODE_TYPE_CELL:
		if (f->cell.has_id)
			fprintf(stderr, " id=\"%d\"", f->cell.id);
		if (f->cell.has_material)
			fprintf(stderr, " material=\"%d\"", f->cell.material);
		if (f->cell.has_fill)
			fprintf(stderr, " fill=\"%d\"", f->cell.fill);
		if (f->cell.has_universe)
			fprintf(stderr, " universe=\"%d\"", f->cell.universe);
		if (f->cell.surfaces != NULL) {
			int i;

			fprintf(stderr, " surfaces=\"");
			for (i = 0; i < f->cell.surfaces_count; i++)
				fprintf(stderr, " %d", f->cell.surfaces[i]);
			fprintf(stderr, "\"");
		}

		break;
	case NODE_TYPE_LATTICE:
		if (f->lattice.has_id)
			fprintf(stderr, " id=\"%d\"", f->lattice.id);
		if (f->lattice.type != NULL)
			fprintf(stderr, " type=\"%s\"", f->lattice.type);
		if (f->lattice.dimmensions != NULL) {
			int i;

			fprintf(stderr, " dimmensions=\"");
			for (i = 0; i < f->lattice.dimmensions_count; i++) {
				fprintf(stderr, " %d",
					f->lattice.dimmensions[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->lattice.universes != NULL) {
			int i;

			fprintf(stderr, " universes=\"");
			for (i = 0; i < f->lattice.universes_count; i++) {
				fprintf(stderr, " %d",
					f->lattice.universes[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->lattice.width != NULL) {
			int i;

			fprintf(stderr, " width=\"");
			for (i = 0; i < f->lattice.width_count; i++) {
				fprintf(stderr, " %f",
					f->lattice.width[i]);
			}
			fprintf(stderr, "\"");
		}
		if (f->lattice.origin != NULL) {
			int i;

			fprintf(stderr, " origin=\"");
			for (i = 0; i < f->lattice.origin_count; i++) {
				fprintf(stderr, " %f",
					f->lattice.origin[i]);
			}
			fprintf(stderr, "\"");
		}
		break;
	case NODE_TYPE_TYPE:
		break;
	case NODE_TYPE_DIMENSION:
		break;
	case NODE_TYPE_ORIGIN:
		break;
	case NODE_TYPE_WIDTH:
		break;
	case NODE_TYPE_UNIVERSES:
		break;
	case NODE_TYPE_SURFACE:
		if (f->surface.has_id)
			fprintf(stderr, " id=\"%d\"", f->surface.id);
		if (f->surface.type != NULL)
			fprintf(stderr, " type=\"%s\"", f->surface.type);
		if (f->surface.coeffs != NULL) {
			int i;

			fprintf(stderr, " coeffs=\"");
			for (i = 0; i < f->surface.coeffs_count; i++) {
				fprintf(stderr, " %f",
					f->surface.coeffs[i]);
			}
			fprintf(stderr, "\"");
		}
		break;
	}

	fprintf(stderr, "\n");
}

int *strtok_int(const char *str, int *count) {
	int *arr;
	int cnt;
	char *duplicated;
	char *st_tmp;
	char *tok;
	int i;
			
	cnt = 0;
	duplicated = strdup(str);
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		cnt++;
	}
	free(duplicated);
	
	arr = (int *) malloc(sizeof(*arr) * cnt);
	
	duplicated = strdup(str);
	i = 0;
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		assert(i < cnt);
		arr[i] = atoi(tok);
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		i++;
	}
	free(duplicated);

	*count = cnt;
	return arr;
}

double *strtok_double(const char *str, int *count) {
	double *arr;
	int cnt;
	char *duplicated;
	char *st_tmp;
	char *tok;
	int i;
			
	cnt = 0;
	duplicated = strdup(str);
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		cnt++;
	}
	free(duplicated);
	
	arr = (double *) malloc(sizeof(*arr) * cnt);
	
	duplicated = strdup(str);
	i = 0;
	tok = strtok_r(duplicated, " \n\t", &st_tmp);
	while (tok != NULL) {
		assert(i < cnt);
		arr[i] = atof(tok);
		tok = strtok_r(NULL, " \n\t", &st_tmp);
		i++;
	}
	free(duplicated);

	*count = cnt;
	return arr;
}

char *astrncat(char *orig, char *str, int len) {
	char *new_data;

	if (orig == NULL) {
		new_data = strndup(str, len);
	} else {
		int olen;

		olen = strlen(orig);
		new_data = (char *) realloc(orig, olen + len + 1);
		new_data[olen+1] = '\0';
		strncat(new_data, str, len);
	}

	return new_data;
}
