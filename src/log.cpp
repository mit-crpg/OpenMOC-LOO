/*
 * log.cpp
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "log.h"

/* Level-based logging */

void log_setlevel(int newlevel) {
    log_level = newlevel;
}

void log_print(int level, const char *format, ...) {
    if (level >= log_level) {
    	va_list args;

		va_start(args, format);
		vprintf(format, args);
		va_end(args);
    }
}
