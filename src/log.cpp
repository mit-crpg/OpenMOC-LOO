/*
 * log.cpp
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#define LOG_C

#include "log.h"

/* Level-based logging */

logLevel log_level = NORMAL;

void log_setlevel(logLevel newlevel) {
    log_level = newlevel;
}

void log_printf(logLevel level, const char *format, ...) {
    if (level >= log_level) {
    	va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
    }
    if (level == ERROR) {
    	printf("\nExiting program\n");
    	exit(1);
    }
}

