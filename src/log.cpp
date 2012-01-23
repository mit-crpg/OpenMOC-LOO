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

int log_level = 0;

void log_setlevel(int newlevel) {
    log_level = newlevel;
}

void log_printf(logLevels level, const char *format, ...) {
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

