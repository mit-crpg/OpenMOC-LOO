/*
 * log.h
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <iostream>

/**
 * Level-based logging
 * Rough idea of levels:
 * 0 - NORMAL
 * 1 - INFO
 * 2 - WARNING
 * 3 - CRITICAL
 * 4 - ERROR
 * 5 - FATAL
 */

#ifndef LOG_C
	extern int log_level;
#endif


typedef enum logLevels {
	NORMAL,
	INFO,
	WARNING,
	CRITICAL,
	ERROR
} logLevel;

void log_setlevel(int newlevel);
void log_printf(logLevels level, const char *format, ...);

#define LOG(LOG_LEVEL, ...) log_print(LOG_LEVEL, __VA_ARGS__)

#endif /* LOG_H_ */
