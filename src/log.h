/*
 * log.h
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 *  Level-based logging module
 */

#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

typedef enum logLevels {
	DEBUG,
	INFO,
	ACTIVE,
	NORMAL,
	WARNING,
	CRITICAL,
	RESULT,
	ERROR
} logLevel;

void log_setlevel(logLevel newlevel);
void log_setlevel(const char* newlevel);
void log_printf(logLevel level, const char *format, ...);

#ifndef LOG_C
	extern logLevel log_level;
#endif


#endif /* LOG_H_ */
