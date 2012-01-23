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

typedef enum logLevels {
	NORMAL,
	INFO,
	WARNING,
	CRITICAL,
	ERROR,
	DEBUG,
	RESULT
} logLevel;

void log_setlevel(logLevel newlevel);
void log_printf(logLevel level, const char *format, ...);

#ifndef LOG_C
	extern logLevel log_level;
#endif


#endif /* LOG_H_ */
