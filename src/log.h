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
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdexcept>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

typedef enum logLevels {
    DEBUG,
    INFO,
    ACTIVE,
    NORMAL,
    DATA, 
    WARNING,
    CRITICAL,
    RESULT,
    ERROR
} logLevel;

void log_setlevel(logLevel newlevel);
void log_setlevel(const char* newlevel);
void log_printf(logLevel level, const char *format, ...)
    __attribute__((__format__ (__printf__, 2, 3) ));
void log_error(const char *format, ...)
    __attribute__((__noreturn__))
    __attribute__((__format__ (__printf__, 1, 2) ));

void setOutputDirectory(char* directory);
const char* getOutputDirectory();

#ifndef LOG_C
extern logLevel log_level;
#endif


#endif /* LOG_H_ */
