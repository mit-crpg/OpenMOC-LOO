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


/* Default logging level is the lowest (most verbose) level */
logLevel log_level = NORMAL;


/**
 * Set the minimum level which will be printed to the console
 * @param newlevel the logging level
 */
void log_setlevel(logLevel newlevel) {
    log_level = newlevel;

    switch (newlevel) {
    	case NORMAL:
    		log_printf(INFO, "Logging level set to NORMAL\n");
    		break;
    	case INFO:
    		log_printf(INFO, "Logging level set to INFO\n");
    		break;
    	case WARNING:
    		log_printf(INFO, "Logging level set to WARNING\n");
    		break;
    	case CRITICAL:
    		log_printf(INFO, "Logging level set to CRITICAL\n");
    		break;
    	case ERROR:
    		log_printf(INFO, "Logging level set to ERROR\n");
    		break;
    	case DEBUG:
    		log_printf(INFO, "Logging level set to DEBUG\n");
    		break;
    	case RESULT:
    		log_printf(INFO, "Logging level set to RESULT\n");
    		break;
    }
}


/**
 * Set the logging level from a character string which must correspond
 * to one of the logLevel enum types (NORMAL, INFO, CRITICAL, WARNING,
 * ERROR, DEBUG, RESULT)
 * @param newlevel a character string loglevel
 */
void log_setlevel(char* newlevel) {

	if (strcmp("NORMAL", newlevel) == 0) {
		log_level = NORMAL;
		log_printf(INFO, "Logging level set to NORMAL\n");
	}
	else if (strcmp("INFO", newlevel) == 0) {
		log_level = INFO;
		log_printf(INFO, "Logging level set to INFO\n");
	}
	else if (strcmp("WARNING", newlevel) == 0) {
		log_level = WARNING;
		log_printf(INFO, "Logging level set to WARNING\n");
	}
	else if (strcmp("CRITICAL", newlevel) == 0) {
		log_level = CRITICAL;
		log_printf(INFO, "Logging level set to CRITICAL\n");
	}
	else if (strcmp("ERROR", newlevel) == 0) {
		log_level = ERROR;
		log_printf(INFO, "Logging level set to ERROR\n");
	}
	else if (strcmp("DEBUG", newlevel) == 0) {
		log_level = DEBUG;
		log_printf(INFO, "Logging level set to DEBUG\n");
	}
	else if (strcmp("RESULT", newlevel) == 0) {
		log_level = RESULT;
		log_printf(INFO, "Logging level set to RESULT\n");
	}

	return;
}



/**
 * Print a formatted message to the console. If logging level is
 * ERROR, this function will end the program
 * @param level the logging level for this message
 * @param *format variable list of C++ formatted i/o
 */
void log_printf(logLevel level, const char *format, ...) {
    if (level >= log_level) {
    	va_list args;

    	/* Append the log level to the message */
    	switch (level) {
    		case (NORMAL):
    			printf("[  NORMAL ]  ");
    			break;
    		case (INFO):
    			printf("[  INFO   ]  ");
    			break;
    		case (WARNING):
    			printf("[ WARNING ]  ");
    			break;
    		case (CRITICAL):
    			printf("[CRITICAL]   ");
    			break;
    		case (ERROR):
    			printf("[  ERROR  ]  ");
    			break;
    		case (DEBUG):
    			printf("[  DEBUG  ]  ");
    			break;
    		case (RESULT):
    			printf("[  RESULT ]  ");
    			break;
    	}

		va_start(args, format);
		vprintf(format, args);
		va_end(args);
    }
    if (level == ERROR) {
    	printf("\n[  EXIT   ]  Exiting program...\n");
    	abort();
    }
}
