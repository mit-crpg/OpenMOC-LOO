/*
 * utility.cpp
 *
 *  Created on: Jan 13, 2012
 *      Author: will
 */

#include "utility.h"

// The formatting logic for the severity level
template< typename CharT, typename TraitsT >
inline std::basic_ostream< CharT, TraitsT >& operator<< (
    std::basic_ostream< CharT, TraitsT >& strm, severity_level lvl)
{
    static const char* const str[] =
    {
        "normal",
        "notification",
        "warning",
        "result",
        "error",
    };
    if (static_cast< std::size_t >(lvl) < (sizeof(str) / sizeof(*str)))
        strm << str[lvl];
    else
        strm << static_cast< int >(lvl);
    return strm;
}

void initLogger(char* logging_level) {
	severity_level level;
	// Second argument must be the logging level for the console
	if (strcmp(logging_level, "normal") == 0)
		level = normal;
	else if (strcmp(logging_level, "notification") == 0)
		level = notification;
	else if (strcmp(logging_level, "warning") == 0)
		level = warning;
	else if (strcmp(logging_level, "result") == 0)
		level = result;
	else if (strcmp(logging_level, "error") == 0)
		level = error;
	// Default logging level
	else
		level = result;


	// Initialize logging to the console
    logging::init_log_to_console(std::clog, keywords::filter = flt::attr< severity_level >("Severity", std::nothrow) >= level);

    // Initialize logging to the file
    logging::init_log_to_file
    (
        "sample.log",
        keywords::filter = flt::attr< severity_level >("Severity", std::nothrow) >= normal,
        keywords::format = fmt::format("[%1%] %2%")
    		% fmt::attr< severity_level >("Severity", std::nothrow)
            % fmt::message()
    );

    return;
}
