/*
 * utility.h
 *
 *  Created on: Jan 13, 2012
 *      Author: will
 */

#ifndef UTILITY_
#define UTILITY_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/common.hpp>
#include <boost/log/formatters.hpp>
#include <boost/log/filters.hpp>
#include <boost/log/utility/init/to_file.hpp>
#include <boost/log/utility/init/to_console.hpp>
#include <boost/log/sources/global_logger_storage.hpp>

#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>


namespace logging = boost::log;
namespace fmt = boost::log::formatters;
namespace flt = boost::log::filters;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;

enum severity_level
{
    normal,
    notification,
    warning,
    result,
    error
};

BOOST_LOG_DECLARE_GLOBAL_LOGGER(test_lg, src::severity_logger_mt< severity_level >)

template< typename CharT, typename TraitsT >
inline std::basic_ostream< CharT, TraitsT >& operator<< (
    std::basic_ostream< CharT, TraitsT >& strm, severity_level lvl);

void initLogger(char* logging_level);


#endif /* UTILITY_ */
