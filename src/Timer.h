/*
 * Timer.h
 *
 *  Created on: Jan 2, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <time.h>
#include <sys/time.h>


class Timer {
protected:
	#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	 	 timeval start_time, end_time;
	#else
		timespec start_time, end_time;
	#endif

	double elapsed_time;
	bool running;
public:
	Timer();
	virtual ~Timer();
	void start();
	void stop();
	void restart();
	void reset();
	double getTime();
	#ifdef __MACH__
		double diff(timeval start, timeval end);
	#else
		double diff(timespec start, timespec end);
	#endif
};

#endif /* TIMER_H_ */
