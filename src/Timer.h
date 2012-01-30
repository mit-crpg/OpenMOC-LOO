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

#ifdef __MACH__		/* For OSX */
#define timespec timeval
#endif


class Timer {
protected:
	timespec start_time, end_time;
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
	double diff(timespec start, timespec end);
};

#endif /* TIMER_H_ */
