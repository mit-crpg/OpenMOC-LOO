/*
 * Timer.h
 *
 *  Created on: Jan 2, 2012
 *      Author: wbinventor
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

using namespace std;


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
