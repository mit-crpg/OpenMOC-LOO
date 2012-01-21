/*
 * Timer.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: wbinventor
 */

#include "Timer.h"

Timer::Timer() {
	this->running = false;
	this->elapsed_time = 0;
}

Timer::~Timer() {
}

void Timer::start() {

	if (!this->running) {
		#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
			gettimeofday(&this->start_time, NULL);
		#else
			clock_gettime(CLOCK_MONOTONIC, &this->start_time);
		#endif
		this->running = true;
	}
	return;
}

void Timer::stop() {
	if (this->running) {
		#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
			gettimeofday(&this->end_time, NULL);
		#else
		  clock_gettime(CLOCK_MONOTONIC, &this->end_time);
		#endif
		this->running = false;
		this->elapsed_time += this->diff(this->start_time, this->end_time);
	}

	return;
}

void Timer::reset() {
	this->elapsed_time = 0;
	this->running = false;
}

void Timer::restart() {
	if (!this->running) {
		this->elapsed_time += this->diff(this->start_time, this->end_time);
		this->start();
	}
}

double Timer::getTime() {
	if (!this->running) {
		#ifdef __MACH__
			return this->elapsed_time * 1.0E-6;
		#else
			return this->elapsed_time * 1.0E-9;
		#endif
	}
	else {
		#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
			timeval temp;
			gettimeofday(&this->start_time, NULL);
		#else
		  timespec temp;
		  clock_gettime(CLOCK_MONOTONIC, &this->start_time);
		#endif

		this->elapsed_time += this->diff(this->start_time, temp);

		#ifdef __MACH__
			return this->elapsed_time * 1.0E-6;
		#else
			return this->elapsed_time * 1.0E-9;
		#endif
	}
}

#ifdef __MACH__
double Timer::diff(timeval start, timeval end) {
	timeval temp;

	if ((end.tv_usec - start.tv_usec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_usec = 1.0E6 + end.tv_usec - start.tv_usec;
	} else {
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_usec = end.tv_usec - start.tv_usec;
	}

	return (temp.tv_sec * 1.0E6 + temp.tv_usec);
}

#else

double Timer::diff(timespec start, timespec end) {
	timespec temp;

	if ((end.tv_nsec - start.tv_nsec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1.0E9 + end.tv_nsec - start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}

	return (temp.tv_sec * 1.0E9 + temp.tv_nsec);
}
#endif
