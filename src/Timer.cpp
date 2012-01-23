/*
 * Timer.cpp
 *
 *  Created on: Jan 2, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 *  The timer class is for profiling code. It outputs running time in
 *  seconds but has resolution of microseconds on OSX and nanoseconds
 *  on Linux machines.
 */

#include "Timer.h"

/**
 * Timer class constructor
 */
Timer::Timer() {
	this->running = false;
	this->elapsed_time = 0;
}


/**
 * Default Timer destructor
 */
Timer::~Timer() { }


/**
 * Starts the Timer - similar to starting a stopwatch
 */
void Timer::start() {

	if (!this->running) {
		#ifdef __MACH__ 	/* For OSX */
			gettimeofday(&this->start_time, NULL);
		#else				/* For linux */
			clock_gettime(CLOCK_MONOTONIC, &this->start_time);
		#endif
		this->running = true;
	}
	return;
}


/**
 * Stops the timer - similar to stopping a stopwatch
 */
void Timer::stop() {
	if (this->running) {
		#ifdef __MACH__ /* For OSX */
			gettimeofday(&this->end_time, NULL);
		#else			/* For linux */
		  clock_gettime(CLOCK_MONOTONIC, &this->end_time);
		#endif
		this->running = false;
		this->elapsed_time += this->diff(this->start_time, this->end_time);
	}
	return;
}


/**
 * Resets the timer - similar to resetting a stopwatch.
 */
void Timer::reset() {
	this->elapsed_time = 0;
	this->running = false;
}


/**
 * Restarts the timer. The elapsed time will accumulate along with the
 * previous time(s) the timer was running. If the timer was already running
 * this function does nothing
 */
void Timer::restart() {
	if (!this->running) {
		this->elapsed_time += this->diff(this->start_time, this->end_time);
		this->start();
	}
}


/**
 * Returns the amount of time elapsed from start to stop of the timer. If the
 * timer is currently runnning, returns the time from the timer start to the present
 * time.
 * @return the elapsed time
 */
double Timer::getTime() {
	/* If the timer is not running */
	if (!this->running) {
		#ifdef __MACH__		/* For OSX */
			return this->elapsed_time * 1.0E-6;
		#else				/* For Linux */
			return this->elapsed_time * 1.0E-9;
		#endif
	}

	/* If the timer is currently running */
	else {
		#ifdef __MACH__ 	/* For OSX */
			timeval temp;
			gettimeofday(&this->start_time, NULL);
		#else				/* For Linux */
		  timespec temp;
		  clock_gettime(CLOCK_MONOTONIC, &this->start_time);
		#endif

		this->elapsed_time += this->diff(this->start_time, temp);

		#ifdef __MACH__		/* For OSX */
			return this->elapsed_time * 1.0E-6;
		#else				/* For Linux */
			return this->elapsed_time * 1.0E-9;
		#endif
	}
}

#ifdef __MACH__		/* For OSX */
/**
 * Helper function which computes the time between the values of
 * two timeval structs.
 * @param start timeval representing the start time
 * @param end timeval representing the end time
 */
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

#else			/* For Linux */

/**
 * Helper function which computes the time between the values of
 * two timespec structs.
 * @param start timespec representing the start time
 * @param end timespec representing the end time
 */
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
