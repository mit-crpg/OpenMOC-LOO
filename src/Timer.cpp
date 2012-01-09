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

Timer::~Timer() { }

void Timer::start() {

	if (!this->running) {
		clock_gettime(CLOCK_MONOTONIC, &this->start_time);
		this->running = true;
	}
	return;
}

void Timer::stop() {
	if (this->running) {
		clock_gettime(CLOCK_MONOTONIC, &this->end_time);
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
	if (!this->running)
		return this->elapsed_time*1.0E-9;
	else {
		timespec temp;
		clock_gettime(CLOCK_MONOTONIC, &temp);
		this->elapsed_time += this->diff(this->start_time, temp);
		return this->elapsed_time*1.0E-9;
	}
}


double Timer::diff(timespec start, timespec end) {
	timespec temp;

	if ((end.tv_nsec - start.tv_nsec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1.0E9 + end.tv_nsec - start.tv_nsec;
	}
	else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}

	return (temp.tv_sec*1.0E9 + temp.tv_nsec);
}
