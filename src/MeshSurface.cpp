/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){
	_current = new double[NUM_ENERGY_GROUPS];

	for (int i = 0; i < NUM_ENERGY_GROUPS; i++){
		_current[i] = 0;
	}
}

MeshSurface::~MeshSurface(){}

meshSurfaceType MeshSurface::getType(){
	return _type;
}

void MeshSurface::setType(meshSurfaceType type){
	_type = type;
}


void MeshSurface::setCurrent(double current, int group){
	_current[group] = current;
}

double MeshSurface::getCurrent(int group){
	return _current[group];
}


void MeshSurface::incrementCurrent(double current, int group){
	_current[group] += current;
}


