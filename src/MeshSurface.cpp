/*
 * MeshSurface.cpp
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#include "MeshSurface.h"

MeshSurface::MeshSurface(){
}

MeshSurface::~MeshSurface(){}

meshSurfaceType MeshSurface::getType(){
	return _type;
}

void MeshSurface::setType(meshSurfaceType type){
	_type = type;
}



