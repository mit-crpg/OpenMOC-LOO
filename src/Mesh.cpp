/*
 * Mesh.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "Mesh.h"

Mesh::Mesh(int width, int height){
	_width  = width;
	_height = height;
	_cells = new MeshCell[width * height];
}

Mesh::~Mesh(){
	delete [] _cells;
}

int Mesh::getWidth(){
	return _width;
}

int Mesh::getHeight(){
	return _height;
}

void Mesh::setWidth(int width){
	_width = width;
}

void Mesh::setHeight(int height){
	_width = height;
}

MeshCell* Mesh::getCells(){
	return _cells;
}


