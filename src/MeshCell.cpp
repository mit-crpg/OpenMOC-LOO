/*
 * MeshCell.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "MeshCell.h"

MeshCell::MeshCell(){}

MeshCell::~MeshCell(){}

double MeshCell::getWidth(){
	return _width;
}

double MeshCell::getHeight(){
	return _height;
}

void MeshCell::setWidth(double width){
	_width = width;
}

void MeshCell::setHeight(double height){
	_height = height;
}

std::list<int>& MeshCell::getFSRs(){
	return _FSRs;
}

void MeshCell::addFSR(int fsr){
	_FSRs.push_back(fsr);
}


