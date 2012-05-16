/*
 * MeshCell.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "MeshCell.h"

MeshCell::MeshCell(){
	std::vector<int> _FSRs;
	_abs_rate = 0;
	_mesh_surfaces = new MeshSurface[8];
	_bounds = new double[4];

	makeSurfaces();
}

MeshCell::~MeshCell(){}

void MeshCell::makeSurfaces(){

	/* set sides */
	_mesh_surfaces[0].setType(SIDEX);
	_mesh_surfaces[1].setType(SIDEY);
	_mesh_surfaces[2].setType(SIDEX);
	_mesh_surfaces[3].setType(SIDEY);
	_mesh_surfaces[4].setType(CORNER);
	_mesh_surfaces[5].setType(CORNER);
	_mesh_surfaces[6].setType(CORNER);
	_mesh_surfaces[7].setType(CORNER);

	/* set normals to sides */
	_mesh_surfaces[0].setNormal(PI);
	_mesh_surfaces[1].setNormal(3.0 * PI / 2.0);
	_mesh_surfaces[2].setNormal(0.0);
	_mesh_surfaces[3].setNormal(PI / 2.0);
	_mesh_surfaces[4].setNormal(PI);
	_mesh_surfaces[5].setNormal(3.0 * PI / 2.0);
	_mesh_surfaces[6].setNormal(0.0);
	_mesh_surfaces[7].setNormal(PI / 2.0);
}

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

void MeshCell::setAbsRate(double absRate){
	_abs_rate = absRate;
}

double MeshCell::getAbsRate(){
	return _abs_rate;
}

int MeshCell::getFSRStart(){
	return _fsr_start;
}

int MeshCell::getFSREnd(){
	return _fsr_end;
}

void MeshCell::setFSRStart(int fsr){
	_fsr_start = fsr;
}

void MeshCell::setFSREnd(int fsr){
	_fsr_end = fsr;
}


std::vector<int>* MeshCell::getFSRs(){
	return &_FSRs;

}

void MeshCell::addFSR(int fsr){
	_FSRs.push_back(fsr);
}

void MeshCell::setBounds(double x, double y){
	_bounds[0] = x;
	_bounds[1] = y;
	_bounds[2] = x + _width;
	_bounds[3] = y + _height;
}

double* MeshCell::getBounds(){
	return _bounds;
}


MeshSurface* MeshCell::findSurface(LocalCoords* coord){
	MeshSurface* meshSurface = NULL;
	int surface = -1;
	double x = coord->getX();
	double y = coord->getY();

	/* find which surface coord is on */
	/* left */
	if (fabs(x - _bounds[0]) < 1e-8){
		meshSurface = &_mesh_surfaces[0];
		surface = 0;
	}
	/* bottom */
	if (fabs(y - _bounds[1]) < 1e-8){
		if (meshSurface == &_mesh_surfaces[0]){
			meshSurface = &_mesh_surfaces[4];
			surface = 4;
		}
		else{
			meshSurface = &_mesh_surfaces[1];
			surface = 1;
		}
	}
	/* right */
	if (fabs(x - _bounds[2]) < 1e-8){
		if (meshSurface == &_mesh_surfaces[1]){
			meshSurface = &_mesh_surfaces[5];
			surface = 5;
		}
		else{
			meshSurface = &_mesh_surfaces[2];
			surface = 2;
		}
	}
	/* top */
	if (fabs(y - _bounds[3]) < 1e-8){
		if (meshSurface == &_mesh_surfaces[2]){
			meshSurface = &_mesh_surfaces[6];
			surface = 6;
		}
		else if(meshSurface == &_mesh_surfaces[0]){
			meshSurface = &_mesh_surfaces[7];
			surface = 7;
		}
		else{
			meshSurface = &_mesh_surfaces[3];
			surface = 3;
		}
	}

	if (surface != -1){
		log_printf(DEBUG, "coord (%f, %f) on surface: %i -> bounds[%f,%f,%f,%f]", x, y, surface, _bounds[0],_bounds[1],_bounds[2],_bounds[3]);
	}

	return meshSurface;
}

MeshSurface* MeshCell::getMeshSurfaces(int surface){
	return &_mesh_surfaces[surface];
}

void MeshCell::setSigmaA(double sigmaA){
	_sigma_a = sigmaA;
}

double MeshCell::getSigmaA(){
	return _sigma_a;
}

void MeshCell::setSigmaT(double sigmaT){
	_sigma_t = sigmaT;
}

double MeshCell::getSigmaT(){
	return _sigma_t;
}

double MeshCell::getSigmaF(){
	return _sigma_f;
}

void MeshCell::setSigmaF(double sigmaF){
	_sigma_f = sigmaF;
}

double MeshCell::getNuSigmaF(){
	return _nu_sigma_f;
}

void MeshCell::setNuSigmaF(double nuSigmaF){
	_nu_sigma_f = nuSigmaF;
}


