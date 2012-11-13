/*
 * MeshCell.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "MeshCell.h"

MeshCell::MeshCell(){
	std::vector<int> _FSRs;
	_volume = 0;
	_bounds = new double[4];
	_cell_id = 0;
}

MeshCell::~MeshCell(){}

void MeshCell::makeSurfaces(int numAzim){

	/* allocate memory for mesh surfaces */
	_mesh_surfaces = new MeshSurface[8];

	for (int i = 0; i < 8; i++){
		/* make array of currents */
		_mesh_surfaces[i].makeCurrents(numAzim);
	}

	/* set side type */
	_mesh_surfaces[0].setType(SIDEX);		/* left */
	_mesh_surfaces[1].setType(SIDEY);		/* bottom */
	_mesh_surfaces[2].setType(SIDEX);		/* right */
	_mesh_surfaces[3].setType(SIDEY);		/* top */
	_mesh_surfaces[4].setType(CORNER);		/* bottom left */
	_mesh_surfaces[5].setType(CORNER);		/* bottom right */
	_mesh_surfaces[6].setType(CORNER);		/* top right */
	_mesh_surfaces[7].setType(CORNER);		/* top left */

	/* set outward normals to sides */
	/* for corners, normal is set to the normal of surface - 4 */
	_mesh_surfaces[0].setNormal(0.0);
	_mesh_surfaces[1].setNormal(PI / 2.0);
	_mesh_surfaces[2].setNormal(0.0);
	_mesh_surfaces[3].setNormal(PI / 2.0);
	_mesh_surfaces[4].setNormal(PI);
	_mesh_surfaces[5].setNormal(PI);
	_mesh_surfaces[6].setNormal(0.0);
	_mesh_surfaces[7].setNormal(0.0);

	/* set surface id */
	_mesh_surfaces[0].setId(0);
	_mesh_surfaces[1].setId(1);
	_mesh_surfaces[2].setId(2);
	_mesh_surfaces[3].setId(3);
	_mesh_surfaces[4].setId(4);
	_mesh_surfaces[5].setId(5);
	_mesh_surfaces[6].setId(6);
	_mesh_surfaces[7].setId(7);


	for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
		_chi[e] = 0.0;
	}

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


MeshSurface* MeshCell::findSurface(LocalCoords* coord, int i){
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
		log_printf(DEBUG, "coord (%f, %f) on surface: %i, cell: %i -> bounds[%f,%f,%f,%f]", x, y, surface, i, _bounds[0],_bounds[1],_bounds[2],_bounds[3]);
	}

	return meshSurface;
}

MeshSurface* MeshCell::getMeshSurfaces(int surface){
	return &_mesh_surfaces[surface];
}

void MeshCell::setChi(double chi, int e){
	_chi[e] = chi;
}

double* MeshCell::getChi(){
	return _chi;
}

void MeshCell::setSigmaA(double sigmaA, int e){
	_sigma_a[e] = sigmaA;
}

double* MeshCell::getSigmaA(){
	return _sigma_a;
}

double* MeshCell::getSigmaS(){
	return _sigma_s;
}

void MeshCell::setSigmaS(double sigmaS, int e, int ep){
	_sigma_s[e*NUM_ENERGY_GROUPS + ep] = sigmaS;
}

double* MeshCell::getNuSigmaF(){
	return _nu_sigma_f;
}

void MeshCell::setNuSigmaF(double nuSigmaF, int e){
	_nu_sigma_f[e] = nuSigmaF;
}

double* MeshCell::getDiffusivity(){
	return _diffusivity;
}

void MeshCell::setDiffusivity(double diffusivity, int e){
	_diffusivity[e] = diffusivity;
}

double* MeshCell::getOldFlux(){
	return _old_flux;
}

void MeshCell::setOldFlux(double flux, int e){
	_old_flux[e] = flux;
}

double* MeshCell::getNewFlux(){
	return _new_flux;
}

void MeshCell::setNewFlux(double flux, int e){
	_new_flux[e] = flux;
}

double MeshCell::getVolume(){
	return _volume;
}

void MeshCell::setVolume(double volume){
	_volume = volume;
}


