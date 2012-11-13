/*
 * MeshCell.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESHCELL_H_
#define MESHCELL_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <list>
#include "log.h"
#include <vector>
#include "LocalCoords.h"
#include "MeshSurface.h"

class MeshCell {
private:
	double _width;
	double _height;
	std::vector<int> _FSRs;
	MeshSurface* _mesh_surfaces;
	double* _bounds;
	int _fsr_start;
	int _fsr_end;
	double _chi[NUM_ENERGY_GROUPS];
	double _sigma_a[NUM_ENERGY_GROUPS];
	double _sigma_s[NUM_ENERGY_GROUPS*NUM_ENERGY_GROUPS];
	double _nu_sigma_f[NUM_ENERGY_GROUPS];
	double _diffusivity[NUM_ENERGY_GROUPS];
	double _old_flux[NUM_ENERGY_GROUPS];
	double _new_flux[NUM_ENERGY_GROUPS];
	double _volume;
	int _cell_id;

public:
	MeshCell();
	virtual ~MeshCell();
	void makeSurfaces(int numAzim);
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	std::vector<int>* getFSRs();
	void addFSR(int fsr);
	void setChi(double chi, int e);
	double* getChi();
	void setSigmaA(double sigmaA, int e);
	double* getSigmaA();
	void setSigmaS(double sigmaS, int e, int ep);
	double* getSigmaS();
	void setNuSigmaF(double nuSigmaF, int e);
	double* getNuSigmaF();
	void setDiffusivity(double diffusivity, int e);
	double* getDiffusivity();
	MeshSurface* findSurface(LocalCoords* coord, int i);
	void setBounds(double x, double y);
	double* getBounds();
	int getFSRStart();
	int getFSREnd();
	void setFSRStart(int fsr);
	void setFSREnd(int fsr);
	MeshSurface* getMeshSurfaces(int surface);
	void setNewFlux(double flux, int e);
	double* getNewFlux();
	void setOldFlux(double flux, int e);
	double* getOldFlux();
	void setVolume(double volume);
	double getVolume();
};


#endif /* MESHCELL_H_ */
