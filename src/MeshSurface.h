/*
 * MeshSurface.h
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#ifndef MESHSURFACE_H_
#define MESHSURFACE_H_

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
#include "configurations.h"
#include "Surface.h"



class MeshSurface {
private:
    double *_current;
    double **_quad_current;
    double **_quad_flux;
    double **_old_quad_flux;
    double *_d_hat;
    double *_d_tilde;
    double *_d_dif;

    /* Surface ID could be 0,1,2,3 and is later set in MeshCell.cpp */
    int _id; 
    int _cell_id;
    boundaryType _boundary_type;

public:
    MeshSurface();
    virtual ~MeshSurface();

    /* LOO Only */
    void setQuadCurrent(double quad_current, int group, int index);
    void incrementQuadCurrent(double quad_current, int group, int index);
    double getQuadCurrent(int group, int index);
    void setQuadFlux(double quad_current, int group, int index);
    double getQuadFlux(int group, int index);
    void setOldQuadFlux(double quad_current, int group, int index);
    double getOldQuadFlux(int group, int index);

    /* CMFD Only */
    void makeCurrents();
    void setCurrent(double current, int group);
    double getCurrent(int group);
    void incrementCurrent(double *current);
    void setDHat(double dHat, int e);
    double* getDHat();
    void setDTilde(double dTilde, int e);
    double* getDTilde();
    void setDDif(double dTilde, int e);
    double* getDDif();

    /* General Purpose */
    int getId();
    void setId(int id);
    int getCellId();
    void setCellId(int id);
    void setBoundary(boundaryType boundary);
    boundaryType getBoundary();
};


#endif /* MESHSURFACE_H_ */
