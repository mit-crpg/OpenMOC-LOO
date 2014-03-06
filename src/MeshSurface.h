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
    double *_total_wt; 
    double _as_tracked_length;
    double *_current;
    double **_quad_current;
    double **_quad_flux;
    double **_old_quad_flux;
    double *_d_hat;
    double *_d_tilde;

    /* Surface ID could be 0,1,2,3 and is later set in MeshCell.cpp */
    int _id; 
    int _cell_id;
    boundaryType _boundary_type;

public:
    MeshSurface();
    virtual ~MeshSurface();

    /* LOO Only */
    double getQuadCurrent(int group, int index);
    void setQuadCurrent(double quad_current, int group, int index);
    void incrementQuadCurrent(double quad_current, int group, int index);
    void incrementTotalWt(double quad_current, int index);
    void setTotalWt(double wt, int index);
    __inline__ double getTotalWt(int index){
        return _total_wt[index];
    }
    void updateQuadCurrent(double factor, int group, int index);
    double getQuadFlux(int group, int index);
    void setQuadFlux(double quad_current, int group, int index);
    void updateQuadFlux(double ratio, int group, int index);
    double getOldQuadFlux(int group, int index);
    void setOldQuadFlux(double quad_current, int group, int index);

    /* CMFD Only */
    void makeCurrents();
    void setCurrent(double current, int group);
    double getCurrent(int group);
    void updateCurrent(double factor, int group);
    void incrementCurrents(double *current);
    void incrementCurrent(double current, int group);
    double getAsTrackedLength();
    void incrementAsTrackedLength(double length);
    void setDHat(double dHat, int e);
    double* getDHat();
    void setDTilde(double dTilde, int e);
    double* getDTilde();

    /* General Purpose */
    int getId();
    void setId(int id);
    int getCellId();
    void setCellId(int id);
    void setBoundary(boundaryType boundary);
    boundaryType getBoundary();
};


#endif /* MESHSURFACE_H_ */
