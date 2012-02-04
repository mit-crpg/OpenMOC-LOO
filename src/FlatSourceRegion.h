/*
 * FlatSourceRegion.h
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef FLATSOURCEREGION_H_
#define FLATSOURCEREGION_H_

#include "configurations.h"


class FlatSourceRegion {
private:
	int _id;
	int _cell_id;
	int _material_id;
	double _volume;
	double _flux[NUM_ENERGY_GROUPS];
	double _old_flux[NUM_ENERGY_GROUPS];
	double _source[NUM_ENERGY_GROUPS];
	double _old_source[NUM_ENERGY_GROUPS];
public:
	FlatSourceRegion();
	virtual ~FlatSourceRegion();
    int getId() const;
    int getCellId() const;
    int getMaterialId() const;
    double getVolume() const;
    double* getFlux();
    double* getOldFlux();
    double* getOldSource();
    double* getSource();
    void setId(int id);
    void setCellId(int cellId);
    void setMaterialId(int material_id);
    void setVolume(double volume);
};


#endif /* FLATSOURCEREGION_H_ */
