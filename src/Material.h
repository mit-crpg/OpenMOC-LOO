/*
 * Material.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "configurations.h"

class Material {
private:
	int _id;
	double _sigma_t[NUM_ENERGY_GROUPS];
	double _nu_sigma_f[NUM_ENERGY_GROUPS];
	double _sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS];
	double _chi[NUM_ENERGY_GROUPS];
public:
	Material(int id);
	virtual ~Material();
    double* getChi();
    int getId() const;
    double* getNuSigmaF();
    double* getSigmaS();
    double* getSigmaT();
    void setChi(double chi[NUM_ENERGY_GROUPS]);
    void setId(int id);
    void setNuSigmaF(double nu_sigma_f[NUM_ENERGY_GROUPS]);
    void setSigmaS(double sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]);
    void setSigmaT(double sigma_t[NUM_ENERGY_GROUPS]);
};

#endif /* MATERIAL_H_ */
