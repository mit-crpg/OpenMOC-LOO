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

#include <sstream>
#include "configurations.h"
#include "log.h"


class Material {
private:
	static int _n;				/* Counts the number of materials */
	int _uid;					/* monotonically increasing id based on n */
	int _id;
	double _sigma_t[NUM_ENERGY_GROUPS];
	double _nu_sigma_f[NUM_ENERGY_GROUPS];
	double _sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS];
	double _chi[NUM_ENERGY_GROUPS];
public:
	Material(const int id);
	virtual ~Material();
	int getUid() const;
    int getId() const;
    double* getNuSigmaF();
    double* getSigmaS();
    double* getSigmaT();
    double* getChi();
    void setChi(double chi[NUM_ENERGY_GROUPS]);
    void setNuSigmaF(double nu_sigma_f[NUM_ENERGY_GROUPS]);
    void setSigmaS(double sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]);
    void setSigmaT(double sigma_t[NUM_ENERGY_GROUPS]);
    const char* toString();
};

#endif /* MATERIAL_H_ */
