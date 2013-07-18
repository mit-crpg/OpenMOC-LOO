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
#include <string>
#include "configurations.h"
#include "log.h"


class Material {
private:
    static int _n; /* Counts the number of materials */
    int _uid;      /* monotonically increasing id based on n */
    int _id;
    double _sigma_t[NUM_ENERGY_GROUPS];
    double _sigma_a[NUM_ENERGY_GROUPS];
    double _sigma_f[NUM_ENERGY_GROUPS];
    double _nu_sigma_f[NUM_ENERGY_GROUPS];
    double _chi[NUM_ENERGY_GROUPS];
    /* first index is row number; second index is column number */
    double _sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]; 

    /* Indices for the start and end of nonzero elements */
    int _sigma_t_start, _sigma_t_end;
    int _sigma_a_start, _sigma_a_end;
    int _sigma_f_start, _sigma_f_end;
    int _nu_sigma_f_start, _nu_sigma_f_end;
    int _chi_start, _chi_end;
    int _sigma_s_start[NUM_ENERGY_GROUPS];
    int _sigma_s_end[NUM_ENERGY_GROUPS];
public:
    Material(int id,
             double *sigma_a, int sigma_a_cnt,
             double *sigma_t, int sigma_t_cnt,
             double *sigma_f, int sigma_f_cnt,
             double *nu_sigma_f, int nu_sigma_f_cnt,
             double *chi, int chi_cnt,
             double *sigma_s, int sigma_s_cnt);
    virtual ~Material();

    int getUid() const;
    int getId() const;
    double* getSigmaF();
    double* getNuSigmaF();
    double* getSigmaS();
    double* getSigmaT();
    double* getSigmaA();
    double* getChi();

    int getSigmaTStart();
    int getSigmaTEnd();
    int getSigmaAStart();
    int getSigmaAEnd();
    int getSigmaFStart();
    int getSigmaFEnd();
    int getNuSigmaFStart();
    int getNuSigmaFEnd();
    int getChiStart();
    int getChiEnd();
    int getSigmaSStart(int group);
    int getSigmaSEnd(int group);

    void setChi(double chi[NUM_ENERGY_GROUPS]);
    void setSigmaF(double sigma_f[NUM_ENERGY_GROUPS]);
    void setNuSigmaF(double nu_sigma_f[NUM_ENERGY_GROUPS]);
    void setSigmaS(double sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]);
    void setSigmaT(double sigma_t[NUM_ENERGY_GROUPS]);
    void setSigmaA(double sigma_a[NUM_ENERGY_GROUPS]);

    void checkSigmaT();
    std::string toString();
    void compressCrossSections();
};

#endif /* MATERIAL_H_ */
