/*
 * Material.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Material.h"

/* _n keeps track of the number of materials instantiated */
int Material::_n = 0;

/**
 * Material constructor
 * @param id the material's id
 */
Material::Material(int id,
				   double *sigma_a, int sigma_a_cnt,
				   double *sigma_t, int sigma_t_cnt,
				   double *nu_sigma_f, int nu_sigma_f_cnt,
				   double *chi, int chi_cnt,
				   double *sigma_s, int sigma_s_cnt) {
	_uid = _n;
	_id = id;
	_n++;

	if (sigma_a_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_a");
	memcpy(_sigma_a, sigma_a, NUM_ENERGY_GROUPS*sizeof(*_sigma_a));

	if (sigma_t_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_t");
	memcpy(_sigma_t, sigma_t, NUM_ENERGY_GROUPS*sizeof(*_sigma_t));

	if (nu_sigma_f_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of nu_sigma_f");
	memcpy(_nu_sigma_f, nu_sigma_f, NUM_ENERGY_GROUPS*sizeof(*_nu_sigma_f));

	if (chi_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of chi");
	memcpy(_chi, chi, NUM_ENERGY_GROUPS*sizeof(*_chi));

	if (sigma_s_cnt != NUM_ENERGY_GROUPS*NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_s");
	
	for (int i=0; i<NUM_ENERGY_GROUPS; i++) {
		for (int j=0; j<NUM_ENERGY_GROUPS; j++) {
			_sigma_s[i][j] = sigma_s[i*NUM_ENERGY_GROUPS+j];
		}
	}
}

/**
 * Destructor
 */
Material::~Material() { }



/**
 * Get the material's chi array
 * @return the material's chi array
 */
double* Material::getChi() {
    return _chi;
}



/**
 * Return the material's uid
 * @return the material's uid
 */
int Material::getUid() const {
	return _uid;
}

/**
 * Return the material's id
 * @return the material's id
 */
int Material::getId() const {
    return _id;
}


/**
 * Return the material's nu*sigma_f array
 * @return the material's nu*sigma_f array
 */
double* Material::getNuSigmaF() {
    return _nu_sigma_f;
}


/**
 * Return the material's scattering matrix
 * @return the material's scattering matrix
 */
double* Material::getSigmaS() {
    return *_sigma_s;
}


/**
 * Return the material's total cross-section array
 * @return the material's total cross-section array
 */
double* Material::getSigmaT() {
    return _sigma_t;
}


/**
 * Return the material's absorption cross section array
 * @return the material's absorption cross-section array
 */
double* Material::getSigmaA() {
	return _sigma_a;
}


/**
 * Set the material's chi array
 * @param chi the chi array
 */
void Material::setChi(double chi[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_chi[i] = chi[i];
}


/**
 * Set the material's nu*sigma_f array
 * @param nu_sigma_f the nu*sigma_f array
 */
void Material::setNuSigmaF(double nu_sigma_f[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_nu_sigma_f[i] = nu_sigma_f[i];
}


/**
 * Set the material's scattering matrix
 * @param sigma_s the material's scattering matrix
 */
void Material::setSigmaS(double sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
		for (int j=0; i < NUM_ENERGY_GROUPS; i++)
			_sigma_s[i][j] = sigma_s[i][j];
	}
}


/**
 * Set the material's total scattering cross-section array
 * @param sigma_t the material's total scattering cross-section
 */
void Material::setSigmaT(double sigma_t[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_t[i] = sigma_t[i];
}


/**
 * Set the material's absorption scattering cross-section array
 * @param sigma_a the material's absorption scattering cross-section
 */
void Material::setSigmaA(double sigma_a[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_a[i] = sigma_a[i];
}


/**
 * Converts this material's attributes to a character array representation
 * @param a character array of this member's attributes
 */
std::string Material::toString() {
	std::stringstream string;

	string << "Material id = " << _id;

	string << "\n\t\tSigma_a = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _sigma_a[e] << ", ";

	string << "\n\t\tSigma_t = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _sigma_t[e] << ", ";

	string << "\n\t\tnu_sigma_f = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _nu_sigma_f[e] << ", ";

	string << "\n\t\tSigma_s = \n\t\t";
	for (int G = 0; G < NUM_ENERGY_GROUPS; G++) {
		for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
			string << _sigma_s[G][g] << "\t\t ";
		string << "\n\t\t";
	}

	string << "Chi = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _chi[e] << ", ";

	return string.str();
}
