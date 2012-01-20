/*
 * Material.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Material.h"


/**
 * Material constructor
 * @param id the material's id
 */
Material::Material(int id) {
	_id = id;
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
