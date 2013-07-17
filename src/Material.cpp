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
				   double *sigma_f, int sigma_f_cnt,
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
		log_printf(ERROR, "Wrong number of sigma_f");
	memcpy(_sigma_f, sigma_f, NUM_ENERGY_GROUPS*sizeof(*_sigma_f));

	if (nu_sigma_f_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of nu_sigma_f");
	memcpy(_nu_sigma_f, nu_sigma_f, NUM_ENERGY_GROUPS*sizeof(*_nu_sigma_f));

	if (chi_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of chi");

	/* Normalize chi: each material, the sum of chi over all energy groups
	 * should add up to 1 */
	double chi_tot = 0;
	for (int i = 0; i < NUM_ENERGY_GROUPS; i++)
		chi_tot += chi[i];
	if (chi_tot > 0.0)
	{
		for (int i = 0; i < NUM_ENERGY_GROUPS; i++)
			_chi[i] = chi[i] / chi_tot;
	}
	else
		memcpy(_chi, chi, NUM_ENERGY_GROUPS*sizeof(*_chi));

	if (sigma_s_cnt != NUM_ENERGY_GROUPS*NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_s");
	
	/* Set the material's scattering matrix. This assumes that the scattering
	 * matrix passed in has the notation: the ij element is for scattering
	 * from group i to j. For efficient caching of the elements of this matrix
	 * during fixed source iteration, the matrix transpose is what is actually
	 * stored in the material
     */
	for (int i = 0; i < NUM_ENERGY_GROUPS; i++) 
	{
		for (int j = 0; j < NUM_ENERGY_GROUPS; j++) 
			_sigma_s[i][j] = sigma_s[j * NUM_ENERGY_GROUPS + i];
	}

	/* Uncompressed indices for the start and end of nonzero elements */
	_sigma_t_start = 0;
	_sigma_t_end = NUM_ENERGY_GROUPS;
	_sigma_a_start = 0;
	_sigma_a_end = NUM_ENERGY_GROUPS;
	_sigma_f_start = 0;
	_sigma_f_end = NUM_ENERGY_GROUPS;
	_nu_sigma_f_start = 0;
	_nu_sigma_f_end = NUM_ENERGY_GROUPS;
	_chi_start = 0;
	_chi_end = NUM_ENERGY_GROUPS;

	for (int e=0; e < NUM_ENERGY_GROUPS; e++) 
	{
		_sigma_s_start[e] = 0;
		_sigma_s_end[e] = NUM_ENERGY_GROUPS;
	}

	//FIXME
	compressCrossSections();
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
 * Return the material's fission cross-section array
 * @return the material's fission cross-section array
 */
double* Material::getSigmaF() {
    return _sigma_f;
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
 * Return the index of the first non-zero total cross-section in the array
 * of total cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param starting index for total cross-sections
 */
int Material::getSigmaTStart() {
	return _sigma_t_start;
}


/**
 * Return the index of the final non-zero total cross-section in the array
 * of total cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param ending index for total cross-sections
 */
int Material::getSigmaTEnd() {
	return _sigma_t_end;
}


/**
 * Return the index of the first non-zero absorption cross-section in the array
 * of absorption cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param starting index for absorption cross-sections
 */
int Material::getSigmaAStart() {
	return _sigma_a_start;
}


/**
 * Return the index of the final non-zero absorption cross-section in the array
 * of total cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param ending index for absorption cross-sections
 */
int Material::getSigmaAEnd() {
	return _sigma_a_end;
}


/**
 * Return the index of the first non-zero fission cross-section in the array
 * of fission cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param starting index for fission cross-sections
 */
int Material::getSigmaFStart() {
	return _sigma_f_start;
}


/**
 * Return the index of the final non-zero fission cross-section in the array
 * of total cross-sections if the cross-sections have been compressed.
 * Otherwise, returns 0.
 * @param ending index for fission cross-sections
 */
int Material::getSigmaFEnd() {
	return _sigma_f_end;
}


/**
 * Return the index of the first non-zero nu times fission cross-section in
 * the array of fission cross-sections if the cross-sections have been
 * compressed. Otherwise, returns 0.
 * @param starting index for num times fission cross-sections
 */
int Material::getNuSigmaFStart() {
	return _nu_sigma_f_start;
}



/**
 * Return the index of the final non-zero num times fission cross-section
 * in the array of total cross-sections if the cross-sections have been
 * compressed. Otherwise, returns 0.
 * @param ending index for nu times fission cross-sections
 */
int Material::getNuSigmaFEnd() {
	return _nu_sigma_f_end;
}


/**
 * Return the index of the first non-zero value of chi in the array
 * of chi if the cross-sections have been compressed. Otherwise, returns 0.
 * @param starting index for chi
 */
int Material::getChiStart() {
	return _chi_start;
}


/**
 * Return the index of the final non-zero value of chi in the array
 * of chi if the cross-sections have been compressed. Otherwise, returns 0.
 * @param ending index for chi
 */
int Material::getChiEnd() {
	return _chi_end;
}


/**
 * Return the index of the first non-zero scattering cross-section in the
 * scattering matrix for an energy group (row) if the cross-sections
 * have been compressed. Otherwise, returns 0.
 * @param starting index for a row of scattering cross-sections
 */
int Material::getSigmaSStart(int group) {

	if (group < 0 || group >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Unable to return the starting index for sigma_s "
				"start for group %d since it is out of the energy group "
				"bounds", group);

	return _sigma_s_start[group];
}


/**
 * Return the index of the final non-zero scattering cross-section in the
 * scattering matrix for an energy group (row) if the cross-sections
 * have been compressed. Otherwise, returns 0.
 * @param ending index for a row of scattering cross-sections
 */
int Material::getSigmaSEnd(int group) {

	if (group < 0 || group >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Unable to return the ending index for sigma_s for "
				"group %d since it is out of the energy group bounds", group);

	return _sigma_s_end[group];
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
 * Set the material's fission cross-section array
 * @param nu_sigma_f the fission cross-section array
 */
void Material::setSigmaF(double sigma_f[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_f[i] = sigma_f[i];
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
 * Set the material's scattering matrix. This assumes that the scattering
 * matrix passed in has the standard notation: the ij element is for scattering
 * from group i to j. For efficient caching of the elements of this matrix
 * during fixed source iteration, the matrix transpose is what is actually
 * stored in the material
 * @param sigma_s the material's scattering matrix
 */
void Material::setSigmaS(double sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			for (int j=0; j < NUM_ENERGY_GROUPS; j++)
			_sigma_s[i][j] = sigma_s[j][i];
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
 * Checks if the total cross-section for this material is equal to the
 * absorption plus scattering cross-sections for all energy groups
 */
void Material::checkSigmaT() 
{
	double calc_sigma_t;
	for (int i = 0; i < NUM_ENERGY_GROUPS; i++) 
	{
		/* Initialize the calculated total xs to the absorption xs */
		calc_sigma_t = _sigma_a[i];

		/* Increment calculated total xs by scatter xs for each energy group */
		for (int j=0; j < NUM_ENERGY_GROUPS; j++)
			calc_sigma_t += _sigma_s[j][i];

		log_printf(DEBUG, " material %d energy %d calculated sigma_t = %.10f,"
				   " given sigma_t = %.10f", _id, i, calc_sigma_t, _sigma_t[i]);
		_sigma_t[i] = calc_sigma_t;
	}

	return;
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




void Material::compressCrossSections() {

	log_printf(INFO, "Compressing cross-sections for material id = %d", _id);

	bool total_set = false;
	bool absorb_set = false;
	bool nu_fission_set = false;
	bool fission_set = false;
	bool chi_set = false;
	bool scatter_set = false;

	/* Initialize completely compressed indices */
	_sigma_t_start = 0;
	_sigma_t_end = 0;
	_sigma_a_start = 0;
	_sigma_a_end = 0;
	_sigma_f_start = 0;
	_sigma_f_end = 0;
	_nu_sigma_f_start = 0;
	_nu_sigma_f_end = 0;
	_chi_start = 0;
	_chi_end = 0;


	/* Find the indices of the first non-zero cross-section values */
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
		if (!total_set && _sigma_t[i] != 0) {
			_sigma_t_start = i;
			total_set = true;
		}
		if (!absorb_set && _sigma_a[i] != 0) {
			_sigma_a_start = i;
			absorb_set = true;
		}
		if (!nu_fission_set && _nu_sigma_f[i] != 0) {
			_nu_sigma_f_start = i;
			nu_fission_set = true;
		}
		if (!fission_set && _sigma_f[i] != 0) {
			_sigma_f_start = i;
			fission_set = true;
		}
		if (!chi_set && _chi[i] != 0) {
			_chi_start = i;
			chi_set = true;
		}
	}

//	log_printf(NORMAL, "sigma_t_start = %d, sigma_a_start = %d, "
//			"nu_sigma_f_start = %d, sigma_f_start = %d, chi_start = %d",
//			_sigma_t_start, _sigma_a_start, _nu_sigma_f_start, _sigma_f_start,
//			_chi_start);

	total_set = false;
	absorb_set = false;
	nu_fission_set = false;
	fission_set = false;
	chi_set = false;

	/* Find the indices of the final non-zero cross-section values */
	for (int i=NUM_ENERGY_GROUPS-1; i > -1; i--) {
		if (!total_set && _sigma_t[i] != 0) {
			_sigma_t_end = i+1;
			total_set = true;
		}
		if (!absorb_set && _sigma_a[i] != 0) {
			_sigma_a_end = i+1;
			absorb_set = true;
		}
		if (!nu_fission_set && _nu_sigma_f[i] != 0) {
			_nu_sigma_f_end = i+1;
			nu_fission_set = true;
		}
		if (!fission_set && _sigma_f[i] != 0) {
			_sigma_f_end = i+1;
			fission_set = true;
		}
		if (!chi_set && _chi[i] != 0) {
			_chi_end = i+1;
			chi_set = true;
		}
	}

//	log_printf(NORMAL, "sigma_t_end = %d, sigma_a_end = %d, "
//			"nu_sigma_f_end = %d, sigma_f_end = %d, chi_end = %d",
//			_sigma_t_end, _sigma_a_end, _nu_sigma_f_end, _sigma_f_end,
//			_chi_end);

	/* Scattering cross-section matrix starting indices */
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {

		_sigma_s_start[i] = 0;
		_sigma_s_end[i] = 0;
		scatter_set = false;

		for (int j=0; j < NUM_ENERGY_GROUPS; j++) {
			if (!scatter_set && _sigma_s[i][j] != 0) {
				_sigma_s_start[i] = j;
				scatter_set = true;
				break;
			}
		}

		scatter_set = false;

		for (int j=NUM_ENERGY_GROUPS-1; j > -1; j--) {
			if (!scatter_set && _sigma_s[i][j] != 0) {
				_sigma_s_end[i] = j+1;
				scatter_set = true;
				break;
			}
		}

//		log_printf(NORMAL, "i = %d, sigma_s_start = %d, sigma_s_end = %d",
//								i, _sigma_s_start[i], _sigma_s_end[i]);
	}

	return;
}

