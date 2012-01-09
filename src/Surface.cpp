/*
 * Surface.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Surface.h"


/**
 * Constructor sets each coefficient
 * to 0 by default and the sense to true.
 */
Surface::Surface(int uid, surf_type type, int num_coeffs, float* coeffs){

	if (type == CIRCLE) {
		if (num_coeffs != 6) {
			std::cout << "Circle surface needs 3 coefficients but " << num_coeffs << " given. Exiting program." << std::endl;
			exit(1);
		}
		else {
			this->A = 1;
			this->B = 1;
			this->C = 0;
			this->D = -2*coeffs[0];
			this->E = -2*coeffs[1];
			this->F = (coeffs[0]*coeffs[0] + coeffs[1]*coeffs[1] - coeffs[2]*coeffs[2]);
		}

	}

	else if (type == PLANE) {
		if (num_coeffs != 6) {
			std::cout << "Plane type surface needs 6 coefficients but " << num_coeffs << " given. Exiting program." << std::endl;
			exit(1);
		}
		else {
			this->A = coeffs[0];
			this->B = coeffs[1];
			this->C = coeffs[2];
			this->D = coeffs[3];
			this->E = coeffs[4];
			this->F = coeffs[5];
		}
	}

	else if (type == XPLANE) {
		if (num_coeffs != 6) {
			std::cout << "XPlane type surface needs 1 coefficient but " << num_coeffs << " given. Exiting program." << std::endl;
			exit(1);
		}
		this->A = 0;
		this->B = 0;
		this->C = 0;
		this->D = 0;
		this->E = 1;
		this->F = -coeffs[0];

	}

	else if (type == YPLANE) {
		if (num_coeffs != 6) {
			std::cout << "YPlane type surface needs 1 coefficient but " << num_coeffs << " given. Exiting program." << std::endl;
			exit(1);
		}
		else {
			this->A = 0;
			this->B = 0;
			this->C = 0;
			this->D = 1;
			this->E = 0;
			this->F = -coeffs[0];
		}
	}

	else if (type == QUADRATIC) {
		if (num_coeffs != 6) {
			std::cout << "Quadratic type surface needs 6 coefficients but " << num_coeffs << " given. Exiting program." << std::endl;
			exit(1);
		}
		else {
			this->A = coeffs[0];
			this->B = coeffs[1];
			this->C = coeffs[2];
			this->D = coeffs[3];
			this->E = coeffs[4];
			this->F = coeffs[5];
		}
	}

	else {
		std::cout << "Surface type " << type << " not recognized. Exiting program" << std::endl;
		exit(1);
	}

	this->uid = uid;
	this->type = type;
}


Surface::~Surface() { }


float Surface::evaluate(Point* point) {
	float x = point->getX();
	float y = point->getY();
	return (this->A*x*x + this->B*y*y + this->C*x*y + this->D*x + this->E*y + this->F);
}
