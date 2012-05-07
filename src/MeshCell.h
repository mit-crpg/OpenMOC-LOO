/*
 * MeshCell.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESHCELL_H_
#define MESHCELL_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <list>

class MeshCell {
private:
	double _width;
	double _height;
	std::list<int> _FSRs;

public:
	MeshCell();
	virtual ~MeshCell();
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	std::list<int>& getFSRs();
	void addFSR(int fsr);

};


#endif /* MESHCELL_H_ */
