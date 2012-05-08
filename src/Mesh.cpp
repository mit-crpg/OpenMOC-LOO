/*
 * Mesh.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "Mesh.h"

Mesh::Mesh(int cellWidth, int cellHeight){
	_cell_width  = cellWidth;
	_cell_height = cellHeight;
	_cells = new MeshCell[cellWidth * cellHeight];
}

Mesh::~Mesh(){
	delete [] _cells;
}

int Mesh::getCellWidth(){
	return _cell_width;
}

int Mesh::getCellHeight(){
	return _cell_height;
}

void Mesh::setCellWidth(int cellWidth){
	_cell_width = cellWidth;
}

void Mesh::setCellHeight(int cellHeight){
	_cell_height = cellHeight;
}

MeshCell* Mesh::getCells(){
	return _cells;
}

int Mesh::findMeshCell(double pointX, double pointY){

	double left;
	double top = _height / 2.0;
	bool flag = false;
	int cell = 0;
	//log_printf(NORMAL, "left: %f top: %f", - _width / 2.0, top);

	for (int y = 0; y < _cell_height; y++){
		left = - _width / 2.0;
		for (int x = 0; x < _cell_width; x++){
			cell = y * _cell_width + x;
			//log_printf(NORMAL, "cell width: %f, cell height: %f",getCells()[cell].getWidth(), getCells()[cell].getHeight());
			if (pointX >= left && pointX <= left + getCells()[cell].getWidth()
					&& pointY <= top && pointY >= top - getCells()[cell].getHeight()){
				flag = true;
				break;
			}
			left = left + getCells()[cell].getWidth();
		}
		if (flag == true){
			break;
		}
		top = top - getCells()[cell].getHeight();
	}

	return cell;
}

void Mesh::setWidth(double width){
	_width = width;
}

void Mesh::setHeight(double height){
	_height = height;
}

double Mesh::getWidth(){
	return _width;
}

double Mesh::getHeight(){
	return _height;
}



