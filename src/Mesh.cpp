/*
 * Mesh.cpp
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#include "Mesh.h"

Mesh::Mesh(){}

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

MeshCell* Mesh::getCells(int cell){
	return &_cells[cell];
}

void Mesh::makeMeshCells(){
	_cells = new MeshCell[_cell_width * _cell_height];

	/* give each surface an ID to its MeshCell */
	for (int y = 0; y < _cell_height; y++){
		for (int x = 0; x < _cell_width; x++){
			for (int s = 0; s < 8; s++){
				_cells[y*_cell_width + x].getMeshSurfaces(s)->setMeshCell(y*_cell_width + x);
				_cells[y*_cell_width + x].getMeshSurfaces(s)->setSurfaceNum(s);
			}
		}
	}
}

int Mesh::findMeshCell(double pointX, double pointY){

	double left;
	double top = _height / 2.0;
	bool flag = false;
	int cell = 0;

	for (int y = 0; y < _cell_height; y++){
		left = - _width / 2.0;
		cell = y * _cell_width;
		if (pointY <= top && pointY >= top - getCells(cell)->getHeight()){
			for (int x = 0; x < _cell_width; x++){
				cell = y * _cell_width + x;
				if (pointX >= left && pointX <= left + getCells(cell)->getWidth()){
					flag = true;
					break;
				}
				left = left + getCells(cell)->getWidth();
			}
			if (flag == true){
				break;
			}
		}
		top = top - getCells(cell)->getHeight();
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

void Mesh::setCellBounds(){

	double x = -_width / 2.0;
	double y = _height / 2.0;

	/* loop over MeshCells and set bounds */
	for (int i = 0; i < _cell_height; i++){
		x = -_width / 2.0;
		y = y - _cells[i * _cell_width].getHeight();
		for (int j = 0; j < _cell_width; j++){
			_cells[i * _cell_width + j].setBounds(x,y);
			x = x + _cells[i * _cell_width + j].getWidth();
		}
	}

}

void Mesh::setFSRBounds(){

	int min;
	int max;
	int fsr;

	for (int i = 0; i < _cell_height * _cell_width; i++){
		min = _cells[i].getFSRs()->front();
		max = _cells[i].getFSRs()->front();

		std::vector<int>::iterator iter;
		for (iter = _cells[i].getFSRs()->begin(); iter != _cells[i].getFSRs()->end(); ++iter) {
			fsr = *iter;
			min = std::min(fsr, min);
			max = std::max(fsr, max);
		}

		_cells[i].setFSRStart(min);
		_cells[i].setFSREnd(max);
	}
}



MeshSurface* Mesh::findMeshSurface(int fsr_id, LocalCoords* coord){
	MeshSurface* meshSurface = NULL;

	/* find which MeshCell fsr_id is in -> get meshSuface that coord is on*/
	for (int i = 0; i < _cell_width * _cell_height; i++){
		if (fsr_id >= _cells[i].getFSRStart() && fsr_id <= _cells[i].getFSREnd()){
			meshSurface = _cells[i].findSurface(coord);
			break;
		}

	}

	return meshSurface;
}

void Mesh::printBounds(){

	double* bounds;

	for (int i = 0; i < _cell_width * _cell_height; i++){
		bounds = _cells[i].getBounds();
		log_printf(NORMAL, "cell: %i bounds [%f, %f, %f, %f]", i, bounds[0], bounds[1], bounds[2], bounds[3]);
	}

}


void Mesh::printCurrents(){

	double current;

	for (int i = 0; i < _cell_width * _cell_height; i++){
		for (int surface = 0; surface < 8; surface++){
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
				current = _cells[i].getMeshSurfaces(surface)->getCurrent(group);
				log_printf(NORMAL, "cell: %i, surface: %i, group: %i, current: %f", i, surface, group, current);
			}
		}
	}
}


void Mesh::splitCorners(){

	MeshSurface* surfaceSide;
	MeshSurface* surfaceCorner1;
	MeshSurface* surfaceCorner2;

	for (int i = 0; i < _cell_width * _cell_height; i++){

/* side 0 */
		surfaceSide = _cells[i].getMeshSurfaces(0);
		surfaceCorner1 = _cells[i].getMeshSurfaces(4);
		surfaceCorner2 = _cells[i].getMeshSurfaces(7);

		for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
			surfaceSide->incrementCurrent(0.5 * surfaceCorner1->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementCurrent(0.5 * surfaceCorner2->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner1->getFlux(group), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner2->getFlux(group), group);
		}

/* side 1 */
		surfaceSide = _cells[i].getMeshSurfaces(1);
		surfaceCorner2 = _cells[i].getMeshSurfaces(5);

		for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
			surfaceSide->incrementCurrent(0.5 * surfaceCorner1->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementCurrent(0.5 * surfaceCorner2->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner1->getFlux(group), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner2->getFlux(group), group);
		}

/* side 2 */
		surfaceSide = _cells[i].getMeshSurfaces(2);
		surfaceCorner1 = _cells[i].getMeshSurfaces(6);

		for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
			surfaceSide->incrementCurrent(0.5 * surfaceCorner1->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementCurrent(0.5 * surfaceCorner2->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner1->getFlux(group), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner2->getFlux(group), group);
		}

/* side 3 */
		surfaceSide = _cells[i].getMeshSurfaces(3);
		surfaceCorner2 = _cells[i].getMeshSurfaces(7);

		for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
			surfaceSide->incrementCurrent(0.5 * surfaceCorner1->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementCurrent(0.5 * surfaceCorner2->getCurrent(group), surfaceSide->getNormal(), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner1->getFlux(group), group);
			surfaceSide->incrementFlux(0.5 * surfaceCorner2->getFlux(group), group);
		}

/* zero corner currents and surface fluxes */
		for (int side = 4; side < 8; side++){
			surfaceCorner1 = _cells[i].getMeshSurfaces(side);
			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
				surfaceCorner1->setCurrent(0, group);
				surfaceCorner1->setFlux(0, group);
			}
		}
	}
}

}


void Mesh::normalizeSurfaces(){

	MeshSurface* meshSurface;
	MeshCell* meshCell;

	for (int i = 0; i < _cell_width * _cell_height; i++){
		meshCell = getCells(i);
		log_printf(NORMAL, "Normalizing surface %i", i);
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
			meshSurface = meshCell->getMeshSurfaces(0);
			meshSurface->setCurrent(meshSurface->getCurrent(e) / meshSurface->getCrossings(),e);
			meshSurface->setFlux(meshSurface->getFlux(e) / meshSurface->getCrossings(),e);
			meshSurface = meshCell->getMeshSurfaces(1);
			meshSurface->setCurrent(meshSurface->getCurrent(e) / meshSurface->getCrossings(),e);
			meshSurface->setFlux(meshSurface->getFlux(e) / meshSurface->getCrossings(),e);
			meshSurface = meshCell->getMeshSurfaces(2);
			meshSurface->setCurrent(meshSurface->getCurrent(e) / meshSurface->getCrossings(),e);
			meshSurface->setFlux(meshSurface->getFlux(e) / meshSurface->getCrossings(),e);
			meshSurface = meshCell->getMeshSurfaces(3);
			meshSurface->setCurrent(meshSurface->getCurrent(e) / meshSurface->getCrossings(),e);
			meshSurface->setFlux(meshSurface->getFlux(e) / meshSurface->getCrossings(),e);
		}
	}

}









