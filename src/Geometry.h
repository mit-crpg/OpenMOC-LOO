/*
 * Geometry.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <map>
#include <list>
#include <utility>
#include <sstream>
#include <string>
#include <math.h>
#include <limits.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include "Parser.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Universe.h"
#include "Lattice.h"
#include "LocalCoords.h"
#include "Track.h"
#include "log.h"
#include "configurations.h"
#include "Point.h"
#include "silo.h"
#include "quickplot.h"
#include "Mesh.h"
#include "MeshCell.h"


class Geometry {
private:
	double _x_min, _y_min, _x_max, _y_max; 		/* the corners */
	int _base_universe;
	int _num_FSRs;
	int* _FSRs_to_cells;
	int* _FSRs_to_materials;
	double _max_seg_length;
	double _min_seg_length;
	std::map<int, Material*> _materials;
	std::map<int, Surface*> _surfaces;
	std::map<int, Cell*> _cells;
	std::map<int, Universe*> _universes;
	std::map<int, Lattice*> _lattices;

	std::vector<int> _surf_flags;
	std::vector<double> _surf_coeffs;
	std::vector<int> _oper_flags;
	std::vector<int> _left_ids;
	std::vector<int> _right_ids;
	std::vector<int> _zones;

public:
	Geometry(Parser* parser);
	virtual ~Geometry();
	double getWidth() const;
	double getHeight() const;
	int getNumRings() const;
	int getNumSectors() const;
	double getSectorOffset() const;
	int getNumFSRs() const;
	double getMaxSegmentLength() const;
	double getMinSegmentLength() const;
	int* getFSRtoCellMap() const;
	int* getFSRtoMaterialMap() const;

	void addMaterial(Material* material);
	Material* getMaterial(int id);
	void addSurface(Surface* surface);
	Surface* getSurface(int id);
	void addCell(Cell *cell);
	Cell* getCell(int id);
	void addUniverse(Universe* universe);
	Universe* getUniverse(int id);
	void addLattice(Lattice* lattice);
	Lattice* getLattice(int id);
	std::string toString();
	void printString();

	void adjustKeys();
	void buildNeighborsLists();
	Cell* findCell(LocalCoords* coords);
	Cell* findFirstCell(LocalCoords* coords, double angle);
	Cell* findCell(int fsr_id);
	Cell* findCell(Universe* univ, int fsr_id);
	Cell* findNextCell(LocalCoords* coords, double angle);
	int findFSRId(LocalCoords* coords);
	void segmentize(Track* track);

	void compressCrossSections();
	void computePinPowers(double* FSRs_to_powers, double* FSRs_to_pin_powers);
	double computePinPowers(Universe* univ, char* output_file_prefix,
			int FSR_id, double* FSRs_to_powers, double* FSRs_to_pin_powers);

	void generateCSG();

	template <class K, class V>
	bool mapContainsKey(std::map<K, V> map, K key);

	void makeCMFDMesh();
	void findNumLattices(Universe* univ,  int* numLattices);
	void findMeshWidth(Universe* univ, int* width, int depth);
	void findMeshHeight(Universe* univ, int* height, int depth);
	void defineMesh(Mesh* mesh, Universe* univ, int depth, int* meshCellNum, int row, bool base, int fsr_id);
	void findFSRs(Universe* univ, MeshCell meshCell, int* fsr_id);
	int nextLatticeHeight(Universe* curr);
	void plotCMFDMesh(Mesh* mesh);

};

#endif /* GEOMETRY_H_ */
