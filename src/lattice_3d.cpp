#include <random>
#include "ranxoshi256.h"
#include "lattice_3d.h"
#include "linalg.h"
//#include <immintrin.h> 
#include <vector>
#include <iostream>

//#define ZORDERINDEXING

using namespace std;

typedef Lattice3d::LatticePoint LatticePoint;
typedef Lattice3d::Point Point;

Lattice3d::Lattice3d(int dimension): _dimension(dimension) {
    _actValues = new int[dimension * dimension * dimension]();
    _cellIds = new unsigned int[dimension * dimension * dimension]();
    _field = new double[3*dimension * dimension * dimension]();
}

Lattice3d::~Lattice3d() {
    delete[] _actValues;
    delete[] _cellIds;
    delete[] _field;
}

void Lattice3d::setPoint(int cellId, int x, int y, int z, int time, int type) {
    unsigned int id = cellId + (type << 24);
    _cellIds[index(x,y,z)] = id;
    _actValues[index(x,y,z)] = time;


    updateBorderTrackingAround(x, y, z);
}

unsigned int* Lattice3d::getCellIds() {
    return _cellIds;
}

void Lattice3d::updateBorderTrackingAround(int x, int y, int z) {
    if (isPartOfBorder(x,y,z)) {
        _borderIndices.add(index(x,y,z));
    } else {
        _borderIndices.remove(index(x,y,z));
    }
    for (int i = 0; i < getNeighborCount(); i++) {
        LatticePoint neighbor = getNeighbor(x,y,z, i);
        int nX = neighbor.x;
        int nY = neighbor.y;
        int nZ = neighbor.z;
        if (isPartOfBorder(nX,nY,nZ)) {
            _borderIndices.add(index(nX,nY,nZ));
        } else {
            _borderIndices.remove(index(nX,nY,nZ));
        }
    }
}

bool Lattice3d::isPartOfBorder(int x, int y, int z) {
    int cellId = _cellIds[index(x,y,z)];
    for (int i = 0; i < getNeighborCount(); i++) {
        auto neighbor = getNeighbor(x, y, z, i);
        if (neighbor.cellId != cellId)
            return true;
    }
    return false;
}


LatticePoint Lattice3d::getPoint(int x, int y, int z) {
    auto id = _cellIds[index(x,y,z)] & 16777215U;
    auto act = 0;
    if (_actToggle)
        act = _actValues[index(x,y,z)];
    char type = _cellIds[index(x,y,z)] >> 24;
    return {id, type, x, y, z, act};
}


Point Lattice3d::getFieldPoint(LatticePoint& point) {
    auto i = index(point.x, point.y, point.z);
    double x = _field[i];
    double y = _field[i + 1 * size()];
    double z = _field[i + 2 * size()];
    return {x, y, z};
}

LatticePoint Lattice3d::getPoint(int i) {
    auto id = _cellIds[i] & 16777215U;
    auto act = 0;
    if (_actToggle)
        act = _actValues[i];
    char type = _cellIds[i] >> 24;
#ifdef ZORDERINDEXING
    int x = _pext_u64(i, 0x9249249249249249);
    int y = _pext_u64(i, 0x2492492492492492);
    int z = _pext_u64(i, 0x4924924924924924);
#else
    int x = i % _dimension;
    int y = (i %(_dimension * _dimension)) / _dimension;
    int z = i / (_dimension*_dimension);
#endif
    return {id, type, x, y, z, act};
}

int Lattice3d::size() {
    return _dimension * _dimension * _dimension;
}



LatticePoint Lattice3d::getRandomBorderLocation(ranxoshi256& xoshi) {
    int r = ranxoshi256Next(&xoshi) % _borderIndices.size();
    int borderIndex = _borderIndices.get(r);

    
#ifdef ZORDERINDEXING
    int x = _pext_u64(borderIndex, 0x9249249249249249);
    int y = _pext_u64(borderIndex, 0x2492492492492492);
    int z = _pext_u64(borderIndex, 0x4924924924924924);
#else
    int x = borderIndex % _dimension;
    int y = (borderIndex %(_dimension * _dimension)) / _dimension;
    int z = borderIndex / (_dimension*_dimension);
#endif
    

    return getPoint(x,y,z);
}

LatticePoint Lattice3d::getRandomNeighbor(LatticePoint& point, 
        ranxoshi256& xoshi) {
    const int n = ranxoshi256Next(&xoshi) % getNeighborCount();
    return getNeighbor(point, n);
}

LatticePoint Lattice3d::getNeighbor(LatticePoint& point, int i) {
    return getNeighbor(point.x, point.y, point.z, i);
}

LatticePoint Lattice3d::getNeighbor(int x, int y, int z, int i) {
    const int directions[26][3] = {{-1,1,1}, {0,1,1}, {1,1,1}, {1,0,1}, 
        {1,-1,1}, {0,-1,1}, {-1,-1,1}, {-1,0,1}, {0,0,1}, {-1,1,0}, {0,1,0}, 
        {1,1,0}, {1,0,0}, {1,-1,0}, {0,-1,0}, {-1,-1,0}, {-1,0,0}, {-1,1,-1}, 
        {0,1,-1}, {1,1,-1}, {1,0,-1}, {1,-1,-1}, {0,-1,-1}, {-1,-1,-1}, 
        {-1,0,-1}, {0,0,-1}};
    auto offset = directions[i];
    int nX = (x + offset[0] ) & (_dimension-1);
    int nY = (y + offset[1] ) & (_dimension-1);
    int nZ = (z + offset[2] ) & (_dimension-1);
    return getPoint(nX,nY,nZ);
}


void Lattice3d::copy(LatticePoint& source, LatticePoint& target, int time) {
    setPoint(source.cellId, target.x, target.y, target.z, time, source.type);
}

int Lattice3d::index(unsigned int x, unsigned int y, unsigned int z) {
#ifdef ZORDERINDEXING
    return _pdep_u64(x, 0x9249249249249249) | 
        _pdep_u64(y, 0x2492492492492492) | 
        _pdep_u64(z, 0x4924924924924924);
#else
    return z * _dimension * _dimension + y * _dimension + x;
#endif
}

vector<vec3> Lattice3d::getPoints(int cellId) {
    vector<vec3> points;
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            for (int z = 0; z < _dimension; z++) {
                auto p = getPoint(x, y, z);
                if (p.cellId == cellId) {
                    points.push_back({float(x), float(y), float(z)});
                }
            }
        }
    }
    return points;
}


void Lattice3d::resetType(int cellId, int type) {
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            for (int z = 0; z < _dimension; z++) {
                auto p = getPoint(x, y, z);
                if (p.cellId == cellId) {
                    setPoint(cellId, x, y, z, p.act, type);
                }
            }
        }
    }
}

void Lattice3d::setPoints(int id, const vector<vec3>& points, int type) {
    for (auto point: points) {
        auto act = _actValues[index(point.x, point.y, point.z)];
        setPoint(id, point.x, point.y, point.z, act, type);
    }
    for (auto point: points) {
        updateBorderTrackingAround(point.x, point.y, point.z);
    }
}

void Lattice3d::remove(int id) {
    //FIX: border tracking
    for (int i = 0; i < _dimension * _dimension * _dimension; i++) {
        auto currentId = _cellIds[i] % (1<<24);
        if (currentId == id) {
            _cellIds[i] = 0;
        }
    }
}


Point Lattice3d::getCenterOfMass(int id) {
    double comX = 0;
    double comY = 0;
    double comZ = 0;
    int count = 0;
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            for (int z = 0; z < _dimension; z++) {
                auto p = getPoint(x, y, z);
                if (p.cellId == id) {
                    comX += x;
                    comY += y;
                    comZ += z;
                    count++;
                }
            }
        }
    }

    if (count > 0) {
        comX /= count;
        comY /= count;
        comZ /= count;
    }

    return {comX, comY, comZ};
}



int Lattice3d::getNeighborCount() {
    return 26;
}

void Lattice3d::setAct(bool actToggle) {
    _actToggle = actToggle;
}
