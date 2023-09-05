#include <random>
#include "ranxoshi256.h"
#include "lattice_2d.h"
#include "linalg.h"
//#include <immintrin.h>
#include <vector>
#include <iostream>

//#define ZORDERINDEXING

using namespace std;

typedef Lattice2d::LatticePoint LatticePoint;
typedef Lattice2d::Point Point;

Lattice2d::Lattice2d(int dimension): _dimension(dimension) {
    _actValues = new int[dimension * dimension]();
    _cellIds = new unsigned int[dimension * dimension]();
    _field = new double[2*dimension * dimension]();
#ifdef ZORDERINDEXING
    _cellIdsReshaped = new unsigned int[dimension * dimension]();
#endif
}

Lattice2d::~Lattice2d() {
    delete[] _actValues;
    delete[] _cellIds;
    delete[] _field;
}
unsigned int* Lattice2d::getCellIds() {
#ifdef ZORDERINDEXING
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            _cellIdsReshaped[y*_dimension+x] = _cellIds[index(x,y)];
        }
    }
    return _cellIdsReshaped;
#else
    return _cellIds;
#endif
}

void Lattice2d::setPoint(int cellId, int x, int y, int time, int type) {
    unsigned int id = cellId + (type << 24);
    _cellIds[index(x,y)] = id;
    _actValues[index(x,y)] = time;
    updateBorderTrackingAround(x, y);
}

void Lattice2d::updateBorderTrackingAround(int x, int y) {
    if (isPartOfBorder(x,y)) {
        _borderIndices.add(index(x,y));
    } else {
        _borderIndices.remove(index(x,y));
    }
    for (int i = 0; i < 8; i++) {
        LatticePoint neighbor = getNeighbor(x,y, i);
        int nX = neighbor.x;
        int nY = neighbor.y;
        if (isPartOfBorder(nX,nY)) {
            _borderIndices.add(index(nX,nY));
        } else {
            _borderIndices.remove(index(nX,nY));
        }
    }
}

bool Lattice2d::isPartOfBorder(int x, int y) {
    int cellId = _cellIds[index(x,y)];
    for (int i = 0; i < 8; i++) {
        auto neighbor = getNeighbor(x, y, i);
        if (neighbor.cellId != cellId)
            return true;
    }
    return false;
}

LatticePoint Lattice2d::getPoint(int i) {
    auto id = _cellIds[i] & 16777215U;
    auto act = 0;
    if (_actToggle)
        act = _actValues[i];
    char type = _cellIds[i] >> 24;
#ifdef ZORDERINDEXING
    int x = _pext_u64(i, 0x5555555555555555);
    int y = _pext_u64(i, 0xaaaaaaaaaaaaaaaa);
#else
    int x = i % _dimension;
    int y = i / _dimension;
#endif
    return {id, type, x, y, act};
}


Point Lattice2d::getFieldPoint(LatticePoint& point) {
    auto i = index(point.x, point.y);
    double x = _field[i];
    double y = _field[i + 1 * size()];
    return {x, y};
}

LatticePoint Lattice2d::getPoint(int x, int y) {
    auto id = _cellIds[index(x,y)] & 16777215U;
    auto act = 0;
    if (_actToggle)
        act = _actValues[index(x,y)];
    char type = _cellIds[index(x,y)] >> 24;
    return {id, type, x, y, act};
}

int Lattice2d::size() {
    return _dimension * _dimension;
}



LatticePoint Lattice2d::getRandomBorderLocation(ranxoshi256& xoshi) {
    int r = ranxoshi256Next(&xoshi) % _borderIndices.size();
    int borderIndex = _borderIndices.get(r);

#ifdef ZORDERINDEXING
    int x = _pext_u64(borderIndex, 0x5555555555555555);
    int y = _pext_u64(borderIndex, 0xaaaaaaaaaaaaaaaa);
#else
    int x = borderIndex % _dimension;
    int y = borderIndex / _dimension;
#endif
    

    return getPoint(x,y);
}

LatticePoint Lattice2d::getRandomNeighbor(LatticePoint& point, 
        ranxoshi256& xoshi) {
    const int n = ranxoshi256Next(&xoshi) % 8;
    return getNeighbor(point, n);
}

LatticePoint Lattice2d::getNeighbor(LatticePoint& point, int i) {
    return getNeighbor(point.x, point.y, i);
}

LatticePoint Lattice2d::getNeighbor(int x, int y, int i) {
    const int directions[8][2] = {{-1,1}, {0,1}, {1,1}, {1,0}, {1,-1},
        {0,-1}, {-1,-1}, {-1,0}};
    auto offset = directions[i];
    int nX = (x + offset[0] ) & (_dimension-1);
    int nY = (y + offset[1] ) & (_dimension-1);
    return getPoint(nX,nY);
}


void Lattice2d::copy(LatticePoint& source, LatticePoint& target, int time) {
    setPoint(source.cellId, target.x, target.y, time, source.type);
}

int Lattice2d::index(unsigned int x, unsigned int y) {
#ifdef ZORDERINDEXING
    return _pdep_u64(x, 0x5555555555555555) | _pdep_u64(y,0xaaaaaaaaaaaaaaaa);
#else
    return y * _dimension + x;
#endif
}

vector<vec2> Lattice2d::getPoints(int cellId) {
    vector<vec2> points;
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            auto p = getPoint(x, y);
            if (p.cellId == cellId) {
                points.push_back({float(x), float(y)});
            }
        }
    }
    return points;
}


void Lattice2d::resetType(int cellId, int type) {
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            auto p = getPoint(x, y);
            if (p.cellId == cellId) {
                setPoint(cellId, x, y, p.act, type);
            }
        }
    }
}

void Lattice2d::setPoints(int id, const vector<vec2>& points, int type) {
    for (auto point: points) {
        auto act = _actValues[index(point.x, point.y)];
        setPoint(id, point.x, point.y, act, type);
    }
    for (auto point: points) {
        updateBorderTrackingAround(point.x, point.y);
    }
}

void Lattice2d::remove(int id) {
    //FIX: border tracking
    for (int i = 0; i < _dimension * _dimension; i++) {
        auto currentId = _cellIds[i] % (1<<24);
        if (currentId == id) {
            _cellIds[i] = 0;
        }
    }
}


Point Lattice2d::getCenterOfMass(int id) {
    double comX = 0;
    double comY = 0;
    int count = 0;
    for (int x = 0; x < _dimension; x++) {
        for (int y = 0; y < _dimension; y++) {
            auto p = getPoint(x, y);
            if (p.cellId == id) {
                comX += x;
                comY += y;
                count++;
            }
        }
    }

    if (count > 0) {
        comX /= count;
        comY /= count;
    }

    return {comX, comY};
}



int Lattice2d::getNeighborCount() {
    return 8;
}

void Lattice2d::setAct(bool actToggle) {
    _actToggle = actToggle;
}
