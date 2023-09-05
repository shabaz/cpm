#include <iostream>
#include <set>
#include "cell_states.h"
#include "lattice.h"

using namespace std;


template <typename L>
void CellStates<L>::addCell(int area, int perimeter, int type) {
    _areas.push_back(area);
    _perimeters.push_back(perimeter);
    _types.push_back(type);
}

template <typename L>
int CellStates<L>::getArea(int cellId) {
    return _areas[cellId - 1];
}

template <typename L>
void CellStates<L>::removeArea(int cellId, int area) {
    _areas[cellId - 1] -= area;
}

template <typename L>
int CellStates<L>::getPerimeter(int cellId) {
    return _perimeters[cellId - 1];
}

template <typename L>
int CellStates<L>::getType(int cellId) {
    return _types[cellId - 1];
}
template <typename L>
void CellStates<L>::updateType(int cellId, int type) {
    _types[cellId - 1] = type;
}

template <typename L>
void CellStates<L>::updateAreas(LatticePoint& source, LatticePoint& target) {
    if (source.cellId != 0)
        _areas[source.cellId - 1] += 1;
    if (target.cellId != 0)
        _areas[target.cellId - 1] -= 1;
}

template <typename L>
void CellStates<L>::updatePerimeters(LatticePoint& source, LatticePoint& target, L &lattice) {
    int oldPerimeterSource = 0;
    int newPerimeterSource = 0;

    int oldPerimeterTarget = 0;
    int newPerimeterTarget = 0;

    for (int i = 0; i < lattice.getNeighborCount(); i++) {
        LatticePoint neighbor = lattice.getNeighbor(target, i);
        if (neighbor.cellId == source.cellId)
            oldPerimeterSource++;
        if (neighbor.cellId != source.cellId)
            newPerimeterSource++;
        if (neighbor.cellId != target.cellId)
            oldPerimeterTarget++;
        if (neighbor.cellId == target.cellId)
            newPerimeterTarget++;
    }


    if (source.cellId != 0) {
        int difference = -oldPerimeterSource + newPerimeterSource;
        _perimeters[source.cellId-1] += difference;
    }

    if (target.cellId != 0) {
        int difference = -oldPerimeterTarget + newPerimeterTarget;
        _perimeters[target.cellId-1] += difference;
    }

}

template <typename L>
void CellStates<L>::print() {
    for (int i = 0; i < _areas.size(); i++) {
        cout << " cell " << i+1 << " has an area of " << _areas[i] <<
            " and perimeter of " << _perimeters[i] << endl;
    }
}

template <typename L>
int CellStates<L>::nextId() {
    return _areas.size() + 1;
}

template <typename L>
void CellStates<L>::initializeFromGrid(L& lattice, int nrOfCells) {
    for (int i = 0; i < nrOfCells; i++) {
        _areas.push_back(0);
        _perimeters.push_back(0);
        //TODO: don't assume all cells in grid are cell type 1
        _types.push_back(1);
    }

    auto size = lattice.size();
    for (int i = 0; i < size; i++) {
        auto c = lattice.getPoint(i);
        if (c.cellId != 0) {
            _areas[c.cellId-1]++;
        } 

        std::vector<int> neighbors;
        for (int l = 0; l < lattice.getNeighborCount(); l++) {
            auto n = lattice.getNeighbor(c, l);
            if (n.cellId != c.cellId && n.cellId != 0) {
                neighbors.push_back(n.cellId);
            }
        }

        for (auto n: neighbors) {
            _perimeters[n-1]++;
        }
    }

}

template <typename L>
int CellStates<L>::countType(int type) {
    int count = 0;
    for (int i = 0; i < _areas.size(); i++) {
        if (_types[i] == type && _areas[i] > 0) {
            count++;
        }
    }
    return count;
}

template <typename L>
std::vector<int> CellStates<L>::getCellIds(int type) {
    std::vector<int> ids;
    for (int i = 0; i < _areas.size(); i++) {
        if (_types[i] == type && _areas[i] > 0) {
            ids.push_back(i+1);
        }
    }
    return ids;
}

template <typename L>
void CellStates<L>::recalcPerimeter(L& lattice, int id) {
    _perimeters[id-1] = 0;

    auto size = lattice.size();
    for (int i = 0; i < size; i++) {
        auto c = lattice.getPoint(i);
        for (int l = 0; l < lattice.getNeighborCount(); l++) {
            auto n = lattice.getNeighbor(c, l);
            if (n.cellId != c.cellId && c.cellId == id) {
                _perimeters[id-1]++;
            }
        } 
    }
}

template <typename L>
void CellStates<L>::kill(int id) {
    _areas[id-1] = 0;
    _perimeters[id-1] = 0;
}

template class CellStates<Lattice2d>;
template class CellStates<Lattice3d>;
