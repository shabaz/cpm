#include <iostream>
#include <cstdlib>

#include "centroids.h"
#include "lattice_2d.h"
#include "lattice_3d.h"
#include "cell_states.h"

using namespace std;

template <typename L>
void Centroids<L>::initializeFromGrid(L& lattice, int nrOfCells) {
    for (int i = 0; i < nrOfCells; i++) {
        _centers.push_back(IntPoint());
        _counts.push_back(0);
        auto p = Point();
        p.unitRandomize();
        _preferredDirections.push_back(p);
        _history.push_back(list<Point>());
    }

    auto size = lattice.size();
    for (int i = 0; i < size; i++) {
        auto c = lattice.getPoint(i);
        if (c.cellId != 0) {
            _counts[c.cellId-1]++;
            _centers[c.cellId-1] = _centers[c.cellId-1].add(c.times(1));
        } 

    }
}

template <typename L>
Centroids<L>::Centroids(int dimension, int numberOfTypes, CellStates<L>& cellStates): 
    _dimension(dimension), _cellStates(cellStates) {
        _persistenceValues = new double[numberOfTypes]();
        _historyLengths = new int[numberOfTypes];
        for (int i = 0; i < numberOfTypes; i++) {
            _historyLengths[i] = 1;
        }
}

template <typename L>
Centroids<L>::~Centroids() {
    delete[] _persistenceValues;
    delete[] _historyLengths;
}

template <typename L>
void Centroids<L>::setHistoryLength(int type, int historyLength) {
    _historyLengths[type] = historyLength;
}

template <typename L>
void Centroids<L>::setPersistence(int type, double persistence) {
    _persistenceValues[type] = persistence;
}

template <typename L>
typename L::Point Centroids<L>::getPrefDir(int cellId) {
    return _preferredDirections[cellId-1];
}

template <typename L>
void Centroids<L>::update(LatticePoint& source, LatticePoint& target) {
    if (source.cellId) {

        auto offset = _centers[source.cellId-1].subtract( 
            target.times(_counts[source.cellId-1]));
        if (_counts[source.cellId-1] > 0)
            offset = offset.intDiv(_dimension/2 * _counts[source.cellId-1]).mul(_dimension);
        _counts[source.cellId-1]++;
        _centers[source.cellId-1] = _centers[source.cellId-1].add(target);
        _centers[source.cellId-1] = _centers[source.cellId-1].add(offset);
        _centers[source.cellId-1] = _centers[source.cellId-1].modulo(
                _counts[source.cellId-1] * _dimension);


    }

    if (target.cellId) {
        auto offset = _centers[target.cellId-1].subtract( 
            target.times(_counts[target.cellId-1]));
        if (_counts[target.cellId-1] > 0)
            offset = offset.intDiv(_dimension/2 * _counts[target.cellId-1]).mul(_dimension);
        _counts[target.cellId-1]--;
        _centers[target.cellId-1] = _centers[target.cellId-1].subtract(target);
        _centers[target.cellId-1] = _centers[target.cellId-1].subtract(offset);
        if (_counts[target.cellId-1] > 0)
            _centers[target.cellId-1] = _centers[target.cellId-1].modulo(
                    _counts[target.cellId-1] * _dimension);
    }

}

template <typename L>
void Centroids<L>::addCentroid(IntPoint center, int count) {
    _centers.push_back(center);
    _counts.push_back(count);
    auto p = Point();
    p.unitRandomize();
    _preferredDirections.push_back(p);
    _history.push_back(list<Point>());
}

template <typename L>
void Centroids<L>::print() {
    for (int i = 0; i < _counts.size(); i++) {
        cout << _counts[i] << endl;
    }

}



template <typename L>
std::vector<typename L::Point> Centroids<L>::getCentroids() {
    std::vector<Point> points;
    for (int i = 0; i < _centers.size(); i++) {
        points.push_back(_centers[i].divide(_counts[i]));
    }
    return points;
}


template <typename L>
void Centroids<L>::addCheckpoint() {
    auto current = getCentroids();
    for (int i = 0; i < current.size(); i++) {
        auto type = _cellStates.getType(i+1);
        while (_history[i].size() >= _historyLengths[type]) {
            _history[i].pop_front();
        }
        _history[i].push_back(current[i]);
    }
}

template <typename L>
void Centroids<L>::updatePreferentialDirection() {
    auto currentCentroids = getCentroids();
    for (int i = 0; i < _preferredDirections.size(); i++) {
        auto& currentPrefDir = _preferredDirections[i];
        auto currentDir = currentCentroids[i].subtract(_history[i].front());

        if (currentDir.length() == 0) {
            continue;
        }

        currentDir = currentDir.wrap(_dimension);

        currentDir = currentDir.normalize();

        const auto type = _cellStates.getType(i+1);
        const auto persistence = _persistenceValues[type];

        currentPrefDir = currentDir.times(1-persistence).add(currentPrefDir.times(persistence)).normalize();

    }
    
}

template class Centroids<Lattice2d>;
template class Centroids<Lattice3d>;
