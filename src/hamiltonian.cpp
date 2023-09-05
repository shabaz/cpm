#include <cmath>
#include <iostream>
#include "hamiltonian.h"
#include "lattice.h"
#include "cell_states.h"
#include "centroids.h"

using namespace std;

template <typename L>
Hamiltonian<L>::Hamiltonian(int numberOfTypes, double temperature): 
    _numberOfTypes(numberOfTypes), _temperature(temperature) {
    _adhesionMatrix = new int[numberOfTypes * numberOfTypes]();
    _areaLambdas = new double[numberOfTypes]();
    _areaTargets = new double[numberOfTypes]();
    _perimeterLambdas = new double[numberOfTypes]();
    _perimeterTargets = new double[numberOfTypes]();
    _actLambdas = new double[numberOfTypes]();
    _actMaxima = new double[numberOfTypes]();
    _connectedLambdas = new double[numberOfTypes]();
    _chemotaxisLambdas = new double[numberOfTypes]();
    _fixedCelltype = new bool[numberOfTypes]();
    _persistenceLambdas = new double[numberOfTypes]();
}

template <typename L>
Hamiltonian<L>::~Hamiltonian() {
    delete[] _adhesionMatrix;
    delete[] _areaLambdas;
    delete[] _areaTargets;
    delete[] _perimeterLambdas;
    delete[] _perimeterTargets;
    delete[] _actLambdas;
    delete[] _actMaxima;
    delete[] _connectedLambdas;
    delete[] _chemotaxisLambdas;
    delete[] _fixedCelltype;
    delete[] _persistenceLambdas;
}

template <typename L>
void Hamiltonian<L>::setFixedCelltype(int type, bool fixed) {
    _fixedCelltype[type] = fixed;
}

template <typename L>
void Hamiltonian<L>::setPersistence(int type, double lambda) {
    _persistenceLambdas[type] = lambda;
}

template <typename L>
bool Hamiltonian<L>::getFixedCelltype(int type) {
    return _fixedCelltype[type];
}

template <typename L>
void Hamiltonian<L>::setConnectedConstraints(int type, double lambda) {
    _connectedLambdas[type] = lambda;
}

template <typename L>
void Hamiltonian<L>::setAdhesionBetweenTypes(int typeA, int typeB, int type) {
    _adhesionMatrix[typeA * _numberOfTypes + typeB] = type;
    _adhesionMatrix[typeB * _numberOfTypes + typeA] = type;
}

template <typename L>
void Hamiltonian<L>::setAreaConstraints(int type, double lambda, int target) {
    _areaLambdas[type] = lambda;
    _areaTargets[type] = target;
}

template <typename L>
void Hamiltonian<L>::setPerimeterConstraints(int type, double lambda, int target) {
    _perimeterLambdas[type] = lambda;
    _perimeterTargets[type] = target;
}

template <typename L>
void Hamiltonian<L>::setActConstraints(int type, double lambda, int max) {
    _actLambdas[type] = lambda;
    _actMaxima[type] = max;
}

template <typename L>
void Hamiltonian<L>::setChemotaxisConstraints(int type, double lambda) {
    _chemotaxisLambdas[type] = lambda;
}

template <typename L>
int Hamiltonian<L>::getAdhesionBetween(LatticePoint& a, LatticePoint& b) {
    return _adhesionMatrix[b.type * _numberOfTypes + a.type];
}

template <typename L>
int Hamiltonian<L>::adhesionDelta(LatticePoint& source, LatticePoint& target, 
        L& lattice) {
    int previousAdhesion = 0;
    int nextAdhesion = 0;
    for (int i = 0; i < lattice.getNeighborCount(); i++) {
        LatticePoint neighbor = lattice.getNeighbor(target, i);
        if (neighbor.cellId != target.cellId)
            previousAdhesion += getAdhesionBetween(target, neighbor);
        if (neighbor.cellId != source.cellId)
            nextAdhesion += getAdhesionBetween(source, neighbor);
    }
    return nextAdhesion - previousAdhesion;
}

template <typename L>
double Hamiltonian<L>::energyAreaDelta(int area, int newArea, LatticePoint& point) {
    int targetArea = _areaTargets[point.type];
    int lambda = _areaLambdas[point.type];
    return lambda *((newArea - targetArea)*(newArea - targetArea) - 
                (area - targetArea)*(area - targetArea));
}

template <typename L>
double Hamiltonian<L>::areasDelta(LatticePoint& source, LatticePoint& target, 
        CellStates<L>& cellStates) {
    double delta = 0;
    if (source.cellId != 0) {
        int area = cellStates.getArea(source.cellId);
        delta += energyAreaDelta(area, area+1, source); 
    }
    if (target.cellId != 0) {
        int area = cellStates.getArea(target.cellId);
        delta += energyAreaDelta(area, area-1, target); 
    }
    return delta;
}

template <typename L>
double Hamiltonian<L>::perimeterDelta(LatticePoint& source, LatticePoint& target, 
        L& lattice, CellStates<L>& cellStates) {
    double delta = 0;

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
        int currentPerimeter = cellStates.getPerimeter(source.cellId);
        int newPerimeter = currentPerimeter - oldPerimeterSource+
            newPerimeterSource;
        int targetPerimeter = _perimeterTargets[source.type];
        double lambda = _perimeterLambdas[source.type];

        delta += 
            lambda *((newPerimeter - targetPerimeter)*(newPerimeter - targetPerimeter) - 
            (currentPerimeter - targetPerimeter)*(currentPerimeter - targetPerimeter));
    }
    if (target.cellId != 0) {
        int currentPerimeter = cellStates.getPerimeter(target.cellId);
        int newPerimeter = currentPerimeter - oldPerimeterTarget +
            newPerimeterTarget;
        int targetPerimeter = _perimeterTargets[target.type];
        double lambda = _perimeterLambdas[target.type];
        delta += 
            lambda *((newPerimeter - targetPerimeter)*(newPerimeter - targetPerimeter) - 
            (currentPerimeter - targetPerimeter)*(currentPerimeter - targetPerimeter));
    }

    return delta;
}

template <typename L>
double Hamiltonian<L>::actDelta(LatticePoint& source, LatticePoint& target, 
        L& lattice, int time) {

    auto lambda = _actLambdas[source.type];
    auto maxAct = _actMaxima[source.type];

    if (source.type == 0) {
        lambda = _actLambdas[target.type];
        maxAct = _actMaxima[target.type];
    }
    if (lambda == 0) 
        return 0;


    double sourceAct = max(_actMaxima[source.type]+source.act-time, 0.0);
    if (source.cellId != 0) {
        const auto maxAct = _actMaxima[source.type];
        int count = 1;
        for (int i = 0; i < lattice.getNeighborCount(); i++) {
            auto neighbor = lattice.getNeighbor(source, i);
            if (neighbor.cellId == source.cellId) {
                auto act = max(maxAct+neighbor.act-time, 0.0);
                sourceAct *= act;
                count++;
            }
        }
        sourceAct = pow(sourceAct, 1.0/count);

    } else {
        sourceAct = 0;
    }


    double targetAct = max(_actMaxima[target.type]+target.act-time, 0.0);
    if (target.cellId != 0) {
        const auto maxAct = _actMaxima[target.type];
        int count = 1;
        for (int i = 0; i < lattice.getNeighborCount(); i++) {
            auto neighbor = lattice.getNeighbor(target, i);
            if (neighbor.cellId == target.cellId) {
                auto act = max(maxAct+neighbor.act-time, 0.0);
                targetAct *= act;
                count++;
            }
        }
        targetAct = pow(targetAct, 1.0/count);
    } else {
        targetAct = 0;
    }


    if (maxAct == 0)
        return 0;
    else
        return -lambda/maxAct * (sourceAct - targetAct);

}

template <typename L>
double Hamiltonian<L>::persistenceDelta(LatticePoint& source, LatticePoint& target,
        L& lattice, Centroids<L>& centroids) {
    if  (_persistenceLambdas[source.type] == 0)
        return 0;

    auto currentDir = target.subtract(source);
    currentDir = currentDir.wrap(lattice._dimension);
    currentDir = currentDir.normalize();
    return -currentDir.dot(centroids.getPrefDir(source.cellId)) * _persistenceLambdas[source.type];


}

template <typename L>
double Hamiltonian<L>::connectedDelta(LatticePoint& source, LatticePoint& target,
        L& lattice) {
    if  (_connectedLambdas[target.type] == 0)
        return 0;
    int transitions = 0;
    for (int i = 0; i < lattice.getNeighborCount(); i++) {
        auto prev = lattice.getNeighbor(target, (i+lattice.getNeighborCount()-1) % lattice.getNeighborCount());
        auto current = lattice.getNeighbor(target, i);
        if ((prev.cellId != target.cellId && current.cellId == target.cellId) ||
            (prev.cellId == target.cellId && current.cellId != target.cellId)) {
            transitions++;
        }
    }
    if (transitions < 3)
        return 0;
    else {
        return _connectedLambdas[target.type];
    }
}

template <typename L>
double Hamiltonian<L>::chemotaxisDelta(LatticePoint& source, LatticePoint& target,
        L& lattice) {

    auto lambda = _chemotaxisLambdas[source.type];
    if (source.type == 0)
        lambda = _chemotaxisLambdas[target.type];
    if (lambda != 0.0) {
        auto fieldVector = lattice.getFieldPoint(source);
        auto delta = target.subtract(source).normalize();
        return -lambda * delta.dot(fieldVector);
    }
    return 0;
}

template <typename L>
double Hamiltonian<L>::directionDelta(LatticePoint& source, LatticePoint& target,
        L& lattice) {
    vec2 v;
    v.x = target.x - source.x;
    v.y = target.y - source.y;
    v = v.normalize();

    vec2 c;
    c.y =  - (lattice._dimension/2 - source.x);
    c.x =  lattice._dimension/2 - source.y;
    c = c.normalize();

    return -c.dot(v) * 10;
}

template <typename L>
double Hamiltonian<L>::energyDelta(LatticePoint& source, LatticePoint& target, 
        L& lattice, CellStates<L>& cellStates, Centroids<L>& centroids,
        ChemokineField* field, int time) {
    double energyDelta = 0;

    energyDelta += areasDelta(source, target, cellStates);
    energyDelta += adhesionDelta(source, target, lattice);
    //energyDelta += directionDelta(source, target, lattice);
    if(_perimeterEnabled)
        energyDelta += perimeterDelta(source, target, lattice, cellStates);
    if(_actEnabled)
        energyDelta += actDelta(source, target, lattice, time);
    if(_connectedEnabled)
        energyDelta += connectedDelta(source, target, lattice);
    if(_persistenceEnabled)
        energyDelta += persistenceDelta(source, target, lattice, centroids);
    //if (field)
    if(_chemotaxisEnabled) 
        energyDelta += chemotaxisDelta(source, target, lattice);
    return energyDelta;
}

template <typename L>
double Hamiltonian<L>::boltzmannProbability(double energyDelta) {
    return exp(-energyDelta/_temperature);

}
template <typename L>
bool Hamiltonian<L>::getActEnabled() {
    return _actEnabled;
}

template <typename L>
void Hamiltonian<L>::updateConstraintToggles() {
    _perimeterEnabled = false;
    _actEnabled = false;
    _chemotaxisEnabled = false;
    _persistenceEnabled = false;
    _connectedEnabled = false;
    for (int i = 0; i < _numberOfTypes; i++) {
        if (_perimeterLambdas[i] != 0) _perimeterEnabled = true;
        if (_actLambdas[i] != 0) _actEnabled = true;
        if (_chemotaxisLambdas[i] != 0) _chemotaxisEnabled = true;
        if (_persistenceLambdas[i] != 0) _persistenceEnabled = true;
        if (_connectedLambdas[i] != 0) _connectedEnabled = true;
    }

}

template class Hamiltonian<Lattice2d>;
template class Hamiltonian<Lattice3d>;
