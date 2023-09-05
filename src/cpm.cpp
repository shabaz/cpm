#include "cpm.h"

using namespace std;

template <typename L>
Cpm<L>::Cpm(int dimension, int numberOfTypes, double temperature):
    _lattice(dimension), _hamiltonian(numberOfTypes, temperature), 
    _centroids(dimension, numberOfTypes, _cellStates), _simulation(_lattice, _hamiltonian, _cellStates, 
    _centroids, nullptr), nrOfCells(0)
{
    _thread = nullptr;
}

template <typename L>
Cpm<L>::~Cpm() {
}

template <typename L>
void Cpm<L>::updateType(int id, int type) { 
    _cellStates.updateType(id, type);
}

template <typename L>
void Cpm<L>::setAdhesionBetweenTypes(int type, int other, int adhesion) { 
    _hamiltonian.setAdhesionBetweenTypes(type, other, adhesion);
}

template <typename L>
void Cpm<L>::setPersistenceConstraints(int type, double lambda, int history, 
        double persistence) { 
    _hamiltonian.setPersistence(type, lambda);
    _centroids.setHistoryLength(type, history);
    _centroids.setPersistence(type, persistence);
}

template <typename L>
void Cpm<L>::setAreaConstraints(int type, double lambda, int target) { 
    _hamiltonian.setAreaConstraints(type, lambda, target);
}

template <typename L>
void Cpm<L>::setPerimeterConstraints(int type, double lambda, int target) { 
    _hamiltonian.setPerimeterConstraints(type, lambda, target);
}

template <typename L>
void Cpm<L>::setActConstraints(int type, double lambda, int max) { 
    _hamiltonian.setActConstraints(type, lambda, max);
}

template <typename L>
void Cpm<L>::setFixedConstraint(int type, bool fixed) { 
    _hamiltonian.setFixedCelltype(type, fixed);
}

template <typename L>
void Cpm<L>::setConnectedConstraints(int type, double lambda) { 
    _hamiltonian.setConnectedConstraints(type, lambda);
}

template <typename L>
void Cpm<L>::setChemotaxisConstraints(int type, double lambda) { 
    _hamiltonian.setChemotaxisConstraints(type, lambda);
}

template <typename L>
void Cpm<L>::updateCellProps(int nrOfCells) {
    _cellStates.initializeFromGrid(_lattice, nrOfCells);
    _centroids.initializeFromGrid(_lattice, nrOfCells);
}

template <typename L>
void Cpm<L>::run(int ticks) {
    _hamiltonian.updateConstraintToggles();
    _lattice.setAct(_hamiltonian.getActEnabled());
    for (int i = 0; i < ticks; i++) {
        _simulation.monteCarloStep();
        
    }
}

template <typename L>
void Cpm<L>::runAsync(int ticks) {
    _thread = new thread(&Cpm::run, this, ticks);
}

template <typename L>
void Cpm<L>::join() {
    _thread->join();
    _thread = nullptr;
}

template <typename L>
unsigned int* Cpm<L>::getData() {
    return _lattice.getCellIds();
}

template <typename L>
double* Cpm<L>::getField() {
    return _lattice._field;
}

template <typename L>
int* Cpm<L>::getActData() {
    return _lattice._actValues;
}

template <typename L>
int Cpm<L>::getDimension() {
    return _lattice._dimension;
}

template <typename L>
std::vector<typename L::Point> Cpm<L>::getCentroids() {
    return _centroids.getCentroids();
}

template class Cpm<Lattice2d>;
template class Cpm<Lattice3d>;
