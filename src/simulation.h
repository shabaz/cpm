#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "ranxoshi256.h"

template <typename L> class Hamiltonian;
template <typename L> class CellStates;
template <typename L> class Centroids;

class ChemokineField;


template <typename L>
class Simulation {
    public:
        typedef typename L::LatticePoint LatticePoint;
        Simulation(L& lattice, Hamiltonian<L>& hamiltonian, CellStates<L>& cellStates, 
                 Centroids<L>& centroids, ChemokineField* field);
        void monteCarloStep();
        void stratifiedMonteCarloStep(int* indices);
        int copyAttempt(LatticePoint& source, LatticePoint& target);
        int _time;
    private:
        L& _lattice;
        Hamiltonian<L>& _hamiltonian;
        CellStates<L>& _cellStates;
        Centroids<L>& _centroids;
        ChemokineField* _field;
        ranxoshi256 xoshi;
}; 

#endif // SIMULATION_H_
