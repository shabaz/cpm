#include <random>
#include <iostream>
#include "simulation.h"
#include "lattice.h"
#include "hamiltonian.h"
#include "cell_states.h"
#include "centroids.h"

using namespace std;


class ChemokineField {
};

template <typename L>
Simulation<L>::Simulation(L& lattice, Hamiltonian<L>& hamiltonian, 
        CellStates<L>& cellStates, Centroids<L>& centroids, ChemokineField* field): 
    _lattice(lattice), _hamiltonian(hamiltonian), _cellStates(cellStates), 
    _centroids(centroids), _field(field) {
    random_device rd;
    _time = 0;
    unsigned int seed[8];
    for (int i = 0; i < 8; i++)
        seed[i] = rd();
    ranxoshi256Seed(&xoshi, (unsigned char*)&seed);
}

template <typename L>
void Simulation<L>::monteCarloStep() {
    _time += 1;
    float timestep = 0;
    int succesful = 0;
    int attempts = 0;
    while (timestep < 1) {
        if (_lattice._borderIndices.size() == 0)
            return;
        timestep += 1.0/_lattice._borderIndices.size();
        auto source = _lattice.getRandomBorderLocation(xoshi);
        auto target = _lattice.getRandomNeighbor(source, xoshi);
        if(_hamiltonian.getFixedCelltype(source.type)) {
            source.type = 0;
            source.cellId = 0;
        }
        if (source.cellId != target.cellId && 
                !_hamiltonian.getFixedCelltype(target.type)

                ) {
            attempts ++;
            succesful += copyAttempt(source, target);
        }
    }
    _centroids.addCheckpoint();
    _centroids.updatePreferentialDirection();
    //cout << "attempts: " << attempts << " succesful attempts this MCS: " << succesful << endl;
}

template <typename L>
void Simulation<L>::stratifiedMonteCarloStep(int* indices) {
    /*_time += 1;
    int succesful = 0;
    int attempts = 0;

    const int offsets[4][2] = {{0,0}, {0,1}, {1,0}, {1,1}};

    //for (int k = 0; k < 2; k++) {
        //for (int l = 0; l < 2; l++) {
        for (int w = 0; w < 4; w++) {
            int k = offsets[indices[w]][0];
            int l = offsets[indices[w]][1];
            for (int i = 0; i < _lattice._dimension/4; i++) {
                for (int j = 0; j < _lattice._dimension/4; j++) {
                    //for (int a = 0; a < 2; a++) {
                        //for (int b = 0; b < 2; b++) {
                    for (int q = 0; q < 4; q++) {
                        int a = ranxoshi256Next(&xoshi) % 2;
                        int b = ranxoshi256Next(&xoshi) % 2;


                        auto source = _lattice.getPoint(i * 4 + a + k * 2,j * 4 + b + l * 2);
                        auto target = _lattice.getRandomNeighbor(source, xoshi);
                        if (source.cellId != target.cellId) {
                            attempts ++;
                            succesful += copyAttempt(source, target);
                        }
                    }
                }


                        //}
                    //}
                }
            }
        //}
    //}



    //cout << "attempts: " << attempts << " succesful attempts this MCS: " << succesful << endl;
    */
}




template <typename L>
int Simulation<L>::copyAttempt(LatticePoint& source, LatticePoint& target) {
    auto energyDelta = _hamiltonian.energyDelta(source, target, _lattice, 
            _cellStates, _centroids, _field, _time);
    uniform_real_distribution<double> dist(0, 1);
    if (energyDelta < 0 || 
            _hamiltonian.boltzmannProbability(energyDelta) > 
            ranxoshi256DoubleCO(&xoshi)) {
        
        _lattice.copy(source, target, _time);
        _cellStates.updateAreas(source, target);
        _cellStates.updatePerimeters(source, target, _lattice);
        _centroids.update(source, target);
        return 1;
    }
    return 0;
}

template class Simulation<Lattice3d>;
template class Simulation<Lattice2d>;
