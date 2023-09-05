#ifndef CPM_H_
#define CPM_H_

#include <thread>
#include <type_traits>



#include "lattice.h"
#include "cell_states.h"
#include "hamiltonian.h"
#include "centroids.h"
#include "simulation.h"



template <typename L>
class Cpm {
    public:
        typedef typename L::Point Point;
        Cpm(int dimension, int numberOfTypes, double temperature);
        ~Cpm();
        void setFixedConstraint(int type, bool fixed);
        void setAdhesionBetweenTypes(int type, int other, int adhesion);
        void setAreaConstraints(int type, double lambda, int target);
        void setPerimeterConstraints(int type, double lambda, int target);
        void setActConstraints(int type, double lambda, int max);
        void setConnectedConstraints(int type, double lambda);
        void setChemotaxisConstraints(int type, double lambda);
        void setPersistenceConstraints(int type, double lambda, int history, 
                double persistence);
        void updateType(int id, int type);
        void run(int ticks);
        void runAsync(int ticks);
        void join();
        unsigned int* getData();
        double* getField();
        int* getActData();
        void updateCellProps(int nrOfCells);
        int getDimension();

        std::vector<Point> getCentroids();

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice2d>::value, int>::type = 0>
            void updatePoint(int x, int y, int type) {
                auto source = _lattice.getPoint(x, y);
                auto target = source;
                source.type = type;
                source.cellId = nrOfCells;
                _lattice.copy(source, target, source.act);
                _cellStates.updateAreas(source, target);
                _cellStates.updatePerimeters(source, target, _lattice);
                _centroids.update(source, target);
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice3d>::value, int>::type = 0>
            void updatePoint(int x, int y, int z, int type) {
                auto source = _lattice.getPoint(x, y, z);
                auto target = source;
                source.type = type;
                source.cellId = nrOfCells;
                _lattice.copy(source, target, source.act);
                _cellStates.updateAreas(source, target);
                _cellStates.updatePerimeters(source, target, _lattice);
                _centroids.update(source, target);
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice2d>::value, int>::type = 0>
            void setPoint(int x, int y, int cellId, int type) {
                _lattice.setPoint(cellId, x, y, 0, type);
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice3d>::value, int>::type = 0>
            void setPoint(int x, int y, int z, int cellId, int type) {
                _lattice.setPoint(cellId, x, y, z, 0, type);
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice2d>::value, int>::type = 0>
            void addCell(int type) {
                auto cellId = _cellStates.nextId();
                _cellStates.addCell(0, 0, type);
                _centroids.addCentroid({0,0}, 0);
                nrOfCells++;
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice3d>::value, int>::type = 0>
            void addCell(int type) {
                auto cellId = _cellStates.nextId();
                _cellStates.addCell(1, 26, type);
                _centroids.addCentroid({0,0,0}, 0);
                nrOfCells++;
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice2d>::value, int>::type = 0>
            void addCell(int x, int y, int type) {
                auto cellId = _cellStates.nextId();
                _lattice.setPoint(cellId, x, y, 0, type);
                _cellStates.addCell(1, 8, type);
                _centroids.addCentroid({x,y}, 1);
                nrOfCells++;
            }

        template<typename U = L, typename std::enable_if<
            std::is_same<U, Lattice3d>::value, int>::type = 0>
            void addCell(int x, int y, int z, int type) {
                auto cellId = _cellStates.nextId();
                _lattice.setPoint(cellId, x, y, z, 0, type);
                _cellStates.addCell(1, 26, type);
                _centroids.addCentroid({x,y,z}, 1);
                nrOfCells++;
            }

    private:
        int nrOfCells;
        L _lattice;
        CellStates<L> _cellStates;
        Centroids<L> _centroids;
        Hamiltonian<L> _hamiltonian;
        Simulation<L> _simulation;
        std::thread* _thread;
};


#endif // CPM_H_
