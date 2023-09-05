#ifndef CENTROIDS_H
#define CENTROIDS_H

#include <vector>
#include <list>

template <typename L> class CellStates;

template <typename L>
class Centroids {
    public:
        typedef typename L::LatticePoint LatticePoint;
        typedef typename L::Point Point;
        typedef typename L::IntPoint IntPoint;

        Centroids(int dimension, int numberOfTypes, CellStates<L>& cellStates);
        ~Centroids();

        void addCentroid(IntPoint center, int count);
        void update(LatticePoint& source, LatticePoint& target);
        std::vector<Point> getCentroids();
        void addCheckpoint();
        void updatePreferentialDirection();
        void print();
        void setHistoryLength(int type, int historyLength);
        void setPersistence(int type, double persistence);
        Point getPrefDir(int cellId);
        void initializeFromGrid(L& lattice, int nrOfCells);
    private:
        std::vector<IntPoint> _centers;
        std::vector<int> _counts;

        std::vector<std::list<Point>> _history;
        std::vector<Point> _preferredDirections;
        int _dimension;
        int* _historyLengths;
        double* _persistenceValues;
        CellStates<L>& _cellStates;
};

#endif // CENTROIDS_H
