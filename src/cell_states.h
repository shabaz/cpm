#ifndef CELL_STATES_H
#define CELL_STATES_H

#include <vector>

template <typename L>
class CellStates {
    public:
        typedef typename L::LatticePoint LatticePoint;
        void addCell(int area, int perimeter, int type);
        int getArea(int cellId);
        void removeArea(int cellId, int area);
        int getPerimeter(int cellId);
        int getType(int cellId);
        void updateType(int cellId, int type);
        void updateAreas(LatticePoint& source, LatticePoint& target);
        void updatePerimeters(LatticePoint& source, LatticePoint& target, 
                L& lattice);
        void print();
        int nextId();
        void initializeFromGrid(L& lattice, int nrOfCells);
        void recalcPerimeter(L& lattice, int id);
        int countType(int type);
        void kill(int id);
        std::vector<int> getCellIds(int type);
    private:
        std::vector<int> _areas;
        std::vector<int> _perimeters;
        std::vector<int> _types;
};

#endif // CELL_STATES_H
