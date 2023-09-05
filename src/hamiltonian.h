#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

template <typename L> class CellStates;
class ChemokineField;
template <typename L> class Centroids;


template <typename L>
class Hamiltonian {
    public:
        typedef typename L::LatticePoint LatticePoint;
        Hamiltonian(int types, double temperature);
        ~Hamiltonian();
        void setAdhesionBetweenTypes(int typeA, int typeB, int type);
        void setAreaConstraints(int type, double lambda, int target);
        void setPerimeterConstraints(int type, double lambda, int target);
        void setActConstraints(int type, double lambda, int max);
        void setConnectedConstraints(int type, double lambda);
        void setChemotaxisConstraints(int type, double lambda);
        void setFixedCelltype(int type, bool fixed);
        void setPersistence(int type, double lambda);
        int getAdhesionBetween(LatticePoint& a, LatticePoint& b);
        bool getFixedCelltype(int type);
        void updateConstraintToggles();
        double persistenceDelta(LatticePoint& source, LatticePoint& target,
                L& lattice, Centroids<L>& centroids);
        int adhesionDelta(LatticePoint& source, LatticePoint& target, 
                L& lattice);
        double areasDelta(LatticePoint& source, LatticePoint& target, 
                CellStates<L>& cellStates);
        double perimeterDelta(LatticePoint& source, LatticePoint& target, 
                L& lattice, CellStates<L>& cellStates);
        double actDelta(LatticePoint& source, LatticePoint& target, 
                L& lattice, int time);
        double connectedDelta(LatticePoint& source, LatticePoint& target,
                L& lattice);
        double chemotaxisDelta(LatticePoint& source, LatticePoint& target,
                L& lattice);
        double directionDelta(LatticePoint& source, LatticePoint& target, 
                L& lattice);
        double energyDelta(LatticePoint& source, LatticePoint& target, 
                L& lattice, CellStates<L>& cellStates, Centroids<L>& centroids,
                ChemokineField* field, int time);
        double boltzmannProbability(double energyDelta);
        bool getActEnabled();
    private:
        double energyAreaDelta(int area, int newArea, LatticePoint& point);
        int* _adhesionMatrix;
        double* _areaLambdas;
        double* _areaTargets;
        double* _perimeterLambdas;
        double* _perimeterTargets;
        double* _actLambdas;
        double* _actMaxima;
        double* _connectedLambdas;
        double* _chemotaxisLambdas;
        bool* _fixedCelltype;
        double* _persistenceLambdas;
        int _numberOfTypes;
        double _temperature;

        bool _perimeterEnabled;
        bool _actEnabled;
        bool _chemotaxisEnabled;
        bool _persistenceEnabled;
        bool _connectedEnabled;
};

#endif // HAMILTONIAN_H_
