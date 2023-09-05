#ifndef LATTICE2D_H_
#define LATTICE2D_H_

#include <vector>
#include <iostream>
#include <random>
#include <string>
#include "dice_set.h"
#include "linalg.h"

using namespace std;

struct ranxoshi256;

class Lattice2d {
    public:
        struct IntPoint;
        struct Point;
        struct LatticePoint {
            unsigned int cellId;
            char type;
            int x;
            int y;
            int act;

            IntPoint times(int r) {
                IntPoint a;
                a.x = x * r;
                a.y = y * r;
                return a;
            }

            Point subtract(LatticePoint& rhs) {
                Point p;
                p.x = x - rhs.x;
                p.y = y - rhs.y;
                return p;
            }

            void print() {
                cout << x << " " << y << endl;
            }

        };

        struct Point {
            double x;
            double y;
            void print() {
                cout << x << " " << y << endl;
            }
            double length() {
                return sqrt(x*x + y*y);
            }
            Point subtract(Point& rhs) {
                Point p;
                p.x = x - rhs.x;
                p.y = y - rhs.y;
                return p;
            }
            Point add(const Point& rhs) {
                Point p;
                p.x = x + rhs.x;
                p.y = y + rhs.y;
                return p;
            }

            double dot(const Point& rhs) {
                return rhs.x * x + rhs.y * y;
            }

            Point times(double r) {
                Point p;
                p.x = x * r;
                p.y = y * r;
                return p;
            }


            Point wrap(int dimension) {
                Point p = *this;
                if (p.x > dimension/2)
                    p.x -= dimension;
                if (p.x < -dimension/2)
                    p.x += dimension;
                if (p.y > dimension/2)
                    p.y -= dimension;
                if (p.y < -dimension/2)
                    p.y += dimension;
                return p;
            }

            Point normalize() {
                Point p = *this;
                auto l = length();
                p.x = x / l;
                p.y = y / l;
                return p;
            }

            void unitRandomize() {
                random_device rng;
                mt19937 urng(rng());
                normal_distribution<double> distribution(0.0,1.0);
                x = distribution(urng);
                y = distribution(urng);
                *this = normalize();
            }
        }
        ;

        struct IntPoint {
            int x;
            int y;

            IntPoint(int a, int b) {
                x = a;
                y = b;
            }

            IntPoint() {
                x = 0;
                y = 0;
            }

            bool isZero(int dimension) {
                return x == 0 and y == 0;
            }

            void print(string s) {
                cout << x << " " << y <<  endl;
            }

            IntPoint add(const IntPoint& p) {
                IntPoint r;
                r.x = x + p.x;
                r.y = y + p.y;
                return r;
            }

            IntPoint subtract(const IntPoint& p) {
                IntPoint r;
                r.x = x - p.x;
                r.y = y - p.y;
                return r;
            }

            Point subtract(const Point& p) {
                Point r;
                r.x = x - p.x;
                r.y = y - p.y;
                return r;
            }

            IntPoint add(const LatticePoint& p) {
                IntPoint r;
                r.x = x + p.x;
                r.y = y + p.y;
                return r;
            }

            IntPoint subtract(const LatticePoint& p) {
                IntPoint r;
                r.x = x - p.x;
                r.y = y - p.y;
                return r;
            }

            IntPoint modulo(int m) {
                IntPoint r;
                r.x = (x+m)%m;
                r.y = (y+m)%m;
                return r;
            }

            IntPoint mul(int r) {
                IntPoint p;
                p.x = x * r;
                p.y = y * r;
                return p;
            }

            IntPoint intDiv(int r) {
                IntPoint p;
                p.x = x/r;
                p.y = y/r;
                return p;
            }

            Point divide(double r) {
                Point p;
                p.x = x / r;
                p.y = y / r;
                return p;
            }
                
        };

        Lattice2d(int dimension);
        int getNeighborCount();
        void setPoint(int cellId, int x, int y, int time, int type);
        LatticePoint getPoint(int i);
        LatticePoint getPoint(int x, int y);
        int size();
        LatticePoint getRandomBorderLocation(ranxoshi256& xoshi);
        LatticePoint getRandomNeighbor(LatticePoint& point, ranxoshi256& xoshi);
        LatticePoint getNeighbor(LatticePoint& point, int i);
        LatticePoint getNeighbor(int x, int y, int i);
        Point getCenterOfMass(int id);
        bool isPartOfBorder(int x, int y);
        void copy(LatticePoint& source, LatticePoint& target, int time);
        void updateBorderTrackingAround(int x, int y);
        int index(unsigned int x, unsigned int y);
        std::vector<vec2> getPoints(int cellId);
        void setPoints(int id, const std::vector<vec2>& points, int type);
        void resetType(int cellId, int type);
        void remove(int id);
        Point getFieldPoint(LatticePoint& point);
        unsigned int* getCellIds();
        void setAct(bool actToggle);
        ~Lattice2d();

        DiceSet _borderIndices;
        unsigned int* _cellIds;
        unsigned int* _cellIdsReshaped;
        double* _field;
        int* _actValues;
        const int _dimension;
        bool _actToggle;
    private:
}; 

#endif // LATTICE_H_
