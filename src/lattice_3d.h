#ifndef LATTICE3D_H_
#define LATTICE3D_H_

#include <vector>
#include <iostream>
#include <string>
#include "dice_set.h"
#include "linalg.h"

using namespace std;

struct ranxoshi256;

class Lattice3d {
    public:

        struct IntPoint;
        struct Point;
        struct LatticePoint {
            unsigned int cellId;
            char type;
            int x;
            int y;
            int z;
            int act;

            IntPoint times(int r) {
                IntPoint a;
                a.x = x * r;
                a.y = y * r;
                a.z = z * r;
                return a;
            }

            Point subtract(LatticePoint& rhs) {
                Point p;
                p.x = x - rhs.x;
                p.y = y - rhs.y;
                p.z = z - rhs.z;
                return p;
            }




            void print() {
                cout << x << " " << y <<  " " << z << endl;
            }
        };

        struct Point {
            double x;
            double y;
            double z;
            void print() {
                cout << x << " " << y <<  " " << z << endl;
            }

            double length() {
                return sqrt(x*x + y*y + z*z);
            }
            Point subtract(Point& rhs) {
                Point p;
                p.x = x - rhs.x;
                p.y = y - rhs.y;
                p.z = z - rhs.z;
                return p;
            }

            Point add(const Point& rhs) {
                Point p;
                p.x = x + rhs.x;
                p.y = y + rhs.y;
                p.z = z + rhs.z;
                return p;
            }

            double dot(const Point& rhs) {
                return rhs.x * x + rhs.y * y + rhs.z * z;
            }

            Point times(double r) {
                Point p;
                p.x = x * r;
                p.y = y * r;
                p.z = z * r;
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
                if (p.z > dimension/2)
                    p.z -= dimension;
                if (p.z < -dimension/2)
                    p.z += dimension;
                return p;
            }

            Point normalize() {
                Point p = *this;
                auto l = length();
                p.x = x / l;
                p.y = y / l;
                p.z = z / l;
                return p;
            }

            void unitRandomize() {
                random_device rng;
                mt19937 urng(rng());
                normal_distribution<double> distribution(0.0,1.0);
                x = distribution(urng);
                y = distribution(urng);
                z = distribution(urng);
                *this = normalize();
            }
        };

        struct IntPoint {
            int x;
            int y;
            int z;

            IntPoint(int a, int b, int c) {
                x = a;
                y = b;
                z = c;
            }

            IntPoint() {
                x = 0;
                y = 0;
                z = 0;
            }

            bool isZero(int dimension) {
                return x > dimension or x < -dimension or
                    y > dimension or y < -dimension or
                    z > dimension or z < -dimension;
            }

            void print(string s) {
                cout << s << endl;
                cout << x << " " << y <<  " " << z << endl;
            }

            IntPoint add(const IntPoint& p) {
                IntPoint r;
                r.x = x + p.x;
                r.y = y + p.y;
                r.z = z + p.z;
                return r;
            }

            IntPoint subtract(const IntPoint& p) {
                IntPoint r;
                r.x = x - p.x;
                r.y = y - p.y;
                r.z = z - p.z;
                return r;
            }

            IntPoint add(const LatticePoint& p) {
                IntPoint r;
                r.x = x + p.x;
                r.y = y + p.y;
                r.z = z + p.z;
                return r;
            }

            IntPoint subtract(const LatticePoint& p) {
                IntPoint r;
                r.x = x - p.x;
                r.y = y - p.y;
                r.z = z - p.z;
                return r;
            }

            Point subtract(const Point& p) {
                Point r;
                r.x = x - p.x;
                r.y = y - p.y;
                r.z = z - p.z;
                return r;
            }

            IntPoint modulo(int m) {
                IntPoint r;
                r.x = (x+m)%m;
                r.y = (y+m)%m;
                r.z = (z+m)%m;
                return r;
            }

            IntPoint mul(int r) {
                IntPoint p;
                p.x = x * r;
                p.y = y * r;
                p.z = z * r;
                return p;
            }

            IntPoint intDiv(int r) {
                IntPoint p;
                p.x = x/r;
                p.y = y/r;
                p.z = z/r;
                return p;
            }

            Point divide(double r) {
                Point p;
                p.x = x / r;
                p.y = y / r;
                p.z = z / r;
                return p;
            }
        };


        Lattice3d(int dimension);
        int getNeighborCount();
        void setPoint(int cellId, int x, int y, int z, int time, int type);
        LatticePoint getPoint(int x, int y, int z);
        LatticePoint getPoint(int i);
        int size();
        LatticePoint getRandomBorderLocation(ranxoshi256& xoshi);
        LatticePoint getRandomNeighbor(LatticePoint& point, ranxoshi256& xoshi);
        LatticePoint getNeighbor(LatticePoint& point, int i);
        LatticePoint getNeighbor(int x, int y, int z, int i);
        Point getCenterOfMass(int id);
        bool isPartOfBorder(int x, int y, int z);
        void copy(LatticePoint& source, LatticePoint& target, int time);
        void updateBorderTrackingAround(int x, int y, int z);
        int index(unsigned int x, unsigned int y, unsigned int z);
        std::vector<vec3> getPoints(int cellId);
        void setPoints(int id, const std::vector<vec3>& points, int type);
        unsigned int* getCellIds();
        void resetType(int cellId, int type);
        void remove(int id);
        Point getFieldPoint(LatticePoint& point);
        void setAct(bool actToggle);
        ~Lattice3d();

        DiceSet _borderIndices;
        unsigned int* _cellIds;
        int* _actValues;
        double* _field;
        const int _dimension;

        bool _actToggle;


    private:
}; 

#endif // LATTICE_H_
