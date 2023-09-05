#include "linalg.h"
#include <math.h>

using namespace std;

vec2 calculateAverage(std::vector<vec2> points) {
    vec2 avg = {0,0};
    for (auto point: points) {
        avg += point;
    }
    return avg / points.size();
}

vec3 calculateAverage(std::vector<vec3> points) {
    vec3 avg = {0,0,0};
    for (auto point: points) {
        avg += point;
    }
    return avg / points.size();
}

mat2 calculateCovarianceMatrix(std::vector<vec2> points, vec2 average) {
    float xx = 0;
    float yy = 0;
    float xy = 0;
    for (auto p: points) {
        auto point = p - average;
        xx += point.x * point.x;
        yy += point.y * point.y;
        xy += point.x * point.y;
    }
    return {xx, xy, xy, yy};
}

std::vector<float> calculateEigenValues(mat2 m) {
    float trace = m.v[1] + m.v[2];
    float determinant = m.v[1] * m.v[2] - m.v[0] * m.v[3];

    float val1 = sqrt(0.25 * trace * trace - determinant) + 0.5 * trace;
    float val2 = -sqrt(0.25 * trace * trace - determinant) + 0.5 * trace;
    return {val1, val2};
}

vec2 calculateEigenVector(mat2 m, float eigenvalue) {
    if (m.v[2] != 0) {
        return {eigenvalue - m.v[3],  m.v[2]};
    } else if (m.v[1] != 0) {
        return {m.v[1], eigenvalue - m.v[0]};
    } else {
        return {1, 0};
    }
}
