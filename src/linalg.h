#ifndef LINALG_H_
#define LINALG_H_

#include <vector>
#include <math.h>

struct vec2 {
    float x;
    float y;

    vec2& operator+=(const vec2& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    vec2 operator/(const float& denominator) {
        vec2 v = *this;
        v.x /= denominator;
        v.y /= denominator;
        return v;
    }

    vec2 operator-(const vec2& rhs) {
        vec2 v = *this;
        v.x -= rhs.x;
        v.y -= rhs.y;
        return v;
    }

    vec2 normalize() {
        vec2 v = *this;
        float length = sqrt(x * x + y * y);
        if (length > 0) {
            v.x /= length;
            v.y /= length;
        }
        return v;
    }

    float dot(const vec2& rhs) {
        return x * rhs.x + y * rhs.y;
    }
};


struct vec3 {
    float x;
    float y;
    float z;

    vec3& operator+=(const vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    vec3 operator/(const float& denominator) {
        vec3 v = *this;
        v.x /= denominator;
        v.y /= denominator;
        v.z /= denominator;
        return v;
    }

    vec3 operator-(const vec3& rhs) {
        vec3 v = *this;
        v.x -= rhs.x;
        v.y -= rhs.y;
        v.z -= rhs.z;
        return v;
    }

    vec3 normalize() {
        vec3 v = *this;
        float length = sqrt(x * x + y * y + z * z);
        if (length > 0) {
            v.x /= length;
            v.y /= length;
            v.z /= length;
        }
        return v;
    }

    float dot(const vec3& rhs) {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }
};


struct mat2 {
    float v[4];
};

vec2 calculateAverage(std::vector<vec2> points);
vec3 calculateAverage(std::vector<vec3> points);
mat2 calculateCovarianceMatrix(std::vector<vec2> points, vec2 average);
vec2 calculateEigenVector(mat2 m, float eigenvalue);
std::vector<float> calculateEigenValues(mat2 m);

#endif // LINALG_H_
