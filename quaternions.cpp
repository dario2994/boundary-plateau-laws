#ifndef QUATERNIONS
#define QUATERNIONS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"

template<typename T>
struct Quaternion {
    T q[4];
    Quaternion(): q{T(0),T(0),T(0),T(0)} {}
    Quaternion(T a, T b, T c, T d): q{a,b,c,d} {}
    Quaternion(T r): q{r, T(0), T(0), T(0)} {}

    T& operator [](int i) { return q[i]; }
    const T& operator [](int i) const { return q[i]; }
};

template<typename T>
Quaternion<T> operator +(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] + y[i];
    return z;
}

template<typename T>
Quaternion<T> operator -(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] - y[i];
    return z;
}

template<typename T>
Quaternion<T> operator *(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    z[0] = x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3];
    z[1] = x[0]*y[1] + x[1]*y[0] + x[2]*y[3] - x[3]*y[2];
    z[2] = x[2]*y[0] + x[0]*y[2] + x[3]*y[1] - x[1]*y[3];
    z[3] = x[3]*y[0] + x[0]*y[3] + x[1]*y[2] - x[2]*y[1];
    return z;
}

template<typename T>
Quaternion<T>& operator *=(Quaternion<T>& x, const Quaternion<T>& y) {
    x = x * y;
    return x;
}

template<typename T>
Quaternion<T> conj(const Quaternion<T>& z) {
    return Quaternion<T>(z[0], -z[1], -z[2], -z[3]);
}


template<typename T>
ostream& operator <<(ostream& out, const Quaternion<T>& q) {
    out << q[0] << " + i" << q[1] << " + j" << q[2] << " + k" << q[3];
    return out;
}

#endif
