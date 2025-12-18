// -----------------------#ifndef VEC3_H
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <algorithm>
#include "TRandom3.h"
#include "TMath.h"


struct Vec3
{
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3 &o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3 &o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(double a) const { return Vec3(a * x, a * y, a * z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); } // <- add this
};

static inline double Dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline double Norm(const Vec3 &a) { return std::sqrt(Dot(a, a)); }
static inline Vec3 Unit(const Vec3 &a)
{
    double n = Norm(a);
    if (n <= 0)
        return Vec3(1, 0, 0);
    return a * (1.0 / n);
}

static inline Vec3 ReflectSpecular(const Vec3 &v, const Vec3 &nHat)
{
    return v - nHat * (2.0 * Dot(v, nHat));
}

static inline Vec3 SampleIsotropic(TRandom3 &rng)
{
    double u = rng.Uniform(-1.0, 1.0);
    double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    double s = std::sqrt(std::max(0.0, 1.0 - u * u));
    return Vec3(s * std::cos(phi), s * std::sin(phi), u);
}

// Cosine-weighted (Lambertian) hemisphere around nHat
static inline Vec3 SampleLambert(const Vec3 &nHat, TRandom3 &rng)
{
    Vec3 w = Unit(nHat);
    Vec3 a = (std::fabs(w.x) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    Vec3 u = Unit(Vec3(
        w.y * a.z - w.z * a.y,
        w.z * a.x - w.x * a.z,
        w.x * a.y - w.y * a.x));
    Vec3 v = Vec3(
        w.y * u.z - w.z * u.y,
        w.z * u.x - w.x * u.z,
        w.x * u.y - w.y * u.x);

    double r1 = rng.Uniform();
    double r2 = rng.Uniform();
    double phi = 2.0 * TMath::Pi() * r1;
    double cosTheta = std::sqrt(1.0 - r2);
    double sinTheta = std::sqrt(r2);

    Vec3 dirLocal(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
    Vec3 dir = u * dirLocal.x + v * dirLocal.y + w * dirLocal.z;
    return Unit(dir);
}
#endif
