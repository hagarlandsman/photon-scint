// ToyScintOptics.C
// ROOT macro: toy optical photon transport in a rectangular scintillator bar
// Usage in ROOT:
//   .x ToyScintOptics.C(200000)
// or compiled:
//   .x ToyScintOptics.C+(200000)

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include <cmath>
#include <algorithm>
#include <iostream>

struct Vec3 {
  double x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
  Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
  Vec3 operator*(double a) const { return Vec3(a*x, a*y, a*z); }
};

static inline double Dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double Norm(const Vec3& a) { return std::sqrt(Dot(a,a)); }
static inline Vec3 Unit(const Vec3& a) {
  double n = Norm(a);
  if (n <= 0) return Vec3(1,0,0);
  return a*(1.0/n);
}
static inline Vec3 ReflectSpecular(const Vec3& v, const Vec3& nHat) {
  // v' = v - 2 (v·n) n
  return v - nHat*(2.0*Dot(v, nHat));
}

// Sample Lambertian direction in hemisphere around nHat
static inline Vec3 SampleLambert(const Vec3& nHat, TRandom3& rng) {
  // Build an orthonormal basis (u, v, w=nHat)
  Vec3 w = Unit(nHat);
  Vec3 a = (std::fabs(w.x) < 0.9) ? Vec3(1,0,0) : Vec3(0,1,0);
  Vec3 u = Unit(Vec3(
    w.y*a.z - w.z*a.y,
    w.z*a.x - w.x*a.z,
    w.x*a.y - w.y*a.x
  ));
  Vec3 v = Vec3(
    w.y*u.z - w.z*u.y,
    w.z*u.x - w.x*u.z,
    w.x*u.y - w.y*u.x
  );

  // Cosine-weighted hemisphere sampling
  double r1 = rng.Uniform();
  double r2 = rng.Uniform();
  double phi = 2.0*TMath::Pi()*r1;
  double cosTheta = std::sqrt(1.0 - r2);
  double sinTheta = std::sqrt(r2);

  Vec3 dirLocal(
    std::cos(phi)*sinTheta,
    std::sin(phi)*sinTheta,
    cosTheta
  );

  // Transform to global
  Vec3 dir = u*dirLocal.x + v*dirLocal.y + w*dirLocal.z;
  return Unit(dir);
}

// Uniform isotropic direction on sphere
static inline Vec3 SampleIsotropic(TRandom3& rng) {
  double u = rng.Uniform(-1.0, 1.0);
  double phi = rng.Uniform(0.0, 2.0*TMath::Pi());
  double s = std::sqrt(std::max(0.0, 1.0 - u*u));
  return Vec3(s*std::cos(phi), s*std::sin(phi), u);
}

double ToyScintOptics(
  long long N = 200000,
  double L = 90.0,   // cm (x direction) (distance between PMT)
  double W = 90.0,   // cm (y direction)
  double T = 1,    // cm (z direction)
  double nScint = 1.58,
  double nOut = 1.0,
  double absLen = 300.0, //300.0,  // cm bulk attenuation length
  double Rwrap = 0.9,    // wrapping reflectivity when TIR fails (0 means no wrap help)
  int maxSteps = 2000
){
  TRandom3 rng(0);

  // Critical condition for TIR:
  // TIR if incidence angle alpha > theta_c
  // Equivalent: |v·n|/|v| < cos(theta_c)
  double sinThetaC = nOut / nScint;
  if (sinThetaC >= 1.0) sinThetaC = 0.999999;
  double thetaC = std::asin(sinThetaC);
  double cosThetaC = std::cos(thetaC);

  std::cout << "ToyScintOptics\n";
  std::cout << "L,W,T = " << L << " " << W << " " << T << " cm\n";
  std::cout << "nScint = " << nScint << ", nOut = " << nOut << "\n";
  std::cout << "theta_c = " << thetaC*180.0/TMath::Pi() << " deg, cos(theta_c) = " << cosThetaC << "\n";
  std::cout << "absLen = " << absLen << " cm, Rwrap = " << Rwrap << "\n";

  TH1D* hEnd = new TH1D("hEnd", "Which end reached;end (0=none,1=x=0,2=x=L)", 3, -0.5, 2.5);
  TH1D* hBounces = new TH1D("hBounces", "Number of reflections;reflections;photons", 60, -0.5, 59.5);
  TH1D* hPath = new TH1D("hPath", "Total path length traveled;path length (cm);photons", 120, 0.0, 600.0);
  TH2D* hXYhit = new TH2D("hXYhit", "Hit position at exit end; y (cm); z (cm)", 120, -W/2.0, W/2.0, 80, -T/2.0, T/2.0);

  long long nReach0 = 0;
  long long nReachL = 0;
  long long nAbsorb = 0;
  long long nEscape = 0;

  for (long long i = 0; i < N; i++) {

    // Random emission point inside bar
    Vec3 pos(
      rng.Uniform(0.0, L),
      rng.Uniform(-W/2.0, W/2.0),
      rng.Uniform(-T/2.0, T/2.0)
    );

    // Isotropic direction
    Vec3 dir = Unit(SampleIsotropic(rng));

    int bounces = 0;
    double path = 0.0;
    int endReached = 0;

    for (int step = 0; step < maxSteps; step++) {

      // Compute distance to each boundary plane along dir
      // Planes: x=0, x=L, y=-W/2, y=+W/2, z=-T/2, z=+T/2
      double tx = 1e99, ty = 1e99, tz = 1e99;
      int hitPlane = -1;

      // x planes
      if (dir.x > 1e-15) {
        tx = (L - pos.x) / dir.x; // to x=L
        hitPlane = 1;
      } else if (dir.x < -1e-15) {
        tx = (0.0 - pos.x) / dir.x; // to x=0 (dir.x negative makes positive tx)
        hitPlane = 0;
      }

      // y planes
      double ty1 = 1e99, ty2 = 1e99;
      if (dir.y > 1e-15) ty1 = ( W/2.0 - pos.y) / dir.y;
      if (dir.y < -1e-15) ty2 = (-W/2.0 - pos.y) / dir.y;

      // z planes
      double tz1 = 1e99, tz2 = 1e99;
      if (dir.z > 1e-15) tz1 = ( T/2.0 - pos.z) / dir.z;
      if (dir.z < -1e-15) tz2 = (-T/2.0 - pos.z) / dir.z;

      // Pick nearest positive intersection
      double tmin = tx;
      int plane = -1;
      if (tx > 1e-12) plane = (hitPlane == 0) ? 0 : 1;

      if (ty1 > 1e-12 && ty1 < tmin) { tmin = ty1; plane = 3; } // y=+W/2
      if (ty2 > 1e-12 && ty2 < tmin) { tmin = ty2; plane = 2; } // y=-W/2
      if (tz1 > 1e-12 && tz1 < tmin) { tmin = tz1; plane = 5; } // z=+T/2
      if (tz2 > 1e-12 && tz2 < tmin) { tmin = tz2; plane = 4; } // z=-T/2

      if (plane < 0 || tmin > 1e98) {
        // numerical issue
        break;
      }

      // Bulk absorption over distance tmin
      // survival probability exp(-tmin/absLen)
      if (absLen > 0) {
        double survive = std::exp(-tmin / absLen);
        if (rng.Uniform() > survive) {
          path += tmin;
          nAbsorb++;
          endReached = 0;
          break;
        }
      }

      // Move to boundary
      pos = pos + dir * tmin;
      path += tmin;

      // Check if we reached an end (x plane). If yes, count as detected at that end.
      if (plane == 0) { // x=0
        endReached = 1;
        nReach0++;
        hXYhit->Fill(pos.y, pos.z);
        break;
      }
      if (plane == 1) { // x=L
        endReached = 2;
        nReachL++;
        hXYhit->Fill(pos.y, pos.z);
        break;
      }

      // Otherwise, we hit a side face. Decide reflect or escape.
      Vec3 nHat;
      if (plane == 2) nHat = Vec3(0, -1, 0); // y=-W/2 outward normal
      if (plane == 3) nHat = Vec3(0, +1, 0); // y=+W/2
      if (plane == 4) nHat = Vec3(0, 0, -1); // z=-T/2
      if (plane == 5) nHat = Vec3(0, 0, +1); // z=+T/2
      nHat = Unit(nHat);

      double cosInc = std::fabs(Dot(dir, nHat)); // since dir is unit, this is cos(angle to normal)
      bool isTIR = (cosInc < cosThetaC);

      if (isTIR) {
        dir = Unit(ReflectSpecular(dir, nHat));
        bounces++;
        continue;
      } else {
        // Not TIR. Either escape, or reflect off wrapping with probability Rwrap.
        if (rng.Uniform() < Rwrap) {
          // Diffuse Lambertian reflection back into the scintillator
          // Use inward normal: -nHat
          dir = SampleLambert(nHat * (-1.0), rng);
          bounces++;
          continue;
        } else {
          nEscape++;
          endReached = 0;
          break;
        }
      }
    }

    hEnd->Fill(endReached);
    hBounces->Fill(std::min(bounces, 59));
    hPath->Fill(std::min(path, 599.0));
  }

  double f0 = double(nReach0) / double(N);
  double fL = double(nReachL) / double(N);
  double fTot = double(nReach0 + nReachL) / double(N);

  std::cout << "\nResults (toy):\n";
  std::cout << "  Reached x=0: " << nReach0 << " (" << f0 << ")\n";
  std::cout << "  Reached x=L: " << nReachL << " (" << fL << ")\n";
  std::cout << "  Total reached an end: " << (nReach0 + nReachL) << " (" << fTot << ")\n";
  std::cout << "  Absorbed in bulk: " << nAbsorb << "\n";
  std::cout << "  Escaped (no reflect): " << nEscape << "\n";
  std::cout << "  Effective area: " << (fTot)*L*W << "\n";

  // Draw quick diagnostics
  TCanvas* c1 = new TCanvas("c1", "ToyScintOptics", 1200, 900);
  c1->Divide(2,2);
  c1->cd(1); hEnd->Draw();
  c1->cd(2); hBounces->Draw();
  c1->cd(3); hPath->Draw();
  c1->cd(4); hXYhit->Draw("COLZ");
  c1->Update();

  return((fTot)* L * W );
}

