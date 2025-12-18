// ToyOpticsWedge.C
// ROOT macro for toy optical photon transport in a rectangular scintillator bar
// with optional end light-guide wedges (taper in y only, no taper in z).
//
// Compile in ROOT:
//   .L ToyOpticsWedge.C+
//
// Example run (writes a tree):
//   OpticsConfig cfg;
//   cfg.savePath = true;      // store xPath/yPath/zPath
//   cfg.useWedge = true;      // enable wedges at both ends
//   cfg.wedgeLen = 20.0;      // cm
//   cfg.wedgeTipW = 5.0;      // cm (width at x=0 and x=L)
//   RunManySites(50, 200, 120, 20, 1, 1, cfg, "toyOptics.root");
//
// Draw one event from the file (if path was saved):
//   DrawEventFromTree("toyOptics.root", 0, true);


#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TPaveText.h"
#include "TString.h"
#include "TROOT.h"
#include "TPolyMarker3D.h"
#include "TPad.h"
#include "TLegend.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// -----------------------------
// Basic vector math
// -----------------------------
struct Vec3 {
  double x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
  Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
  Vec3 operator*(double a) const { return Vec3(a * x, a * y, a * z); }
  Vec3 operator-() const { return Vec3(-x, -y, -z); } // <- add this

};

static inline double Dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double Norm(const Vec3& a) { return std::sqrt(Dot(a,a)); }
static inline Vec3 Unit(const Vec3& a) {
  double n = Norm(a);
  if (n <= 0) return Vec3(1,0,0);
  return a*(1.0/n);
}

static inline Vec3 ReflectSpecular(const Vec3& v, const Vec3& nHat) {
  return v - nHat*(2.0*Dot(v, nHat));
}

static inline Vec3 SampleIsotropic(TRandom3& rng) {
  double u = rng.Uniform(-1.0, 1.0);
  double phi = rng.Uniform(0.0, 2.0*TMath::Pi());
  double s = std::sqrt(std::max(0.0, 1.0 - u*u));
  return Vec3(s*std::cos(phi), s*std::sin(phi), u);
}

// Cosine-weighted (Lambertian) hemisphere around nHat
static inline Vec3 SampleLambert(const Vec3& nHat, TRandom3& rng) {
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

  double r1 = rng.Uniform();
  double r2 = rng.Uniform();
  double phi = 2.0*TMath::Pi()*r1;
  double cosTheta = std::sqrt(1.0 - r2);
  double sinTheta = std::sqrt(r2);

  Vec3 dirLocal(std::cos(phi)*sinTheta, std::sin(phi)*sinTheta, cosTheta);
  Vec3 dir = u*dirLocal.x + v*dirLocal.y + w*dirLocal.z;
  return Unit(dir);
}

// -----------------------------
// Optics + geometry helpers
// -----------------------------
static inline bool InsidePMTCircle(double y, double z, double rPMT) {
  return (y*y + z*z) <= (rPMT*rPMT);
}

static inline double epsCoupleExp(double y, double z, double rPMT, double eps0, double lambdaC) {
  double d = std::sqrt(y*y + z*z);
  if (d <= rPMT) return 1.0;
  return eps0 * std::exp(-d / lambdaC);
}

// Wedge logic: taper in y only near both ends.
// Left wedge active for x in [0, Lg], right wedge for x in [L-Lg, L].
// Width shrinks linearly from W at x=Lg (or x=L-Lg) to wTip at x=0 (or x=L).
static inline bool InLeftWedge(double x, double Lg) { return x >= 0.0 && x <= Lg; }
static inline bool InRightWedge(double x, double L, double Lg) { return x >= (L - Lg) && x <= L; }

// -----------------------------
// Config and per-photon output
// -----------------------------
struct OpticsConfig {
  double nScint = 1.58;
  double nOut   = 1.0;

  double absLen = 300.0;
  double Rwrap  = 0.95;
  int maxSteps  = 2000;

  double rPMT     = 2.5;   // cm
  double epsCouple= 0.90;
  double pde      = 0.20;

  // coupling position dependence
  double eps0     = 0.90;
  double lambdaC  = 120.0;

  // toggles
  bool savePath = false;

  // wedge parameters
  bool useWedge   = false;
  double wedgeLen = 20.0;  // cm
  double wedgeTipW= 5.0;   // cm (width at end plane)
};

struct PhotonResult {
  // stamped inputs
  double L=0, W=0, T=0;
  int wrap=1;            // 1=PTFE, 2=Mylar
  int site_number=0;

  // site + end point
  double x0=0,y0=0,z0=0;
  double xf=0,yf=0,zf=0;

  // outcome
  int endPlane=-1;       // 0..5 box, 6/7 wedge walls
  int pmt_side=0;        // 1=x=0, 2=x=L
  int inPMT=0;
  int detected=0;
  int absorbed=0;
  int escaped=0;
  int reachedEnd=0;

  int nBounces=0;
  double path=0;
  double epsRaysPMT=0;

  // optional path
  std::vector<float> xPath, yPath, zPath;
};

// -----------------------------
// Propagate one photon from a given site position
// -----------------------------
PhotonResult PropagateOnePhoton(
  TRandom3& rng,
  const Vec3& sitePos,
  int site_number,
  double L, double W, double T,
  int wrap,
  const OpticsConfig& cfg
){
  PhotonResult r;
  r.L=L; r.W=W; r.T=T; r.wrap=wrap; r.site_number=site_number;
  r.x0=sitePos.x; r.y0=sitePos.y; r.z0=sitePos.z;

  const int PTFE  = (wrap == 1);
  const int mylar = (wrap == 2);

  double sinThetaC = cfg.nOut / cfg.nScint;
  if (sinThetaC >= 1.0) sinThetaC = 0.999999;
  double cosThetaC = std::cos(std::asin(sinThetaC));

  Vec3 pos = sitePos;
  Vec3 dir = Unit(SampleIsotropic(rng));

  if (cfg.savePath) {
    r.xPath.clear(); r.yPath.clear(); r.zPath.clear();
    r.xPath.push_back((float)pos.x);
    r.yPath.push_back((float)pos.y);
    r.zPath.push_back((float)pos.z);
  }

  for (int step=0; step<cfg.maxSteps; step++) {

    // Candidate intersections with standard box
    double tx=1e99; int xPlane=-1;
    if (dir.x > 1e-15) { tx = (L - pos.x)/dir.x; xPlane=1; }
    if (dir.x < -1e-15) { tx = (0.0 - pos.x)/dir.x; xPlane=0; }

    double tyP=1e99, tyM=1e99;
    if (dir.y > 1e-15) tyP = ( W/2.0 - pos.y)/dir.y;
    if (dir.y < -1e-15) tyM = (-W/2.0 - pos.y)/dir.y;

    double tzP=1e99, tzM=1e99;
    if (dir.z > 1e-15) tzP = ( T/2.0 - pos.z)/dir.z;
    if (dir.z < -1e-15) tzM = (-T/2.0 - pos.z)/dir.z;

    // Optional wedge planes (replace y=+-W/2 inside wedge regions)
    double tWedgeP=1e99, tWedgeM=1e99;
    int wedgeSide=0; // 0 none, 1 left, 2 right

    if (cfg.useWedge) {
      if (InLeftWedge(pos.x, cfg.wedgeLen)) wedgeSide = 1;
      else if (InRightWedge(pos.x, L, cfg.wedgeLen)) wedgeSide = 2;

      if (wedgeSide != 0 && cfg.wedgeLen > 0 && cfg.wedgeTipW > 0 && cfg.wedgeTipW <= W) {

        // slope b for yMax(x) = wTip/2 + b*x (left side)
        // for right side it is mirrored.
        double b = (0.5*W - 0.5*cfg.wedgeTipW) / cfg.wedgeLen;

        // We represent the two wedge planes in a linear implicit form and solve for t.
        // Left wedge:  y - (a + b*x) = 0   and   y + (a + b*x) = 0
        // Right wedge: y - (a + b*(L-x)) = 0 and symmetric
        double a = 0.5*cfg.wedgeTipW;

        // Encode both left and right with y - a0 + bx*x = 0 for the + wall
        double a0 = a;
        double bx = -b; // left: y - a - b*x = 0
        if (wedgeSide == 2) {
          a0 = a + b*L; // y - (a + b*L - b*x) = 0 -> y - a0 + b*x = 0
          bx = +b;
        }

        // + wedge surface
        {
          double denom = (dir.y + bx*dir.x);
          if (std::fabs(denom) > 1e-15) {
            double t = (a0 - pos.y - bx*pos.x)/denom;
            if (t > 1e-12) {
              double xh = pos.x + t*dir.x;
              double zh = pos.z + t*dir.z;
              bool okx = (wedgeSide==1) ? (xh>=0.0 && xh<=cfg.wedgeLen) : (xh>=L-cfg.wedgeLen && xh<=L);
              if (okx && std::fabs(zh) <= T/2.0) tWedgeP = t;
            }
          }
        }

        // - wedge surface
        {
          double denom = (dir.y - bx*dir.x);
          if (std::fabs(denom) > 1e-15) {
            double t = (-a0 - pos.y + bx*pos.x)/denom;
            if (t > 1e-12) {
              double xh = pos.x + t*dir.x;
              double zh = pos.z + t*dir.z;
              bool okx = (wedgeSide==1) ? (xh>=0.0 && xh<=cfg.wedgeLen) : (xh>=L-cfg.wedgeLen && xh<=L);
              if (okx && std::fabs(zh) <= T/2.0) tWedgeM = t;
            }
          }
        }

        // Inside wedge region, the wedge surfaces replace the full y=+-W/2 walls
        tyP = 1e99;
        tyM = 1e99;
      }
    }

    // Pick the nearest intersection
    double tmin=1e99; int plane=-1;
    if (tx  > 1e-12 && tx  < tmin) { tmin=tx;  plane=xPlane; }
    if (tyM > 1e-12 && tyM < tmin) { tmin=tyM; plane=2; }
    if (tyP > 1e-12 && tyP < tmin) { tmin=tyP; plane=3; }
    if (tzM > 1e-12 && tzM < tmin) { tmin=tzM; plane=4; }
    if (tzP > 1e-12 && tzP < tmin) { tmin=tzP; plane=5; }
    if (tWedgeP > 1e-12 && tWedgeP < tmin) { tmin=tWedgeP; plane=6; }
    if (tWedgeM > 1e-12 && tWedgeM < tmin) { tmin=tWedgeM; plane=7; }

    if (plane < 0 || tmin > 1e98) break;

    // Bulk absorption
    if (cfg.absLen > 0) {
      double survive = std::exp(-tmin / cfg.absLen);
      if (rng.Uniform() > survive) {
        r.path += tmin;
        r.absorbed = 1;
        r.endPlane = plane;

        r.xf = pos.x + dir.x*tmin;
        r.yf = pos.y + dir.y*tmin;
        r.zf = pos.z + dir.z*tmin;

        if (cfg.savePath) {
          r.xPath.push_back((float)r.xf);
          r.yPath.push_back((float)r.yf);
          r.zPath.push_back((float)r.zf);
        }
        return r;
      }
    }

    // Move to boundary
    pos = pos + dir*tmin;
    r.path += tmin;
    r.endPlane = plane;
    r.xf = pos.x; r.yf = pos.y; r.zf = pos.z;

    if (cfg.savePath) {
      r.xPath.push_back((float)pos.x);
      r.yPath.push_back((float)pos.y);
      r.zPath.push_back((float)pos.z);
    }

    // End planes: detection
    if (plane==0 || plane==1) {
      r.reachedEnd = 1;
      r.pmt_side = (plane==0) ? 1 : 2;
      r.inPMT = InsidePMTCircle(pos.y, pos.z, cfg.rPMT) ? 1 : 0;

      r.epsRaysPMT = epsCoupleExp(pos.y, pos.z, cfg.rPMT, cfg.eps0, cfg.lambdaC);
      double pEnd = cfg.epsCouple * cfg.pde * r.epsRaysPMT;
      if (rng.Uniform() < pEnd) r.detected = 1;

      return r;
    }

    // Surface normal
    Vec3 nHat;
    if (plane==2) nHat = Vec3(0,-1,0);
    if (plane==3) nHat = Vec3(0,+1,0);
    if (plane==4) nHat = Vec3(0,0,-1);
    if (plane==5) nHat = Vec3(0,0,+1);

    if (plane==6 || plane==7) {
      // wedge normals, approximate: slope dy/dx = +/-b
      double b = (cfg.wedgeLen > 0) ? (0.5*W - 0.5*cfg.wedgeTipW) / cfg.wedgeLen : 0.0;
      double sgn = (plane==6) ? +1.0 : -1.0;
      nHat = Unit(Vec3(-b, sgn, 0));
    }

    nHat = Unit(nHat);

    double cosInc = std::fabs(Dot(dir, nHat));
    bool isTIR = (cosInc < cosThetaC);

    if (isTIR) {
      dir = Unit(ReflectSpecular(dir, nHat));
      r.nBounces++;
      continue;
    }

    // Non-TIR: reflect with probability Rwrap, else escape
    if (rng.Uniform() >= cfg.Rwrap) {
      r.escaped = 1;
      return r;
    }

    // Wrap model: PTFE -> diffuse, Mylar -> specular
    if (PTFE && !mylar) dir = SampleLambert(-nHat, rng);
    else dir = Unit(ReflectSpecular(dir, nHat));

    r.nBounces++;
  }

  return r;
}

// -----------------------------
// Tree writer
// -----------------------------
class TreeWriter {
public:
  TreeWriter(const char* outFile, bool savePath)
  : savePath_(savePath)
  {
    fout_ = new TFile(outFile, "RECREATE");
    t_ = new TTree("tPhot","Toy optical photon transport");

    t_->Branch("L",&r_.L,"L/D");
    t_->Branch("W",&r_.W,"W/D");
    t_->Branch("T",&r_.T,"T/D");
    t_->Branch("wrap",&r_.wrap,"wrap/I");
    t_->Branch("site_number",&r_.site_number,"site_number/I");

    t_->Branch("x0",&r_.x0,"x0/D");
    t_->Branch("y0",&r_.y0,"y0/D");
    t_->Branch("z0",&r_.z0,"z0/D");

    t_->Branch("xf",&r_.xf,"xf/D");
    t_->Branch("yf",&r_.yf,"yf/D");
    t_->Branch("zf",&r_.zf,"zf/D");

    t_->Branch("epsRaysPMT",&r_.epsRaysPMT,"epsRaysPMT/D");

    t_->Branch("endPlane",&r_.endPlane,"endPlane/I");
    t_->Branch("pmt_side",&r_.pmt_side,"pmt_side/I");
    t_->Branch("inPMT",&r_.inPMT,"inPMT/I");
    t_->Branch("detected",&r_.detected,"detected/I");
    t_->Branch("absorbed",&r_.absorbed,"absorbed/I");
    t_->Branch("escaped",&r_.escaped,"escaped/I");
    t_->Branch("reachedEnd",&r_.reachedEnd,"reachedEnd/I");

    t_->Branch("nBounces",&r_.nBounces,"nBounces/I");
    t_->Branch("path",&r_.path,"path/D");

    if (savePath_) {
      t_->Branch("xPath",&r_.xPath);
      t_->Branch("yPath",&r_.yPath);
      t_->Branch("zPath",&r_.zPath);
    }
  }

  void Fill(const PhotonResult& in) {
    r_ = in;
    if (!savePath_) {
      r_.xPath.clear();
      r_.yPath.clear();
      r_.zPath.clear();
    }
    t_->Fill();
  }

  void Close() {
    fout_->cd();
    t_->Write();
    fout_->Close();
    delete fout_;
    fout_ = nullptr;
    t_ = nullptr;
  }

private:
  TFile* fout_ = nullptr;
  TTree* t_ = nullptr;
  PhotonResult r_;
  bool savePath_ = false;
};

// -----------------------------
// Wrapper: many sites and many photons per site
// -----------------------------
void RunManySites(
  long long Nsites,
  long long NphotPerSite,
  double L, double W, double T,
  int wrap,
  const OpticsConfig& cfg,
  const char* outFile="toyOptics.root"
){
  TRandom3 rng(0);
  TreeWriter wr(outFile, cfg.savePath);

  long long nDet = 0;
  for (long long j=0; j<Nsites; j++) {

    Vec3 site(
      rng.Uniform(0.0, L),
      rng.Uniform(-W/2.0, W/2.0),
      rng.Uniform(-T/2.0, T/2.0)
    );

    for (long long i=0; i<NphotPerSite; i++) {
      PhotonResult res = PropagateOnePhoton(rng, site, (int)j, L, W, T, wrap, cfg);
      if (res.detected) nDet++;
      wr.Fill(res);
    }
  }

  wr.Close();

  double frac = (Nsites>0 && NphotPerSite>0) ? double(nDet)/double(Nsites*NphotPerSite) : 0.0;
  std::cout << "RunManySites: detected fraction = " << frac
            << " (Nsites=" << Nsites << ", NphotPerSite=" << NphotPerSite << ")\n";
}

// -----------------------------
// Drawing helpers
// -----------------------------
static TPolyLine3D* Edge(double x1,double y1,double z1,double x2,double y2,double z2){
  auto e = new TPolyLine3D(2);
  e->SetPoint(0,x1,y1,z1);
  e->SetPoint(1,x2,y2,z2);
  e->SetLineWidth(2);
  return e;
}

static void DrawBoxWireframe(double L,double W,double T){
  double x0=0, x1=L;
  double y0=-W/2, y1p=+W/2;
  double z0=-T/2, z1p=+T/2;

  // x=0
  Edge(x0,y0,z0, x0,y1p,z0)->Draw("same");
  Edge(x0,y1p,z0, x0,y1p,z1p)->Draw("same");
  Edge(x0,y1p,z1p, x0,y0,z1p)->Draw("same");
  Edge(x0,y0,z1p, x0,y0,z0)->Draw("same");

  // x=L
  Edge(x1,y0,z0, x1,y1p,z0)->Draw("same");
  Edge(x1,y1p,z0, x1,y1p,z1p)->Draw("same");
  Edge(x1,y1p,z1p, x1,y0,z1p)->Draw("same");
  Edge(x1,y0,z1p, x1,y0,z0)->Draw("same");

  // connect
  Edge(x0,y0,z0, x1,y0,z0)->Draw("same");
  Edge(x0,y1p,z0, x1,y1p,z0)->Draw("same");
  Edge(x0,y1p,z1p, x1,y1p,z1p)->Draw("same");
  Edge(x0,y0,z1p, x1,y0,z1p)->Draw("same");
}

// Wedge outline (visual only): drawn in the z=0 plane as an x-y outline
static void DrawWedgeOutlineXY(double L,double W,double wTip,double Lg){
  double z=0.0;

  // left end
  Edge(0, -wTip/2, z, 0, +wTip/2, z)->Draw("same");
  Edge(Lg, -W/2, z, Lg, +W/2, z)->Draw("same");
  Edge(0, +wTip/2, z, Lg, +W/2, z)->Draw("same");
  Edge(0, -wTip/2, z, Lg, -W/2, z)->Draw("same");

  // right end
  Edge(L, -wTip/2, z, L, +wTip/2, z)->Draw("same");
  Edge(L-Lg, -W/2, z, L-Lg, +W/2, z)->Draw("same");
  Edge(L, +wTip/2, z, L-Lg, +W/2, z)->Draw("same");
  Edge(L, -wTip/2, z, L-Lg, -W/2, z)->Draw("same");
}
static void DrawEventInfoBox(
  Long64_t ievt,
  double L, double W, double T,
  int wrap,
  int reachedEnd, int pmt_side, int inPMT, int detected,
  int absorbed, int escaped,
  int endPlane, int nBounces,
  double pathLen
){
  const char* wrapName = "wrap=?";
  if (wrap == 1) wrapName = "PTFE";
  else if (wrap == 2) wrapName = "Mylar";

  auto p = new TPaveText(0.12, 0.78, 0.88, 0.93, "NDC"); // x1,y1,x2,y2
  p->SetFillColor(0);
  p->SetFillStyle(0);
  p->SetBorderSize(1);
  p->SetTextAlign(12);   // left, vertically centered
  p->SetTextSize(0.028); // adjust if needed

  p->AddText(Form("event %lld | geometry: L=%.1f cm, W=%.1f cm, T=%.1f cm", (long long)ievt, L, W, T));
  p->AddText(Form("wrap: %s | endPlane=%d | bounces=%d | path=%.2f cm", wrapName, endPlane, nBounces, pathLen));
  p->AddText(Form("end status: reachedEnd=%d (pmtSide=%d, inPMT=%d, detected=%d) | absorbed=%d | escaped=%d",
                  reachedEnd, pmt_side, inPMT, detected, absorbed, escaped));

  p->Draw();
}
struct PathMarkers {
  TPolyMarker3D* start = nullptr;
  TPolyMarker3D* end   = nullptr;
};

static PathMarkers DrawStartEndMarkers3D(
  double xs, double ys, double zs,
  double xe, double ye, double ze
){
  PathMarkers pm;

  pm.start = new TPolyMarker3D(1);
  pm.start->SetPoint(0, xs, ys, zs);
  pm.start->SetMarkerStyle(20);
  pm.start->SetMarkerSize(1.6);
  pm.start->SetMarkerColor(kGreen+2);
  pm.start->Draw("same");

  pm.end = new TPolyMarker3D(1);
  pm.end->SetPoint(0, xe, ye, ze);
  pm.end->SetMarkerStyle(29);
  pm.end->SetMarkerSize(2.0);
  pm.end->SetMarkerColor(kRed+1);
  pm.end->Draw("same");

  return pm;
}


void DrawEventFromTree(
  const char* fn="toyOptics.root",
  Long64_t ievt=0,
  bool drawWedge=true,
  double wedgeLen=20.0,
  double wedgeTipW=5.0,
  bool zoom=true,
  double zoomMargin=0.15 // 15% padding around the path
){
  TFile f(fn);
  auto t = (TTree*)f.Get("tPhot");
  if (!t) {
    std::cout << "DrawEventFromTree: cannot find tPhot\n";
    return;
  }

  // Always present in the tree (per our writer)
  double L=0,W=0,T=0;
  int wrap=0;

  int reachedEnd=0, pmt_side=0, inPMT=0, detected=0;
  int absorbed=0, escaped=0, endPlane=0, nBounces=0;
  double pathLen=0;

  t->SetBranchAddress("L",&L);
  t->SetBranchAddress("W",&W);
  t->SetBranchAddress("T",&T);
  t->SetBranchAddress("wrap",&wrap);

  // These should exist if you used the TreeWriter block I gave
  if (t->GetBranch("reachedEnd")) t->SetBranchAddress("reachedEnd",&reachedEnd);
  if (t->GetBranch("pmt_side"))   t->SetBranchAddress("pmt_side",&pmt_side);
  if (t->GetBranch("inPMT"))      t->SetBranchAddress("inPMT",&inPMT);
  if (t->GetBranch("detected"))   t->SetBranchAddress("detected",&detected);

  if (t->GetBranch("absorbed"))   t->SetBranchAddress("absorbed",&absorbed);
  if (t->GetBranch("escaped"))    t->SetBranchAddress("escaped",&escaped);
  if (t->GetBranch("endPlane"))   t->SetBranchAddress("endPlane",&endPlane);
  if (t->GetBranch("nBounces"))   t->SetBranchAddress("nBounces",&nBounces);
  if (t->GetBranch("path"))       t->SetBranchAddress("path",&pathLen);

  // Path branches exist only if you ran with cfg.savePath=true
  std::vector<float>* xP=nullptr;
  std::vector<float>* yP=nullptr;
  std::vector<float>* zP=nullptr;

  bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
  if (hasPath) {
    t->SetBranchAddress("xPath",&xP);
    t->SetBranchAddress("yPath",&yP);
    t->SetBranchAddress("zPath",&zP);
  }

  t->GetEntry(ievt);

  const char* wrapName = "wrap=?";
  if (wrap==1) wrapName = "PTFE";
  else if (wrap==2) wrapName = "Mylar";

  // Build a descriptive title

  TString title;
  title.Form("evt %lld",(long long) ievt);
 /* title.Form(
    "evt %lld | geom L=%.1f W=%.1f T=%.1f cm | wrap=%s | reachedEnd=%d pmtSide=%d inPMT=%d detected=%d | absorbed=%d escaped=%d | endPlane=%d bounces=%d path=%.2f cm; x (cm); y (cm); z (cm)",
    (long long)ievt, L, W, T, wrapName,
    reachedEnd, pmt_side, inPMT, detected,
    absorbed, escaped,
    endPlane, nBounces, pathLen
  );
*/

  auto c = new TCanvas("cEvt","event",1100,800);

  // Default view ranges
  double xmin=0, xmax=L;
  double ymin=-W/2, ymax=+W/2;
  double zmin=-T/2, zmax=+T/2;

  // Zoom range around the path if requested and available
  if (zoom && hasPath && xP && xP->size() >= 2) {
    double xlo=(*xP)[0], xhi=(*xP)[0];
    double ylo=(*yP)[0], yhi=(*yP)[0];
    double zlo=(*zP)[0], zhi=(*zP)[0];

    for (size_t i=1;i<xP->size();i++) {
      xlo = std::min(xlo, (double)(*xP)[i]); xhi = std::max(xhi, (double)(*xP)[i]);
      ylo = std::min(ylo, (double)(*yP)[i]); yhi = std::max(yhi, (double)(*yP)[i]);
      zlo = std::min(zlo, (double)(*zP)[i]); zhi = std::max(zhi, (double)(*zP)[i]);
    }

    // Add margin
    double dx = (xhi-xlo); if (dx<=0) dx = 1.0;
    double dy = (yhi-ylo); if (dy<=0) dy = 1.0;
    double dz = (zhi-zlo); if (dz<=0) dz = 1.0;

    xlo -= zoomMargin*dx; xhi += zoomMargin*dx;
    ylo -= zoomMargin*dy; yhi += zoomMargin*dy;
    zlo -= zoomMargin*dz; zhi += zoomMargin*dz;

    // Keep inside detector bounds (optional but usually nice)
    xlo = std::max(xlo, 0.0);  xhi = std::min(xhi, L);
    ylo = std::max(ylo, -W/2); yhi = std::min(yhi, W/2);
    zlo = std::max(zlo, -T/2); zhi = std::min(zhi, T/2);

    xmin=xlo; xmax=xhi;
    ymin=ylo; ymax=yhi;
    zmin=zlo; zmax=zhi;
  } else {
    // Keep a 1:1:1 feel by expanding to a cube around the detector
    double maxDim = std::max(L, std::max(W, T));
    double xPad = 0.5*(maxDim - L);
    xmin = -xPad; xmax = L + xPad;
    ymin = -0.5*maxDim; ymax = 0.5*maxDim;
    zmin = -0.5*maxDim; zmax = 0.5*maxDim;
  }

  auto frame = new TH3D("fr", title,
                        10, xmin, xmax,
                        10, ymin, ymax,
                        10, zmin, zmax);
  frame->SetStats(0);
  frame->Draw();

  DrawBoxWireframe(L,W,T);
  if (drawWedge) DrawWedgeOutlineXY(L,W,wedgeTipW,wedgeLen);


  if (hasPath && xP && xP->size() >= 2) {
    auto pl = new TPolyLine3D((int)xP->size());
    for (int i=0;i<(int)xP->size();i++) pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
    pl->SetLineWidth(3);
    pl->Draw("same");
    const int n = (int)xP->size();
auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                (*xP)[n-1], (*yP)[n-1], (*zP)[n-1]);

auto leg = new TLegend(0.62, 0.22, 0.99, 0.36); // x1,y1,x2,y2 in NDC
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.05);

leg->AddEntry(pl, "Photon path", "l");
leg->AddEntry(pm.start, "Start point", "p");
leg->AddEntry(pm.end, "End point", "p");
leg->Draw();

  } else {
    std::cout << "DrawEventFromTree: no saved path in file (run with cfg.savePath=true)\n";
  }

  DrawEventInfoBox(ievt, L, W, T, wrap,
                 reachedEnd, pmt_side, inPMT, detected,
                 absorbed, escaped,
                 endPlane, nBounces, pathLen);

  c->Modified();
  c->Update();
}

void DrawEventFromTreeOld(
  const char* fn="toyOptics.root",
  Long64_t ievt=0,
  bool drawWedge=true,
  double wedgeLen=20.0,
  double wedgeTipW=5.0
){
  TFile f(fn);
  auto t = (TTree*)f.Get("tPhot");
  if (!t) {
    std::cout << "DrawEventFromTree: cannot find tPhot\n";
    return;
  }

  double L=0,W=0,T=0;
  int wrap=0;

  std::vector<float>* xP=nullptr;
  std::vector<float>* yP=nullptr;
  std::vector<float>* zP=nullptr;

  t->SetBranchAddress("L",&L);
  t->SetBranchAddress("W",&W);
  t->SetBranchAddress("T",&T);
  t->SetBranchAddress("wrap",&wrap);

  // These branches exist only if you ran with cfg.savePath=true
  bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
  if (hasPath) {
    t->SetBranchAddress("xPath",&xP);
    t->SetBranchAddress("yPath",&yP);
    t->SetBranchAddress("zPath",&zP);
  }

  t->GetEntry(ievt);

  auto c = new TCanvas("cEvt","event",900,700);

  // Force a 1:1:1 feel by drawing a cubic frame (equal axis numeric spans)
  double maxDim = std::max(L, std::max(W, T));
  double xPad = 0.5*(maxDim - L);

  auto frame = new TH3D("fr","Scintillator and one photon; x (cm); y (cm); z (cm)",
                        10, -xPad, L + xPad,
                        10, -0.5*maxDim, 0.5*maxDim,
                        10, -0.5*maxDim, 0.5*maxDim);
  frame->SetStats(0);
  frame->Draw();

  DrawBoxWireframe(L,W,T);
  if (drawWedge) DrawWedgeOutlineXY(L,W,wedgeTipW,wedgeLen);

  if (hasPath && xP && xP->size() >= 2) {
    auto pl = new TPolyLine3D((int)xP->size());
    for (int i=0;i<(int)xP->size();i++) pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
    pl->SetLineWidth(3);
    pl->Draw("same");
    const int n = (int)xP->size();
auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                (*xP)[n-1], (*yP)[n-1], (*zP)[n-1]);

auto leg = new TLegend(0.62, 0.22, 0.99, 0.36); // x1,y1,x2,y2 in NDC
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.05);

leg->AddEntry(pl, "Photon path", "l");
leg->AddEntry(pm.start, "Start point", "p");
leg->AddEntry(pm.end, "End point", "p");
leg->Draw();


  } else {
    std::cout << "DrawEventFromTree: no saved path in file (run with cfg.savePath=true)\n";
  }

  c->Modified();
  c->Update();
}


void DrawEventSplitViewFromTree(
  const char* fn="toyOptics.root",
  Long64_t ievt=0,
  bool drawWedge=true,
  double wedgeLen=20.0,
  double wedgeTipW=5.0,
  double zoomMargin=0.12   // padding for the zoomed-in view
){
  TFile f(fn);
  auto t = (TTree*)f.Get("tPhot");
  if (!t) { std::cout << "DrawEventSplitViewFromTree: cannot find tPhot\n"; return; }

  double L=0,W=0,T=0;
  int wrap=0;
  int reachedEnd=0, pmt_side=0, inPMT=0, detected=0;
  int absorbed=0, escaped=0, endPlane=0, nBounces=0;
  double pathLen=0;

  t->SetBranchAddress("L",&L);
  t->SetBranchAddress("W",&W);
  t->SetBranchAddress("T",&T);
  t->SetBranchAddress("wrap",&wrap);

  if (t->GetBranch("reachedEnd")) t->SetBranchAddress("reachedEnd",&reachedEnd);
  if (t->GetBranch("pmt_side"))   t->SetBranchAddress("pmt_side",&pmt_side);
  if (t->GetBranch("inPMT"))      t->SetBranchAddress("inPMT",&inPMT);
  if (t->GetBranch("detected"))   t->SetBranchAddress("detected",&detected);

  if (t->GetBranch("absorbed"))   t->SetBranchAddress("absorbed",&absorbed);
  if (t->GetBranch("escaped"))    t->SetBranchAddress("escaped",&escaped);
  if (t->GetBranch("endPlane"))   t->SetBranchAddress("endPlane",&endPlane);
  if (t->GetBranch("nBounces"))   t->SetBranchAddress("nBounces",&nBounces);
  if (t->GetBranch("path"))       t->SetBranchAddress("path",&pathLen);

  std::vector<float>* xP=nullptr;
  std::vector<float>* yP=nullptr;
  std::vector<float>* zP=nullptr;

  bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
  if (!hasPath) { std::cout << "No xPath/yPath/zPath in file. Run with cfg.savePath=true\n"; return; }

  t->SetBranchAddress("xPath",&xP);
  t->SetBranchAddress("yPath",&yP);
  t->SetBranchAddress("zPath",&zP);

  t->GetEntry(ievt);
  if (!xP || xP->size() < 2) { std::cout << "Event has no stored path points\n"; return; }

  // Compute zoomed-in bounds from path
  double xlo=(*xP)[0], xhi=(*xP)[0];
  double ylo=(*yP)[0], yhi=(*yP)[0];
  double zlo=(*zP)[0], zhi=(*zP)[0];

  for (size_t i=1;i<xP->size();i++) {
    xlo = std::min(xlo, (double)(*xP)[i]); xhi = std::max(xhi, (double)(*xP)[i]);
    ylo = std::min(ylo, (double)(*yP)[i]); yhi = std::max(yhi, (double)(*yP)[i]);
    zlo = std::min(zlo, (double)(*zP)[i]); zhi = std::max(zhi, (double)(*zP)[i]);
  }

  double dx = xhi-xlo; if (dx<=0) dx=1;
  double dy = yhi-ylo; if (dy<=0) dy=1;
  double dz = zhi-zlo; if (dz<=0) dz=1;

  xlo -= zoomMargin*dx; xhi += zoomMargin*dx;
  ylo -= zoomMargin*dy; yhi += zoomMargin*dy;
  zlo -= zoomMargin*dz; zhi += zoomMargin*dz;

  // Clamp zoom-in to detector bounds
  xlo = std::max(xlo, 0.0);  xhi = std::min(xhi, L);
  ylo = std::max(ylo, -W/2); yhi = std::min(yhi,  W/2);
  zlo = std::max(zlo, -T/2); zhi = std::min(zhi,  T/2);

  // Zoomed-out cube bounds (for 1:1:1 feel)
  double maxDim = std::max(L, std::max(W, T));
  double xPad = 0.5*(maxDim - L);

  double xo_min = -xPad, xo_max = L + xPad;
  double yo_min = -0.5*maxDim, yo_max = 0.5*maxDim;
  double zo_min = -0.5*maxDim, zo_max = 0.5*maxDim;

  // Unique canvas name
  TString cname = Form("cSplit_%lld", (long long)ievt);
  if (auto old = gROOT->FindObject(cname)) delete old;

  auto c = new TCanvas(cname, "event split view", 1400, 1000);
  c->Divide(1,2);

  auto drawOnePad = [&](int ipad,
                        double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,
                        const char* padTitle)
  {
    c->cd(ipad);

    TString frname = Form("fr_%lld_%d", (long long)ievt, ipad);
    auto frame = new TH3D(frname, Form("%s; x (cm); y (cm); z (cm)", padTitle),
                          10, xmin, xmax,
                          10, ymin, ymax,
                          10, zmin, zmax);
    frame->SetDirectory(nullptr);
    frame->SetStats(0);
    frame->Draw();

    DrawBoxWireframe(L,W,T);
    if (drawWedge) DrawWedgeOutlineXY(L,W,wedgeTipW,wedgeLen);

    auto pl = new TPolyLine3D((int)xP->size());
//    pl->SetName(Form("pl_%lld_%d", (long long)ievt, ipad));
    for (int i=0;i<(int)xP->size();i++) pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
    pl->SetLineWidth(3);
    pl->Draw("same");
const int n = (int)xP->size();
auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                (*xP)[n-1], (*yP)[n-1], (*zP)[n-1]);

    // Put info box only once (left pad) or on both if you want
    if (ipad == 1) {
      DrawEventInfoBox(ievt, L, W, T, wrap,
                       reachedEnd, pmt_side, inPMT, detected,
                       absorbed, escaped,
                       endPlane, nBounces, pathLen);
                       auto leg = new TLegend(0.62, 0.22, 0.99, 0.36); // x1,y1,x2,y2 in NDC
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.05);

leg->AddEntry(pl, "Photon path", "l");
leg->AddEntry(pm.start, "Start point", "p");
leg->AddEntry(pm.end, "End point", "p");
leg->Draw();

    }
  };

  drawOnePad(1, xo_min,xo_max, yo_min,yo_max, zo_min,zo_max, "Zoom out");
  drawOnePad(2, xlo,xhi,     ylo,yhi,     zlo,zhi,     "Zoom in");

  c->Modified();
  c->Update();
}
