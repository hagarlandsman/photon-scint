// ScanToyOptics_withTree.C
// ROOT macro: toy optical photon transport in a rectangular scintillator bar
// Includes optional TTree output with per-photon details.
// Usage example in ROOT:
//
//   .L ScanToyOptics_withTree.C+
//   TH2D *hFrac=0, *hC1=0, *hC2=0;
//   ScanToyOptics(hFrac, hC1, hC2, 20000);  // scan map (no tree output)
//
//   // Single run with tree output (writes toyOptics.root):
//   TRandom3 rng(0);
//   double f = ToyScintOptics(rng, 200000, 90.0, 30.0, 1, 1.0, 1.58, 1.0, 300.0, 0.95, 2000,
//                            2.5, 0.90, 0.20,
//                            true, "toyOptics.root", 90.0, 30.0, 1.0);
//
// Then inspect:
//   TFile f("toyOptics.root");
//   tPhot->Scan("x0:y0:z0:pmt_side:inPMT:detected:nBounces:path:xf:yf:zf:endPlane","","colsize=12");
/*
TRandom3 rng(0);

// Geometry 1
.L ScanToyOptics_withTree.C+
//   TH2D *hFrac=0, *hC1=0, *hC2=0;
ToyScintOptics(rng, 100000, 90, 30, 1, 1.58, 1.0, 300, 0.95, 2000, 2.5, 0.9, 0.2, true, "toyOptics.root");

// Geometry 2
ToyScintOptics(rng, 100000, 20, 20, 1,
            1.58, 1.0, 300, 0.95, 2000,
            2.5, 0.9, 0.2,
            true, "toyOptics.root");
TFile f("toyOptics.root");
TTree* tPhot = (TTree*)f.Get("tPhot");
tPhot->Draw("W");
tPhot->Draw("L");
*/

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH3D.h"

#include <vector>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>

struct Vec3
{
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3 &o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3 &o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(double a) const { return Vec3(a * x, a * y, a * z); }
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
static void DrawOnePhotonPath(
    const std::vector<double> &xs,
    const std::vector<double> &ys,
    const std::vector<double> &zs,
    double L, double W, double T,
    const char *outPng = "one_photon_path.png")
{
    if (xs.size() < 2)
        return;

    auto c = new TCanvas("cPhotonPath", "Photon path", 900, 700);

    // Frame only for axes, title, and ranges
    auto frame = new TH3D(
        "framePhotonPath",
        "One photon path; x (cm); y (cm); z (cm)",
        10, 0.0, L,
        10, -W / 2.0, W / 2.0,
        10, -T / 2.0, T / 2.0);
    frame->SetStats(0);
    frame->Draw();

    auto pl = new TPolyLine3D((int)xs.size());
    for (int i = 0; i < (int)xs.size(); i++)
        pl->SetPoint(i, xs[i], ys[i], zs[i]);
    pl->Draw("same"); // draw on top of the frame

    c->Modified();
    c->Update();
    c->SaveAs(outPng);
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

    Vec3 dirLocal(
        std::cos(phi) * sinTheta,
        std::sin(phi) * sinTheta,
        cosTheta);

    Vec3 dir = u * dirLocal.x + v * dirLocal.y + w * dirLocal.z;
    return Unit(dir);
}

static inline double epsCoupleExp(double x, double y, double rPMT, double eps0, double lambdaC)
{
    double d = sqrt(x * x + y * y);
    if (d <= rPMT)
        return 1;
    else
        return eps0 * std::exp(-d / lambdaC);
}

static inline bool InsidePMTCircle(double y, double z, double rPMT)
{
    return (y * y + z * z) <= (rPMT * rPMT);
}

double ToyScintOptics(
    TRandom3 &rng,
    long long NperPoint = 200000,
    double L = 90.0,
    double W = 30.0,
    double T = 1.0,
    int wrap = 1,
    double nScint = 1.58,
    double nOut = 1.0,
    double absLen = 300.0,
    double Rwrap = 0.95,
    int maxSteps = 2000,
    double rPMT = 0.1,
    double epsCouple = 0.90,
    double pde = 0.20,
    bool writeTree = false,
    const char *outFile = "toyOptics.root")
{
    double sinThetaC = nOut / nScint;
    if (sinThetaC >= 1.0)
        sinThetaC = 0.999999;
    double thetaC = std::asin(sinThetaC);
    double cosThetaC = std::cos(thetaC);


    double bL = L, bW = W, bT = T;
    double x0 = 0, y0 = 0, z0 = 0;
    double xf = 0, yf = 0, zf = 0;
    int endPlane = -1;
    int pmt_side = 0;
    int inPMT = 0;
    int detected = 0;
    int absorbed = 0;
    int escaped = 0;
    int nBounces = 0;
    int reachedEnd = 0;
    double pathTot = 0;
    double epsRaysPMT = 0;
    int PTFE = 0;  // wrap = 1
    int mylar = 0; // wrap = 2
    //  std::vector<float> xPath, yPath, zPath;

    // double geoTagL = bL, geoTagW = bW, geoTagT = bT;
    std::vector<double> pathX, pathY, pathZ;
    bool recordThisPhoton = false;
    bool drawThisPhoton = false;
    TFile* fout = new TFile(outFile, "RECREATE");
    TTree* tPhot = new TTree("tPhot", "Toy optical photon transport");

    printf("Creating branches \n");

    tPhot->Branch("L", &bL, "L/D");
    tPhot->Branch("W", &bW, "W/D");
    tPhot->Branch("T", &bT, "T/D");

    tPhot->Branch("wrap", &wrap, "wrap/I");

    tPhot->Branch("x0", &x0, "x0/D");
    tPhot->Branch("y0", &y0, "y0/D");
    tPhot->Branch("z0", &z0, "z0/D");

    tPhot->Branch("xf", &xf, "xf/D");
    tPhot->Branch("yf", &yf, "yf/D");
    tPhot->Branch("zf", &zf, "zf/D");

    tPhot->Branch("epsRaysPMT", &epsRaysPMT, "epsRaysPMT/D");

    tPhot->Branch("endPlane", &endPlane, "endPlane/I");
    tPhot->Branch("pmt_side", &pmt_side, "pmt_side/I");
    tPhot->Branch("inPMT", &inPMT, "inPMT/I");
    tPhot->Branch("detected", &detected, "detected/I");
    tPhot->Branch("absorbed", &absorbed, "absorbed/I");
    tPhot->Branch("escaped", &escaped, "escaped/I");
    tPhot->Branch("reachedEnd", &reachedEnd, "reached_end/I");

    tPhot->Branch("nBounces", &nBounces, "nBounces/I");
    tPhot->Branch("path", &pathTot, "path/D");
    tPhot->Branch("xPath", &pathX);
    tPhot->Branch("yPath", &pathY);
    tPhot->Branch("zPath", &pathZ);

    if (wrap == 1)
    {
        PTFE = 1;
        mylar = 0;
    }
    else if (wrap == 2)
    {
        PTFE = 0;
        mylar = 1;
    }

    long long nDet0 = 0;
    long long nDetL = 0;
    long long nAbsorb = 0;
    long long nEscape = 0;
    long long nReachedEnds = 0;
    for (long long i = 0; i < NperPoint; i++)
    {
        recordThisPhoton = 1; //(i == 0); // record the first photon only
        if (recordThisPhoton)
        {
            pathX.clear();
            pathY.clear();
            pathZ.clear();
        }

        Vec3 pos(
            rng.Uniform(0.0, L),
            rng.Uniform(-W / 2.0, W / 2.0),
            rng.Uniform(-T / 2.0, T / 2.0));

        if (recordThisPhoton)
        {
            pathX.push_back(pos.x);
            pathY.push_back(pos.y);
            pathZ.push_back(pos.z);
        }
        // xPath.push_back((float)pos.x);
        // yPath.push_back((float)pos.y);
        // zPath.push_back((float)pos.z);

        x0 = pos.x;
        y0 = pos.y;
        z0 = pos.z;

        Vec3 dir = Unit(SampleIsotropic(rng));

        int bounces = 0;
        double path = 0.0;

        endPlane = -1;
        pmt_side = 0;
        inPMT = 0;
        detected = 0;
        absorbed = 0;
        escaped = 0;
        reachedEnd = 0;
        xf = pos.x;
        yf = pos.y;
        zf = pos.z;
        for (int step = 0; step < maxSteps; step++)
        {

            double tx = 1e99;
            int xPlane = -1;
            if (dir.x > 1e-15)
            {
                tx = (L - pos.x) / dir.x;
                xPlane = 1;
            }
            if (dir.x < -1e-15)
            {
                tx = (0.0 - pos.x) / dir.x;
                xPlane = 0;
            }

            double ty1 = 1e99, ty2 = 1e99;
            if (dir.y > 1e-15)
                ty1 = (W / 2.0 - pos.y) / dir.y;
            if (dir.y < -1e-15)
                ty2 = (-W / 2.0 - pos.y) / dir.y;

            double tz1 = 1e99, tz2 = 1e99;
            if (dir.z > 1e-15)
                tz1 = (T / 2.0 - pos.z) / dir.z;
            if (dir.z < -1e-15)
                tz2 = (-T / 2.0 - pos.z) / dir.z;

            double tmin = 1e99;
            int plane = -1;

            if (tx > 1e-12 && tx < tmin)
            {
                tmin = tx;
                plane = xPlane;
            }
            if (ty2 > 1e-12 && ty2 < tmin)
            {
                tmin = ty2;
                plane = 2;
            }
            if (ty1 > 1e-12 && ty1 < tmin)
            {
                tmin = ty1;
                plane = 3;
            }
            if (tz2 > 1e-12 && tz2 < tmin)
            {
                tmin = tz2;
                plane = 4;
            }
            if (tz1 > 1e-12 && tz1 < tmin)
            {
                tmin = tz1;
                plane = 5;
            }

            if (plane < 0 || tmin > 1e98)
                break;

            if (absLen > 0)
            {
                double survive = std::exp(-tmin / absLen);
                if (rng.Uniform() > survive)
                {
                    path += tmin;
                    nAbsorb++;
                    absorbed = 1;

                    xf = pos.x + dir.x * tmin;
                    yf = pos.y + dir.y * tmin;
                    zf = pos.z + dir.z * tmin;

                    endPlane = plane;
                    break;
                }
            }

            pos = pos + dir * tmin;
            path += tmin;
            if (recordThisPhoton)
            {
                pathX.push_back(pos.x);
                pathY.push_back(pos.y);
                pathZ.push_back(pos.z);
            }
            //  xPath.push_back((float)pos.x);
            //   yPath.push_back((float)pos.y);
            //   zPath.push_back((float)pos.z);

            xf = pos.x;
            yf = pos.y;
            zf = pos.z;
            endPlane = plane;
            double eps0 = 0.90;
            double lambdaC = 120.0;

            if (plane == 0)
            {
                reachedEnd++;
                pmt_side = 1;
                inPMT = InsidePMTCircle(pos.y, pos.z, rPMT) ? 1 : 0;
                epsRaysPMT = epsCoupleExp(pos.y, pos.z, rPMT, eps0, lambdaC);
                double pEnd = epsCouple * pde * epsRaysPMT;
                if (rng.Uniform() < pEnd)
                {
                    detected = 1;
                    nDet0++;
                }

                if (recordThisPhoton)
                {
                    pathX.push_back(xf);
                    pathY.push_back(yf);
                    pathZ.push_back(zf);
                }

                break;
            }
            if (plane == 1)
            {
                reachedEnd++;
                pmt_side = 2;
                inPMT = InsidePMTCircle(pos.y, pos.z, rPMT) ? 1 : 0;
                epsRaysPMT = epsCoupleExp(pos.y, pos.z, rPMT, eps0, lambdaC);
                double pEnd = epsCouple * pde * epsRaysPMT;
                if (rng.Uniform() < pEnd)
                {
                    detected = 1;
                    nDet0++;
                }
                if (recordThisPhoton)
                {
                    pathX.push_back(xf);
                    pathY.push_back(yf);
                    pathZ.push_back(zf);
                }

                break;
            }

            Vec3 nHat;
            if (plane == 2)
                nHat = Vec3(0, -1, 0);
            if (plane == 3)
                nHat = Vec3(0, +1, 0);
            if (plane == 4)
                nHat = Vec3(0, 0, -1);
            if (plane == 5)
                nHat = Vec3(0, 0, +1);
            nHat = Unit(nHat);

            double cosInc = std::fabs(Dot(dir, nHat));
            bool isTIR = (cosInc < cosThetaC);

            if (isTIR)
            {
                dir = Unit(ReflectSpecular(dir, nHat));
                bounces++;
            }
            else
            {
                if (rng.Uniform() < Rwrap)
                {
                    if (PTFE && !mylar)
                        dir = SampleLambert(nHat * (-1.0), rng); // Tyvek/PTFE
                    else if (mylar && !PTFE)
                        dir = Unit(ReflectSpecular(dir, nHat)); // Mylar
                    else
                    {
                        printf("no input regarding mylar/tyvek assumin specular reflection! \n");
                        dir = Unit(ReflectSpecular(dir, nHat)); // Mylar
                    }
                    bounces++;
                }
                else
                {
                    nEscape++;
                    escaped = 1;
                    break;
                }
            }
        }

        nBounces = bounces;
        pathTot = path;

            tPhot->Fill();
        if (drawThisPhoton)
        {
            DrawOnePhotonPath(pathX, pathY, pathZ, L, W, T, "one_photon_path.png");
            // only do it once
            drawThisPhoton = false;
        }
    }

    double fTot = double(nDet0 + nDetL) / double(NperPoint);


        fout->cd();
        tPhot->Write();
        fout->Close();
        delete fout;


    std::cout << "ToyScintOptics: detected fraction = " << fTot
              << " (nDet0=" << nDet0 << ", nDetL=" << nDetL
              << ", absorbed=" << nAbsorb << ", escaped=" << nEscape << ")\n";

    return fTot;
}

void ScanToyOptics(
    TH2D *&hFrac,
    TH2D *&hFracCost1,
    TH2D *&hFracCost2,
    long long NperPoint = 50000,
    double T = 1.0,
    double nScint = 1.58,
    double nOut = 1.0,
    double absLen = 300.0,
    double Rwrap = 0.95,
    int maxSteps = 2000)
{
    TRandom3 rng(0);

    const double cost_pmt = 1269.0;
    const double cost_cm2 = 0.4;

    const double minCM = 30.0;
    const double maxCM = 200.0;
    const double stepCM = 10.0;

    const int nBins = int((maxCM - minCM) / stepCM) + 1;

    hFrac = new TH2D(
        "hFrac",
        "Toy detected fraction;Width W (cm);Length L (cm)",
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0,
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0);

    hFracCost1 = new TH2D(
        "hFracCost1",
        "cost1/(fraction*W*L) with 1 PMT;Width W (cm);Length L (cm)",
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0,
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0);

    hFracCost2 = new TH2D(
        "hFracCost2",
        "cost2/(fraction*W*L) with 2 PMTs;Width W (cm);Length L (cm)",
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0,
        nBins, minCM - stepCM / 2.0, maxCM + stepCM / 2.0);
    int wrap = 1;
    for (int ix = 1; ix <= nBins; ix++)
    {
        double W = minCM + (ix - 1) * stepCM;

        for (int iy = 1; iy <= nBins; iy++)
        {
            double L = minCM + (iy - 1) * stepCM;

            double cost1 = W * L * cost_cm2 + cost_pmt;
            double cost2 = W * L * cost_cm2 + 2.0 * cost_pmt;

            double frac = ToyScintOptics(rng, NperPoint, L, W, T, wrap, nScint, nOut, absLen, Rwrap, maxSteps,
                                         2.5, 0.90, 0.20, false, "toyOptics.root");

            hFrac->SetBinContent(ix, iy, frac);

            double denom1 = (frac > 0) ? ((frac / 2.0) * W * L) : 1e-12;
            double denom2 = (frac > 0) ? (frac * W * L) : 1e-12;

            hFracCost1->SetBinContent(ix, iy, cost1 / denom1);
            hFracCost2->SetBinContent(ix, iy, cost2 / denom2);

            std::cout << "W=" << W << " cm, L=" << L << " cm -> frac=" << frac << "\n";
        }
    }
}