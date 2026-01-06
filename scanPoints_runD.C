// scanPoints.C
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStopwatch.h"

#include "Geometry.h"
//#include "Vec3.h"
#include "OpticsConfig.h"
#include "TreeWriter.h"
#include "Transport.h"
#include "DrawHelper.h"
#include <vector> // make sure this exists
// ...

#include "Math/Vector3D.h"

using V3 = ROOT::Math::XYZVector;
V3 SampleIsotropic(TRandom3 &rng);

// TODO: add info if was hit left or right. Run again without saving path information 22/12/2025
// ------------------------------
// Helper: count outcomes inside a per-point ROOT file
struct Counts
{
    Long64_t total = 0;
    Long64_t absorbed = 0;
    Long64_t escaped = 0;
    Long64_t hitPMT = 0; // uses inPMT (or hitPMT if you have that)
    Long64_t detected = 0;
};

// ------------------------------
// Helper: mean and RMS of "path" for detected photons per PMT side
struct PathStats
{
    Long64_t n = 0;
    double mean = 0.0;
    double rms = 0.0; // width (standard deviation)
};

static Long64_t CountDetectedSideFromFile(const char *fn, int pmtSide, const char *tn = "tPhot")
{
    TFile f(fn, "READ");
    if (f.IsZombie())
    {
        std::cout << "CountDetectedSideFromFile: cannot open " << fn << "\n";
        return 0;
    }
    auto t = (TTree *)f.Get(tn);
    if (!t)
    {
        std::cout << "CountDetectedSideFromFile: cannot find tree " << tn << " in " << fn << "\n";
        return 0;
    }
    if (!t->GetBranch("detected") || !t->GetBranch("pmt_side"))
    {
        std::cout << "CountDetectedSideFromFile: missing branches in " << fn
                  << " (need detected, pmt_side)\n";
        return 0;
    }

    TString sel = Form("detected==1 && pmt_side==%d", pmtSide);
    return t->GetEntries(sel);
}

static PathStats PathStatsFromFile(const char *fn, int pmtSide, const char *tn = "tPhot")
{
    PathStats s;

    TFile f(fn, "READ");
    if (f.IsZombie())
    {
        std::cout << "PathStatsFromFile: cannot open " << fn << "\n";
        return s;
    }

    auto t = (TTree *)f.Get(tn);
    if (!t)
    {
        std::cout << "PathStatsFromFile: cannot find tree " << tn << " in " << fn << "\n";
        return s;
    }

    // Require needed branches
    if (!t->GetBranch("path") || !t->GetBranch("detected") || !t->GetBranch("pmt_side"))
    {
        std::cout << "PathStatsFromFile: missing branches in " << fn
                  << " (need path, detected, pmt_side)\n";
        return s;
    }

    double path = 0.0;
    int detected = 0;
    int side = 0;

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("path", 1);
    t->SetBranchStatus("detected", 1);
    t->SetBranchStatus("pmt_side", 1);

    t->SetBranchAddress("path", &path);
    t->SetBranchAddress("detected", &detected);
    t->SetBranchAddress("pmt_side", &side);

    long double sum = 0.0L;
    long double sum2 = 0.0L;

    const Long64_t nent = t->GetEntries();
    for (Long64_t i = 0; i < nent; i++)
    {
        t->GetEntry(i);
        if (detected != 1)
            continue;
        if (side != pmtSide)
            continue;

        s.n++;
        sum += (long double)path;
        sum2 += (long double)path * (long double)path;
    }

    if (s.n > 0)
    {
        const long double mean = sum / (long double)s.n;
        const long double var = sum2 / (long double)s.n - mean * mean;
        s.mean = (double)mean;
        s.rms = (var > 0.0L) ? std::sqrt((double)var) : 0.0;
    }

    return s;
}

static Counts CountFromFile(const char *fn, const char *tn = "tPhot")
{
    Counts c;

    TFile f(fn, "READ");
    if (f.IsZombie())
    {
        std::cout << "CountFromFile: cannot open " << fn << "\n";
        return c;
    }

    auto t = (TTree *)f.Get(tn);
    if (!t)
    {
        std::cout << "CountFromFile: cannot find tree " << tn << " in " << fn << "\n";
        return c;
    }

    auto has = [&](const char *b)
    { return t->GetBranch(b) != nullptr; };

    const char *bAbs = has("absorbed") ? "absorbed" : nullptr;
    const char *bEsc = has("escaped") ? "escaped" : nullptr;
    const char *bHit = has("inPMT") ? "inPMT" : (has("hitPMT") ? "hitPMT" : nullptr);
    const char *bDet = has("detected") ? "detected" : nullptr;

    auto countOnes = [&](const char *bname) -> Long64_t
    {
        if (!bname)
            return 0;
        TString sel = Form("%s!=0", bname);
        return t->GetEntries(sel);
    };

    c.total = t->GetEntries();
    c.absorbed = countOnes(bAbs);
    c.escaped = countOnes(bEsc);
    c.hitPMT = countOnes(bHit);
    c.detected = countOnes(bDet);
    return c;
}

auto mirrorBin = [](int b, int nb)
{ return nb - b + 1; };

static inline V3 SampleIsotropic(TRandom3 &rng)
{
    double u = rng.Uniform(-1.0, 1.0);
    double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    double s = std::sqrt(std::max(0.0, 1.0 - u * u));
    return V3(s * std::cos(phi), s * std::sin(phi), u);
}

static inline V3 SampleInScint(TRandom3 &rng, double L, double W, double T, double wedgeLen, double wedgeTipW)
{

    double x = rng.Uniform(-L * 0.5, +L * 0.5);
    double y = rng.Uniform(-W * 0.5, +W * 0.5);
    double z = rng.Uniform(-T * 0.5, +T * 0.5);
    return V3(x, y, z);
}

// ------------------------------
// Main scan
void scanPoints_runD()
{
    // If you run this as a ROOT macro, you can still compile-load your project files like this
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
    gROOT->ProcessLine(".L DrawHelper.cxx+");

    OpticsConfig cfg;
    cfg.savePath = false;
    cfg.mirrorx = true;
    cfg.mirrory = true;
    cfg.dScan = 2.0; // cm

    // geometry
    cfg.L = 60.0; // was 90
    cfg.W = 30.0; // was 30
    cfg.T = 1.0;
    cfg.wrap = 1; // 1=PTFE, 2=Mylar
    const int Nphot = 100000;

    // optics
    cfg.nScint = 1.58;
    cfg.nOut = 1.0;
    cfg.absLen = 300.0;
    cfg.Rwrap = 0.95;
    cfg.maxSteps = 2000;

    // PMT / coupling
    //    cfg.rPMT = 2.5;
    cfg.rPMT = 1.27; // 1 inch diameter

    cfg.epsCouple = 0.90;
    cfg.pde = 0.20;
    cfg.eps0 = 0.0;
    cfg.lambdaC = 120.0;

    // wedges
    cfg.useWedge = true;
    cfg.wedgeLen = 20;                // 20.0;
    cfg.wedgeTipW = 10;               // 5.0;
    cfg.L = cfg.L + 2 * cfg.wedgeLen; // was 90

    double x0, x1, y0, y1;
    double epsilon = 0.01;

    if (cfg.useWedge)
    {
        x0 = cfg.wedgeLen + epsilon;
        x1 = cfg.L - cfg.wedgeLen - epsilon;
    }
    else
    {
        x0 = 0.0 + epsilon;
        x1 = cfg.L - epsilon;
    }

    TGraph2D *g2 = nullptr;
    std::vector<V3> normals, pointPlane;
    cout << "Making faces...\n";
    makeFaces(cfg, g2, normals, pointPlane);
    cout << "Making faces doe...\n";

    g2->SetMarkerStyle(20);

    g2->Draw("AP LINE");
    TRandom3 rng(0);
    TGraph2D *ghits=new TGraph2D();
    for (int j = 0; j < 1; ++j)
    {
        V3 p0 = SampleInScint(rng, cfg.L, cfg.W, cfg.T, cfg.wedgeLen, cfg.wedgeTipW);
        V3 dir = SampleIsotropic(rng);
        //  printf("Starting point (%f, %f, %f), direction=(%f, %f, %f)\n", p0.x(), p0.y(), p0.z(), dir.x(), dir.y(), dir.z());
        double tout;

        V3 hitPoint;

        int n = 0;
        int pre = 0;
        V3 prehitpoint;
        ghits->SetPoint(ghits->GetN(), p0.x(), p0.y(), p0.z());
        int r = getIntersectionPlan(normals,
                                    pointPlane,
                                    p0,
                                    dir,
                                    tout,
                                    hitPoint);
        cout << r << endl;
    }

    if (false)
    {
        y0 = -cfg.W * 0.5 + epsilon;
        y1 = +cfg.W * 0.5 - epsilon;

        const int NstepsX_full = int((x1 - x0) / cfg.dScan); // scan every dscan cm . Only generate events in the active volume "the box part, not the wedge"
        const int NstepsY_full = int((cfg.W) / cfg.dScan);   // scan every dscan cm
        const int NstepsX_half = (NstepsX_full + 1) / 2;     // include middle if odd
        const int NstepsY_half = (NstepsY_full + 1) / 2;     // include y=0 if odd
        int NstepsX = NstepsX_full;
        int NstepsY = NstepsY_full;

        if (cfg.mirrorx)
            NstepsX = NstepsX_half;
        if (cfg.mirrory)
            NstepsY = NstepsY_half;

        TString name = Form("out/scint_%d_%d_%d_", int(cfg.L), int(cfg.W), int(cfg.T));
        if (cfg.useWedge)
        {
            name += Form("wedge_%d_%d_", int(cfg.wedgeLen), int(cfg.wedgeTipW));
        }
        else
        {
            name += "nowedge_";
        }
        if (cfg.wrap == 1)
        {
            name += "PTFE_";
        }
        else if (cfg.wrap == 2)
        {
            name += "Mylar_";
        }
        else
        {
            name += "nowrap_";
        }
        name += Form("N_%d", Nphot);
        // One summary file with one tree that collects results for all (x,y)
        TFile fsum(name + "_summary.root", "RECREATE");
        TTree tsum("tScan", "Per-point scan summary");

        double x = 0, y = 0, z = 0;
        int ix = 0, iy = 0;
        Long64_t total = 0, absorbed = 0, escaped = 0, hitPMT = 0, detected = 0;
        double frac_absorbed = 0, frac_escaped = 0, frac_hitPMT = 0, frac_detected = 0;
        char outFileName[512];

        tsum.Branch("ix", &ix, "ix/I");
        tsum.Branch("iy", &iy, "iy/I");
        tsum.Branch("x", &x, "x/D");
        tsum.Branch("y", &y, "y/D");
        tsum.Branch("z", &z, "z/D");
        tsum.Branch("outFile", outFileName, "outFile/C");

        tsum.Branch("N", &total, "N/L");
        tsum.Branch("absorbed", &absorbed, "absorbed/L");
        tsum.Branch("escaped", &escaped, "escaped/L");
        tsum.Branch("hitPMT", &hitPMT, "hitPMT/L");
        tsum.Branch("detected", &detected, "detected/L");

        tsum.Branch("frac_absorbed", &frac_absorbed, "frac_absorbed/D");
        tsum.Branch("frac_escaped", &frac_escaped, "frac_escaped/D");
        tsum.Branch("frac_hitPMT", &frac_hitPMT, "frac_hitPMT/D");
        tsum.Branch("frac_detected", &frac_detected, "frac_detected/D");
        TH2D *h_absorbed = new TH2D("h_absorbed", "frac absorbed; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_escaped = (TH2D *)h_absorbed->Clone("h_escaped");
        h_escaped->SetTitle("frac escaped");
        TH2D *h_detected = (TH2D *)h_absorbed->Clone("h_detected");
        h_escaped->SetTitle("frac detected");
        TH2D *h_hitPMT = (TH2D *)h_absorbed->Clone("h_hitPMT");
        h_escaped->SetTitle("frac hit PMT");
        // Per-point path stats for detected photons by PMT side
        TH2D *h_pathMean_det_p2 = new TH2D("h_pathMean_det_p2", "mean(path) detected, pmt_side=2; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_pathMean_det_p1 = new TH2D("h_pathMean_det_p1", "mean(path) detected, pmt_side=1; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_pathRMS_det_p2 = new TH2D("h_pathRMS_det_p2", "RMS(path) detected, pmt_side=2; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_pathRMS_det_p1 = new TH2D("h_pathRMS_det_p1", "RMS(path) detected, pmt_side=1; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_detected_p1 = new TH2D("h_detected_p1", "N detected, pmt_side=1; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        TH2D *h_detected_p2 = new TH2D("h_detected_p2", "N detected, pmt_side=2; x; y", NstepsX, x0, x1, NstepsY, y0, y1);

        double a0, a1, b0, b1;
        int N;
        if (x1 - x0 > y1 - y0)
        {
            double diff = (x1 - x0) - (y1 - y0);
            a0 = x0;
            a1 = x1;
            b0 = y0 - diff / 2;
            b1 = y1 + diff / 2;
        }
        else
        {
            double diff = (y1 - y0) - (x1 - x0);
            a0 = x0 - diff / 2;
            a1 = x1 + diff / 2;
            b0 = y0;
            b1 = y1;
        }
        TH2D *h_absorbed_s = new TH2D("h_absorbed_s", "frac absorbed scale; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
        h_absorbed_s->GetXaxis()->SetRangeUser(a0, a1);
        h_absorbed_s->GetYaxis()->SetRangeUser(b0, b1);
        TH2D *h_escaped_s = (TH2D *)h_absorbed_s->Clone("h_escaped_s");
        h_escaped_s->SetTitle("frac escaped scale");
        TH2D *h_detected_s = (TH2D *)h_absorbed_s->Clone("h_detected_s");
        h_escaped_s->SetTitle("frac detected sclae");
        TH2D *h_hitPMT_s = (TH2D *)h_absorbed_s->Clone("h_hitPMT_s");
        h_escaped_s->SetTitle("frac hit PMT scale");
        TRandom3 rng(0);
        bool saveTree = false;
        int ic = 0;
        TCanvas *c1 = new TCanvas();
        TStopwatch swAll;
        swAll.Start();
        TStopwatch swRow;

        for (ix = 0; ix < NstepsX; ix++)
        {
            swRow.Start(kTRUE); // reset + start

            x = x0 + (x1 - x0) * ix / (NstepsX - 1.0);

            for (iy = 0; iy < NstepsY; iy++)
            {
                rng.SetSeed(12345 + ix * 100000 + iy); // optional

                y = y0 + (y1 - y0) * iy / (NstepsY - 1.0);
                z = 0.0;

                TString outFile = Form(name + "_%03d_%03d.root", ix, iy);
                snprintf(outFileName, sizeof(outFileName), "%s", outFile.Data());
                // Make the per-point file

                // RAII: TreeWriter is destroyed each loop, so no leaks as long as ~TreeWriter is defined
                // TreeWriter wr(outFile.Data(), cfg);
                int Ndetected = 0;
                int Nescaped = 0;
                int Nabsorbed = 0;
                int NinPMT = 0;
                int Nreached = 0;

                /*   for (int i = 0; i < Nphot; i++)
                   {
                       PhotonResult res = PropagateOnePhoton(rng, site, 0, cfg);
                       wr.Fill(res);
                       Ndetected = Ndetected + res.detected;
                       Nescaped = Nescaped +res.escaped;
                       Nabsorbed= Nabsorbed+res.absorbed;
                       NinPMT=NinPMT+res.inPMT;
                       Nreached = Nreached +res.reachedEnd;
                   }
                   wr.Close();
    */

                // accumulators per point
                Long64_t Ntot = 0, Nabs = 0, Nesc = 0, NhPMT = 0, Ndet = 0;
                Long64_t Ndet_p1 = 0, Ndet_p2 = 0;
                long double sum_p1 = 0.0L, sum2_p1 = 0.0L;
                long double sum_p2 = 0.0L, sum2_p2 = 0.0L;

                TRandom3 rng(0);
                V3 site(x, y, z);
                ic++;
                if (ic % int(NstepsX * NstepsY / 5) == 0)
                {
                    saveTree = true;
                }
                else
                {
                    saveTree = false;
                }
                TreeWriter *wr = nullptr;

                if (saveTree)
                    wr = new TreeWriter(outFile.Data(), cfg);

                for (int i = 0; i < Nphot; i++)
                {
                    PhotonResult res = PropagateOnePhoton(rng, site, 0, cfg);
                    if (saveTree)
                        wr->Fill(res);
                    Ntot++;
                    Nabs += (res.absorbed != 0);
                    Nesc += (res.escaped != 0);
                    NhPMT += (res.inPMT != 0);
                    Ndet += (res.detected != 0);

                    if (res.detected)
                    {
                        if (res.pmt_side == 1)
                        {
                            Ndet_p1++;
                            sum_p1 += (long double)res.path;
                            sum2_p1 += (long double)res.path * (long double)res.path;
                        }
                        else if (res.pmt_side == 2)
                        {
                            Ndet_p2++;
                            sum_p2 += (long double)res.path;
                            sum2_p2 += (long double)res.path * (long double)res.path;
                        }
                    }
                }
                if (saveTree)
                    wr->Close();
                // means + rms
                double mean_p1 = 0.0, rms_p1 = 0.0;
                double mean_p2 = 0.0, rms_p2 = 0.0;

                if (Ndet_p1 > 0)
                {
                    long double m = sum_p1 / (long double)Ndet_p1;
                    long double v = sum2_p1 / (long double)Ndet_p1 - m * m;
                    mean_p1 = (double)m;
                    rms_p1 = (v > 0.0L) ? std::sqrt((double)v) : 0.0;
                }
                if (Ndet_p2 > 0)
                {
                    long double m = sum_p2 / (long double)Ndet_p2;
                    long double v = sum2_p2 / (long double)Ndet_p2 - m * m;
                    mean_p2 = (double)m;
                    rms_p2 = (v > 0.0L) ? std::sqrt((double)v) : 0.0;
                }

                total = Ntot;
                absorbed = Nabs;
                escaped = Nesc;
                hitPMT = NhPMT;
                detected = Ndet;

                frac_absorbed = (Ntot > 0) ? double(Nabs) / double(Ntot) : 0.0;
                frac_escaped = (Ntot > 0) ? double(Nesc) / double(Ntot) : 0.0;
                frac_hitPMT = (Ntot > 0) ? double(NhPMT) / double(Ntot) : 0.0;
                frac_detected = (Ntot > 0) ? double(Ndet) / double(Ntot) : 0.0;
                cout << "Point (" << ix << "," << iy << ") : detected % = " << frac_detected * 100 << "\t";
                cout << ": escaped % = " << frac_escaped * 100 << "\t";
                cout << ": absorbed % = " << frac_absorbed * 100 << "\t";
                cout << ": hitpmt % = " << frac_hitPMT * 100 << endl;
                tsum.Fill();
                int bx = ix + 1;
                int by = iy + 1;
                h_pathMean_det_p1->SetBinContent(bx, by, mean_p1);
                h_pathRMS_det_p1->SetBinContent(bx, by, rms_p1);
                h_pathMean_det_p2->SetBinContent(bx, by, mean_p2);
                h_pathRMS_det_p2->SetBinContent(bx, by, rms_p2);
                h_detected->SetBinContent(bx, by, Ndetected);
                h_detected_s->SetBinContent(bx, by, frac_detected);
                h_absorbed->SetBinContent(bx, by, Nabsorbed);
                h_absorbed_s->SetBinContent(bx, by, frac_absorbed);
                h_escaped->SetBinContent(bx, by, Nescaped);
                h_escaped_s->SetBinContent(bx, by, frac_escaped);
                h_hitPMT->SetBinContent(bx, by, NinPMT);
                h_hitPMT_s->SetBinContent(bx, by, frac_hitPMT);
                if (cfg.mirrory)
                {
                    int byM = NstepsY - by + 1;
                    h_pathMean_det_p1->SetBinContent(bx, byM, mean_p1);
                    h_pathRMS_det_p1->SetBinContent(bx, byM, rms_p1);
                    h_pathMean_det_p2->SetBinContent(bx, byM, mean_p2);
                    h_pathRMS_det_p2->SetBinContent(bx, byM, rms_p2);
                    h_detected->SetBinContent(bx, byM, Ndetected);
                    h_detected_s->SetBinContent(bx, byM, frac_detected);
                    h_absorbed->SetBinContent(bx, byM, Nabsorbed);
                    h_absorbed_s->SetBinContent(bx, byM, frac_absorbed);
                    h_escaped->SetBinContent(bx, byM, Nescaped);
                    h_escaped_s->SetBinContent(bx, byM, frac_escaped);
                    h_hitPMT->SetBinContent(bx, byM, NinPMT);
                    h_hitPMT_s->SetBinContent(bx, byM, frac_hitPMT);
                }
                if (cfg.mirrorx)
                {
                    int bxM = NstepsX - bx + 1;
                    h_pathMean_det_p1->SetBinContent(bxM, by, mean_p2);
                    h_pathRMS_det_p1->SetBinContent(bxM, by, rms_p2);
                    h_pathMean_det_p2->SetBinContent(bxM, by, mean_p1);
                    h_pathRMS_det_p2->SetBinContent(bxM, by, rms_p1);
                    h_detected->SetBinContent(bxM, by, Ndetected);
                    h_detected_s->SetBinContent(bxM, by, frac_detected);
                    h_absorbed->SetBinContent(bxM, by, Nabsorbed);
                    h_absorbed_s->SetBinContent(bxM, by, frac_absorbed);
                    h_escaped->SetBinContent(bxM, by, Nescaped);
                    h_escaped_s->SetBinContent(bxM, by, frac_escaped);
                    h_hitPMT->SetBinContent(bxM, by, NinPMT);
                    h_hitPMT_s->SetBinContent(bxM, by, frac_hitPMT);
                }
                if (cfg.mirrorx && cfg.mirrory)
                {
                    int bxM = NstepsX - bx + 1;
                    int byM = NstepsY - by + 1;
                    h_pathMean_det_p1->SetBinContent(bxM, byM, mean_p2);
                    h_pathRMS_det_p1->SetBinContent(bxM, byM, rms_p2);
                    h_pathMean_det_p2->SetBinContent(bxM, byM, mean_p1);
                    h_pathRMS_det_p2->SetBinContent(bxM, byM, rms_p1);
                    h_detected->SetBinContent(bxM, byM, Ndetected);
                    h_detected_s->SetBinContent(bxM, byM, frac_detected);
                    h_absorbed->SetBinContent(bxM, byM, Nabsorbed);
                    h_absorbed_s->SetBinContent(bxM, byM, frac_absorbed);
                    h_escaped->SetBinContent(bxM, byM, Nescaped);
                    h_escaped_s->SetBinContent(bxM, byM, frac_escaped);
                    h_hitPMT->SetBinContent(bxM, byM, NinPMT);
                    h_hitPMT_s->SetBinContent(bxM, by, frac_hitPMT);
                }

                bool verbose = false;
                if (verbose)
                {
                    printf("%d/%d %d/%d\t", ix, iy, NstepsX, NstepsY);
                    printf("frac detected=%f, absorbed=%f, mean_p1=%f, rms_p1=%f, mean_p2=%f, rms_p2=%f, tot=%f\n",
                           frac_detected, frac_absorbed,
                           mean_p1, rms_p1,
                           mean_p2, rms_p2,
                           frac_detected + frac_absorbed + frac_escaped + frac_hitPMT);
                }
                h_detected_p1->SetBinContent(ix + 1, iy + 1, (double)Ndet_p1);
                h_detected_p2->SetBinContent(ix + 1, iy + 1, (double)Ndet_p2);

                // Optional: draw one event from this point
                // DrawEventSplitViewFromTree(outFile.Data(), 12, true, 0.12);
            }

            swRow.Stop();
            double secPerRow = swRow.RealTime();
            double rowsLeft = (NstepsX - 1 - ix);

            std::cout << "Row ix=" << ix << "/" << NstepsX - 1 << " x" << NstepsY
                      << " took real " << swRow.RealTime() << " s"
                      << ", cpu " << swRow.CpuTime() << " s"
                      << " (" << NstepsY << " points)"
                      << "\t ETA ~ " << (secPerRow * rowsLeft / 60.0) << " min\n";
        }
        fsum.cd();
        h_detected->Write("h_detected");
        h_absorbed->Write("h_absorbed");
        h_escaped->Write("h_escaped");
        h_hitPMT->Write("h_hitPMT");
        h_detected_s->Write("h_detected_s");
        h_absorbed_s->Write("h_absorbed_s");
        h_escaped_s->Write("h_escaped_s");
        h_hitPMT_s->Write("h_hitPMT_s");
        h_pathMean_det_p2->Write("h_pathMean_det_p2");
        h_pathMean_det_p1->Write("h_pathMean_det_p1");
        h_pathRMS_det_p2->Write("h_pathRMS_det_p2");
        h_pathRMS_det_p1->Write("h_pathRMS_det_p1");
        h_detected_p1->Write("h_detected_p1");
        h_detected_p2->Write("h_detected_p2");

        fsum.cd();
        tsum.Write();
        fsum.Close();
        swAll.Stop();

        std::cout << "Total scan time: real " << swAll.RealTime()
                  << " s, cpu " << swAll.CpuTime() << " s\n";

        std::cout << "Wrote summary to " << name + "_summary.root \n";
    }
}
