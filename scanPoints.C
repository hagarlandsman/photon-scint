// scanPoints.C
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>

#include "Geometry.h"
#include "Vec3.h"
#include "OpticsConfig.h"
#include "TreeWriter.h"
#include "Transport.h"
#include "DrawHelper.h"

// TODO: add info if was hit left or right. Run again without saving path information 22/12/2025
// ------------------------------
// Helper: count outcomes inside a per-point ROOT file
struct Counts {
    Long64_t total = 0;
    Long64_t absorbed = 0;
    Long64_t escaped = 0;
    Long64_t hitPMT = 0;     // uses inPMT (or hitPMT if you have that)
    Long64_t detected = 0;
};

static Counts CountFromFile(const char* fn, const char* tn="tPhot")
{
    Counts c;

    TFile f(fn, "READ");
    if (f.IsZombie()) {
        std::cout << "CountFromFile: cannot open " << fn << "\n";
        return c;
    }

    auto t = (TTree*)f.Get(tn);
    if (!t) {
        std::cout << "CountFromFile: cannot find tree " << tn << " in " << fn << "\n";
        return c;
    }

    auto has = [&](const char* b){ return t->GetBranch(b) != nullptr; };

    const char* bAbs = has("absorbed") ? "absorbed" : nullptr;
    const char* bEsc = has("escaped")  ? "escaped"  : nullptr;
    const char* bHit = has("inPMT")    ? "inPMT"    : (has("hitPMT") ? "hitPMT" : nullptr);
    const char* bDet = has("detected") ? "detected" : nullptr;

    auto countOnes = [&](const char* bname)->Long64_t {
        if (!bname) return 0;
        TString sel = Form("%s!=0", bname);
        return t->GetEntries(sel);
    };

    c.total    = t->GetEntries();
    c.absorbed = countOnes(bAbs);
    c.escaped  = countOnes(bEsc);
    c.hitPMT   = countOnes(bHit);
    c.detected = countOnes(bDet);
    return c;
}

// ------------------------------
// Main scan
void scanPoints()
{
    // If you run this as a ROOT macro, you can still compile-load your project files like this
    gROOT->ProcessLine(".L TreeWriter.cxx+");
    gROOT->ProcessLine(".L Geometry.cxx+");
    gROOT->ProcessLine(".L Transport.cxx+");
    gROOT->ProcessLine(".L DrawHelper.cxx+");

    OpticsConfig cfg;
    cfg.savePath = true;

    // geometry
    cfg.L = 150.0;
    cfg.W = 30.0;
    cfg.T = 1.0;
    cfg.wrap = 1; // 1=PTFE, 2=Mylar

    // optics
    cfg.nScint = 1.58;
    cfg.nOut = 1.0;
    cfg.absLen = 300.0;
    cfg.Rwrap = 0.95;
    cfg.maxSteps = 2000;

    // PMT / coupling
    cfg.rPMT = 2.5;
    cfg.epsCouple = 0.90;
    cfg.pde = 0.20;
    cfg.eps0 = 0.0;
    cfg.lambdaC = 120.0;

    // wedges
    cfg.useWedge = false;
    cfg.wedgeLen = 20.0;
    cfg.wedgeTipW = 5.0;

    double x0, x1, y0, y1;
    double epsilon = 0.01;

    if (cfg.useWedge) {
        x0 = cfg.wedgeLen + epsilon;
        x1 = cfg.L - cfg.wedgeLen - epsilon;
    } else {
        x0 = 0.0 + epsilon;
        x1 = cfg.L - epsilon;
    }
    y0 = -cfg.W * 0.5 + epsilon;
    y1 = +cfg.W * 0.5 - epsilon;


    const int NstepsX = int ((cfg.L-epsilon*2) / 1); // scan every 1 cm
    const int NstepsY = int ((cfg.W-epsilon*2) / 1); // scan every 1 cm
    const int Nphot   = 10000;

    // One summary file with one tree that collects results for all (x,y)
    const char* summaryFile = "scan_summary.root";
    TFile fsum(summaryFile, "RECREATE");
    TTree tsum("tScan", "Per-point scan summary");

    double x = 0, y = 0, z = 0;
    int ix = 0, iy = 0;
    Long64_t total=0, absorbed=0, escaped=0, hitPMT=0, detected=0;
    double frac_absorbed=0, frac_escaped=0, frac_hitPMT=0, frac_detected=0;
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
    TH2D* h_absorbed = new TH2D("h_absorbed","frac absorbed; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
    TH2D* h_escaped = (TH2D*) h_absorbed->Clone("h_escaped"); h_escaped->SetTitle("frac escaped");
    TH2D* h_detected = (TH2D*) h_absorbed->Clone("h_detected"); h_escaped->SetTitle("frac detected");
    TH2D* h_hitPMT = (TH2D*) h_absorbed->Clone("h_hitPMT"); h_escaped->SetTitle("frac hit PMT");

    double a0,a1,b0,b1;
    int N;
    if (x1 - x0 > y1 - y0) {
        double diff = (x1-x0) - (y1-y0);
        a0 = x0;
        a1 = x1;
        b0 = y0 - diff / 2;
        b1 = y1 + diff / 2;
    }
    else {
        double diff = (y1-y0) - (x1-x0) ;
        a0 = x0 - diff / 2;
        a1 = x1 + diff / 2;
        b0 = y0;
        b1 = y1;
    }
    TH2D* h_absorbed_s = new TH2D("h_absorbed_s","frac absorbed scale; x; y", NstepsX, x0, x1, NstepsY, y0, y1);
    h_absorbed_s->GetXaxis()->SetRangeUser(a0,a1);
    h_absorbed_s->GetYaxis()->SetRangeUser(b0,b1);
    TH2D* h_escaped_s = (TH2D*) h_absorbed_s->Clone("h_escaped_s"); h_escaped_s->SetTitle("frac escaped scale");
    TH2D* h_detected_s = (TH2D*) h_absorbed_s->Clone("h_detected_s"); h_escaped_s->SetTitle("frac detected sclae");
    TH2D* h_hitPMT_s = (TH2D*) h_absorbed_s->Clone("h_hitPMT_s"); h_escaped_s->SetTitle("frac hit PMT scale");

    TCanvas *c1= new TCanvas();
    for (ix = 0; ix < NstepsX; ix++)
    {
        x = x0 + (x1 - x0) * ix / (NstepsX - 1.0);

        for (iy = 0; iy < NstepsY; iy++)
        {
            y = y0 + (y1 - y0) * iy / (NstepsY - 1.0);
            z = 0.0;

            TString outFile = Form("out/out_%03d_%03d.root", ix, iy);
            snprintf(outFileName, sizeof(outFileName), "%s", outFile.Data());

            // Make the per-point file
            TRandom3 rng(0);
            Vec3 site(x, y, z);


                // RAII: TreeWriter is destroyed each loop, so no leaks as long as ~TreeWriter is defined
                TreeWriter wr(outFile.Data(), cfg);

                for (int i = 0; i < Nphot; i++)
                {
                    PhotonResult res = PropagateOnePhoton(rng, site, 0, cfg);
                    wr.Fill(res);
                }
                wr.Close();


            // Count outcomes in that file
            printf ("%d/%d %d/%d \t ",ix,NstepsX,iy,NstepsY);
            Counts c = CountFromFile(outFile.Data(), "tPhot");
            total    = c.total;
            absorbed = c.absorbed;
            escaped  = c.escaped;
            hitPMT   = c.hitPMT;
            detected = c.detected;

            frac_absorbed = (total > 0) ? double(absorbed)/double(total) : 0.0;
            frac_escaped  = (total > 0) ? double(escaped )/double(total) : 0.0;
            frac_hitPMT   = (total > 0) ? double(hitPMT  )/double(total) : 0.0;
            frac_detected = (total > 0) ? double(detected)/double(total) : 0.0;

            std::cout << ix << " " << iy << " x=" << x << " y=" << y
                      << " detFrac=" << frac_detected << "\n";

            // Fill summary tree
            tsum.Fill();
            h_detected->SetBinContent(ix+1,iy+1,detected);
            h_absorbed->SetBinContent(ix+1,iy+1,absorbed);
            h_escaped->SetBinContent(ix+1,iy+1,escaped);
            h_hitPMT->SetBinContent(ix+1,iy+1,hitPMT);
            h_detected_s->SetBinContent(ix+1,iy+1,frac_detected);
            h_absorbed_s->SetBinContent(ix+1,iy+1,frac_absorbed);
            h_escaped_s->SetBinContent(ix+1,iy+1,frac_escaped);
            h_hitPMT_s->SetBinContent(ix+1,iy+1,frac_hitPMT);

            // Optional: draw one event from this point
            // DrawEventSplitViewFromTree(outFile.Data(), 12, true, 0.12);
        }
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

    fsum.cd();
    tsum.Write();
    fsum.Close();

    std::cout << "Wrote summary to " << summaryFile << " (tree tScan)\n";
}

