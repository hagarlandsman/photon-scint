// -----------------------------
// Drawing helpers#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TPolyLine3D.h"
#include "TPolyLine.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TPaveText.h"
#include "TString.h"
#include "TROOT.h"
#include "TPolyMarker3D.h"
#include "TPad.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "OpticsConfig.h"
static TPolyLine3D *Edge(double x1, double y1, double z1, double x2, double y2, double z2)
{
    auto e = new TPolyLine3D(2);
    e->SetPoint(0, x1, y1, z1);
    e->SetPoint(1, x2, y2, z2);
    e->SetLineWidth(2);
    return e;
}

static void DrawBoxWireframe(double L, double W, double T)
{
    double x0 = 0, x1 = L;
    double y0 = -W / 2, y1p = +W / 2;
    double z0 = -T / 2, z1p = +T / 2;

    // x=0
    Edge(x0, y0, z0, x0, y1p, z0)->Draw("same");
    Edge(x0, y1p, z0, x0, y1p, z1p)->Draw("same");
    Edge(x0, y1p, z1p, x0, y0, z1p)->Draw("same");
    Edge(x0, y0, z1p, x0, y0, z0)->Draw("same");

    // x=L
    Edge(x1, y0, z0, x1, y1p, z0)->Draw("same");
    Edge(x1, y1p, z0, x1, y1p, z1p)->Draw("same");
    Edge(x1, y1p, z1p, x1, y0, z1p)->Draw("same");
    Edge(x1, y0, z1p, x1, y0, z0)->Draw("same");

    // connect
    Edge(x0, y0, z0, x1, y0, z0)->Draw("same");
    Edge(x0, y1p, z0, x1, y1p, z0)->Draw("same");
    Edge(x0, y1p, z1p, x1, y1p, z1p)->Draw("same");
    Edge(x0, y0, z1p, x1, y0, z1p)->Draw("same");
}

// Wedge outline (visual only): drawn in the z=0 plane as an x-y outline
static void DrawWedgeOutlineXY(double L, double W, double wTip, double Lg)
{
    double z = 0.0;

    // left end
    Edge(0, -wTip / 2, z, 0, +wTip / 2, z)->Draw("same");
    Edge(Lg, -W / 2, z, Lg, +W / 2, z)->Draw("same");
    Edge(0, +wTip / 2, z, Lg, +W / 2, z)->Draw("same");
    Edge(0, -wTip / 2, z, Lg, -W / 2, z)->Draw("same");

    // right end
    Edge(L, -wTip / 2, z, L, +wTip / 2, z)->Draw("same");
    Edge(L - Lg, -W / 2, z, L - Lg, +W / 2, z)->Draw("same");
    Edge(L, +wTip / 2, z, L - Lg, +W / 2, z)->Draw("same");
    Edge(L, -wTip / 2, z, L - Lg, -W / 2, z)->Draw("same");
}

static void DrawOutlineXY(double L, double W)
{
    double y1 = -W / 2.0, y2 = +W / 2.0;

    auto pl = new TPolyLine(5);
    pl->SetPoint(0, 0.0, y1);
    pl->SetPoint(1, L, y1);
    pl->SetPoint(2, L, y2);
    pl->SetPoint(3, 0.0, y2);
    pl->SetPoint(4, 0.0, y1);
    pl->SetLineWidth(2);
    pl->Draw("same");
}

static void DrawOutlineXZ(double L, double T)
{
    double z1 = -T / 2.0, z2 = +T / 2.0;

    auto pl = new TPolyLine(5);
    pl->SetPoint(0, 0.0, z1);
    pl->SetPoint(1, L, z1);
    pl->SetPoint(2, L, z2);
    pl->SetPoint(3, 0.0, z2);
    pl->SetPoint(4, 0.0, z1);
    pl->SetLineWidth(2);
    pl->Draw("same");
}

// Wedge outline in XY only (taper in y)
static void DrawWedgeOutlineXY2D(double L, double W, double wedgeLen, double wedgeTipW)
{
    if (wedgeLen <= 0 || wedgeTipW <= 0 || wedgeTipW > W)
        return;

    double yTip = wedgeTipW / 2.0;
    double yFull = W / 2.0;

    // left wedge boundary lines (two diagonals)
    auto l1 = new TPolyLine(2);
    l1->SetPoint(0, 0.0, +yTip);
    l1->SetPoint(1, wedgeLen, +yFull);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->Draw("same");

    auto l2 = new TPolyLine(2);
    l2->SetPoint(0, 0.0, -yTip);
    l2->SetPoint(1, wedgeLen, -yFull);
    l2->SetLineStyle(2);
    l2->SetLineWidth(2);
    l2->Draw("same");

    // right wedge boundary lines
    auto r1 = new TPolyLine(2);
    r1->SetPoint(0, L, +yTip);
    r1->SetPoint(1, L - wedgeLen, +yFull);
    r1->SetLineStyle(2);
    r1->SetLineWidth(2);
    r1->Draw("same");

    auto r2 = new TPolyLine(2);
    r2->SetPoint(0, L, -yTip);
    r2->SetPoint(1, L - wedgeLen, -yFull);
    r2->SetLineStyle(2);
    r2->SetLineWidth(2);
    r2->Draw("same");
}

static void DrawInfoBox2D(
    Long64_t ievt, double L, double W, double T, int wrap,
    int reachedEnd, int pmt_side, int inPMT, int detected,
    int absorbed, int escaped, int endPlane, int nBounces, double pathLen)
{
    const char *wrapName = "wrap=?";
    if (wrap == 1)
        wrapName = "PTFE";
    else if (wrap == 2)
        wrapName = "Mylar";

    auto p = new TPaveText(0.10, 0.78, 0.90, 0.93, "NDC");
    p->SetFillStyle(0);
    p->SetBorderSize(1);
    p->SetTextAlign(12);
    p->SetTextSize(0.03);
    printf("event %lld | geometry: L=%.1f cm, W=%.1f cm, T=%.1f cm | wrap: %s \n",(long long)ievt, L, W, T, wrapName);
    printf("status: reachedEnd=%d (pmtSide=%d, inPMT=%d, detected=%d) | absorbed=%d | escaped=%d \n",reachedEnd, pmt_side, inPMT, detected, absorbed, escaped);
    printf("endPlane=%d | bounces=%d | path=%.2f cm \n", endPlane, nBounces, pathLen);

    p->AddText(Form("event %lld | geometry: L=%.1f cm, W=%.1f cm, T=%.1f cm | wrap: %s",
                    (long long)ievt, L, W, T, wrapName));
    p->AddText(Form("status: reachedEnd=%d (pmtSide=%d, inPMT=%d, detected=%d) | absorbed=%d | escaped=%d",
                    reachedEnd, pmt_side, inPMT, detected, absorbed, escaped));
    p->AddText(Form("endPlane=%d | bounces=%d | path=%.2f cm", endPlane, nBounces, pathLen));
    p->Draw();
}

static void DrawEventInfoBox(
    Long64_t ievt,
    double L, double W, double T,
    int wrap,
    int reachedEnd, int pmt_side, int inPMT, int detected,
    int absorbed, int escaped,
    int endPlane, int nBounces,
    double pathLen)
{
    const char *wrapName = "wrap=?";
    if (wrap == 1)
        wrapName = "PTFE";
    else if (wrap == 2)
        wrapName = "Mylar";

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
struct PathMarkers
{
    TPolyMarker3D *start = nullptr;
    TPolyMarker3D *end = nullptr;
};

static PathMarkers DrawStartEndMarkers3D(
    double xs, double ys, double zs,
    double xe, double ye, double ze)
{
    PathMarkers pm;

    pm.start = new TPolyMarker3D(1);
    pm.start->SetPoint(0, xs, ys, zs);
    pm.start->SetMarkerStyle(20);
    pm.start->SetMarkerSize(1.6);
    pm.start->SetMarkerColor(kGreen + 2);
    pm.start->Draw("same");

    pm.end = new TPolyMarker3D(1);
    pm.end->SetPoint(0, xe, ye, ze);
    pm.end->SetMarkerStyle(29);
    pm.end->SetMarkerSize(2.0);
    pm.end->SetMarkerColor(kRed + 1);
    pm.end->Draw("same");

    return pm;
}

void DrawEventFromTree(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    bool drawWedge = true,
    bool zoom = true,
    double zoomMargin = 0.15 // 15% padding around the path
)
{
    TFile f(fn);
    auto t = (TTree *)f.Get("tPhot");
    if (!t)
    {
        std::cout << "DrawEventFromTree: cannot find tPhot\n";
        return;
    }

    // Always present in the tree (per our writer)
    double L = 0, W = 0, T = 0;
    int wrap = 0;

    int reachedEnd = 0, pmt_side = 0, inPMT = 0, detected = 0;
    int absorbed = 0, escaped = 0, endPlane = 0, nBounces = 0;
    double pathLen = 0;

    double wedgeLen = 20.0;
    double wedgeTipW = 5.0;

    t->SetBranchAddress("cfg_wedgeLen", &wedgeLen);
    t->SetBranchAddress("cfg_wedgeTipW", &wedgeTipW);
    t->SetBranchAddress("cfg_L", &L);
    t->SetBranchAddress("cfg_W", &W);
    t->SetBranchAddress("cfg_T", &T);


    // These should exist if you used the TreeWriter block I gave
    if (t->GetBranch("reachedEnd"))
        t->SetBranchAddress("reachedEnd", &reachedEnd);
    if (t->GetBranch("pmt_side"))
        t->SetBranchAddress("pmt_side", &pmt_side);
    if (t->GetBranch("inPMT"))
        t->SetBranchAddress("inPMT", &inPMT);
    if (t->GetBranch("detected"))
        t->SetBranchAddress("detected", &detected);

    if (t->GetBranch("absorbed"))
        t->SetBranchAddress("absorbed", &absorbed);
    if (t->GetBranch("escaped"))
        t->SetBranchAddress("escaped", &escaped);
    if (t->GetBranch("endPlane"))
        t->SetBranchAddress("endPlane", &endPlane);
    if (t->GetBranch("nBounces"))
        t->SetBranchAddress("nBounces", &nBounces);
    if (t->GetBranch("path"))
        t->SetBranchAddress("path", &pathLen);

    // Path branches exist only if you ran with cfg.savePath=true
    std::vector<float> *xP = nullptr;
    std::vector<float> *yP = nullptr;
    std::vector<float> *zP = nullptr;

    bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
    if (hasPath)
    {
        t->SetBranchAddress("xPath", &xP);
        t->SetBranchAddress("yPath", &yP);
        t->SetBranchAddress("zPath", &zP);
    }

    t->GetEntry(ievt);

    const char *wrapName = "wrap=?";
    if (wrap == 1)
        wrapName = "PTFE";
    else if (wrap == 2)
        wrapName = "Mylar";

    // Build a descriptive title

    TString title;
    title.Form("evt %lld", (long long)ievt);
    /* title.Form(
       "evt %lld | geom L=%.1f W=%.1f T=%.1f cm | wrap=%s | reachedEnd=%d pmtSide=%d inPMT=%d detected=%d | absorbed=%d escaped=%d | endPlane=%d bounces=%d path=%.2f cm; x (cm); y (cm); z (cm)",
       (long long)ievt, L, W, T, wrapName,
       reachedEnd, pmt_side, inPMT, detected,
       absorbed, escaped,
       endPlane, nBounces, pathLen
     );
   */

     TCanvas*   c = new TCanvas("cEvt", "event", 1100, 800);

    // Default view ranges
    double xmin = 0, xmax = L;
    double ymin = -W / 2, ymax = +W / 2;
    double zmin = -T / 2, zmax = +T / 2;

    // Zoom range around the path if requested and available
    if (zoom && hasPath && xP && xP->size() >= 2)
    {
        double xlo = (*xP)[0], xhi = (*xP)[0];
        double ylo = (*yP)[0], yhi = (*yP)[0];
        double zlo = (*zP)[0], zhi = (*zP)[0];

        for (size_t i = 1; i < xP->size(); i++)
        {
            xlo = std::min(xlo, (double)(*xP)[i]);
            xhi = std::max(xhi, (double)(*xP)[i]);
            ylo = std::min(ylo, (double)(*yP)[i]);
            yhi = std::max(yhi, (double)(*yP)[i]);
            zlo = std::min(zlo, (double)(*zP)[i]);
            zhi = std::max(zhi, (double)(*zP)[i]);
        }

        // Add margin
        double dx = (xhi - xlo);
        if (dx <= 0)
            dx = 1.0;
        double dy = (yhi - ylo);
        if (dy <= 0)
            dy = 1.0;
        double dz = (zhi - zlo);
        if (dz <= 0)
            dz = 1.0;

        xlo -= zoomMargin * dx;
        xhi += zoomMargin * dx;
        ylo -= zoomMargin * dy;
        yhi += zoomMargin * dy;
        zlo -= zoomMargin * dz;
        zhi += zoomMargin * dz;

        // Keep inside detector bounds (optional but usually nice)
        xlo = std::max(xlo, 0.0);
        xhi = std::min(xhi, L);
        ylo = std::max(ylo, -W / 2);
        yhi = std::min(yhi, W / 2);
        zlo = std::max(zlo, -T / 2);
        zhi = std::min(zhi, T / 2);

        xmin = xlo;
        xmax = xhi;
        ymin = ylo;
        ymax = yhi;
        zmin = zlo;
        zmax = zhi;
    }
    else
    {
        // Keep a 1:1:1 feel by expanding to a cube around the detector
        double maxDim = std::max(L, std::max(W, T));
        double xPad = 0.5 * (maxDim - L);
        xmin = -xPad;
        xmax = L + xPad;
        ymin = -0.5 * maxDim;
        ymax = 0.5 * maxDim;
        zmin = -0.5 * maxDim;
        zmax = 0.5 * maxDim;
    }

    auto frame = new TH3D("fr", title,
                          10, xmin, xmax,
                          10, ymin, ymax,
                          10, zmin, zmax);
    frame->SetStats(0);
    frame->Draw();

    DrawBoxWireframe(L, W, T);
    if (drawWedge)
        DrawWedgeOutlineXY(L, W, wedgeTipW, wedgeLen);

    if (hasPath && xP && xP->size() >= 2)
    {
        auto pl = new TPolyLine3D((int)xP->size());
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");
        const int n = (int)xP->size();
        auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                        (*xP)[n - 1], (*yP)[n - 1], (*zP)[n - 1]);

        auto leg = new TLegend(0.62, 0.22, 0.99, 0.36); // x1,y1,x2,y2 in NDC
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);

        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(pm.start, "Start point", "p");
        leg->AddEntry(pm.end, "End point", "p");
        leg->Draw();
    }
    else
    {
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
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    bool drawWedge = true
    )
{
    TFile f(fn);
    auto t = (TTree *)f.Get("tPhot");
    if (!t)
    {
        std::cout << "DrawEventFromTree: cannot find tPhot\n";
        return;
    }

    double L = 0, W = 0, T = 0;
    int wrap = 0;
    std::vector<float> *xP = nullptr;
    std::vector<float> *yP = nullptr;
    std::vector<float> *zP = nullptr;
    double wedgeLen = 20.0;
    double wedgeTipW = 5.0;

    t->SetBranchAddress("cfg_wedgeLen", &wedgeLen);
    t->SetBranchAddress("cfg_wedgeTipW", &wedgeTipW);

    t->SetBranchAddress("cfg_L", &L);
    t->SetBranchAddress("cfg_W", &W);
    t->SetBranchAddress("cfg_T", &T);
    t->SetBranchAddress("cfg_wrap", &wrap);

    // These branches exist only if you ran with cfg.savePath=true
    bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
    if (hasPath)
    {
        t->SetBranchAddress("xPath", &xP);
        t->SetBranchAddress("yPath", &yP);
        t->SetBranchAddress("zPath", &zP);
    }

    t->GetEntry(ievt);

     TCanvas*  c = new TCanvas("cEvt", "event", 900, 700);

    // Force a 1:1:1 feel by drawing a cubic frame (equal axis numeric spans)
    double maxDim = std::max(L, std::max(W, T));
    double xPad = 0.5 * (maxDim - L);

    auto frame = new TH3D("fr", "Scintillator and one photon; x (cm); y (cm); z (cm)",
                          10, -xPad, L + xPad,
                          10, -0.5 * maxDim, 0.5 * maxDim,
                          10, -0.5 * maxDim, 0.5 * maxDim);
    frame->SetStats(0);
    frame->Draw();

    DrawBoxWireframe(L, W, T);
    if (drawWedge)
        DrawWedgeOutlineXY(L, W, wedgeTipW, wedgeLen);

    if (hasPath && xP && xP->size() >= 2)
    {
        auto pl = new TPolyLine3D((int)xP->size());
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");
        const int n = (int)xP->size();
        auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                        (*xP)[n - 1], (*yP)[n - 1], (*zP)[n - 1]);

        auto leg = new TLegend(0.62, 0.22, 0.99, 0.36); // x1,y1,x2,y2 in NDC
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);

        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(pm.start, "Start point", "p");
        leg->AddEntry(pm.end, "End point", "p");
        leg->Draw();
    }
    else
    {
        std::cout << "DrawEventFromTree: no saved path in file (run with cfg.savePath=true)\n";
    }

    c->Modified();
    c->Update();
}

void DrawEventSplitViewFromTree(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    bool drawWedge = true,
    double zoomMargin = 0.12 // padding for the zoomed-in view
)
{
    TFile f(fn);
    auto t = (TTree *)f.Get("tPhot");
    if (!t)
    {
        std::cout << "DrawEventSplitViewFromTree: cannot find tPhot\n";
        return;
    }

    double L = 0, W = 0, T = 0;
    int wrap = 0;
    int reachedEnd = 0, pmt_side = 0, inPMT = 0, detected = 0;
    int absorbed = 0, escaped = 0, endPlane = 0, nBounces = 0;
    double pathLen = 0;

    double wedgeLen = 20.0;
    double wedgeTipW = 5.0;

    t->SetBranchAddress("cfg_wedgeLen", &wedgeLen);
    t->SetBranchAddress("cfg_wedgeTipW", &wedgeTipW);


    t->SetBranchAddress("cfg_L", &L);
    t->SetBranchAddress("cfg_W", &W);
    t->SetBranchAddress("cfg_T", &T);
    t->SetBranchAddress("cfg_wrap", &wrap);

    if (t->GetBranch("reachedEnd"))
        t->SetBranchAddress("reachedEnd", &reachedEnd);
    if (t->GetBranch("pmt_side"))
        t->SetBranchAddress("pmt_side", &pmt_side);
    if (t->GetBranch("inPMT"))
        t->SetBranchAddress("inPMT", &inPMT);
    if (t->GetBranch("detected"))
        t->SetBranchAddress("detected", &detected);

    if (t->GetBranch("absorbed"))
        t->SetBranchAddress("absorbed", &absorbed);
    if (t->GetBranch("escaped"))
        t->SetBranchAddress("escaped", &escaped);
    if (t->GetBranch("endPlane"))
        t->SetBranchAddress("endPlane", &endPlane);
    if (t->GetBranch("nBounces"))
        t->SetBranchAddress("nBounces", &nBounces);
    if (t->GetBranch("path"))
        t->SetBranchAddress("path", &pathLen);

    std::vector<float> *xP = nullptr;
    std::vector<float> *yP = nullptr;
    std::vector<float> *zP = nullptr;

    bool hasPath = (t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath"));
    if (!hasPath)
    {
        std::cout << "No xPath/yPath/zPath in file. Run with cfg.savePath=true\n";
        return;
    }

    t->SetBranchAddress("xPath", &xP);
    t->SetBranchAddress("yPath", &yP);
    t->SetBranchAddress("zPath", &zP);

    t->GetEntry(ievt);
    if (!xP || xP->size() < 2)
    {
        std::cout << "Event has no stored path points\n";
        return;
    }

    // Compute zoomed-in bounds from path
    double xlo = (*xP)[0], xhi = (*xP)[0];
    double ylo = (*yP)[0], yhi = (*yP)[0];
    double zlo = (*zP)[0], zhi = (*zP)[0];

    for (size_t i = 1; i < xP->size(); i++)
    {
        xlo = std::min(xlo, (double)(*xP)[i]);
        xhi = std::max(xhi, (double)(*xP)[i]);
        ylo = std::min(ylo, (double)(*yP)[i]);
        yhi = std::max(yhi, (double)(*yP)[i]);
        zlo = std::min(zlo, (double)(*zP)[i]);
        zhi = std::max(zhi, (double)(*zP)[i]);
    }

    double dx = xhi - xlo;
    if (dx <= 0)
        dx = 1;
    double dy = yhi - ylo;
    if (dy <= 0)
        dy = 1;
    double dz = zhi - zlo;
    if (dz <= 0)
        dz = 1;

    xlo -= zoomMargin * dx;
    xhi += zoomMargin * dx;
    ylo -= zoomMargin * dy;
    yhi += zoomMargin * dy;
    zlo -= zoomMargin * dz;
    zhi += zoomMargin * dz;

    // Clamp zoom-in to detector bounds
    xlo = std::max(xlo, 0.0);
    xhi = std::min(xhi, L);
    ylo = std::max(ylo, -W / 2);
    yhi = std::min(yhi, W / 2);
    zlo = std::max(zlo, -T / 2);
    zhi = std::min(zhi, T / 2);

    // Zoomed-out cube bounds (for 1:1:1 feel)
    double maxDim = std::max(L, std::max(W, T));
    double xPad = 0.5 * (maxDim - L);

    double xo_min = -xPad, xo_max = L + xPad;
    double yo_min = -0.5 * maxDim, yo_max = 0.5 * maxDim;
    double zo_min = -0.5 * maxDim, zo_max = 0.5 * maxDim;

    // Unique canvas name
    TString cname = Form("cSplit_%lld", (long long)ievt);
    if (auto old = gROOT->FindObject(cname))
        delete old;

    TCanvas* c = new TCanvas(cname, "event split view", 1400, 1000);
    c->Divide(1, 2);

    auto drawOnePad = [&](int ipad,
                          double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                          const char *padTitle)
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

        DrawBoxWireframe(L, W, T);
        if (drawWedge)
            DrawWedgeOutlineXY(L, W, wedgeTipW, wedgeLen);

        auto pl = new TPolyLine3D((int)xP->size());
        //    pl->SetName(Form("pl_%lld_%d", (long long)ievt, ipad));
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");
        const int n = (int)xP->size();
        auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                        (*xP)[n - 1], (*yP)[n - 1], (*zP)[n - 1]);

        // Put info box only once (left pad) or on both if you want
        if (ipad == 1)
        {
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

    drawOnePad(1, xo_min, xo_max, yo_min, yo_max, zo_min, zo_max, "Zoom out");
    drawOnePad(2, xlo, xhi, ylo, yhi, zlo, zhi, "Zoom in");

    c->Modified();
    c->Update();
}
void DrawEvent4ViewFromTree(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    bool drawWedge = true,
    double zoomMargin = 0.12)
{
    TFile f(fn);
    auto t = (TTree *)f.Get("tPhot");
    if (!t)
    {
        std::cout << "DrawEvent4ViewFromTree: cannot find tPhot\n";
        return;
    }

    double L = 0, W = 0, T = 0;
    int wrap = 0;
    int reachedEnd = 0, pmt_side = 0, inPMT = 0, detected = 0;
    int absorbed = 0, escaped = 0, endPlane = 0, nBounces = 0;
    double pathLen = 0;
    double wedgeLen = 20.0;
    double wedgeTipW = 5.0;

    t->SetBranchAddress("cfg_wedgeLen", &wedgeLen);
    t->SetBranchAddress("cfg_wedgeTipW", &wedgeTipW);

    t->SetBranchAddress("cfg_L", &L);
    t->SetBranchAddress("cfg_W", &W);
    t->SetBranchAddress("cfg_T", &T);
    t->SetBranchAddress("cfg_wrap", &wrap);

    if (t->GetBranch("reachedEnd"))
        t->SetBranchAddress("reachedEnd", &reachedEnd);
    if (t->GetBranch("pmt_side"))
        t->SetBranchAddress("pmt_side", &pmt_side);
    if (t->GetBranch("inPMT"))
        t->SetBranchAddress("inPMT", &inPMT);
    if (t->GetBranch("detected"))
        t->SetBranchAddress("detected", &detected);
    if (t->GetBranch("absorbed"))
        t->SetBranchAddress("absorbed", &absorbed);
    if (t->GetBranch("escaped"))
        t->SetBranchAddress("escaped", &escaped);
    if (t->GetBranch("endPlane"))
        t->SetBranchAddress("endPlane", &endPlane);
    if (t->GetBranch("nBounces"))
        t->SetBranchAddress("nBounces", &nBounces);
    if (t->GetBranch("path"))
        t->SetBranchAddress("path", &pathLen);

    std::vector<float> *xP = nullptr;
    std::vector<float> *yP = nullptr;
    std::vector<float> *zP = nullptr;

    if (!t->GetBranch("xPath") || !t->GetBranch("yPath") || !t->GetBranch("zPath"))
    {
        std::cout << "DrawEvent4ViewFromTree: no xPath/yPath/zPath. Run with cfg.savePath=true\n";
        return;
    }
    t->SetBranchAddress("xPath", &xP);
    t->SetBranchAddress("yPath", &yP);
    t->SetBranchAddress("zPath", &zP);

    t->GetEntry(ievt);
    if (!xP || xP->size() < 2)
    {
        std::cout << "Event has no stored path points\n";
        return;
    }

    // Compute path bounds
    double xlo = (*xP)[0], xhi = (*xP)[0];
    double ylo = (*yP)[0], yhi = (*yP)[0];
    double zlo = (*zP)[0], zhi = (*zP)[0];

    for (size_t i = 1; i < xP->size(); i++)
    {
        xlo = std::min(xlo, (double)(*xP)[i]);
        xhi = std::max(xhi, (double)(*xP)[i]);
        ylo = std::min(ylo, (double)(*yP)[i]);
        yhi = std::max(yhi, (double)(*yP)[i]);
        zlo = std::min(zlo, (double)(*zP)[i]);
        zhi = std::max(zhi, (double)(*zP)[i]);
    }

    // Zoom bounds with margin, clamped to detector
    double dx = xhi - xlo;
    if (dx <= 0)
        dx = 1;
    double dy = yhi - ylo;
    if (dy <= 0)
        dy = 1;
    double dz = zhi - zlo;
    if (dz <= 0)
        dz = 1;

    double zx1 = std::max(0.0, xlo - zoomMargin * dx);
    double zx2 = std::min(L, xhi + zoomMargin * dx);
    double zy1 = std::max(-W / 2.0, ylo - zoomMargin * dy);
    double zy2 = std::min(+W / 2.0, yhi + zoomMargin * dy);
    double zz1 = std::max(-T / 2.0, zlo - zoomMargin * dz);
    double zz2 = std::min(+T / 2.0, zhi + zoomMargin * dz);

    // Zoom-out bounds in each projection
    double xo1 = 0.0, xo2 = L;
    double yo1 = -W / 2.0, yo2 = +W / 2.0;
    double zo1 = -T / 2.0, zo2 = +T / 2.0;

    // Unique canvas name
    //TString cname = Form("c4_%lld", (long long)ievt);
    TString cname = "c4_scint";

    if (auto old = gROOT->FindObject(cname))
        delete old;

     TCanvas*  c = new TCanvas(cname, "XY and XZ, zoomed and unzoomed", 1400, 1000);
    c->Divide(2, 2, 0.01, 0.01);

    auto drawXY = [&](int pad, bool zoomed)
    {
        c->cd(pad);
        gPad->SetGrid(0, 0);
        gPad->SetFixedAspectRatio();

        double xmin = zoomed ? zx1 : xo1;
        double xmax = zoomed ? zx2 : xo2;
        double ymin = zoomed ? zy1 : yo1;
        double ymax = zoomed ? zy2 : yo2;

        TString hname = Form("hxy_%lld_%d", (long long)ievt, pad);
        auto h = new TH2D(hname, zoomed ? "XY zoom in; x (cm); y (cm)" : "XY zoom out; x (cm); y (cm)",
                          10, xmin, xmax, 10, ymin, ymax);
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->Draw();

        DrawOutlineXY(L, W);
        if (drawWedge)
            DrawWedgeOutlineXY2D(L, W, wedgeLen, wedgeTipW);

        // Path polyline in XY
        auto pl = new TPolyLine((int)xP->size());
        pl->SetLineWidth(3);
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*yP)[i]);
        pl->Draw("same");

        // Start (green) and end (red)
        const int n = (int)xP->size();
        auto mS = new TMarker((*xP)[0], (*yP)[0], 20);
        mS->SetMarkerColor(kGreen + 2);
        mS->SetMarkerSize(1.5);
        mS->Draw("same");

        auto mE = new TMarker((*xP)[n - 1], (*yP)[n - 1], 29);
        mE->SetMarkerColor(kRed + 1);
        mE->SetMarkerSize(1.7);
        mE->Draw("same");

        // Legend
        auto leg = new TLegend(0.12, 0.12, 0.42, 0.27);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);
        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(mS, "Start point", "p");
        leg->AddEntry(mE, "End point", "p");
        leg->Draw();

        // Info box only on the top-left pad
        if (pad == 1)
        {
            DrawInfoBox2D(ievt, L, W, T, wrap,
                          reachedEnd, pmt_side, inPMT, detected,
                          absorbed, escaped, endPlane, nBounces, pathLen);
        }
    };

    auto drawXZ = [&](int pad, bool zoomed)
    {
        c->cd(pad);
        gPad->SetGrid(0, 0);
        gPad->SetFixedAspectRatio();

        double xmin = zoomed ? zx1 : xo1;
        double xmax = zoomed ? zx2 : xo2;
        double zmin = zoomed ? zz1 : zo1;
        double zmax = zoomed ? zz2 : zo2;

        TString hname = Form("hxz_%lld_%d", (long long)ievt, pad);
        auto h = new TH2D(hname, zoomed ? "XZ zoom in; x (cm); z (cm)" : "XZ zoom out; x (cm); z (cm)",
                          10, xmin, xmax, 10, zmin, zmax);
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->Draw();

        DrawOutlineXZ(L, T);

        // Path polyline in XZ
        auto pl = new TPolyLine((int)xP->size());
        pl->SetLineWidth(3);
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*zP)[i]);
        pl->Draw("same");

        // Start (green) and end (red)
        const int n = (int)xP->size();
        auto mS = new TMarker((*xP)[0], (*zP)[0], 20);
        mS->SetMarkerColor(kGreen + 2);
        mS->SetMarkerSize(1.5);
        mS->Draw("same");

        auto mE = new TMarker((*xP)[n - 1], (*zP)[n - 1], 29);
        mE->SetMarkerColor(kRed + 1);
        mE->SetMarkerSize(1.7);
        mE->Draw("same");

        // Legend
        auto leg = new TLegend(0.12, 0.12, 0.42, 0.27);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);
        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(mS, "Start point", "p");
        leg->AddEntry(mE, "End point", "p");
        leg->Draw();
    };

    // Pad mapping:
    // 1: XY zoom out, 2: XY zoom in, 3: XZ zoom out, 4: XZ zoom in
    drawXY(1, false);
    drawXY(2, true);
    drawXZ(3, false);
    drawXZ(4, true);

    ////
 // Compute zoomed-in bounds from path
     xlo = (*xP)[0], xhi = (*xP)[0];
     ylo = (*yP)[0], yhi = (*yP)[0];
     zlo = (*zP)[0], zhi = (*zP)[0];

    for (size_t i = 1; i < xP->size(); i++)
    {
        xlo = std::min(xlo, (double)(*xP)[i]);
        xhi = std::max(xhi, (double)(*xP)[i]);
        ylo = std::min(ylo, (double)(*yP)[i]);
        yhi = std::max(yhi, (double)(*yP)[i]);
        zlo = std::min(zlo, (double)(*zP)[i]);
        zhi = std::max(zhi, (double)(*zP)[i]);
    }

     dx = xhi - xlo;
    if (dx <= 0)
        dx = 1;
     dy = yhi - ylo;
    if (dy <= 0)
        dy = 1;
     dz = zhi - zlo;
    if (dz <= 0)
        dz = 1;

    xlo -= zoomMargin * dx;
    xhi += zoomMargin * dx;
    ylo -= zoomMargin * dy;
    yhi += zoomMargin * dy;
    zlo -= zoomMargin * dz;
    zhi += zoomMargin * dz;

    // Clamp zoom-in to detector bounds
    xlo = std::max(xlo, 0.0);
    xhi = std::min(xhi, L);
    ylo = std::max(ylo, -W / 2);
    yhi = std::min(yhi, W / 2);
    zlo = std::max(zlo, -T / 2);
    zhi = std::min(zhi, T / 2);

    // Zoomed-out cube bounds (for 1:1:1 feel)
    double maxDim = std::max(L, std::max(W, T));
    double xPad = 0.5 * (maxDim - L);

    double xo_min = -xPad, xo_max = L + xPad;
    double yo_min = -0.5 * maxDim, yo_max = 0.5 * maxDim;
    double zo_min = -0.5 * maxDim, zo_max = 0.5 * maxDim;

    // Unique canvas name




    auto drawOnePad = [&](int ipad,
                          double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                          const char *padTitle)
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

        DrawBoxWireframe(L, W, T);
        if (drawWedge)
            DrawWedgeOutlineXY(L, W, wedgeTipW, wedgeLen);

        auto pl = new TPolyLine3D((int)xP->size());
        //    pl->SetName(Form("pl_%lld_%d", (long long)ievt, ipad));
        for (int i = 0; i < (int)xP->size(); i++)
            pl->SetPoint(i, (*xP)[i], (*yP)[i], (*zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");
        const int n = (int)xP->size();
        auto pm = DrawStartEndMarkers3D((*xP)[0], (*yP)[0], (*zP)[0],
                                        (*xP)[n - 1], (*yP)[n - 1], (*zP)[n - 1]);

        // Put info box only once (left pad) or on both if you want
        if (ipad == 1)
        {
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

    drawOnePad(1, xo_min, xo_max, yo_min, yo_max, zo_min, zo_max, "Zoom out");
    drawOnePad(3, xlo, xhi, ylo, yhi, zlo, zhi, "Zoom in");

    c->Modified();
    c->Update();

    /////
    c->Modified();
    c->Update();
}

void DrawFracColz(const char* fn="scan_summary.root", const char* tname = "tScan")
{
  TFile f(fn);
  auto t = (TTree*)f.Get(tname);
  if(!t){ std::cout<<"no tScan\n"; return; }
    else {printf ("N=%lld \n",t->GetEntries());}
  // choose binning (match your grid)
  int nx = 3, ny = 3;

  // detected fraction
  t->Draw(Form("y:x>>hDet(%d,0,90,%d,-15,15)", nx, ny), "frac_detected", "colz");

  // if you want, set titles after the draw:
  auto hDet = (TH2*)gROOT->FindObject("hDet");
  if(hDet){
    hDet->SetTitle("Detected fraction; x (cm); y (cm)");
    gPad->Modified(); gPad->Update();
  }
}

void PrintCountsFromTree(const char* fn="toyOptics.root", const char* tn="tPhot")
{
    TFile f(fn, "READ");
    if (f.IsZombie()) {
        std::cout << "Cannot open file: " << fn << "\n";
        return;
    }

    auto t = (TTree*)f.Get(tn);
    if (!t) {
        std::cout << "Cannot find tree: " << tn << "\n";
        return;
    }

    auto has = [&](const char* b){ return t->GetBranch(b) != nullptr; };

    // You might have one of these depending on how you wrote the tree
    // absorbed: "absorbed"
    // escaped:  "escaped"
    // hitPMT:   "inPMT" or "hitPMT"
    // detected: "detected"
    const char* bAbs = has("absorbed") ? "absorbed" : nullptr;
    const char* bEsc = has("escaped")  ? "escaped"  : nullptr;
    const char* bHit = has("inPMT")    ? "inPMT"    : (has("hitPMT") ? "hitPMT" : nullptr);
    const char* bDet = has("detected") ? "detected" : nullptr;

    if (!bAbs) std::cout << "Missing branch: absorbed\n";
    if (!bEsc) std::cout << "Missing branch: escaped\n";
    if (!bHit) std::cout << "Missing branch: inPMT or hitPMT\n";
    if (!bDet) std::cout << "Missing branch: detected\n";

    auto countOnes = [&](const char* bname)->Long64_t {
        if (!bname) return 0;
        // assumes 0/1 per entry
        TString sel = Form("%s!=0", bname);
        return t->GetEntries(sel);
    };

    Long64_t nTot = t->GetEntries();

    Long64_t nAbs = countOnes(bAbs);
    Long64_t nEsc = countOnes(bEsc);
    Long64_t nHit = countOnes(bHit);
    Long64_t nDet = countOnes(bDet);

    std::cout << "File: " << fn << " | Tree: " << tn << "\n";
    std::cout << "Total entries: " << nTot << "\n";
    std::cout << "absorbed: "  << nAbs << "\n";
    std::cout << "escaped: "   << nEsc << "\n";
    std::cout << "hitPMT: "    << nHit << "  (branch=" << (bHit ? bHit : "none") << ")\n";
    std::cout << "detected: "  << nDet << "\n";

    OpticsConfig* cfg = nullptr;

    TTree *tCfg = (TTree*)f.Get("tCfg");
    tCfg->SetBranchAddress("cfg",&cfg);
   // tCfg->GetEntry(0);
   // tCfg->Print();
   printf ("W=%f \n",cfg->W);
    tCfg->Scan("L:W:T:nScint:nOut");
    //tCfg->Show(0);

}
