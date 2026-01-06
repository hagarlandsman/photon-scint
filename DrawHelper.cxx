// DrawHelper.C (or DrawHelper.cxx)
//
// Full helper file:
// - DrawEventFromTree3D: single 3D view (geometry + path)
// - DrawEvent4ViewFromTree: 4 pads (XY/XZ zoom-out + zoom-in) with 2D projections of BOTH geometry and path
// - DrawEventSplitViewFromTree: 2 pads (3D zoom-out + 3D zoom-in) with geometry + path
// - DrawFracColz: quick colz plot
// - PrintCountsFromTree: quick counters
//
// Key fixes:
// - Geometry is read once from TGraph2D "scintGeometry" and converted into
//   (A) a 3D polyline (for 3D pads) and
//   (B) 2D polylines for XY and XZ projections (for 2D pads).
// - No drawing 3D objects in 2D pads.
// - No passing TPolyLine3D by value.
// - Zoom logic uses correct clamping for your [0,L] and [-W/2,W/2], [-T/2,T/2] conventions.

#include "TRandom3.h"
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
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "OpticsConfig.h"

using std::cout;
using std::endl;

// -----------------------------
// Small drawing helpers

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


static void DrawInfoBox2D(
    Long64_t ievt, double L, double W, double T, int wrap,
    int reachedEnd, int pmt_side, int inPMT, int detected,
    int absorbed, int escaped, int endPlane, int nBounces, double pathLen)
{
    const char *wrapName = "wrap=?";
    if (wrap == 1) wrapName = "PTFE";
    else if (wrap == 2) wrapName = "Mylar";

    auto p = new TPaveText(0.10, 0.78, 0.90, 0.93, "NDC");
    p->SetFillStyle(0);
    p->SetBorderSize(1);
    p->SetTextAlign(12);
    p->SetTextSize(0.03);

    p->AddText(Form("event %lld | L=%.1f W=%.1f T=%.1f cm | wrap: %s",
                    (long long)ievt, L, W, T, wrapName));
    p->AddText(Form("reachedEnd=%d (side=%d, inPMT=%d, detected=%d) | absorbed=%d | escaped=%d",
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
    if (wrap == 1) wrapName = "PTFE";
    else if (wrap == 2) wrapName = "Mylar";

    auto p = new TPaveText(0.12, 0.78, 0.88, 0.93, "NDC");
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetBorderSize(1);
    p->SetTextAlign(12);
    p->SetTextSize(0.028);

    p->AddText(Form("event %lld | geometry: L=%.1f cm, W=%.1f cm, T=%.1f cm",
                    (long long)ievt, L, W, T));
    p->AddText(Form("wrap: %s | endPlane=%d | bounces=%d | path=%.2f cm",
                    wrapName, endPlane, nBounces, pathLen));
    p->AddText(Form("end: reachedEnd=%d (side=%d, inPMT=%d, detected=%d) | absorbed=%d | escaped=%d",
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

// -----------------------------
// Geometry conversion helpers

static TPolyLine3D* BuildFrameGeometry3D(TGraph2D* g2, bool closeLoop=true)
{
    const int n = g2->GetN();
    if (n <= 0) return nullptr;

    auto pl = new TPolyLine3D(closeLoop ? (n + 1) : n);

    double x,y,z;
    for (int i=0;i<n;++i) {
        g2->GetPoint(i,x,y,z);
        pl->SetPoint(i,x,y,z);
    }
    if (closeLoop) {
        g2->GetPoint(0,x,y,z);
        pl->SetPoint(n,x,y,z);
    }

    pl->SetLineWidth(2);
    pl->SetLineColor(kBlack);
    return pl;
}

static TPolyLine* BuildFrameGeometryXY(TGraph2D* g2, bool closeLoop=true)
{
    const int n = g2->GetN();
    if (n <= 0) return nullptr;

    auto pl = new TPolyLine(closeLoop ? (n + 1) : n);

    double x,y,z;
    for (int i=0;i<n;++i) {
        g2->GetPoint(i,x,y,z);
        pl->SetPoint(i, x, y);
    }
    if (closeLoop) {
        g2->GetPoint(0,x,y,z);
        pl->SetPoint(n, x, y);
    }

    pl->SetLineWidth(2);
    pl->SetLineColor(kBlack);
    return pl;
}

static TPolyLine* BuildFrameGeometryXZ(TGraph2D* g2, bool closeLoop=true)
{
    const int n = g2->GetN();
    if (n <= 0) return nullptr;

    auto pl = new TPolyLine(closeLoop ? (n + 1) : n);

    double x,y,z;
    for (int i=0;i<n;++i) {
        g2->GetPoint(i,x,y,z);
        pl->SetPoint(i, x, z);
    }
    if (closeLoop) {
        g2->GetPoint(0,x,y,z);
        pl->SetPoint(n, x, z);
    }

    pl->SetLineWidth(2);
    pl->SetLineColor(kBlack);
    return pl;
}

// -----------------------------
// Read common event branches

struct EventData
{
    double L=0, W=0, T=0;
    int wrap=0;

    int reachedEnd=0, pmt_side=0, inPMT=0, detected=0;
    int absorbed=0, escaped=0, endPlane=0, nBounces=0;
    double pathLen=0;

    double wedgeLen=0, wedgeTipW=0;

    std::vector<float>* xP=nullptr;
    std::vector<float>* yP=nullptr;
    std::vector<float>* zP=nullptr;
    bool hasPath=false;
};

static bool LoadEvent(TTree* t, Long64_t ievt, EventData& ev)
{
    if (!t) return false;

    // cfg
    t->SetBranchAddress("cfg_L",&ev.L);
    t->SetBranchAddress("cfg_W",&ev.W);
    t->SetBranchAddress("cfg_T",&ev.T);
    if (t->GetBranch("cfg_wrap")) t->SetBranchAddress("cfg_wrap",&ev.wrap);
    if (t->GetBranch("cfg_wedgeLen")) t->SetBranchAddress("cfg_wedgeLen",&ev.wedgeLen);
    if (t->GetBranch("cfg_wedgeTipW")) t->SetBranchAddress("cfg_wedgeTipW",&ev.wedgeTipW);

    // results (optional)
    if (t->GetBranch("reachedEnd")) t->SetBranchAddress("reachedEnd",&ev.reachedEnd);
    if (t->GetBranch("pmt_side"))  t->SetBranchAddress("pmt_side",&ev.pmt_side);
    if (t->GetBranch("inPMT"))     t->SetBranchAddress("inPMT",&ev.inPMT);
    if (t->GetBranch("detected"))  t->SetBranchAddress("detected",&ev.detected);
    if (t->GetBranch("absorbed"))  t->SetBranchAddress("absorbed",&ev.absorbed);
    if (t->GetBranch("escaped"))   t->SetBranchAddress("escaped",&ev.escaped);
    if (t->GetBranch("endPlane"))  t->SetBranchAddress("endPlane",&ev.endPlane);
    if (t->GetBranch("nBounces"))  t->SetBranchAddress("nBounces",&ev.nBounces);
    if (t->GetBranch("path"))      t->SetBranchAddress("path",&ev.pathLen);

    // paths
    ev.hasPath = t->GetBranch("xPath") && t->GetBranch("yPath") && t->GetBranch("zPath");
    if (ev.hasPath) {
        t->SetBranchAddress("xPath",&ev.xP);
        t->SetBranchAddress("yPath",&ev.yP);
        t->SetBranchAddress("zPath",&ev.zP);
        }

    t->GetEntry(ievt);
    cout<<"Loaded event "<<ievt<<"\n";
    cout<<"  L="<<ev.L<<" W="<<ev.W<<" T="<<ev.T<<" wrap="<<ev.wrap<<"\n";
    cout<<" ev length="<<ev.xP->size()<<"\n";
    return true;
}

// -----------------------------
// 3D single view

void DrawEventFromTree3D(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0
    )
{
    TFile f(fn);
    if (f.IsZombie()) return;

    auto t = (TTree*)f.Get("tPhot");
    if (!t) return;

    auto g2 = (TGraph2D*)f.Get("scintGeometry");
    if (!g2) return;

    EventData ev;
    if (!LoadEvent(t, ievt, ev)) return;

    auto frame3D = BuildFrameGeometry3D(g2, true);

    TString cname = Form("c3D_%lld", (long long)ievt);
    if (auto old = gROOT->FindObject(cname)) delete old;
    auto c = new TCanvas(cname, "3D event", 1100, 800);

    double maxDim = std::max(ev.L, std::max(ev.W, ev.T));
    auto frame = new TH3D("fr3d", "3D view; x (cm); y (cm); z (cm)",
                          10, -0.5*maxDim, 0.5*maxDim,
                          10, -0.5*maxDim, 0.5*maxDim,
                          10, -0.5*maxDim, 0.5*maxDim);
    frame->SetStats(0);
    frame->Draw();

    if (frame3D) frame3D->Draw("same");

    // optional wedge wireframe overlay

    if (ev.hasPath && ev.xP && ev.xP->size()>1) {
        auto pl = new TPolyLine3D((int)ev.xP->size());
        for (int i=0;i<(int)ev.xP->size();++i)
            pl->SetPoint(i, (*ev.xP)[i], (*ev.yP)[i], (*ev.zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");

        const int n = (int)ev.xP->size();
        auto pm = DrawStartEndMarkers3D((*ev.xP)[0],(*ev.yP)[0],(*ev.zP)[0],
                                        (*ev.xP)[n-1],(*ev.yP)[n-1],(*ev.zP)[n-1]);

        auto leg = new TLegend(0.62, 0.22, 0.99, 0.36);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(pm.start, "Start", "p");
        leg->AddEntry(pm.end, "End", "p");
        leg->Draw();
    }

    DrawEventInfoBox(ievt, ev.L, ev.W, ev.T, ev.wrap,
                     ev.reachedEnd, ev.pmt_side, ev.inPMT, ev.detected,
                     ev.absorbed, ev.escaped,
                     ev.endPlane, ev.nBounces, ev.pathLen);

    c->Modified();
    c->Update();
}

// -----------------------------
// 4 pads XY/XZ zoom-out + zoom-in (2D)

void DrawEvent4ViewFromTree(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    double zoomMargin = 0.12)
{
    TFile f(fn);
    if (f.IsZombie()) return;

    auto t = (TTree*)f.Get("tPhot");
    if (!t) return;

    auto g2 = (TGraph2D*)f.Get("scintGeometry");
    if (!g2) return;

    EventData ev;
    if (!LoadEvent(t, ievt, ev)) return;

    if (!ev.hasPath || !ev.xP || ev.xP->size() < 2) {
        cout << ev.hasPath << " " << ev.xP << " "
             << (ev.xP ? ev.xP->size() : 0) << "\n";
        cout << "DrawEvent4ViewFromTree: no stored path (cfg.savePath=true)\n";
        return;
    }

    // Build 2D geometry projections
    auto geomXY = BuildFrameGeometryXY(g2, true);
    auto geomXZ = BuildFrameGeometryXZ(g2, true);

    // Compute path bounds
    double xlo=(*ev.xP)[0], xhi=(*ev.xP)[0];
    double ylo=(*ev.yP)[0], yhi=(*ev.yP)[0];
    double zlo=(*ev.zP)[0], zhi=(*ev.zP)[0];

    for (size_t i=1;i<ev.xP->size();++i) {
        xlo = std::min(xlo, (double)(*ev.xP)[i]);
        xhi = std::max(xhi, (double)(*ev.xP)[i]);
        ylo = std::min(ylo, (double)(*ev.yP)[i]);
        yhi = std::max(yhi, (double)(*ev.yP)[i]);
        zlo = std::min(zlo, (double)(*ev.zP)[i]);
        zhi = std::max(zhi, (double)(*ev.zP)[i]);
    }

    double dx = xhi-xlo; if (dx<=0) dx=1;
    double dy = yhi-ylo; if (dy<=0) dy=1;
    double dz = zhi-zlo; if (dz<=0) dz=1;

    // Zoom-in bounds (correct clamping for your conventions)
    double zx1 = std::max(-(ev.L/2. + ev.wedgeLen), xlo - zoomMargin*dx);
    double zx2 = std::min(ev.L, xhi + zoomMargin*dx);
    double zy1 = std::max(-ev.W/2.0, ylo - zoomMargin*dy);
    double zy2 = std::min(+ev.W/2.0, yhi + zoomMargin*dy);
    double zz1 = std::max(-ev.T/2.0, zlo - zoomMargin*dz);
    double zz2 = std::min(+ev.T/2.0, zhi + zoomMargin*dz);

    // Zoom-out bounds
    double xo1 = -(ev.L/2.0+ev.wedgeLen) * 1.3, xo2 = (ev.L/2.0+ev.wedgeLen) * 1.3;
    double yo1 = -ev.W/2.0* 1.3, yo2 = +ev.W/2.0* 1.3;
    double zo1 = -ev.T/2.0* 1.3, zo2 = +ev.T/2.0* 1.3;
    cout<<"Zoom-out X:["<<xo1<<","<<xo2<<"] Y:["<<yo1<<","<<yo2<<"]\n";
    TString cname = Form("c4_%lld", (long long)ievt);
    if (auto old = gROOT->FindObject(cname)) delete old;

    auto c = new TCanvas(cname, "4-view (XY/XZ, zoom-out/zoom-in)", 1400, 1000);
    c->Divide(2,2,0.01,0.01);

    auto drawXY = [&](int pad, bool zoomed)
    {
        c->cd(pad);
        gPad->SetGrid(0,0);
        gPad->SetFixedAspectRatio();

        double xmin = zoomed ? zx1 : xo1;
        double xmax = zoomed ? zx2 : xo2;
        double ymin = zoomed ? zy1 : yo1;
        double ymax = zoomed ? zy2 : yo2;
        printf ("drawXY pad=%d zoomed=%d x:[%.2f,%.2f] y:[%.2f,%.2f]\n",
                pad, (int)zoomed, xmin, xmax, ymin, ymax);
        TString hname = Form("hxy_%lld_%d", (long long)ievt, pad);
        auto h = new TH2D(hname,
                          zoomed ? "XY zoom in; x (cm); y (cm)" : "XY zoom out; x (cm); y (cm)",
                          10, xmin, xmax, 10, ymin, ymax);
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->Draw();

        // Detector outline (optional)
        DrawOutlineXY(ev.L, ev.W);

        // Geometry projection
        if (geomXY) geomXY->Draw("same");

        // Path polyline
        auto pl = new TPolyLine((int)ev.xP->size());
        pl->SetLineWidth(3);
        for (int i=0;i<(int)ev.xP->size();++i)
            pl->SetPoint(i, (*ev.xP)[i], (*ev.yP)[i]);
        pl->Draw("same");

        // Start/end markers
        const int n = (int)ev.xP->size();
        auto mS = new TMarker((*ev.xP)[0], (*ev.yP)[0], 20);
        mS->SetMarkerColor(kGreen+2);
        mS->SetMarkerSize(1.3);
        mS->Draw("same");

        auto mE = new TMarker((*ev.xP)[n-1], (*ev.yP)[n-1], 29);
        mE->SetMarkerColor(kRed+1);
        mE->SetMarkerSize(1.4);
        mE->Draw("same");

        auto leg = new TLegend(0.12, 0.12, 0.45, 0.28);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);
        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(mS, "Start", "p");
        leg->AddEntry(mE, "End", "p");
        leg->Draw();

        if (pad == 1) {
            DrawInfoBox2D(ievt, ev.L, ev.W, ev.T, ev.wrap,
                          ev.reachedEnd, ev.pmt_side, ev.inPMT, ev.detected,
                          ev.absorbed, ev.escaped, ev.endPlane, ev.nBounces, ev.pathLen);
        }
    };

    auto drawXZ = [&](int pad, bool zoomed)
    {
        c->cd(pad);
        gPad->SetGrid(0,0);
        gPad->SetFixedAspectRatio();

        double xmin = zoomed ? zx1 : xo1;
        double xmax = zoomed ? zx2 : xo2;
        double zmin = zoomed ? zz1 : zo1;
        double zmax = zoomed ? zz2 : zo2;

        TString hname = Form("hxz_%lld_%d", (long long)ievt, pad);
        auto h = new TH2D(hname,
                          zoomed ? "XZ zoom in; x (cm); z (cm)" : "XZ zoom out; x (cm); z (cm)",
                          10, xmin, xmax, 10, zmin, zmax);
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->Draw();

        DrawOutlineXZ(ev.L, ev.T);

        // Geometry projection
        if (geomXZ) geomXZ->Draw("same");

        auto pl = new TPolyLine((int)ev.xP->size());
        pl->SetLineWidth(3);
        for (int i=0;i<(int)ev.xP->size();++i)
            pl->SetPoint(i, (*ev.xP)[i], (*ev.zP)[i]);
        pl->Draw("same");

        const int n = (int)ev.xP->size();
        auto mS = new TMarker((*ev.xP)[0], (*ev.zP)[0], 20);
        mS->SetMarkerColor(kGreen+2);
        mS->SetMarkerSize(1.3);
        mS->Draw("same");

        auto mE = new TMarker((*ev.xP)[n-1], (*ev.zP)[n-1], 29);
        mE->SetMarkerColor(kRed+1);
        mE->SetMarkerSize(1.4);
        mE->Draw("same");

        auto leg = new TLegend(0.12, 0.12, 0.45, 0.28);
        leg->SetFillStyle(0);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.03);
        leg->AddEntry(pl, "Photon path", "l");
        leg->AddEntry(mS, "Start", "p");
        leg->AddEntry(mE, "End", "p");
        leg->Draw();
    };

    // Pads:
    // 1: XY zoom out, 2: XY zoom in, 3: XZ zoom out, 4: XZ zoom in
    drawXY(1, false);
    drawXY(2, true);
    drawXZ(3, false);
    drawXZ(4, true);

    c->Modified();
    c->Update();
}

// -----------------------------
// 2 pads 3D zoom-out + zoom-in

void DrawEventSplitViewFromTree(
    const char *fn = "toyOptics.root",
    Long64_t ievt = 0,
    double zoomMargin = 0.12)
{
    TFile f(fn);
    if (f.IsZombie()) return;

    auto t = (TTree*)f.Get("tPhot");
    if (!t) return;

    auto g2 = (TGraph2D*)f.Get("scintGeometry");
    if (!g2) return;

    EventData ev;
    if (!LoadEvent(t, ievt, ev)) return;

    if (!ev.hasPath || !ev.xP || ev.xP->size() < 2) {
        cout << "DrawEventSplitViewFromTree: no stored path (cfg.savePath=true)\n";
        return;
    }

    auto frame3D = BuildFrameGeometry3D(g2, true);

    // Compute zoom bounds from path
    double xlo=(*ev.xP)[0], xhi=(*ev.xP)[0];
    double ylo=(*ev.yP)[0], yhi=(*ev.yP)[0];
    double zlo=(*ev.zP)[0], zhi=(*ev.zP)[0];

    for (size_t i=1;i<ev.xP->size();++i) {
        xlo = std::min(xlo, (double)(*ev.xP)[i]);
        xhi = std::max(xhi, (double)(*ev.xP)[i]);
        ylo = std::min(ylo, (double)(*ev.yP)[i]);
        yhi = std::max(yhi, (double)(*ev.yP)[i]);
        zlo = std::min(zlo, (double)(*ev.zP)[i]);
        zhi = std::max(zhi, (double)(*ev.zP)[i]);
    }

    double dx=xhi-xlo; if (dx<=0) dx=1;
    double dy=yhi-ylo; if (dy<=0) dy=1;
    double dz=zhi-zlo; if (dz<=0) dz=1;

    xlo = std::max(0.0, xlo - zoomMargin*dx);
    xhi = std::min(ev.L, xhi + zoomMargin*dx);
    ylo = std::max(-ev.W/2.0, ylo - zoomMargin*dy);
    yhi = std::min(+ev.W/2.0, yhi + zoomMargin*dy);
    zlo = std::max(-ev.T/2.0, zlo - zoomMargin*dz);
    zhi = std::min(+ev.T/2.0, zhi + zoomMargin*dz);

    // Zoomed-out cube bounds for nice 1:1:1
    double maxDim = std::max(ev.L, std::max(ev.W, ev.T));
    double xo_min = -0.5*maxDim, xo_max = 0.5*maxDim;
    double yo_min = -0.5*maxDim, yo_max = 0.5*maxDim;
    double zo_min = -0.5*maxDim, zo_max = 0.5*maxDim;

    TString cname = Form("cSplit_%lld", (long long)ievt);
    if (auto old = gROOT->FindObject(cname)) delete old;

    auto c = new TCanvas(cname, "3D split view", 1400, 1000);
    c->Divide(1,2);

    auto drawOnePad3D = [&](int ipad,
                            double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                            const char* padTitle)
    {
        c->cd(ipad);
        TString frname = Form("fr3_%lld_%d", (long long)ievt, ipad);
        auto frame = new TH3D(frname, Form("%s; x (cm); y (cm); z (cm)", padTitle),
                              10, xmin, xmax,
                              10, ymin, ymax,
                              10, zmin, zmax);
        frame->SetDirectory(nullptr);
        frame->SetStats(0);
        frame->Draw();

        if (frame3D) frame3D->Draw("same");

        auto pl = new TPolyLine3D((int)ev.xP->size());
        for (int i=0;i<(int)ev.xP->size();++i)
            pl->SetPoint(i, (*ev.xP)[i], (*ev.yP)[i], (*ev.zP)[i]);
        pl->SetLineWidth(3);
        pl->Draw("same");

        const int n = (int)ev.xP->size();
        auto pm = DrawStartEndMarkers3D((*ev.xP)[0],(*ev.yP)[0],(*ev.zP)[0],
                                        (*ev.xP)[n-1],(*ev.yP)[n-1],(*ev.zP)[n-1]);

        if (ipad == 1) {
            DrawEventInfoBox(ievt, ev.L, ev.W, ev.T, ev.wrap,
                             ev.reachedEnd, ev.pmt_side, ev.inPMT, ev.detected,
                             ev.absorbed, ev.escaped,
                             ev.endPlane, ev.nBounces, ev.pathLen);

            auto leg = new TLegend(0.62, 0.22, 0.99, 0.36);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.05);
            leg->AddEntry(pl, "Photon path", "l");
            leg->AddEntry(pm.start, "Start", "p");
            leg->AddEntry(pm.end, "End", "p");
            leg->Draw();
        }
    };

    drawOnePad3D(1, xo_min, xo_max, yo_min, yo_max, zo_min, zo_max, "Zoom out");
    drawOnePad3D(2, xlo, xhi, ylo, yhi, zlo, zhi, "Zoom in");

    c->Modified();
    c->Update();
}

// -----------------------------
// colz quick view

void DrawFracColz(const char *fn = "scan_summary.root", const char *tname = "tScan")
{
    TFile f(fn);
    auto t = (TTree *)f.Get(tname);
    if (!t) {
        std::cout << "no " << tname << "\n";
        return;
    }

    printf("N=%lld\n", t->GetEntries());

    int nx = 3, ny = 3;
    t->Draw(Form("y:x>>hDet(%d,0,90,%d,-15,15)", nx, ny), "frac_detected", "colz");

    auto hDet = (TH2 *)gROOT->FindObject("hDet");
    if (hDet) {
        hDet->SetTitle("Detected fraction; x (cm); y (cm)");
        gPad->Modified();
        gPad->Update();
    }
}

// -----------------------------
// quick counters + cfg scan

void PrintCountsFromTree(const char *fn = "toyOptics.root", const char *tn = "tPhot")
{
    TFile f(fn, "READ");
    if (f.IsZombie()) {
        std::cout << "Cannot open file: " << fn << "\n";
        return;
    }

    auto t = (TTree *)f.Get(tn);
    if (!t) {
        std::cout << "Cannot find tree: " << tn << "\n";
        return;
    }

    auto has = [&](const char *b) { return t->GetBranch(b) != nullptr; };

    const char *bAbs = has("absorbed") ? "absorbed" : nullptr;
    const char *bEsc = has("escaped") ? "escaped" : nullptr;
    const char *bHit = has("inPMT") ? "inPMT" : (has("hitPMT") ? "hitPMT" : nullptr);
    const char *bDet = has("detected") ? "detected" : nullptr;

    auto countOnes = [&](const char *bname) -> Long64_t
    {
        if (!bname) return 0;
        TString sel = Form("%s!=0", bname);
        return t->GetEntries(sel);
    };

    Long64_t nTot = t->GetEntries();
    Long64_t nAbs = countOnes(bAbs);
    Long64_t nEsc = countOnes(bEsc);
    Long64_t nHit = countOnes(bHit);
    Long64_t nDet = countOnes(bDet);

    std::cout << "File: " << fn << " | Tree: " << tn << "\n";
    std::cout << "Total: " << nTot << "\n";
    std::cout << "absorbed: " << nAbs << "\n";
    std::cout << "escaped:  " << nEsc << "\n";
    std::cout << "hitPMT:   " << nHit << "  (branch=" << (bHit ? bHit : "none") << ")\n";
    std::cout << "detected: " << nDet << "\n";

    // Optional cfg tree scan if it exists
    TTree *tCfg = (TTree *)f.Get("tCfg");
    if (!tCfg) return;

    OpticsConfig *cfg = nullptr;
    tCfg->SetBranchAddress("cfg", &cfg);
    tCfg->GetEntry(0);
    if (cfg) {
        printf("cfg: L=%f W=%f T=%f nScint=%f nOut=%f\n", cfg->L, cfg->W, cfg->T, cfg->nScint, cfg->nOut);
    }
    tCfg->Scan("cfg.L:cfg.W:cfg.T:cfg.nScint:cfg.nOut");
}
