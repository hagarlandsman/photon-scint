// runFracMap.C
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TString.h"

#include <iostream>

// Forward declaration so this file compiles
double ToyScintOptics(
  TRandom3& rng,
  long long NperPoint,
  double L,
  double W,
  double T,
  int wrap,
  double nScint,
  double nOut,
  double absLen,
  double Rwrap,
  int maxSteps,
  double rPMT,
  double epsCouple,
  double pde,
  bool writeTree,
  const char* outFile
);
/*
TFile f("toyOptics.root");

TTree* t = (TTree*)f.Get("tPhot");
t->Draw("yf:zf","inPMT==0 && W==30 && L==90")
*/
void anaFile() {
  // Build frac(W,L) = N(detected)/N(total)
// Sanity: does the file exist and is it non-trivial size?
  TString outFile = TString(gSystem->WorkingDirectory()) + "/toyOptics.root";
Long_t id=0, flags=0, modtime=0;
Long64_t sz=0;
int st = gSystem->GetPathInfo(outFile.Data(), &id, &sz, &flags, &modtime);
if (st != 0) {
  std::cout << "ERROR: output file not found: " << outFile << "\n";
  return;
}
std::cout << "File size [bytes]: " << sz << "\n";


  // Open and read tree
  TFile f(outFile, "READ");
  if (f.IsZombie()) {
    std::cout << "ERROR: cannot open " << outFile << "\n";
    return;
  }

  TTree* t = (TTree*)f.Get("tPhot");
  if (!t) {
    std::cout << "ERROR: cannot find tPhot in " << outFile << "\n";
    f.ls();
    return;
  }

  std::cout << "tPhot entries: " << t->GetEntries() << "\n";


  // Build frac map: N(detected)/N(total)
  gROOT->cd();
  TH2D* hNum = new TH2D("hNum","N(detected);W (cm);L (cm)", 20, 10, 200, 20, 10, 200);
  TH2D* hDen = new TH2D("hDen","N(total);W (cm);L (cm)",    20, 10, 200, 20, 10, 200);

  // IMPORTANT: use the existing histograms (no name collisions)
  t->Draw("L:W>>+hNum", "detected==1 && inPMT==1", "goff");
  t->Draw("L:W>>+hDen", "",            "goff");

  std::cout << "hDen integral: " << hDen->Integral() << "\n";
  std::cout << "hNum integral: " << hNum->Integral() << "\n";

  TH2D* hEff = (TH2D*)hNum->Clone("hEff");
  hEff->SetTitle("frac(W,L) = N(detected)/N(total);W (cm);L (cm)");
  hEff->Divide(hDen);

  // Draw
  TCanvas* c = new TCanvas("cFrac","frac map", 900, 700);
  hEff->Draw("COLZ");
  c->Update();


}
#include "TPolyLine3D.h"

static TPolyLine3D* MakeEdge(double x1,double y1,double z1, double x2,double y2,double z2)
{
  auto e = new TPolyLine3D(2);
  e->SetPoint(0, x1,y1,z1);
  e->SetPoint(1, x2,y2,z2);
  e->SetLineWidth(2);
  return e;
}

static void DrawBarWireframe(double L, double W, double T)
{
  const double y0 = -W/2.0, y1 = +W/2.0;
  const double z0 = -T/2.0, z1 = +T/2.0;
  const double x0 = 0.0,   x1p = L;

  // 4 edges at x = 0
  MakeEdge(x0,y0,z0, x0,y1,z0)->Draw("same");
  MakeEdge(x0,y1,z0, x0,y1,z1)->Draw("same");
  MakeEdge(x0,y1,z1, x0,y0,z1)->Draw("same");
  MakeEdge(x0,y0,z1, x0,y0,z0)->Draw("same");

  // 4 edges at x = L
  MakeEdge(x1p,y0,z0, x1p,y1,z0)->Draw("same");
  MakeEdge(x1p,y1,z0, x1p,y1,z1)->Draw("same");
  MakeEdge(x1p,y1,z1, x1p,y0,z1)->Draw("same");
  MakeEdge(x1p,y0,z1, x1p,y0,z0)->Draw("same");

  // 4 edges connecting x=0 to x=L
  MakeEdge(x0,y0,z0, x1p,y0,z0)->Draw("same");
  MakeEdge(x0,y1,z0, x1p,y1,z0)->Draw("same");
  MakeEdge(x0,y1,z1, x1p,y1,z1)->Draw("same");
  MakeEdge(x0,y0,z1, x1p,y0,z1)->Draw("same");
}

void DrawEvent(const char* fn="toyOptics.root", Long64_t ievt=0)
{
  TFile f(fn);
  auto t = (TTree*)f.Get("tPhot");
gROOT->cd();
  std::vector<float>* x=0;
  std::vector<float>* y=0;
  std::vector<float>* z=0;
  int endPlane=-1;
  int pmt_side=0;
  int inPMT=0;
  int detected=0;
  int absorbed=0;
  int escaped=0;
  int nBounces=0;
  int reachedEnd=0;
  double L=0,W=0,T=0;
  int wrap = 0 ;
  t->SetBranchAddress("L",&L);
  t->SetBranchAddress("W",&W);
  t->SetBranchAddress("T",&T);
 t->SetBranchAddress("wrap",&wrap);

  t->SetBranchAddress("xPath",&x);
  t->SetBranchAddress("yPath",&y);
  t->SetBranchAddress("zPath",&z);

  t->SetBranchAddress("endPlane",&endPlane);
  t->SetBranchAddress("pmt_side",&pmt_side);
  t->SetBranchAddress("inPMT",&inPMT);
  t->SetBranchAddress("detected",&detected);
  t->SetBranchAddress("absorbed",&absorbed);
  t->SetBranchAddress("escaped",&escaped);
  t->SetBranchAddress("reachedEnd",&reachedEnd);
  t->SetBranchAddress("nBounces",&nBounces);

  t->GetEntry(ievt);

  auto c = new TCanvas("c","path",900,700);
   double maxDim = std::max(L, std::max(W, T));

  // Center the bar in a cube in y,z, and pad x symmetrically around [0,L]
  double xPad = 0.5 * (maxDim - L);
  double xMin = 0.0 - xPad;
  double xMax = L   + xPad;

  double yMin = -0.5 * maxDim;
  double yMax = +0.5 * maxDim;

  double zMin = -0.5 * maxDim;
  double zMax = +0.5 * maxDim;
  TString Title("");
  if (reachedEnd) {
    Title="Reached end, ";
    if (inPMT) Title+="Hit PMT, ";
        else Title+="Didn't hit PMT, ";
    if (detected) Title+="Detected, ";
        if (!detected) Title+="Not detected, ";


  }
  else Title="Didn't reach end, ";

  if (absorbed) Title+="Absorbed";
  if (escaped) Title += "Escaped";
  if (wrap==1) Title += "(PTFE)";
  else if (wrap==2) Title += "(Mylar)";
  else  Title+="(unknown wrap)";
  auto frame = new TH3D("fr",Form("Photon %lld path: %s, bounces=%d ; x (cm); y (cm); z (cm)",ievt,Title.Data(),nBounces),
                        10, xMin, xMax,
                        10, yMin, yMax,
                        10, zMin, zMax);
  frame->SetStats(0);
  frame->Draw();


  DrawBarWireframe(L, W, T);

  auto pl = new TPolyLine3D((int)x->size());
  for (int i=0;i<(int)x->size();i++) {
    pl->SetPoint(i, (*x)[i], (*y)[i], (*z)[i]);
    printf ("%f\t%f\t%f\n",(*x)[i], (*y)[i], (*z)[i]);
  }
  pl->Draw("same");
}


void run1() {
  gROOT->ProcessLine(".L ScanToyOptics_withTree.C+");
  TRandom3 rng(0);

  const double W = 30.0;
  const double L = 90.0;
  const long long N = 100000;
  const double T = 1.0;
  const double nScint = 1.58;
  const double nOut = 1.0;
  const double absLen = 300.0;
  const double Rwrap = 0.95;
  const int maxSteps = 2000;
  const double rPMT = 0.5;
  const double epsCouple = 0.90;
  const double pde = 0.20;
  TString outFile("one.root");
  int wrap = 1;
  ToyScintOptics(rng, N, L, W, T, wrap,
                     nScint, nOut, absLen, Rwrap, maxSteps,
                     rPMT, epsCouple, pde,
                     true, outFile.Data());

  printf ("TFile f(\"%s\", \"READ\"); \n"
    "TTree* t = (TTree*)f.Get(\"tPhot\"); \n"
    "gROOT->cd(); \n"
    "printf (\"Entries=%%lld\\n\",t->GetEntries()); \n"
    "t->Draw(\"yf:zf\",\"inPMT==0\")\n",
    outFile.Data());

}


void runFracMap() {

  // Load/compile the macro that DEFINES ToyScintOptics
  // Adjust path/filename if needed
  gROOT->ProcessLine(".L ScanToyOptics_withTree.C+");
  TString outFile = TString(gSystem->WorkingDirectory()) + "/toyOptics_WLscan.root";
  // Start clean
  gSystem->Unlink(outFile);
  TRandom3 rng(0);

  const long long N = 10000;
  const double T = 1.0;
  const double nScint = 1.58;
  const double nOut = 1.0;
  const double absLen = 300.0;
  const double Rwrap = 0.95;
  const int maxSteps = 2000;
  const double rPMT = 2.5;
  const double epsCouple = 0.90;
  const double pde = 0.20;

  // Choose a few values
  //double Ws[] = {10,20, 30, 40, 50,60,70,80, 90,100,110,120,130,140,150,160,170,180,190,200};
  //double Ls[] = {10,20, 30, 40, 50,60,70,80, 90,100,110,120,130,140,150,160,170,180,190,200};
  double Ws[] = {10,20, 30, 40, 50,60,70,80, 90,100};
  double Ls[] = {10,20, 30, 40, 50,60,70,80, 90,100};

  /*
   for (int iW = 0; iW < (int)(sizeof(Ws)/sizeof(double)); iW++) {
    for (int iL = 0; iL < (int)(sizeof(Ls)/sizeof(double)); iL++) {
      double W = Ws[iW];
      double L = Ls[iL]; */
 int wrap = 1;  // PTFE  wrap = 1, Mylar wrap = 2

 for (double W=20; W<=200; W=W+10){
 for (double L=20; L<=200; L=L+10){
      TString fname=Form("scintOpt_W%d_L%d_T%d_wrap%d.root",int(W),int(L),int(T),wrap);
      std::cout << "Running L=" << L << " W=" << W << "\n";

      ToyScintOptics(rng, N, L, W, T, wrap,
                     nScint, nOut, absLen, Rwrap, maxSteps,
                     rPMT, epsCouple, pde,
                     true, fname.Data());
    }
  }


}
/*
Check with thicker scintillaotr.
Check with different wrap.
Check with different light collection profile.
Check normalization. maybe also divide by W*L the rate - no need
HERE
*/
/*
// ana1:
TChain ch("tPhot");
ch.Add("scintOpt_W*");
TH2D* hNum = new TH2D("hNum","N(Reached ends);W (cm);L (cm)", 20, 20, 210, 20, 20, 210);
TH2D* hPMT1 = (TH2D*)hNum->Clone("hPMT1"); hPMT1->Reset("ICES");  hPMT1->SetTitle("N(Reached_PMT)")
TH2D* hDen = (TH2D*)hNum->Clone("hDen"); hDen->Reset("ICES");
TH2D* hCost1 = (TH2D*)hNum->Clone("hCost1"); hCost1->Reset("ICES"); hCost1->SetTitle("cost per cm^2 with 1 PMT");
TH2D* hCost2 = (TH2D*)hNum->Clone("hCost2"); hCost2->Reset("ICES"); hCost2->SetTitle("cost per cm^2 with 2 pmts")
TH2D* hFrac = (TH2D*)hNum->Clone("hFrac"); hFrac->Reset("ICES"); hFrac->SetTitle("cost per cm^2 with 1 PMT");
const double cost_pmt = 1269.0;
const double cost_cm2 = 0.4;


for (int ix=1; ix<hCost1->GetNbinsX(); ix++){
    for (int iy=1; iy<hCost1->GetNbinsY(); iy++) {
        double X = hCost1->GetXaxis()->GetBinCenter(ix);
        double Y = hCost1->GetYaxis()->GetBinCenter(iy);
        hCost2->Fill(X,Y, (X*Y*cost_cm2 + 2*cost_pmt)/X/Y);

        hCost1->Fill(X,Y, (X*Y*cost_cm2 + cost_pmt)/X/Y);
    }
}


  // IMPORTANT: use the existing histograms (no name collisions)
  ch.Draw("L:W>>hNum", "reachedEnd==1");
   ch.Draw("L:W>>hPMT1", "inPMT==1");
 ch.Draw("L:W>>hDen", "");
  hNum->Divide(hDen);
  hPMT1->Divide(hDen);


  TH2D* hEff2 = (TH2D*)hNum->Clone("hEff2");
  hEff2->SetTitle("frac(W,L) = N(reachedEnd)/N(total);W (cm);L (cm)");
  hEff2->Divide(hDen);
  TH2D* hEff1 = (TH2D*)hEff2->Clone("hEff2");
  hEff1->Scale(0.5);

  TH2D* hEffcost1 = (TH2D*)hCost1->Clone("hEffcost1");
  hEffcost1->Divide(hEff1);

  // Draw
  TCanvas* c = new TCanvas("cFrac","frac map", 900, 700);
  hEff->Draw("COLZ");
  c->Update();

==1", "goff");
  ch.Draw("L:W>>+hDen", "",            "goff");

  std::cout << "hDen integral: " << hDen->Integral() << "\n";
  std::cout << "hNum integral: " << hNum->Integral() << "\n";

  TH2D* hEff = (TH2D*)hNum->Clone("hEff");
  hEff->SetTitle("frac(W,L) = N(detected)/N(total);W (cm);L (cm)");
  hEff->Divide(hDen);

  // Draw
  TCanvas* c = new TCanvas("cFrac","frac map", 900, 700);
  hEff->Draw("COLZ");
  c->Update();


*/