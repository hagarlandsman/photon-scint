// scanDetectedVsArea.C
#include <vector>

{
  cout << "hi" << endl;
  TString dir = "out";
  TSystemDirectory d("outdir", dir);
  TList *files = d.GetListOfFiles();
  if (!files)
  {
    Error("scan_scint_files", "Cannot list files in directory: %s", dir.Data());
    return;
  }
  std::vector<int> Ws;
  std::vector<int> Ls;
  std::vector<int> Ns;

  printf("statrt \n");
  TIter next(files);
  while (auto obj = next())
  {
    auto f = dynamic_cast<TSystemFile *>(obj);
    if (!f)
      continue;

    TString fname = f->GetName();
    if (f->IsDirectory())
      continue;
    if (!fname.EndsWith("_summary.root"))
      continue;
    if (!fname.BeginsWith("scint_"))
      continue;

    int L = 0, W = 0, N = 0;
    int c = 0, d = 0;
    int nmatched = sscanf(fname.Data(),
                          "scint_%d_%d_1_nowedge_PTFE_N_%d_%d_%dsummary.root",
                          &L, &W, &N, &c, &d);
    if (nmatched != 3)
      continue;

    std::cout << "OK: " << dir << "/" << fname
              << " -> L=" << L << " W=" << W << " N=" << N << "\n";
    Ws.push_back(W);
    Ls.push_back(L);
    Ns.push_back(N);
  }
  TGraph *gr_detected_30 = new TGraph();
  gr_detected_30->SetTitle("Integral(h_detected) 30; L; Integral");
  gr_detected_30->SetMarkerStyle(20);
  gr_detected_30->SetMarkerColor(kRed);
  gr_detected_30->SetLineColor(kRed);
  TGraph *gr_dt_30 = new TGraph();
  gr_dt_30->SetTitle("dT( error ) 30; L; dt error");
  gr_dt_30->SetMarkerStyle(20);
  gr_dt_30->SetMarkerColor(kRed);
  gr_dt_30->SetLineColor(kRed);
   TGraph *gr_dtmax_30 = new TGraph();
  gr_dtmax_30->SetTitle("dT( error ) 30; L; dt error");
  gr_dtmax_30->SetMarkerStyle(20);
  gr_dtmax_30->SetMarkerColor(kRed);
  gr_dtmax_30->SetLineColor(kRed);
  TGraph *gr_detected_50 = new TGraph();
  gr_detected_50->SetTitle("Integral(h_detected) 50; L; Integral");
  gr_detected_50->SetMarkerStyle(20);
  gr_detected_50->SetMarkerColor(kBlue);
  gr_detected_50->SetLineColor(kBlue);
  TGraph *gr_dt_50 = new TGraph();
  gr_dt_50->SetTitle("dT( error ) 50; L; dt error");
  gr_dt_50->SetMarkerStyle(20);
  gr_dt_50->SetMarkerColor(kBlue);
  gr_dt_50->SetLineColor(kBlue);
   TGraph *gr_dtmax_50 = new TGraph();
  gr_dtmax_50->SetTitle("dT( error ) 50; L; dt error");
  gr_dtmax_50->SetMarkerStyle(20);
  gr_dtmax_50->SetMarkerColor(kBlue);
  gr_dtmax_50->SetLineColor(kBlue);
  TGraph *gr_detected_20 = new TGraph();
  gr_detected_20->SetTitle("Integral(h_detected) 20; L; Integral");
  gr_detected_20->SetMarkerStyle(20);
  gr_detected_20->SetMarkerColor(kGreen);
  gr_detected_20->SetLineColor(kGreen);
  TGraph *gr_dt_20 = new TGraph();
  gr_dt_20->SetTitle("dT( error ) 20; L; dt error");
  gr_dt_20->SetMarkerStyle(20);
  gr_dt_20->SetMarkerColor(kGreen);
  gr_dt_20->SetLineColor(kGreen);
 TGraph *gr_dtmax_20 = new TGraph();
  gr_dtmax_20->SetTitle("dT( error ) 20; L; dt error");
  gr_dtmax_20->SetMarkerStyle(20);
  gr_dtmax_20->SetMarkerColor(kGreen);
  gr_dtmax_20->SetLineColor(kGreen);

  TGraph *gr_detected_10 = new TGraph();
  gr_detected_10->SetTitle("Integral(h_detected) 10; L; Integral");
  gr_detected_10->SetMarkerStyle(20);
  gr_detected_10->SetMarkerColor(kOrange);
  gr_detected_10->SetLineColor(kOrange);

  TGraph *gr_dt_10 = new TGraph();
  gr_dt_10->SetTitle("dT( error ) 10; L; dt error");
  gr_dt_10->SetMarkerStyle(20);
  gr_dt_10->SetMarkerColor(kOrange);
  gr_dt_10->SetLineColor(kOrange);

 TGraph *gr_dtmax_10 = new TGraph();
  gr_dtmax_10->SetTitle("dT( error ) 10; L; dt error");
  gr_dtmax_10->SetMarkerStyle(20);
  gr_dtmax_10->SetMarkerColor(kOrange);
  gr_dtmax_10->SetLineColor(kOrange);

  TGraph *gr_escaped = new TGraph();
  gr_escaped->SetTitle("Integral(h_escaped); L; Integral");
  gr_escaped->SetMarkerStyle(20);
  gr_escaped->SetMarkerColor(kBlue);
  TGraph *gr_absorbed = new TGraph();
  gr_absorbed->SetTitle("Integral(h_absorbed); L; Integral");
  gr_absorbed->SetMarkerStyle(20);
  gr_absorbed->SetMarkerColor(kGreen);
  TGraph *gr_hitPMT = new TGraph();
  gr_hitPMT->SetTitle("Integral(h_hitPMT); L; Integral");
  gr_hitPMT->SetMarkerStyle(20);
  gr_hitPMT->SetMarkerColor(kBlack);

  // std::vector<int> Ws = {30,30,30 , 30 , 30 ,30 ,50,50 ,50 ,50 ,20,50};//, 60}; // 60}  // put your L values here
  // std::vector<int> Ls = {50,90,110, 120, 150,200,50,100,150,200,90,90};//, 100 }; // 100}  // matching W values here (same length as Ls)
  TH2F *h2 = new TH2F("h2", "Integral(h_detected); L; W", 20, 0, 200, 10, 0, 100);

  int ip = 0;
  double cost_pmt = 1269;
  double cost_cm2 = 0.4;
  for (size_t i = 0; i < Ls.size() + 1; i++)
  {
    int L = Ls[i];
    int W = Ws[i];
    int N = Ns[i];
    int area = L * W;
    TString fn = Form("out/scint_%d_%d_1_nowedge_PTFE_N_%d_summary.root", L, W, N);
    printf ("processing file %s \n", fn.Data());
    TFile f(fn, "READ");
    if (f.IsZombie())
    {
      printf("Cannot open %s\n", fn.Data());
      continue;
    }

    TH2 *h_absorbed = (TH2 *)f.Get("h_absorbed_s");
    TH2 *h_detected = (TH2 *)f.Get("h_detected_s");
    TH2 *h_escaped = (TH2 *)f.Get("h_escaped_s");
    TH2 *h_hitPMT = (TH2 *)f.Get("h_hitPMT_s");
    if (!h_detected)
    {
      printf("Missing h_detected in %s\n", fn.Data());
      continue;
    }

    double cost1 = W * L * cost_cm2 + cost_pmt;

    double cost2 = W * L * cost_cm2 + 2 * cost_pmt;


    int bx = h2->GetXaxis()->FindBin(L);
    int by = h2->GetYaxis()->FindBin(W);
    h2->SetBinContent(bx, by, h_detected->Integral() / cost2);

    TH2 *h_pathMean_p1 = (TH2 *)f.Get("h_pathMean_det_p1");
    TH2 *h_pathMean_p2 = (TH2 *)f.Get("h_pathMean_det_p2");
    TH2 *h_pathRMS_p1 = (TH2 *)f.Get("h_pathRMS_det_p1");
    TH2 *h_pathRMS_p2 = (TH2 *)f.Get("h_pathRMS_det_p2");

    TH2D* h_sigDT = (TH2D *)h_pathRMS_p1->Clone("h_sigDT");
    h_sigDT->Reset();
    h_sigDT->SetTitle("sigma(Delta t);x;y;ns");

    int nx = h_sigDT->GetNbinsX();
    int ny = h_sigDT->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix)
    {
      for (int iy = 1; iy <= ny; ++iy)
      {
        double sL1 = h_pathRMS_p1->GetBinContent(ix, iy);
        double sL2 = h_pathRMS_p2->GetBinContent(ix, iy);

        if (sL1 <= 0 || sL2 <= 0)
        {
          h_sigDT->SetBinContent(ix, iy, 0);
          continue;
        }

        double sDT = std::sqrt(sL1 * sL1 + sL2 * sL2); // in cm, divide by v for ns where v = 29.9792458 / n_group;
        h_sigDT->SetBinContent(ix, iy, sDT);
      }
    }

    double sumW = 0.0;
    double v = 1;  // v= 29.9792458 / 1.58; // cm/ns
    // For E[Var|bin]
    double sumW_var = 0.0;

    // For Var(E[dt|bin]) = E[mu^2] - (E[mu])^2
    double sumW_mu = 0.0;
    double sumW_mu2 = 0.0;

     nx = h_pathMean_p1->GetNbinsX();
     ny = h_pathMean_p1->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix)
    {
      for (int iy = 1; iy <= ny; ++iy)
      {

        double w = 1; //(hW ? hW->GetBinContent(ix, iy) : 1.0);
        if (w <= 0)
          continue;

        double muL1 = h_pathMean_p1->GetBinContent(ix, iy);
        double muL2 = h_pathMean_p2->GetBinContent(ix, iy);
        double sL1 = h_pathRMS_p1->GetBinContent(ix, iy);
        double sL2 = h_pathRMS_p2->GetBinContent(ix, iy);

        if (sL1 <= 0 || sL2 <= 0)
          continue;

        double muDT = (muL1 - muL2) / v;                        // ns
        double varDT_given = (sL1 * sL1 + sL2 * sL2) / (v * v); // ns^2

        sumW += w;
        sumW_var += w * varDT_given;
        sumW_mu += w * muDT;
        sumW_mu2 += w * muDT * muDT;
      }
    }

    double E_var = sumW_var / sumW;
    double E_mu = sumW_mu / sumW;
    double Var_mu = (sumW_mu2 / sumW) - (E_mu * E_mu);

    double sigma_global = std::sqrt(E_var + Var_mu);



    printf("Global sigma(Delta t) over (x,y) = %.4f ns\n", sigma_global);

    h_pathMean_p1->Add(h_pathMean_p2, -1);
    double m=h_pathMean_p1->GetMaximum();

      if (W == 50) {
      gr_dtmax_50->SetPoint(gr_dtmax_50->GetN(), L, m);
      gr_dt_50->SetPoint(gr_dt_50->GetN(), L, sigma_global);
      gr_detected_50->SetPoint(gr_detected_50->GetN(), L, h_detected->Integral("width") / cost2);
    }
    if (W == 30) {
      gr_dtmax_30->SetPoint(gr_dtmax_30->GetN(), L, m);
      gr_dt_30->SetPoint(gr_dt_30->GetN(), L, sigma_global);
      gr_detected_30->SetPoint(gr_detected_30->GetN(), L, h_detected->Integral("width") / cost2);
    }
    if (W == 20) {
      gr_dtmax_20->SetPoint(gr_dtmax_20->GetN(), L, m);
      gr_dt_20->SetPoint(gr_dt_20->GetN(), L, sigma_global);
      gr_detected_20->SetPoint(gr_detected_20->GetN(), L, h_detected->Integral("width") / cost2);
    }
    if (W == 10) {
      gr_dtmax_10->SetPoint(gr_dtmax_10->GetN(), L, m);

      gr_dt_10->SetPoint(gr_dt_10->GetN(), L, sigma_global);
      gr_detected_10->SetPoint(gr_detected_10->GetN(), L, h_detected->Integral("width") / cost2);
    }


  }

  auto mg = new TMultiGraph();
  mg->SetTitle("Cost yield vs L; L; events/cost");
  mg->Add(gr_detected_10);
  mg->Add(gr_detected_20);
  mg->Add(gr_detected_30);
  mg->Add(gr_detected_50);
  gr_detected_30->Sort();
  gr_detected_50->Sort();
  gr_detected_20->Sort();
  gr_detected_10->Sort();

  TCanvas *c = new TCanvas("c", "c", 900, 650);
  mg->Draw("AlP");
  c->Update();

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(gr_detected_10, "detected 10", "p");
  leg->AddEntry(gr_detected_20, "detected 20", "p");
  leg->AddEntry(gr_detected_30, "detected 30", "p");

  leg->AddEntry(gr_detected_50, "detected 50", "p");
  leg->Draw();
  TCanvas *c2 = new TCanvas("c2", "c2", 900, 650);
  h2->Draw("COLZ");

    auto mg2 = new TMultiGraph();
    mg2->SetTitle("overall dt error vs L; L; overall error (cm)");
  mg2->Add(gr_dt_10);
  mg2->Add(gr_dt_20);
  mg2->Add(gr_dt_30);
  mg2->Add(gr_dt_50);
  gr_dt_30->Sort();
  gr_dt_50->Sort();
  gr_dt_20->Sort();
  gr_dt_10->Sort();
  TCanvas *c3 = new TCanvas("c3", "c3", 900, 650);
  mg2->Draw("AlP");
  c3->Update();
  leg->Draw();

  auto mg3 = new TMultiGraph();
  mg3->SetTitle("Max dt vs L; L; max dt (cm)");
  mg3->Add(gr_dtmax_10);
  mg3->Add(gr_dtmax_20);
  mg3->Add(gr_dtmax_30);
  mg3->Add(gr_dtmax_50);
  gr_dtmax_30->Sort();
  gr_dtmax_50->Sort();
  gr_dtmax_20->Sort();
  gr_dtmax_10->Sort();
  TCanvas *c4 = new TCanvas("c4", "c4", 900, 650);
  mg3->Draw("AlP");
  c4->Update();
  leg->Draw();
c2->SaveAs("c.png");
c2->SaveAs("c2.png");
c3->SaveAs("c3.png");
c4->SaveAs("c4.png");

  TFile out("detected_vs_area.root", "RECREATE");
  gr_dtmax_30->Write("gr_dtmax_30");
  gr_dtmax_20->Write("gr_dtmax_20");
  gr_dtmax_10->Write("gr_dtmax_10");
  gr_dtmax_50->Write("gr_dtmax_50");
  gr_dt_10->Write("gr_dt_10");
  gr_dt_20->Write("gr_dt_20");
  gr_dt_30->Write("gr_dt_30");
  gr_dt_50->Write("gr_dt_50");
  gr_detected_10->Write("gr_detected_10");
  gr_detected_20->Write("gr_detected_20");
  gr_detected_30->Write("gr_detected_30");
  gr_detected_50->Write("gr_detected_50");

  h2->Write("h2");
  out.Close();
}
