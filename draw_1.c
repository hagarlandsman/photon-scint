{
    TTree t("t", "from csv");
    t.ReadFile("clean.csv", "Ntot/D:detf:err:absf:escf:hpmtf:cost:fdet_tot:cost_per:L:W:T:useWedge:wedgeL:wedgeT:nscint:absLen:Rwrap:wrap:rPMT", ',');
    /* TString file = "results.root";
     TFile f(file.Data(), "READ");
     TTree *t = (TTree *)f.Get("t");
 */

    TCut cut1 = "useWedge == 0 && T==1 &&  wrap ==1 && rPMT==1.27 ";

    Double_t ws[] = {10, 20, 30, 40, 50, 70, 90, 110};
    Double_t ls[] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130};
    TGraph *temp;
    TMultiGraph *mgW = new TMultiGraph();
    TMultiGraph *mgL = new TMultiGraph();
    TLegend *legW = new TLegend(0.65, 0.60, 0.88, 0.88);
    legW->SetBorderSize(0);
    legW->SetFillStyle(0);
    TLegend *legL = new TLegend(0.65, 0.60, 0.88, 0.88);
    legL->SetBorderSize(0);
    legL->SetFillStyle(0);
    const int nws = sizeof(ws) / sizeof(ws[0]);
    const int nls = sizeof(ls) / sizeof(ls[0]);
    for (int i = 0; i < nls; ++i)
    {
        TCut cutL = cut1 && Form("L==%g", ws[i]);
        Long64_t n = t.Draw("cost_per:W", cutL && cut1, "goff");
        if (n <= 0)
            continue;
        const double *x = t.GetV2(); // L
        const double *y = t.GetV1(); // NdetFrac
        TGraph *g = new TGraph((int)n, x, y);
        g->SetName(Form("g_L%g", ws[i]));

        g->Sort(); // sorts by X (L), good for drawing lines
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetLineColor(i + 1);
        g->SetMarkerColor(i + 1);
        mgL->Add(g, "LP");
        mgL->SetTitle("Cost per effective area vs W (by L); W (cm); Cost per effective area ($/cm^{2})");
        legL->AddEntry(g, Form("L = %g", ls[i]), "lp");
    }
    for (int i = 0; i < nws; ++i)
    {
        TCut cutL = cut1 && Form("W==%g", ws[i]);
        Long64_t n = t.Draw("cost_per:L", cutL && cut1, "goff");
        if (n <= 0)
            continue;
        const double *x = t.GetV2(); // L
        const double *y = t.GetV1(); // NdetFrac
        TGraph *g = new TGraph((int)n, x, y);
        g->SetName(Form("g_W%g", ws[i]));

        g->Sort(); // sorts by X (L), good for drawing lines
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetLineColor(i + 1);
        g->SetMarkerColor(i + 1);
        mgW->Add(g, "LP");
        mgW->SetTitle("Cost per effective area vs L (by W); L (cm); Cost per effective area ($/cm^{2})");
        legW->AddEntry(g, Form("W = %g", ws[i]), "lp");
    }

    TCanvas *c = new TCanvas("c", "NdetFrac vs L (by W)", 900, 700);
    c->Divide(2, 2);
    c->cd(1);
    mgL->Draw("Ac");
    legL->Draw();
    c->cd(2);
    mgW->Draw("Ac");
    legW->Draw();
    c->Modified();
    c->Update();
    ////////////////////////////////////////////////////////////////////
    TCut cut2 = "useWedge == 0 && T==1 &&  wrap ==1 && rPMT==1.27  && L==90";
    TCut cut3 = "useWedge == 1 && wedgeL==20 && T==1 &&  wrap ==1 && rPMT==1.27  && L==90";
    TCut cut4 = "useWedge == 1 && wedgeL==W/2 && T==1 &&  wrap ==1 && rPMT==1.27 && L==90";
    TCut cut22 = "useWedge == 0 && T==1 &&  wrap ==2 && rPMT==1.27  && L==90";
    TCut cut33 = "useWedge == 1 && wedgeL==20 && T==2 &&  wrap ==1 && rPMT==1.27  && L==90";
    TCut cut44 = "useWedge == 1 && wedgeL==W/2 && T==2 &&  wrap ==1 && rPMT==1.27 && L==90";

    TMultiGraph *mgW2 = new TMultiGraph();
    TLegend *legW2 = new TLegend(0.65, 0.60, 0.88, 0.88);
    legW2->SetBorderSize(0);
    legW2->SetFillStyle(0);
    mgW2->SetTitle("Cost per effective area vs W (by wedge option); W (cm); Cost per effective area ($/cm^{2})");
    Long64_t n;
    n = t.Draw("cost_per:W", cut2, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g2 = new TGraph((int)n, x, y);
        g2->SetName("g2");
        g2->Sort(); // sorts by X (L), good for drawing lines
        g2->SetLineWidth(2);
        g2->SetMarkerStyle(20);
        g2->SetLineColor(1);
        g2->SetMarkerColor(1);
        mgW2->Add(g2, "LP");
        legW2->AddEntry(g2, "no wedge, PTFE wrap=1", "lp");
    }
  n = t.Draw("cost_per:W", cut22, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g22 = new TGraph((int)n, x, y);
        g22->SetName("g22");
        g22->Sort(); // sorts by X (L), good for drawing lines
        g22->SetLineWidth(2);
        g22->SetLineStyle(2);
        g22->SetMarkerStyle(20);
        g22->SetLineColor(1);
        g22->SetMarkerColor(1);
        mgW2->Add(g22, "LP");
        legW2->AddEntry(g22, "no wedge,  Mylar wrap=2", "lp");
    }


     n = t.Draw("cost_per:W", cut3, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g3 = new TGraph((int)n, x, y);
        g3->SetName("g3");
        g3->Sort(); // sorts by X (L), good for drawing lines
        g3->SetLineWidth(2);
        g3->SetMarkerStyle(20);
        g3->SetLineColor(2);
        g3->SetMarkerColor(2);
        legW2->AddEntry(g3, "wedgeL=20, PTFE wrap=1", "lp");

        mgW2->Add(g3, "LP");
    }
    n = t.Draw("cost_per:W", cut33, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g33 = new TGraph((int)n, x, y);
        g33->SetName("g33");
        g33->Sort(); // sorts by X (L), good for drawing lines
        g33->SetLineWidth(2);
        g33->SetLineStyle(2);
        g33->SetMarkerStyle(20);
        g33->SetLineColor(2);
        g33->SetMarkerColor(2);
        legW2->AddEntry(g33, "wedgeL=20, Mylar wrap=2 ", "lp");

        mgW2->Add(g33, "LP");
    }
    n = t.Draw("cost_per:W", cut4, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g4 = new TGraph((int)n, x, y);
        g4->SetName("g4");
        g4->Sort(); // sorts by X (L), good for drawing lines
        g4->SetLineWidth(2);
        g4->SetMarkerStyle(20);
        g4->SetLineColor(4);
        g4->SetMarkerColor(4);
        mgW2->Add(g4, "LP");
        legW2->AddEntry(g4, "wedgeL=W/2, PTFE wrap=1", "lp");
    }

n = t.Draw("cost_per:W", cut44, "goff");
    if (n > 0)
    {
        double *x = t.GetV2(); // L
        double *y = t.GetV1(); // NdetFrac
        TGraph *g44 = new TGraph((int)n, x, y);
        g44->SetName("g44");
        g44->Sort(); // sorts by X (L), good for drawing lines
        g44->SetLineWidth(2);
                g44->SetLineStyle(2);

        g44->SetMarkerStyle(20);
        g44->SetLineColor(4);
        g44->SetMarkerColor(4);
        mgW2->Add(g44, "LP");
        legW2->AddEntry(g44 , "wedgeL=W/2, Mylar wrap=2", "lp");
    }
    c->cd(3);
    mgW2->Draw("Ac");
legW2->Draw();
c->Modified();
c->Update();

    TCut cutT1 = "useWedge == 0 && T==1 &&  wrap ==1 && rPMT==1.27 ";
    TCut cutT2 = "useWedge == 0 && T==2 &&  wrap ==1 && rPMT==1.27 ";

    TMultiGraph *mgL2 = new TMultiGraph();
    TLegend *legL2 = new TLegend(0.65, 0.60, 0.88, 0.88);
    legL2->SetBorderSize(0);
    legL2->SetFillStyle(0);
    for (int i = 0; i < nls; ++i)
    {
        TCut cutL = cutT1 && Form("L==%g", ws[i]);
         n = t.Draw("detf:W", cutL , "goff");
        if (n <= 0)
            continue;
        const double *x = t.GetV2(); // L
        const double *y = t.GetV1(); // NdetFrac
        TGraph *g = new TGraph((int)n, x, y);
        g->SetName(Form("g_L%g", ws[i]));

        g->Sort(); // sorts by X (L), good for drawing lines
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetLineColor(i + 1);
        g->SetMarkerColor(i + 1);
        mgL2->Add(g, "LP");
        legL2->AddEntry(g, Form("L = %g, T=1", ls[i]), "lp");

        TCut cutL2 = cutT2 && Form("L==%g", ws[i]);
         n = t.Draw("detf:W", cutL2  , "goff");
        if (n <= 0) {
            cout<<"skip L2="<<ls[i]<<"\n";
            continue;
        }
         cout<<"draw L2="<<ls[i]<<"\n";
        const double *x2= t.GetV2(); // L
        const double *y2 = t.GetV1(); // NdetFrac
        TGraph *g2 = new TGraph((int)n, x2, y2);
        g2->SetName(Form("g2_L%g, T=2", ws[i]));

        g2->Sort(); // sorts by X (L), good for drawing lines
        g2->SetLineWidth(2);
        g2->SetLineStyle(2);
        g2->SetMarkerStyle(20);
        g2->SetLineColor(i + 1);
        g2->SetMarkerColor(i + 1);
        mgL2->Add(g2, "LP");
        legL2->AddEntry(g2, Form("L = %g,  T=2", ls[i]), "lp");

    }

            mgL2->SetTitle("Cost per effective area vs W (by L, dashed T=2,); W (cm); Cost per effective area ($/cm^{2})");

    c->cd(4);
    mgL2->Draw("Ac");
    legL2   ->Draw();
    c->Modified();
    c->Update();
    // mg->Draw("A");
}
