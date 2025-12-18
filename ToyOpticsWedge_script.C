.L ToyOpticsWedge.C+
// No wedges, save path and draw
OpticsConfig cfg;
cfg.savePath = true;
cfg.useWedge = false;      // no wedges
RunManySites(
  1,        // Nsites
  1,        // NphotPerSite (one photon total)
  120, 20, 1,   // L,W,T
  1,        // wrap (1=PTFE)
  cfg,
  "nowedge.root"
);
DrawEventFromTree("nowedge.root", 0, false);
////////////////////////////////////////////////////////////////////////////////
// Yes wedges, save path and draw
.L ToyOpticsWedge.C+

OpticsConfig cfg2;
cfg2.savePath = true;
cfg2.useWedge = true;      // wedges on
cfg2.wedgeLen  = 20.0;
cfg2.wedgeTipW = 2.54;
cfg2.rPMT = 2.54; // cm
cfg2.eps0 = 0; // kills all non-direct PMT hits
RunManySites(
  1,
  1000,
  cfg2,
  "wedge.root"
);
// No zoom
DrawEventFromTree("wedge.root", 0, true, false);
DrawEventSplitViewFromTree("wedge.root", 1, true, 20.0, 5.0, 0.10);
DrawEvent4ViewFromTree("wedge.root",   1, true,  20.0, 5.0, 0.12);
// Zoom in tighter:
DrawEventFromTree("wedge.root", 0, true, 20.0, 5.0, true, 0.05);

////////////////////////////////////////////////////////////////////////////////
// Yes wedges, save path and draw
OpticsConfig cfg2;
cfg2.savePath = true;
cfg2.useWedge = true;      // wedges on
cfg2.wedgeLen  = 20.0;
cfg2.wedgeTipW = 5.0;
RunManySites(
  1,
  1,
  120, 20, 1,
  1,
  cfg2,
  "wedge.root"
);

OpticsConfig cfg;
cfg.savePath = true;
cfg.useWedge = false;      // no wedges
RunManySites(
  1,        // Nsites
  1,        // NphotPerSite (one photon total)
  120, 20, 1,   // L,W,T
  1,        // wrap (1=PTFE)
  cfg,
  "nowedge.root"
);
DrawEventFromTree("nowedge.root", 0, false);
DrawEventFromTree("nowedge.root", 0, true, 20.0, 5.0, false);
