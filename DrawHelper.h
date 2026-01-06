void DrawEventSplitViewFromTree(
    const char *,
    Long64_t,
    double zoomMargin // padding for the zoomed-in view
);
void DrawEvent4ViewFromTree   (
    const char *,
    Long64_t,
    double zoomMargin // padding for the zoomed-in view
);
   void DrawEventFromTree3D(
    const char *,
    Long64_t
);

void DrawFracColz(const char* , const char* );
void PrintCountsFromTree(const char* fname = "onePoint_100.root", const char* treename = "tPhot");
