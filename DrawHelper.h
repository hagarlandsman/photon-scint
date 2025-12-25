void DrawEventSplitViewFromTree(
    const char *,
    Long64_t,
    bool drawWedge ,
    double zoomMargin // padding for the zoomed-in view
);
void DrawEvent4ViewFromTree   (
    const char *,
    Long64_t,
    bool drawWedge ,
    double zoomMargin // padding for the zoomed-in view
);

void DrawFracColz(const char* , const char* );
void PrintCountsFromTree(const char* fname = "onePoint_100.root", const char* treename = "tPhot");
