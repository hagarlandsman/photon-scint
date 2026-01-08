#include <TPolyLine3D.h>
#include <TGraph2D.h>
#include "Rtypes.h"

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

 TPolyLine3D* BuildFrameGeometry3D(TGraph2D* g2, bool closeLoop=true );

void DrawFracColz(const char* , const char* );
void PrintCountsFromTree(const char* fname = "onePoint_100.root", const char* treename = "tPhot");
void SetViewEqualXYZ(double xmin, double xmax,
                     double ymin, double ymax,
                     double zmin, double zmax,
                     bool orthographic = true);
