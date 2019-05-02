/*

alienv load O2/latest
root -l  
gSystem->Load("libO2TPCFastTransformation");
.x IrregularSpline2D3DCalibratorTest.C++

*/

#include "TFile.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "algorithm"
#include <chrono>

#include "IrregularSpline2D3D.h"
#include "IrregularSpline2DCalibrator.h"


using namespace ali_tpc_common::tpc_fast_transformation ;

const int PolynomDegree = 20;

void F( float u, float v , float &x, float &y, float &z )
{  
  static double cu[PolynomDegree+1], cv[PolynomDegree+1];
  static int doInit=1;
  
  if( doInit ){
    gRandom->SetSeed( 1 );  
    for( int i=0; i<=PolynomDegree; i++ ){
      cu[i] = gRandom->Uniform(-1,1);
      cv[i] = gRandom->Uniform(-1,1);
    }
    doInit = 0;
  }
 
  u-=0.5;
  v-=0.5;
  double uu = 1.;  
  double vv = 1.;  
  x = 0;
  for( int i=0; i<=PolynomDegree; i++ ){
    x+=cu[i]*uu;
    x+=cv[i]*vv;
    uu*=u;
    vv*=v;
  }
  y = v;
  z = ((u-.5)*(u-.5));
}



int IrregularSpline2D3DCalibratorTest()
{

  //generating input Spline parameters & Spline & Grid------------------------------------------------------------------------------------

  cout<<"Splinegeneration Start"<<endl;

  IrregularSpline2DCalibrator finder;

  finder.nAxisTicksU = 70;
  finder.nAxisTicksV = 70;
  finder.u_ticks = 20;
  finder.v_ticks = 20;
  finder.nKnotsU = 40;
  finder.nKnotsV = 40;
  finder.max_tolerated_deviation = 0.5;  

  IrregularSpline2D3D spline;
  float * spline_data = finder.calibrateSpline(spline, F);

  int nKnots = spline.getNumberOfKnots();
 
 
  //Drawing Results------------------------------------------------------

  TCanvas *canv = new TCanvas("cQA","2D splines  QA",1000,1000);
  
  canv->Draw();
  //canv->Divide(3,3);
  canv->Update();

  // knots 

  TGraph2D *gknots = new TGraph2D();
  gknots->SetName("gknots");
  gknots->SetTitle("gknots");
  gknots->SetLineColor(kRed);
  gknots->SetMarkerSize(1.);
  gknots->SetMarkerStyle(8);
  gknots->SetMarkerColor(kRed);

  const IrregularSpline1D& gridU = spline.getGridU();
  const IrregularSpline1D& gridV = spline.getGridV();

  int nu = gridU.getNumberOfKnots();
  int gknotsN=0;
  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    double u = gridU.getKnot( i ).u;
    for( int j=0; j<gridV.getNumberOfKnots(); j++ ){
      double v = gridV.getKnot( j ).u;
      float fx, fy, fz;
      F(u,v,fx,fy,fz);
      gknots->SetPoint(gknotsN++,u,v,fx); 
    }
  }
  
  // ntuple with input and spline
  TNtuple *nt = new TNtuple("nt","nt","u:v:f0:fs");

  
  TGraph2D *gf0 = new TGraph2D();
  gf0->SetName("gf0");
  gf0->SetTitle("Input function");
  gf0->SetLineColor(kBlue);
  int gf0N=0;

  TGraph2D *gfs = new TGraph2D();
  gfs->SetName("gfs");
  gfs->SetTitle("Spline");
  gfs->SetLineColor(kRed);
  int gfsN=0;

  TH1F *qaX = new TH1F("qaX","qaX [um]",1000,-0.1,0.1); 
  TH1F *qaY = new TH1F("qaY","qaY [um]",1000,-0.1,0.1);
  TH1F *qaZ = new TH1F("qaZ","qaZ [um]",1000,-0.1,0.1);

  int iter=0;
  float stepu = 1.e-3;
  float stepv = 1.e-3;

  for( float u=0; u<=1; u+=stepu ){
    for( float v=0; v<=1; v+=stepv ){
      float fx0, fy0, fz0;
      F( u, v, fx0, fy0, fz0);
      float fx1, fy1, fz1;
      spline.getSplineVec( spline_data, u,v, fx1, fy1, fz1);
      if( u>=0 && v>=0 && u<=1 && v<=1 ){
	qaX->Fill( (fx1 - fx0) );
	qaY->Fill( (fy1 - fy0) );
	qaZ->Fill( (fz1 - fz0) );
	nt->Fill(u,v,fx0,fx1);
      }
      gf0->SetPoint(gf0N++,u,v,fx0);
      gfs->SetPoint(gfsN++,u,v,fx1);
    }
  }
  gStyle->SetPalette(1);
  gf0->Draw("surf");
  gfs->Draw("surf, same");
  gknots->Draw("P,same");  
  canv->Update();  

  /*
    Specific drawing options can be used to paint a TGraph2D:

    "TRI"  : The Delaunay triangles are drawn using filled area.
    An hidden surface drawing technique is used. The surface is
    painted with the current fill area color. The edges of each
    triangles are painted with the current line color.
    "TRIW" : The Delaunay triangles are drawn as wire frame
    "TRI1" : The Delaunay triangles are painted with color levels. The edges
    of each triangles are painted with the current line color.
    "TRI2" : the Delaunay triangles are painted with color levels.
    "P"    : Draw a marker at each vertex
    "P0"   : Draw a circle at each vertex. Each circle background is white.
    "PCOL" : Draw a marker at each vertex. The color of each marker is
    defined according to its Z position.
    "CONT" : Draw contours
    "LINE" : Draw a 3D polyline
  */
  return 0;
}
