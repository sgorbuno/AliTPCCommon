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
#include <Vc/Vc>
#include <chrono>

#include "IrregularSpline2D3D.h"
#include "IrregularSpline2DCalibrator.h"


using namespace ali_tpc_common::tpc_fast_transformation ;
const int PolynomDegree = 20;

double cu[PolynomDegree+1], cv[PolynomDegree+1];

void F( float u, float v , float &x, float &y, float &z){
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

float Fx( float u, float v )
{
  u-=0.5;
  v-=0.5;
  double uu = 1.;  
  double vv = 1.;  
  double f = 0;
  for( int i=0; i<=PolynomDegree; i++ ){
    f+=cu[i]*uu;
    f+=cv[i]*vv;
    uu*=u;
    vv*=v;
  }
  return (f);
}


float Fy( float u, float v ){  return (v); }
float Fz( float u, float v ){  return ((u-.5)*(u-.5)); }


float * generate_input_spline (int &nKnotsU, int &nKnotsV, int &nAxisTicksU, int &nAxisTicksV, IrregularSpline2D3D &spline){
  
  float du = 1./(nKnotsU-1);
  float dv = 1./(nKnotsV-1);

  float knotsU[nKnotsU], knotsV[nKnotsV];
  knotsU[0] = 0;
  knotsU[nKnotsU-1] = 1;
  knotsV[0] = 0;
  knotsV[nKnotsV-1] = 1;

  for( int i=1; i<nKnotsU-1; i++ ){
    knotsU[i] = i*du + gRandom->Uniform(-du/3,du/3);
  }

  for( int i=1; i<nKnotsV-1; i++ ){
    knotsV[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }
  spline.construct( nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV );
 
  int nKnotsTot = spline.getNumberOfKnots();
 
  const IrregularSpline1D &gridU = spline.getGridU();
  const IrregularSpline1D &gridV = spline.getGridV();      

  float *data0 = new float[ 3*nKnotsTot ];
  float *data = new float[ 3*nKnotsTot ];

  
  int nu = gridU.getNumberOfKnots();
  float neg_pos = 0.5 , pos_neg = 0.5;


  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    double u = gridU.getKnot( i ).u ;
    if(i%5 == 0){
      if(neg_pos < 0){
	neg_pos = abs(neg_pos);
      }
      else{
	neg_pos = -neg_pos;
      }
    }
    for( int j=0; j<gridV.getNumberOfKnots(); j++ ){ 
      double v =  gridV.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      if (i%5 == 0){
	data0[ind+0] = Fx(u,v)+neg_pos;
      }
      else{
	if (j%5 == 0){
	  data0[ind+0] = Fx(u,v)+pos_neg;
	  if(pos_neg < 0){
	    pos_neg = abs(pos_neg);
	  }
	  else {
	    pos_neg = -pos_neg;
	  }
	}
	else{
	  data0[ind+0] = Fx(u,v);
	}
      }
      data0[ind+1] = Fy(u,v);
      data0[ind+2] = Fz(u,v);
    }
  }
  

  for( int i=0; i<3*nKnotsTot; i++ ){
    data[i] = data0[i];
  }

  spline.correctEdges( data );
  return data;
}


float *  generate_compare_spline (int &nKnotsU, int &nKnotsV, int &nAxisTicksU, int &nAxisTicksV, IrregularSpline2D3D &test_spline, IrregularSpline2D3D &spline, float *input_data){


  float du = 1./(nKnotsU-1);
  float dv = 1./(nKnotsV-1);

  float knotsU[nKnotsU], knotsV[nKnotsV];
  knotsU[0] = 0;
  knotsU[nKnotsU-1] = 1;
  knotsV[0] = 0;
  knotsV[nKnotsV-1] = 1;

  for( int i=1; i<nKnotsU-1; i++ ){
    knotsU[i] = i*du + gRandom->Uniform(-du/3,du/3);
  }

  for( int i=1; i<nKnotsV-1; i++ ){
    knotsV[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }

  
  
  test_spline.construct(nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);

  int nKnotsTot = test_spline.getNumberOfKnots();
 
  const IrregularSpline1D &gridU = test_spline.getGridU();
  const IrregularSpline1D &gridV = test_spline.getGridV();      

  float *data0 = new float[ 3*nKnotsTot ];
  float *data = new float[ 3*nKnotsTot ];

  int nu = gridU.getNumberOfKnots();
  
    
  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    double u = gridU.getKnot( i ).u ;
    for( int j=0; j<gridV.getNumberOfKnots(); j++ ){  
      double v =  gridV.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, u,v,x0,y0,z0);
      data0[ind+0] = x0;
      data0[ind+1] = y0;
      data0[ind+2] = z0;
    }
  }
  

  for( int i=0; i<3*nKnotsTot; i++ ){
    data[i] = data0[i];
  }


  test_spline.correctEdges( data );
  return data;

}


float * generate_spline_u (IrregularSpline2D3D &spline, IrregularSpline2D3D &test_spline, IrregularSpline2D3D &correct_spline, float *input_data, float *test_data, double max_tolerated_deviation, int nAxisTicksU, int nAxisTicksV, int u_ticks, int v_ticks){
  const IrregularSpline1D &gridU = test_spline.getGridU();
  const IrregularSpline1D &gridV = test_spline.getGridV();
  
  float knotsU[gridU.getNumberOfKnots()*u_ticks];
  int nKnotsU = 0; 
  double u_small = 0;
  double u_high = 0;
  int nKnotsV = gridV.getNumberOfKnots();
  float knotsV[nKnotsV];
  for( int k=0; k<nKnotsV; k++ ){
    knotsV[k] = gridV.getKnot(k).u ;
  }


  for( int i=0; i<(gridU.getNumberOfKnots()-1); i++ ){
    u_small = gridU.getKnot( i ).u ;
    u_high = gridU.getKnot( (i+1) ).u ;
    double difference_u = u_high - u_small;
    float du = difference_u/(u_ticks+2);
    knotsU[nKnotsU] = u_small;
    nKnotsU++;
    for( int k=1; k<u_ticks; k++ ){
      bool value_changed = false;
      double current_u = du*k + u_small;
      for( int j=0; j<(gridV.getNumberOfKnots()); j++ ){
	if(!value_changed){
	  double v = gridV.getKnot( j ).u;
	  float x0, y0, z0;
	  spline.getSplineVec(input_data, current_u,v,x0,y0,z0);
	  float x, y, z;
	  test_spline.getSplineVec( test_data, current_u,v,x,y,z);
	  if ((abs(x-x0)) >= max_tolerated_deviation){
	    value_changed = true;
	    knotsU[nKnotsU] = current_u;
	    nKnotsU++;
	  }
	}
      }
    }
  }
  knotsU[nKnotsU] = u_high;
  nKnotsU++;
  
  correct_spline.construct(nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);
  
  int nKnotsTot = correct_spline.getNumberOfKnots();
 
  const IrregularSpline1D &correct_gridU = correct_spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_spline.getGridV();
  
  float *correct_data0 = new float[ 3*nKnotsTot ];
  float *correct_data = new float[ 3*nKnotsTot ];

  int nu = correct_gridU.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU.getKnot( i ).u ;
    for( int j=0; j<correct_gridV.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV.getKnot( j ).u ;
      int correct_ind = (nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0[correct_ind+0] = x0;
      correct_data0[correct_ind+1] = y0;
      correct_data0[correct_ind+2] = z0; 
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    correct_data[i] = correct_data0[i];
  }


  correct_spline.correctEdges( correct_data );
  return correct_data;
}

float * generate_spline_u_vec (IrregularSpline2D3D &spline, IrregularSpline2D3D &test_spline, IrregularSpline2D3D &correct_spline, float *input_data, float *test_data, int &nKnotsV, double max_tolerated_deviation, int nAxisTicksU, int nAxisTicksV, int u_ticks){
  const IrregularSpline1D &gridU = test_spline.getGridU();
  const IrregularSpline1D &gridV = test_spline.getGridV();
  
  /*
  
  correct_spline.construct(nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);
  
  int nKnotsTot = correct_spline.getNumberOfKnots();
 
  const IrregularSpline1D &correct_gridU = correct_spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_spline.getGridV();
  
  float *correct_data0 = new float[ 3*nKnotsTot ];
  float *correct_data = new float[ 3*nKnotsTot ];

  int nu = correct_gridU.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU.getKnot( i ).u ;
    for( int j=0; j<correct_gridV.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV.getKnot( j ).u ;
      int correct_ind = (nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0[correct_ind+0] = x0;
      correct_data0[correct_ind+1] = y0;
      correct_data0[correct_ind+2] = z0; 
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    correct_data[i] = correct_data0[i];
  }

  correct_spline.correctEdges( correct_data );
  return correct_data;
  */
}


float * generate_spline_uv(IrregularSpline2D3D &spline, IrregularSpline2D3D &correct_spline, IrregularSpline2D3D &correct_spline_two, float *input_data, float *correct_data, double max_tolerated_deviation,  int nAxisTicksU, int nAxisTicksV, int v_ticks){
  const IrregularSpline1D &correct_gridU = correct_spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_spline.getGridV();
  
  int nKnotsU = correct_gridU.getNumberOfKnots();
  //int v_ticks  = 60;
  float knotsV[correct_gridV.getNumberOfKnots()*v_ticks];
  int nKnotsV = 0;
  double v_high = 0;
  double v_small = 0;
  float knotsU[correct_gridU.getNumberOfKnots()];
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    knotsU[i] = correct_gridU.getKnot( i ).u ;
  }
  
  for( int i=0; i<(correct_gridV.getNumberOfKnots()-1); i++ ){
    v_small = correct_gridV.getKnot( i ).u ;
    v_high = correct_gridV.getKnot( (i+1) ).u ;
    double difference_v = v_high - v_small;
    float dv = difference_v/(v_ticks+2);
    knotsV[nKnotsV] = v_small;
    nKnotsV++;
    for( int k=1; k<v_ticks; k++ ){
      bool value_changed = false;
      double current_v = dv*k + v_small;
      for( int j=0; j<correct_gridU.getNumberOfKnots(); j++ ){
	if(!value_changed){
	  double u = correct_gridU.getKnot( j ).u ;
	  float x0, y0, z0;
	  spline.getSplineVec( input_data, u,current_v,x0,y0,z0);
	  float x, y, z;
	  correct_spline.getSplineVec( correct_data,u,current_v,x,y,z);
	  if ((abs(x-x0)) >= max_tolerated_deviation){
	    value_changed = true;
	    knotsV[nKnotsV] = current_v;
	    nKnotsV++;
	  }
	}
      }
    }
  }
  knotsV[nKnotsV] = v_high;
  nKnotsV++;

  correct_spline_two.construct( nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV );
 
  int correct_nKnotsTot = correct_spline_two.getNumberOfKnots();
  const IrregularSpline1D &correct_gridU_two = correct_spline_two.getGridU();
  const IrregularSpline1D &correct_gridV_two = correct_spline_two.getGridV();

  float *correct_data0_two = new float[ 3*correct_nKnotsTot ];
  float *correct_data_two = new float[ 3*correct_nKnotsTot ];

  int correct_nu = correct_gridU_two.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU_two.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU_two.getKnot( i ).u ;
    for( int j=0; j<correct_gridV_two.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV_two.getKnot( j ).u ;
      int correct_ind = (correct_nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0_two[correct_ind+0] = x0;
      correct_data0_two[correct_ind+1] = y0;
      correct_data0_two[correct_ind+2] = z0; 
    }
  }
  for( int i=0; i<3*correct_nKnotsTot; i++ ){
    correct_data_two[i] = correct_data0_two[i];
  }

  correct_spline_two.correctEdges( correct_data_two );


  return correct_data_two;
}

float * generate_2spline_u (IrregularSpline2D3D &spline, IrregularSpline2D3D &correct_2spline, float *input_data, int &nKnotsU_start, int &nKnotsV, double max_tolerated_deviation, int nAxisTicksU, int nAxisTicksV, int u_ticks){
  const IrregularSpline1D &gridU = spline.getGridU();
  const IrregularSpline1D &gridV = spline.getGridV();


  float du = 1./(nKnotsU_start-1);
  float dv = 1./(nKnotsV-1);

  float knotsU_start[nKnotsU_start], knotsV[nKnotsV];
  knotsU_start[0] = 0;
  knotsU_start[nKnotsU_start-1] = 1;
  knotsV[0] = 0;
  knotsV[nKnotsV-1] = 1;

  for( int i=1; i<nKnotsU_start-1; i++ ){
    knotsU_start[i] = i*du + gRandom->Uniform(-du/3,du/3);
  }

  for( int i=1; i<nKnotsV-1; i++ ){
    knotsV[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }
  

  float knotsU[nKnotsU_start*u_ticks];
  int nKnotsU = 0; 
  double u_small = 0;
  double u_high = 0;
  
  for( int i=0; i<(nKnotsU_start-1); i++ ){
    u_small = knotsU_start[i] ;
    u_high = knotsU_start[i+1] ;
    double difference_u = u_high - u_small;
    du = difference_u/(u_ticks+2);
    knotsU[nKnotsU] = u_small;
    nKnotsU++;
    for( int k=1; k<u_ticks; k++ ){
      bool value_changed = false;
      double current_u = du*k + u_small;
      for( int j=0; j<(nKnotsV); j++ ){
	double v = knotsV[j] ;
	float x, y, z;
	spline.getSplineVec(input_data, u_small,v,x,y,z);
	float x0, y0, z0;
	spline.getSplineVec(input_data, current_u,v,x0,y0,z0);
	if ((abs(x-x0)) >= max_tolerated_deviation && !value_changed){
	  value_changed = true;
	  knotsU[nKnotsU] = current_u;
	  nKnotsU++;
	}
      }
    }
  }
  knotsU[nKnotsU] = u_high;
  nKnotsU++;


  

  correct_2spline.construct(nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);
  
  int nKnotsTot = correct_2spline.getNumberOfKnots();
 
  const IrregularSpline1D &correct_gridU = correct_2spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_2spline.getGridV();
  

  
  float *correct_data0 = new float[ 3*nKnotsTot ];
  float *correct_data = new float[ 3*nKnotsTot ];

  int nu = correct_gridU.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU.getKnot( i ).u ;
    for( int j=0; j<correct_gridV.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV.getKnot( j ).u ;
      int correct_ind = (nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0[correct_ind+0] = x0;
      correct_data0[correct_ind+1] = y0;
      correct_data0[correct_ind+2] = z0; 
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    correct_data[i] = correct_data0[i];
  }

  correct_2spline.correctEdges( correct_data );
  return correct_data;
}


float * generate_2spline_uv(IrregularSpline2D3D &spline,IrregularSpline2D3D &correct_spline, IrregularSpline2D3D &correct_2spline_two, float *input_data, float *correct_data, double &max_tolerated_deviation,  int nAxisTicksU, int nAxisTicksV, int v_ticks, int nKnotsU_start, int nKnotsV_start){
  const IrregularSpline1D &gridU = correct_spline.getGridU();
  const IrregularSpline1D &gridV = spline.getGridV();

  
  float dv = 1./(nKnotsV_start-1);
  float knotsV_start[nKnotsV_start];
  knotsV_start[0] = 0;
  knotsV_start[nKnotsV_start-1] = 1;
  for( int i=1; i<nKnotsV_start-1; i++ ){
    knotsV_start[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }

  float du = 1./(nKnotsU_start-1);
  float knotsU_start[nKnotsU_start];
  knotsU_start[0] = 0;
  knotsU_start[nKnotsU_start-1] = 1;
  for( int i=1; i<nKnotsU_start-1; i++ ){
    knotsU_start[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }

  
  int nKnotsU = gridU.getNumberOfKnots();
  float knotsV[nKnotsV_start*v_ticks];
  int nKnotsV = 0;
  double v_high = 0;
  double v_small = 0;
  float knotsU[gridU.getNumberOfKnots()];
  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    knotsU[i] = gridU.getKnot( i ).u ;
  }


  for( int i=0; i<(nKnotsV_start-1); i++ ){
    v_small = knotsV_start[i] ;
    v_high = knotsV_start[i+1] ;
    double difference_v = v_high - v_small;
    float dv = difference_v/(v_ticks+2);
    knotsV[nKnotsV] = v_small;
    nKnotsV++;
    for( int k=1; k<v_ticks; k++ ){
      bool value_changed = false;
      double current_v = dv*k + v_small;
      for( int j=0; j<(nKnotsU_start); j++ ){
	double u = knotsU_start[i];
 	float x0, y0, z0;
 	spline.getSplineVec( input_data, u,current_v,x0,y0,z0);
 	float x, y, z;
 	correct_spline.getSplineVec( correct_data,u,current_v,x,y,z);
 	if ((abs(x-x0)) >= max_tolerated_deviation && !value_changed){
 	  value_changed = true;
 	  knotsV[nKnotsV] = current_v;
 	  nKnotsV++;
 	}
      }
    }
  }
  knotsV[nKnotsV] = v_high;
  nKnotsV++;

  correct_2spline_two.construct( nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV );
 
  int correct_nKnotsTot = correct_2spline_two.getNumberOfKnots();
  const IrregularSpline1D &correct_gridU = correct_2spline_two.getGridU();
  const IrregularSpline1D &correct_gridV = correct_2spline_two.getGridV();

  float *correct_data0_two = new float[ 3*correct_nKnotsTot ];
  float *correct_data_two = new float[ 3*correct_nKnotsTot ];


  int correct_nu = correct_gridU.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU.getKnot( i ).u ;
    for( int j=0; j<correct_gridV.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV.getKnot( j ).u ;
      int correct_ind = (correct_nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0_two[correct_ind+0] = x0;
      correct_data0_two[correct_ind+1] = y0;
      correct_data0_two[correct_ind+2] = z0; 
    }
  }
  for( int i=0; i<3*correct_nKnotsTot; i++ ){
    correct_data_two[i] = correct_data0_two[i];
  }

  correct_2spline_two.correctEdges( correct_data_two );

  return correct_data_two;
}


float * generate_1spline (IrregularSpline2D3D &spline, IrregularSpline2D3D &correct_spline, float *input_data, int &nKnotsU_start, int &nKnotsV_start, double max_tolerated_deviation, int nAxisTicksU, int nAxisTicksV, int u_ticks, int v_ticks){
  const IrregularSpline1D &gridU = spline.getGridU();
  const IrregularSpline1D &gridV = spline.getGridV();


  float du = 1./(nKnotsU_start-1);
  float dv = 1./(nKnotsV_start-1);

  float knotsU_start[nKnotsU_start], knotsV_start[nKnotsV_start];
  knotsU_start[0] = 0;
  knotsU_start[nKnotsU_start-1] = 1;
  knotsV_start[0] = 0;
  knotsV_start[nKnotsV_start-1] = 1;

  for( int i=1; i<nKnotsU_start-1; i++ ){
    knotsU_start[i] = i*du + gRandom->Uniform(-du/3,du/3);
  }

  for( int i=1; i<nKnotsV_start-1; i++ ){
    knotsV_start[i] = i*dv + gRandom->Uniform(-dv/3,dv/3);
  }
  
  float knotsU[nKnotsU_start*u_ticks], knotsV[nKnotsV_start*v_ticks];
  int nKnotsU = 0, nKnotsV = 0; 
  double u_small = 0, u_high = 0;
  double v_small = 0, v_high = 0;
  
  for( int i=0; i<(nKnotsU_start-1); i++ ){
    u_small = knotsU_start[i] ;
    u_high = knotsU_start[i+1] ;
    double difference_u = u_high - u_small;
    du = difference_u/(u_ticks+2);
    knotsU[nKnotsU] = u_small;
    nKnotsU++;
    for( int k=1; k<u_ticks; k++ ){
      bool value_changed = false;
      double current_u = du*k + u_small;
      for( int j=0; j<(nKnotsV_start-1); j++ ){
	if(!value_changed){
	  double v = knotsV_start[j];
	  float x, y, z;
	  spline.getSplineVec(input_data, u_small,v,x,y,z);
	  float x0, y0, z0;
	  spline.getSplineVec(input_data, current_u,v,x0,y0,z0);
	  if ((abs(x-x0)) >= max_tolerated_deviation){
	    value_changed = true;
	    knotsU[nKnotsU] = current_u;
	    nKnotsU++;
	  }
	}
      }
    }
  }
  knotsU[nKnotsU] = u_high;
  nKnotsU++;


  for( int i=0; i<(nKnotsV_start-1); i++ ){
    v_small = knotsV_start[i] ;
    v_high = knotsV_start[i+1] ;
    double difference_v = v_high - v_small;
    dv = difference_v/(v_ticks+2);
    knotsV[nKnotsV] = v_small;
    nKnotsV++;
    for( int k=1; k<v_ticks; k++ ){
      bool value_changed = false;
      double current_v = dv*k + v_small;
      for( int j=0; j<(nKnotsU_start-1); j++ ){
	double u = knotsU_start[j];
	float x, y, z;
	spline.getSplineVec(input_data, u,v_small,x,y,z);
	float x0, y0, z0;
	spline.getSplineVec(input_data,u,current_v,x0,y0,z0);
	if ((abs(x-x0)) >= max_tolerated_deviation && !value_changed){
	  value_changed = true;
	  knotsV[nKnotsV] = current_v;
	  nKnotsV++;
	}
      }
    }
  }
  knotsV[nKnotsV] = v_high;
  nKnotsV++;


  

  correct_spline.construct(nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);
  
  int nKnotsTot = correct_spline.getNumberOfKnots();
 
  const IrregularSpline1D &correct_gridU = correct_spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_spline.getGridV();
  
  float *correct_data0 = new float[ 3*nKnotsTot ];
  float *correct_data = new float[ 3*nKnotsTot ];

  int nu = correct_gridU.getNumberOfKnots();
    
  for( int i=0; i<correct_gridU.getNumberOfKnots(); i++ ){
    double correct_u = correct_gridU.getKnot( i ).u ;
    for( int j=0; j<correct_gridV.getNumberOfKnots(); j++ ){  
      double correct_v =  correct_gridV.getKnot( j ).u ;
      int correct_ind = (nu*j+i)*3;
      float x0, y0, z0;
      spline.getSplineVec( input_data, correct_u,correct_v,x0,y0,z0);
      correct_data0[correct_ind+0] = x0;
      correct_data0[correct_ind+1] = y0;
      correct_data0[correct_ind+2] = z0; 
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    correct_data[i] = correct_data0[i];
  }

  correct_spline.correctEdges( correct_data );
  return correct_data;
}



int IrregularSpline2D3DCalibratorTest()
{
  //GENERATIC GENERIC PARAMETERS
  //using namespace ali_tpc_common::tpc_fast_transformation ;

  gRandom->SetSeed(0);
  UInt_t seed = 1;//gRandom->Integer(100000); // 605
  gRandom->SetSeed(seed);
  //cout<<"Random seed: "<<seed<<" "<<gRandom->GetSeed()<<endl;
  
  for( int i=0; i<=PolynomDegree; i++ ){
    cu[i] = gRandom->Uniform(-1,1);
    cv[i] = gRandom->Uniform(-1,1);
    //cout<<"Random Number check "<<cu[i]<<"Check2 "<<cv[i]<<endl;
  }

  //generating input Spline parameters & Spline & Grid------------------------------------------------------------------------------------

  cout<<"Splinegeneration Start"<<endl;

  IrregularSpline2DCalibrator finder;
  int nAxisTicksU = 70, nAxisTicksV = 70;
  int u_ticks = 20, v_ticks = 20;
  int nKnotsU = 40, nKnotsV = 40;
  double max_tolerated_deviation = 0.5;
  float max_deviation = 0.5;


  IrregularSpline2D3D input_spline;
float * input_data = finder.generate_input_spline(input_spline, F);
  const IrregularSpline1D &input_gridU = input_spline.getGridU();
  const IrregularSpline1D &input_gridV = input_spline.getGridV();
    
  nKnotsU = 20;
  nKnotsV = 20;

  IrregularSpline2D3D spline;
  float * data = finder.generate_compare_spline(spline, input_spline, input_data);
  const IrregularSpline1D &gridU = spline.getGridU();
  const IrregularSpline1D &gridV = spline.getGridV();

  cout<<"Input Spline done"<<endl;
  

  //generating result with non vectorized code 3Splines------------------------------------------------------------------------------------

  auto start_compare = std::chrono::high_resolution_clock::now();
  cout<<"Starting non-vectorized version 3 Splines"<<endl;
  nKnotsU = 5;
  nKnotsV = 5;

  IrregularSpline2D3D test_spline;
  float * test_data = finder.generate_compare_spline(test_spline, spline, data);
  const IrregularSpline1D &test_gridU = test_spline.getGridU();
  const IrregularSpline1D &test_gridV = test_spline.getGridV();

  auto finish_compare = std::chrono::high_resolution_clock::now();
  float time_compare = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_compare-start_compare).count()*0.00001; 
  //std::cout << "Compare  Spline Done(non vectorized 3 Splines): " << std::chrono::duration_cast<std::chrono::nanoseconds>(compare_spline-start).count()*0.00001 << " Millisekunden"<< std::endl;

  auto start_normal = std::chrono::high_resolution_clock::now();
  IrregularSpline2D3D correct_spline;
  float * correct_data = finder.generate_spline_u(spline, test_spline, correct_spline, data, test_data);
  const IrregularSpline1D &correct_gridU = correct_spline.getGridU();
  const IrregularSpline1D &correct_gridV = correct_spline.getGridV();

  auto finish_normal_u = std::chrono::high_resolution_clock::now();
  //std::cout << "First Spline Done(non vectorized 3 Splines): " << std::chrono::duration_cast<std::chrono::nanoseconds>(first_spline-start).count()*0.00001 << " Millisekunden"<< std::endl;

IrregularSpline2D3D correct_spline_two;
//float * correct_data_two = finder.generate_spline_uv(spline, correct_spline, correct_spline_two, data, correct_data);
float * correct_data_two = finder.calibrateSpline(correct_spline_two, F);


  int totalknots_3 = correct_spline_two.getNumberOfKnots();
 
  auto finish_normal = std::chrono::high_resolution_clock::now();
  float time_normal_u = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_normal_u - start_normal).count()*0.00001; 
  float time_normal = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_normal - start_normal).count()*0.00001; 
  //std::cout << "Second Spline Done(non vectorized 3 Splines): " << time_3_spline << " Millisekunden"<< std::endl;
  // cout<<"Finished non-vectorized code 3 Splines"<<endl;
 

   //generating result with vectorized code 3Splines------------------------------------------------------------------------------------
   /*
  auto start_vector = std::chrono::high_resolution_clock::now();
  float knotsU_final[test_gridU.getNumberOfKnots()*u_ticks];
  int nKnotsU_final = 0;
   float knotsV[nKnotsV];
   Vc::float_m mask;
   Vc::float_m tmp_mask;
   Vc::float_v tmpv_x0;
   Vc::float_v tmpv_x;
   Vc::float_v tmpv_result;
   Vc::float_v deviation;

   deviation[0] = max_deviation;
   deviation[1] = max_deviation;
   deviation[2] = max_deviation;
   deviation[3] = max_deviation;

   for( int k=0; k<nKnotsV; k++ ){
     knotsV[k] = test_gridV.getKnot(k).u ;
   }
   
   for( int i=0; i<(test_gridU.getNumberOfKnots()-1); i++ ){
     double u_small = test_gridU.getKnot(i).u;
     knotsU_final[nKnotsU_final] = u_small;
     nKnotsU_final++;
     double offset = (test_gridU.getKnot(i+1).u - u_small)/(u_ticks+2);
     for( int k=1; k<u_ticks; k= k+4 ){
       double u_1 = u_small+(offset*k);
       double u_2 = u_1+offset;
       double u_3 = u_2+offset;
       double u_4 = u_3+offset;

       mask[0] = false;
       mask[1] = false;
       mask[2] = false;
       mask[3] = false;
       
       for( int j=0; j<(test_gridV.getNumberOfKnots()); j++ ){
	 double v = test_gridV.getKnot(j).u;

	 float x0, y0, z0, x, y, z;

	 spline.getSplineVec(data, u_1,v,x0,y0,z0);
	 test_spline.getSplineVec( test_data, u_1,v,x,y,z);
	 tmpv_x0[0] = x0;
	 tmpv_x[0] = x;
       
	 spline.getSplineVec(data, u_2,v,x0,y0,z0);
	 test_spline.getSplineVec( test_data, u_2,v,x,y,z);
	 tmpv_x0[1] = x0;
	 tmpv_x[1] = x;
       
	 spline.getSplineVec(data, u_3,v,x0,y0,z0);
	 test_spline.getSplineVec( test_data, u_3,v,x,y,z);
	 tmpv_x0[2] = x0;
	 tmpv_x[2] = x;
      
	 spline.getSplineVec(data, u_4,v,x0,y0,z0);
	 test_spline.getSplineVec( test_data, u_4,v,x,y,z);
	 tmpv_x0[3] = x0;
	 tmpv_x[3] = x;

	 tmpv_result = tmpv_x - tmpv_x0;

	 tmpv_result[0] = abs(tmpv_result[0]);
	 tmpv_result[1] = abs(tmpv_result[1]);
	 tmpv_result[2] = abs(tmpv_result[2]);
	 tmpv_result[3] = abs(tmpv_result[3]);


	 tmp_mask = tmpv_result >= deviation;       
	 mask = mask || tmp_mask;
       }
       if(mask[0]){
	 knotsU_final[nKnotsU_final] = u_1;
	 nKnotsU_final++;
       }
       if(mask[1]){
	 knotsU_final[nKnotsU_final] = u_2;
	 nKnotsU_final++;
       }
       if(mask[2]){
	 knotsU_final[nKnotsU_final] = u_3;
	 nKnotsU_final++;
       }
       if(mask[3]){
	 knotsU_final[nKnotsU_final] = u_4;
	 nKnotsU_final++;
       }
     }
   }
   knotsU_final[nKnotsU_final] = 1;
   
   IrregularSpline2D3D vectoru_spline;
   vectoru_spline.construct(nKnotsU_final, knotsU_final, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV);

   int nKnotsTot_vectoru = vectoru_spline.getNumberOfKnots();
 
   const IrregularSpline1D &vectoru_gridU = vectoru_spline.getGridU();
   const IrregularSpline1D &vectoru_gridV = vectoru_spline.getGridV();
  
   float *vectoru_data0 = new float[ 3*nKnotsTot_vectoru ];
   float *vectoru_data = new float[ 3*nKnotsTot_vectoru ];

   int nu_vectoru = vectoru_gridU.getNumberOfKnots();
    
   for( int i=0; i<vectoru_gridU.getNumberOfKnots(); i++ ){
     double vectoru_u = vectoru_gridU.getKnot( i ).u ;
     for( int j=0; j<vectoru_gridV.getNumberOfKnots(); j++ ){  
       double vectoru_v =  vectoru_gridV.getKnot( j ).u ;
       int vectoru_ind = (nu_vectoru*j+i)*3;
       float x0, y0, z0;
       spline.getSplineVec( data, vectoru_u,vectoru_v,x0,y0,z0);
       vectoru_data0[vectoru_ind+0] = x0;
       vectoru_data0[vectoru_ind+1] = y0;
       vectoru_data0[vectoru_ind+2] = z0; 
     }
   }

   for( int i=0; i<3*nKnotsTot_vectoru; i++ ){
     vectoru_data[i] = vectoru_data0[i];
   }


   vectoru_spline.correctEdges( vectoru_data );

   auto finish_vector_u = std::chrono::high_resolution_clock::now();


   float knotsV_final[vectoru_gridV.getNumberOfKnots()*v_ticks];
   int nKnotsV_final = 0;
   
   for( int i=0; i<(vectoru_gridV.getNumberOfKnots()-1); i++ ){
     double v_small = vectoru_gridV.getKnot(i).u;
     knotsV_final[nKnotsV_final] = v_small;
     nKnotsV_final++;
     double offset = (vectoru_gridV.getKnot(i+1).u - v_small)/(v_ticks+2);
     for( int k=1; k<v_ticks; k= k+4 ){
       double v_1 = v_small+(offset*k);
       double v_2 = v_1+offset;
       double v_3 = v_2+offset;
       double v_4 = v_3+offset;

       mask[0] = false;
       mask[1] = false;
       mask[2] = false;
       mask[3] = false;
       
       for( int j=0; j<(vectoru_gridU.getNumberOfKnots()); j++ ){
	 double u = vectoru_gridU.getKnot(j).u;

	 float x0, y0, z0, x, y, z;

	 spline.getSplineVec(data, u,v_1,x0,y0,z0);
	 vectoru_spline.getSplineVec( vectoru_data, u,v_1,x,y,z);
	 tmpv_x0[0] = x0;
	 tmpv_x[0] = x;
       
	 spline.getSplineVec(data, u,v_2,x0,y0,z0);
	 vectoru_spline.getSplineVec( vectoru_data, u,v_2,x,y,z);
	 tmpv_x0[1] = x0;
	 tmpv_x[1] = x;
       
	 spline.getSplineVec(data, u,v_3,x0,y0,z0);
	 vectoru_spline.getSplineVec( vectoru_data, u,v_3,x,y,z);
	 tmpv_x0[2] = x0;
	 tmpv_x[2] = x;
      
	 spline.getSplineVec(data, u,v_4,x0,y0,z0);
	 vectoru_spline.getSplineVec( vectoru_data, u,v_4,x,y,z);
	 tmpv_x0[3] = x0;
	 tmpv_x[3] = x;

	 tmpv_result = tmpv_x - tmpv_x0;

	 tmpv_result[0] = abs(tmpv_result[0]);
	 tmpv_result[1] = abs(tmpv_result[1]);
	 tmpv_result[2] = abs(tmpv_result[2]);
	 tmpv_result[3] = abs(tmpv_result[3]);


	 tmp_mask = tmpv_result >= deviation;       
	 mask = mask || tmp_mask;
       }
       if(mask[0]){
	 knotsV_final[nKnotsV_final] = v_1;
	 nKnotsV_final++;
       }
       if(mask[1]){
	 knotsV_final[nKnotsV_final] = v_2;
	 nKnotsV_final++;
       }
       if(mask[2]){
	 knotsV_final[nKnotsV_final] = v_3;
	 nKnotsV_final++;
       }
       if(mask[3]){
	 knotsV_final[nKnotsV_final] = v_4;
	 nKnotsV_final++;
       }
     }
   }
   knotsV_final[nKnotsV_final] = 1;
   
   IrregularSpline2D3D vector_spline;
   vector_spline.construct(nKnotsU_final, knotsU_final, nAxisTicksU, nKnotsV_final, knotsV_final, nAxisTicksV);

   int nKnotsTot_vector = vector_spline.getNumberOfKnots();
 
   const IrregularSpline1D &vector_gridU = vector_spline.getGridU();
   const IrregularSpline1D &vector_gridV = vector_spline.getGridV();
  
   float *vector_data0 = new float[ 3*nKnotsTot_vector ];
   float *vector_data = new float[ 3*nKnotsTot_vector ];

   int nu_vector = vector_gridU.getNumberOfKnots();
    
   for( int i=0; i<vector_gridU.getNumberOfKnots(); i++ ){
     double vector_u = vector_gridU.getKnot( i ).u ;
     for( int j=0; j<vector_gridV.getNumberOfKnots(); j++ ){  
       double vector_v =  vector_gridV.getKnot( j ).u ;
       int vector_ind = (nu_vector*j+i)*3;
       float x0, y0, z0;
       spline.getSplineVec( data, vector_u,vector_v,x0,y0,z0);
       vector_data0[vector_ind+0] = x0;
       vector_data0[vector_ind+1] = y0;
       vector_data0[vector_ind+2] = z0; 
     }
   }

   for( int i=0; i<3*nKnotsTot_vector; i++ ){
     vector_data[i] = vector_data0[i];
   }


   vector_spline.correctEdges( vector_data );

   auto finish_vector = std::chrono::high_resolution_clock::now();
   float time_vector_u = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_vector_u - start_vector).count()*0.00001;
   float time_vector = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_vector - start_vector).count()*0.00001;


   cout<<"Time Compare Spline= "<<time_compare<<endl;
   cout<<"Time First Normal Spline= "<<time_normal_u<<endl;
   cout<<"Time Second Normal Spline= "<<time_normal<<endl;
   cout<<"Time First Vector Spline= "<<time_vector_u<<endl;
   cout<<"Time Second Vector Spline= "<<time_vector<<endl;
   */
   /*
  //generating result with non vectorized code 2Splines------------------------------------------------------------------------------------
   //max_tolerated_deviation = 0.5;



  auto start_2 = std::chrono::high_resolution_clock::now();
  cout<<"Starting non vectorized version 2 Splines"<<endl;
  //max_tolerated_deviation = 0.01;

  IrregularSpline2D3D correct_2spline;
  float * correct_2data = generate_2spline_u(spline, correct_2spline, data, nKnotsU, nKnotsV, max_tolerated_deviation, nAxisTicksU, nAxisTicksV, u_ticks);

  auto correct_first = std::chrono::high_resolution_clock::now();
  std::cout << "First Spline Done(non vectorized 2 Splines): " << std::chrono::duration_cast<std::chrono::nanoseconds>(correct_first-start_2).count()*0.00001 << " Millisekunden"<< std::endl;

  IrregularSpline2D3D correct_2spline_two;
  float * correct_2data_two = generate_2spline_uv(spline, correct_2spline, correct_2spline_two, data, correct_2data, max_tolerated_deviation, nAxisTicksU, nAxisTicksV, v_ticks, nKnotsU, nKnotsV);
  int totalknots_2 = correct_2spline_two.getNumberOfKnots();
  auto finish_2 = std::chrono::high_resolution_clock::now();
  float time_2_spline = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_2-start_2).count()*0.00001;
  std::cout << "Second Spline Done(non vectorized 2 Splines): " << time_2_spline << " Millisekunden"<< std::endl;
  cout<<"Finished non vectorized version 2 Splines"<<endl;



  
  //generating result with

  auto start_3 = std::chrono::high_resolution_clock::now();
  cout<<"Starting non vectorized version 1 Spline"<<endl;
  IrregularSpline2D3D spline_one;
  float * data_one =  generate_1spline (spline, spline_one, data, nKnotsU, nKnotsV, max_tolerated_deviation, nAxisTicksU, nAxisTicksV, u_ticks, v_ticks);
  int totalknots_1 = spline_one.getNumberOfKnots();
  auto finish_3 = std::chrono::high_resolution_clock::now();
  float time_1_spline = std::chrono::duration_cast<std::chrono::nanoseconds>(finish_3-start_3).count()*0.00001;
  std::cout << "Spline Done(non vectorized 1 Spline1): " << time_1_spline << " Millisekunden"<< std::endl;
  cout<<"Finished non vectorized version 1 Spline1"<<endl;
   */

  //Drawing Results------------------------------------------------------
  TCanvas *canv = new TCanvas("cQA","2D splines  QA",1500,1500);
  canv->Draw();
  //canv->Divide(3,3);
  canv->Update();

  TGraph2D *gknots = new TGraph2D();
  gknots->SetName("gknots");
  gknots->SetTitle("gknots");
  gknots->SetLineColor(kRed);
  gknots->SetMarkerSize(1.);
  gknots->SetMarkerStyle(8);
  gknots->SetMarkerColor(kRed);


  int nu = gridU.getNumberOfKnots();
  int iterate = 0;
  int gknotsN=0;
  TNtuple *knots = new TNtuple("knots","knots", "u:v:f");
  double diff=0;
  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    double u = gridU.getKnot( i ).u;
    for( int j=0; j<gridV.getNumberOfKnots(); j++ ){
      double v = gridV.getKnot( j ).u;
      int ind = (nu*j+i)*3;
      double fx0 = data[ind+0];
      float x, y, z;
      spline.getSpline( data, u,v,x,y,z);
      

      diff+=(fx0-x)*(fx0-x);
      knots->Fill( u, v, fx0 );
      gknots->SetPoint(gknotsN++,u,v,fx0); 
      iterate++;
    }
  }
  knots->SetMarkerSize(1.);
  knots->SetMarkerStyle(8);
  knots->SetMarkerColor(kRed);

  TNtuple *nt = new TNtuple("nt","nt","u:v:f0:fs");


  TGraph2D *gf0 = new TGraph2D();
  gf0->SetName("gf0");
  gf0->SetTitle("Input Spline");
  gf0->SetLineColor(kBlue);
  int gf0N=0;

  TGraph2D *gfs = new TGraph2D();
  gfs->SetName("gfs");
  gfs->SetTitle("3 Spline");
  gfs->SetLineColor(kRed);
  int gfsN=0;

  /*
  
  TGraph2D *gfx = new TGraph2D();
  gfx->SetName("gfx");
  gfx->SetTitle("2 Splines");
  gfx->SetLineColor(kRed);
  int gfxN=0;

  TGraph2D *gfy = new TGraph2D();
  gfy->SetName("gfy");
  gfy->SetTitle("1 Spline");
  gfy->SetLineColor(kRed);
  int gfyN=0;

  */
  
  TH1F *time = new TH1F("Time","Time;Version;Time(MilliSec)",3,0,100);
  time->SetBarWidth(0.5);
  time->SetBarOffset(0.5);
  time->SetStats(0);

  time->SetFillColor(4);
  time->SetBinContent(1, time_normal+time_compare);
  time->GetXaxis()->SetBinLabel(1, "Time");
//time->SetBinContent(2, time_vector+time_compare);
// time->GetXaxis()->SetBinLabel(2, "Time(Vector)");
  
  /*

  TH1F *barknots = new TH1F("Knots","Knots;Version;Count(Knots)",3,0,100);
  barknots->SetBarWidth(0.5);
  barknots->SetBarOffset(0.5);
  barknots->SetStats(0);

  
  barknots->SetFillColor(4);
  barknots->SetBinContent(1, totalknots_3);
  barknots->GetXaxis()->SetBinLabel(1, "Knots 3 Splines");
  barknots->SetBinContent(2, totalknots_2);
  barknots->GetXaxis()->SetBinLabel(2, "Knots 2 Splines");
  barknots->SetBinContent(3, totalknots_1);
  barknots->GetXaxis()->SetBinLabel(3, "Knots 1 Spline");

  
  TH1F *deviation_3 = new TH1F("Deviation 3 Spline","Deviation 3 Splines ;Amount(Deviation);Count(Deviation)",1000,-1.,1.);
  TH1F *deviation_2 = new TH1F("Deviation 3 Spline","Deviation 2 Splines ;Amount(Deviation);Count(Deviation)",1000,-1.,1.);
  TH1F *deviation_1 = new TH1F("Deviation 3 Spline","Deviation 1 Spline ;Amount(Deviation);Count(Deviation)",1000,-1.,1.);
  */

  TH1F *qaX = new TH1F("qaX","qaX [um]",1000,-0.1,0.1); 
  TH1F *qaY = new TH1F("qaY","qaY [um]",1000,-0.1,0.1);
  TH1F *qaZ = new TH1F("qaZ","qaZ [um]",1000,-0.1,0.1);

  int iter=0;
  float stepu = 1.e-3;
  float stepv = 1.e-3;


  for( float u=0; u<=1; u+=stepu ){
    for( float v=0; v<=1; v+=stepv ){
      float x0, y0, z0;
      spline.getSplineVec( data, u,v,x0,y0,z0);
      float x1, y1, z1;
      correct_spline_two.getSplineVec( correct_data_two, u,v,x1,y1,z1);
      //float x2, y2, z2;
      //correct_spline_two.getSplineVec( correct_data_two, u,v,x1,y1,z1);
      //float x3, y3, z3;
      //spline_one.getSplineVec( data_one, u,v,x3,y3,z3);
      //float x4, y4, z4;
      //vector_spline.getSplineVec( vector_data, u,v,x4,y4,z4);
      if( u>=0 && v>=0 && u<=1 && v<=1 ){
	qaX->Fill( (x1 - x0) );
	qaY->Fill( (y1 - y0) );
	qaZ->Fill( (z1 - z0) );
	nt->Fill(u,v,x0,x0);
	//deviation_3->Fill( (x1 - x0)*1.e2 );
	//deviation_2->Fill( (x2 - x0)*1.e2 );
	//deviation_1->Fill( (x3 - x0)*1.e2 );
      }
      gf0->SetPoint(gf0N++,u,v,x1);
//gfs->SetPoint(gfsN++,u,v,x1);
      //gfx->SetPoint(gfxN++,u,v,x2);
      //gfy->SetPoint(gfyN++,u,v,x3);
    }
  }
gStyle->SetPalette(1);
gf0->Draw("surf");
//gfs->Draw("surf, same");
  //gknots->Draw("P,same");
  //time->Draw("b");
  canv->Update();

  /*
  canv->cd(1);
  gf0->Draw("surf");
  gknots->Draw("P,same");
  canv->Update();

  canv->cd(2);
  time->Draw("b");
  canv->Update();

  canv->cd(3);
  barknots->Draw("b");
  canv->Update();

  canv->cd(4);
  gfs->Draw("surf");
  gf0->Draw("surf,same");
  canv->Update();

  canv->cd(5);
  gfx->Draw("surf");
  gf0->Draw("surf,same");
  canv->Update();

  canv->cd(6);
  gfy->Draw("surf");
  gf0->Draw("surf,same");
  canv->Update();
  
  canv->cd(7);
  deviation_3->Draw("L");
  canv->Update();
  
  canv->cd(8);
  deviation_2->Draw("L");
  canv->Update();

  canv->cd(9);
  deviation_1->Draw("L");
  canv->Update();
  */

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
