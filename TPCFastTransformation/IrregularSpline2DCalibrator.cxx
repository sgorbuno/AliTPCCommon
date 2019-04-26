// // Copyright CERN and copyright holders of ALICE O2. This software is
// // distributed under the terms of the GNU General Public License v3 (GPL
// // Version 3), copied verbatim in the file "COPYING".
// //
// // See http://alice-o2.web.cern.ch/license for full licensing information.
// //
// // In applying this license CERN does not waive the privileges and immunities
// // granted to it by virtue of its status as an Intergovernmental Organization
// // or submit itself to any jurisdiction.


// /// \file IrregularSpline2DCalibrator.cxx
// /// \brief Implementation of IrregularSpline2DCalibrator class
// ///
// /// \author  Oscar Lange <langeoscar96@googlemail.com>


// #include "IrregularSpline2DCalibrator.h"
// #include "AliTPCCommonDef.h"

// #include "TRandom.h"
// #include "TNtuple.h"
// #include "Riostream.h"
// #include "algorithm"


#include "IrregularSpline2D3D.h"
#include "IrregularSpline2DCalibrator.h"


namespace ali_tpc_common {
namespace tpc_fast_transformation {

  IrregularSpline2DCalibrator::IrregularSpline2DCalibrator()
{
  nKnotsU = 20;
  nKnotsV = 20;
  minKnotsU = 5;
  minKnotsV = 5;
  u_ticks = 20;
  v_ticks = 20;
  max_tolerated_deviation = 0.5;
  nAxisTicksU = u_ticks*minKnotsU; 
  nAxisTicksV = v_ticks*minKnotsV;
}

  float *  IrregularSpline2DCalibrator::generate_input_spline (IrregularSpline2D3D &input_spline, void (*F)(float, float, float&, float&, float&)){
  
  float du = 1./(nKnotsU-1);
  float dv = 1./(nKnotsV-1);

  float knotsU[nKnotsU], knotsV[nKnotsV];
  knotsU[0] = 0;
  knotsU[nKnotsU-1] = 1;
  knotsV[0] = 0;
  knotsV[nKnotsV-1] = 1;

  for( int i=1; i<nKnotsU-1; i++ ){
    knotsU[i] = i*du;
  }

  for( int i=1; i<nKnotsV-1; i++ ){
    knotsV[i] = i*dv;
  }
  input_spline.construct( nKnotsU, knotsU, nAxisTicksU, nKnotsV, knotsV, nAxisTicksV );
 
  int nKnotsTot = input_spline.getNumberOfKnots();
 
  const IrregularSpline1D &gridU = input_spline.getGridU();
  const IrregularSpline1D &gridV = input_spline.getGridV();      

  float *input_data0 = new float[ 3*nKnotsTot ];
  float *input_data = new float[ 3*nKnotsTot ];

  
  int nu = gridU.getNumberOfKnots();

  for( int i=0; i<gridU.getNumberOfKnots(); i++ ){
    double u = gridU.getKnot( i ).u ;
    for( int j=0; j<gridV.getNumberOfKnots(); j++ ){ 
      double v =  gridV.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      float x,y,z;
      (F)(u,v,x,y,z);
      input_data0[ind+0] = x;
      input_data0[ind+1] = y;
      input_data0[ind+2] = z;
    }
  }  

  for( int i=0; i<3*nKnotsTot; i++ ){
    input_data[i] = input_data0[i];
  }

  input_spline.correctEdges( input_data );
  return input_data;
}


float *  IrregularSpline2DCalibrator::generate_compare_spline (IrregularSpline2D3D &compare_spline, IrregularSpline2D3D &input_spline, float *input_data){
  float du = 1./(minKnotsU-1);
  float dv = 1./(minKnotsV-1);

  float knotsU[minKnotsU], knotsV[minKnotsV];
  knotsU[0] = 0;
  knotsU[minKnotsU-1] = 1;
  knotsV[0] = 0;
  knotsV[minKnotsV-1] = 1;

  for( int i=1; i<minKnotsU-1; i++ ){
    knotsU[i] = i*du;
  }

  for( int i=1; i<minKnotsV-1; i++ ){
    knotsV[i] = i*dv;
  }
  
  compare_spline.construct(minKnotsU, knotsU, nAxisTicksU, minKnotsV, knotsV, nAxisTicksV);

  int nKnotsTot = compare_spline.getNumberOfKnots();
 
  const IrregularSpline1D &compare_gridU = compare_spline.getGridU();
  const IrregularSpline1D &compare_gridV = compare_spline.getGridV();      

  float *compare_data0 = new float[ 3*nKnotsTot ];
  float *compare_data = new float[ 3*nKnotsTot ];

  int nu = compare_gridU.getNumberOfKnots();
  
    
  for( int i=0; i<compare_gridU.getNumberOfKnots(); i++ ){
    double u = compare_gridU.getKnot( i ).u ;
    for( int j=0; j<compare_gridV.getNumberOfKnots(); j++ ){  
      double v =  compare_gridV.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      float x0, y0, z0;
      input_spline.getSplineVec( input_data, u,v,x0,y0,z0);
      compare_data0[ind+0] = x0;
      compare_data0[ind+1] = y0;
      compare_data0[ind+2] = z0;
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    compare_data[i] = compare_data0[i];
  }
  compare_spline.correctEdges( compare_data );
  return compare_data;
}

float * IrregularSpline2DCalibrator::generate_spline_u (IrregularSpline2D3D &input_spline, IrregularSpline2D3D &compare_spline, IrregularSpline2D3D &spline_u, float *input_data, float *compare_data){
  const IrregularSpline1D &compare_gridU = compare_spline.getGridU();
  const IrregularSpline1D &compare_gridV = compare_spline.getGridV();
  
  float knotsU[compare_gridU.getNumberOfKnots()*u_ticks];
  int nKnotsU_new = 0; 
  double u_small = 0;
  double u_high = 0;
  float knotsV[minKnotsV];
  for( int k=0; k<minKnotsV; k++ ){
    knotsV[k] = compare_gridV.getKnot(k).u ;
  }

  for( int i=0; i<(compare_gridU.getNumberOfKnots()-1); i++ ){
    u_small = compare_gridU.getKnot( i ).u ;
    u_high = compare_gridU.getKnot( (i+1) ).u ;
    double difference_u = u_high - u_small;
    float du = difference_u/(u_ticks+2);
    knotsU[nKnotsU_new] = u_small;
    nKnotsU_new++;
    for( int k=1; k<u_ticks; k++ ){
      bool value_changed = false;
      double current_u = du*k + u_small;
      for( int j=0; j<(compare_gridV.getNumberOfKnots()); j++ ){
	if(!value_changed){
	  double v = compare_gridV.getKnot( j ).u;
	  float x0, y0, z0;
	  input_spline.getSplineVec(input_data, current_u,v,x0,y0,z0);
	  float x, y, z;
	  compare_spline.getSplineVec( compare_data, current_u,v,x,y,z);
	  if ((abs(x-x0)) >= max_tolerated_deviation){
	    value_changed = true;
	    knotsU[nKnotsU_new] = current_u;
	    nKnotsU_new++;
	  }
	}
      }
    }
  }
  knotsU[nKnotsU_new] = u_high;
  nKnotsU_new++;
  
  spline_u.construct(nKnotsU_new, knotsU, nAxisTicksU, minKnotsV, knotsV, nAxisTicksV);
  
  int nKnotsTot = spline_u.getNumberOfKnots();
 
  const IrregularSpline1D &gridU_u = spline_u.getGridU();
  const IrregularSpline1D &gridV_u = spline_u.getGridV();
  
  float *data0_u = new float[ 3*nKnotsTot ];
  float *data_u = new float[ 3*nKnotsTot ];

  int nu = gridU_u.getNumberOfKnots();
    
  for( int i=0; i<gridU_u.getNumberOfKnots(); i++ ){
    double u = gridU_u.getKnot( i ).u ;
    for( int j=0; j<gridV_u.getNumberOfKnots(); j++ ){  
      double v =  gridV_u.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      float x0, y0, z0;
      input_spline.getSplineVec( input_data, u, v,x0,y0,z0);
      data0_u[ind+0] = x0;
      data0_u[ind+1] = y0;
      data0_u[ind+2] = z0; 
    }
  }

  for( int i=0; i<3*nKnotsTot; i++ ){
    data_u[i] = data0_u[i];
  }


  spline_u.correctEdges( data_u );
  return data_u;
}

float *  IrregularSpline2DCalibrator::generate_spline_uv(IrregularSpline2D3D &input_spline, IrregularSpline2D3D &spline_u, IrregularSpline2D3D &spline_uv, float *input_data, float *data_u){
  const IrregularSpline1D &gridU_u = spline_u.getGridU();
  const IrregularSpline1D &gridV_u = spline_u.getGridV();
  
  int nKnotsU_new = gridU_u.getNumberOfKnots();
  float knotsV[gridV_u.getNumberOfKnots()*v_ticks];
  int nKnotsV_new = 0;
  double v_high = 0;
  double v_small = 0;
  float knotsU[nKnotsU_new];
  for( int i=0; i<gridU_u.getNumberOfKnots(); i++ ){
    knotsU[i] = gridU_u.getKnot( i ).u ;
  }
  
  for( int i=0; i<(gridV_u.getNumberOfKnots()-1); i++ ){
    v_small = gridV_u.getKnot( i ).u ;
    v_high = gridV_u.getKnot( (i+1) ).u ;
    double difference_v = v_high - v_small;
    float dv = difference_v/(v_ticks+2);
    knotsV[nKnotsV_new] = v_small;
    nKnotsV_new++;
    for( int k=1; k<v_ticks; k++ ){
      bool value_changed = false;
      double v = dv*k + v_small;
      for( int j=0; j<gridU_u.getNumberOfKnots(); j++ ){
	if(!value_changed){
	  double u = gridU_u.getKnot( j ).u ;
	  float x0, y0, z0;
	  input_spline.getSplineVec( input_data, u,v,x0,y0,z0);
	  float x, y, z;
	  spline_u.getSplineVec( data_u,u,v,x,y,z);
	  if ((abs(x-x0)) >= max_tolerated_deviation){
	    value_changed = true;
	    knotsV[nKnotsV_new] = v;
	    nKnotsV_new++;
	  }
	}
      }
    }
  }
  knotsV[nKnotsV_new] = v_high;
  nKnotsV_new++;

  spline_uv.construct( nKnotsU_new, knotsU, nAxisTicksU, nKnotsV_new, knotsV, nAxisTicksV );
 
  int nKnotsTot = spline_uv.getNumberOfKnots();
  const IrregularSpline1D &gridU_uv = spline_uv.getGridU();
  const IrregularSpline1D &gridV_uv = spline_uv.getGridV();

  float *data0_uv = new float[ 3*nKnotsTot ];
  float *data_uv = new float[ 3*nKnotsTot ];

  int nu = gridU_uv.getNumberOfKnots();
    
  for( int i=0; i<gridU_uv.getNumberOfKnots(); i++ ){
    double u = gridU_uv.getKnot( i ).u ;
    for( int j=0; j<gridV_uv.getNumberOfKnots(); j++ ){  
      double v =  gridV_uv.getKnot( j ).u ;
      int ind = (nu*j+i)*3;
      float x0, y0, z0;
      input_spline.getSplineVec( input_data, u,v,x0,y0,z0);
      data0_uv[ind+0] = x0;
      data0_uv[ind+1] = y0;
      data0_uv[ind+2] = z0; 
    }
  }
  for( int i=0; i<3*nKnotsTot; i++ ){
    data_uv[i] = data0_uv[i];
  }

  spline_uv.correctEdges( data_uv );


  return data_uv;
}

float * IrregularSpline2DCalibrator::calibrateSpline (IrregularSpline2D3D &spline_uv, void (*F)(float, float, float&, float&, float&)){

  IrregularSpline2D3D input_spline;
  float * input_data = generate_input_spline(input_spline, F);
  
  IrregularSpline2D3D compare_spline;
  float * compare_data = generate_compare_spline(compare_spline, input_spline, input_data);

  IrregularSpline2D3D spline_u;
  float * data_u = generate_spline_u(input_spline, compare_spline, spline_u, input_data, compare_data);

  float * data_uv = generate_spline_uv(input_spline, spline_u, spline_uv, input_data, data_u);

  return data_uv;
}

 
}// namespace
}// namespace

