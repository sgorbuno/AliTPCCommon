// @(#) $Id$
//*************************************************************************
// This file is property of and copyright by the ALICE HLT Project        * 
// ALICE Experiment at CERN, All rights reserved.                         *
//                                                                        *
// Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>              *
//                  for The ALICE HLT Project.                            *
//                                                                        *
// Permission to use, copy, modify and distribute this software and its   *
// documentation strictly for non-commercial purposes is hereby granted   *
// without fee, provided that the above copyright notice appears in all   *
// copies and that both the copyright notice and this permission notice   *
// appear in the supporting documentation. The authors make no claims     *
// about the suitability of this software for any purpose. It is          *
// provided "as is" without express or implied warranty.                  *
//*************************************************************************

#include "AliHLTTPCCAHit.h"


ClassImp(AliHLTTPCCAHit);

void AliHLTTPCCAHit::Set( Int_t ID, Double_t Y, Double_t Z, 
			  Double_t ErrY, Double_t ErrZ  )
{
  // set parameters
  fID = ID;
  fY = Y;
  fZ = Z;
  fErrY = ErrY;
  fErrZ = ErrZ;
}
