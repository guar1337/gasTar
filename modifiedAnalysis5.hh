#ifndef modifiedAnalysis5_hh
#define modifiedAnalysis5_hh 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "ROOT/RDataFrame.hxx"
#include <fstream>

enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos, sLang1, sLang2, sLang3, sRang};
std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sTarPos", "sLang1", "sLang2", "sLang3", "sRang"};
std::vector<Double_t> parameters;

enum sCuts {mX1, X1, mY1, Y1, mX2, X2, mY2, Y2, mX3, X3, mY3, Y3};
std::vector<std::string> corrNames = {"-X1", "X1", "-Y1", "Y1", "-X2", "X2", "-Y2", "Y2", "-X3", "X3", "-Y3", "Y3"};
std::vector<Int_t> vCut;
Double_t tarPos = 0.0;
Double_t tarAngle = 33.0 * TMath::DegToRad();
Double_t leftAngle = 70.0 * TMath::DegToRad();
Double_t sqlDist = 170.0;
Double_t rightAngle = (9.0) * TMath::DegToRad();
Double_t sqrDist = 132.0;
Bool_t verbosity;

Double_t tarMass1H = cs::mass1H;
Double_t tarMass2H = cs::mass2H;

Double_t proBeamIntegral, deutBeamIntegral, proTargetThick, deutTargetThick;
Double_t proBeamTargetConst, deutBeamTargetConst;

Int_t protium = 1;
Int_t deuterium = 2;


#endif