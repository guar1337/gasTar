#ifndef helloworld5_hxx
#define helloworld5_hxx 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "TRandomGen.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include <fstream>
#include "TF1.h"


std::string smallFile = "/home/zalewski/dataTmp/small5/small5_fnl.root";
//enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos1, sTarPos2, sTarPos3, sLang1, sLang2, sLang3, sRang};
enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sRang};
//std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sTarPos1", "sTarPos2", "sTarPos3", "sLang1", "sLang2", "sLang3", "sRang"};
std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sRang"};

Double_t tarPos, tarAngle;
Double_t sqlDist = 170.0;
Double_t sqrDist = 132.0;
Double_t rightAngle;
Double_t myRightAngleCorrection;
Double_t sqrDistCorrection = 168.0;
Double_t rightAngleCorrection;
Double_t leftAngle = 70.0 * TMath::DegToRad();

Bool_t calculate;

const Int_t numberOfParameters = 5;
Double_t lowerParamBound[numberOfParameters], upperParamBound[numberOfParameters];
Int_t firstRun, lastRun;

std::string outGeneratedParams = "/home/zalewski/Desktop/6He/analysis/geo5/MWPC_rAng_pp_dd/generated5.txt";
std::string outObtainedParams = "/home/zalewski/Desktop/6He/analysis/geo5/MWPC_rAng_pp_dd/obtained5.txt";

Double_t myMemo[4] = {1000.0, 1000.0, 1000.0, 1000.0};
Double_t tarMass1H = cs::mass1H;
Double_t tarMass2H = cs::mass2H;
#endif
