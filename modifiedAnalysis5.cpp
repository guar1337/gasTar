#include "modifiedAnalysis5.hh"

void trimHistogramLeft(TH1F *tmpHisto, Int_t noBins)
{
	Int_t minimumBin = tmpHisto->FindFirstBinAbove(0.1);
	for (int iii = minimumBin; iii < minimumBin + noBins; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}

}

void trimHistogramRight(TH1F *tmpHisto, Int_t noBins)
{
	Int_t maximumBin = tmpHisto->FindLastBinAbove(0.1);
	for (int iii = maximumBin - noBins; iii <= maximumBin; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}
}

void makeGraph()
{
	TString mcDataPath = "/home/zalewski/Desktop/6He/analysis/data5Out/MCDataCMHisto5.root";
	TFile mcDataFile(mcDataPath.Data(), "UPDATE");

	TString histName = "ppMC_geo5_CM";
	TH1F *mcPPgeo5Hist = (TH1F*)mcDataFile.Get(histName.Data());
	
	histName.ReplaceAll("pp", "dd");
	TH1F *mcDDgeo5Hist = (TH1F*)mcDataFile.Get(histName.Data());

	TString realDataPath = "/home/zalewski/Desktop/6He/analysis/data5Out/realDataCMHisto5.root";
	TFile realDataFile(realDataPath.Data(), "READ");
	histName = "ppReal_geo5_CM";
	TH1F *realPPgeo5Hist = (TH1F*)realDataFile.Get(histName.Data());

	histName.ReplaceAll("pp", "dd").ReplaceAll("geo3","geo1");
	TH1F *realDDgeo5Hist = (TH1F*)realDataFile.Get(histName.Data());

	TString ppGraphTitle = "Elastic scattering 6He+p";
	TString ddGraphTitle = "Elastic scattering 6He+d";
	
	TCanvas *canvasPP = new TCanvas("canvasPP","canvasPP",10,10,1200,700);
	THStack *ppxSectStack = new THStack("PP X-section stack", ppGraphTitle.Data());
	TFile theoHistFile("/home/zalewski/Desktop/6He/analysis/dataOut/source/theoHists.root", "READ");
	TH1F *wolskiExp = (TH1F*)theoHistFile.Get("wolskiExp");
	TH1F *wolskiOM = (TH1F*)theoHistFile.Get("wolskiOM");

	TH1F *xSectPP = new TH1F("xSectPP","xSectPP",90,0,180);
	xSectPP->SetBinContent(15, 1.0*165.578293144036);
	xSectPP->SetBinContent(16, 1.0*140.718631211516);
	xSectPP->SetBinContent(17, 1.0*125.660854195774);
	xSectPP->SetBinContent(18, 1.0*149.31528473443);
	xSectPP->SetBinContent(19, 1.0*124.483452014524);
	xSectPP->SetBinContent(20, 1.0*105.283409917899);
	xSectPP->SetBinContent(21, 1.0*80.955529943601);
	xSectPP->SetBinContent(22, 1.0*80.158040438495);
	xSectPP->SetBinContent(23, 1.0*62.622811969959);
	xSectPP->SetBinContent(24, 1.0*43.349015378232);

/*
	xSectPP->SetBinContent(20, 92.55);
	xSectPP->SetBinContent(21, 78.81);
	xSectPP->SetBinContent(22, 80.72);
	xSectPP->SetBinContent(23, 68.28);
	xSectPP->SetBinContent(24, 57.08);
*/
	xSectPP->SetMarkerStyle(3);
	xSectPP->SetMarkerSize(2);
	xSectPP->SetMarkerColor(kBlack);

	canvasPP->SetLogy();
	Int_t normFactor = 1200;
	realPPgeo5Hist->Divide(mcPPgeo5Hist);
	realPPgeo5Hist->Scale(normFactor);
	realPPgeo5Hist->SetMarkerStyle(20);
	realPPgeo5Hist->SetMarkerSize(1);
	realPPgeo5Hist->SetMarkerColor(kRed);

	wolskiExp->SetMarkerStyle(20);
	wolskiExp->SetMarkerSize(1);
	wolskiExp->SetMarkerColor(kBlack);

	wolskiOM->SetLineStyle(1);
	wolskiOM->SetLineColor(kBlack);

	trimHistogramLeft(realPPgeo5Hist, 7);
	trimHistogramRight(realPPgeo5Hist, 0);

	//ppxSectStack->Add(realPPgeo5Hist);
	ppxSectStack->Add(wolskiExp);
	ppxSectStack->Add(xSectPP);
	ppxSectStack->Add(wolskiOM, "L");

	ppxSectStack->SetMinimum(0.01);
	ppxSectStack->SetMaximum(1000);
	ppxSectStack->Draw();
	ppxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ppxSectStack->GetXaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ppxSectStack->Draw("HIST, nostack,p");

	TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
	legend->AddEntry(wolskiExp,"Wolski data","p");
	legend->AddEntry(wolskiOM,"Wolski OM","L");
	legend->AddEntry(xSectPP,"gas target data","p");
	legend->Draw();

	TString outPPFileName = "/home/zalewski/Desktop/6He/analysis/data5Out/pp.png";
	canvasPP->Print(outPPFileName.Data());
	canvasPP->SaveAs("/home/zalewski/Desktop/6He/analysis/data5Out/pp.C");
	delete ppxSectStack;
	delete canvasPP;
/*
	TCanvas *canvasDD = new TCanvas("canvasDD","canvasDD",10,10,1200,700);
	THStack *ddxSectStack = new THStack("DD X-section stack", ddGraphTitle.Data());
	TH1F *pereyPerey = new TH1F("perey","pereyPerey",90,0,180);
	#include "/home/zalewski/Desktop/6He/analysis/dataOut/source/pereyPerey.hh"
	TH1F *daehnick = new TH1F("dae","daehnick",90,0,180);
	#include "/home/zalewski/Desktop/6He/analysis/dataOut/source/daehnick.hh"
	canvasDD->SetLogy();
	normFactor = 5500;
	Double_t deuterMod = 2.0;
	realDDgeo5Hist->Divide(mcDDgeo5Hist);
	realDDgeo5Hist->Scale(normFactor);
	realDDgeo5Hist->SetMarkerStyle(20);
	realDDgeo5Hist->SetMarkerSize(1);
	realDDgeo5Hist->SetMarkerColor(kRed);

	pereyPerey->SetLineStyle(1);
	pereyPerey->SetLineColor(kBlack);

	daehnick->SetLineStyle(10);
	daehnick->SetLineColor(kBlack);

	trimHistogramLeft(realDDgeo5Hist, 0);
	trimHistogramRight(realDDgeo5Hist, 0);

	ddxSectStack->Add(realDDgeo5Hist);
	ddxSectStack->Add(daehnick, "L");
	ddxSectStack->Add(pereyPerey, "L");

	ddxSectStack->SetMinimum(0.01);
	ddxSectStack->SetMaximum(10000);
	ddxSectStack->Draw();
	ddxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ddxSectStack->GetXaxis()->CenterTitle();
	ddxSectStack->GetXaxis()->SetRangeUser(0,140);
	ddxSectStack->GetYaxis()->CenterTitle();
	ddxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ddxSectStack->Draw("HIST, nostack,p");

	TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
	legend->AddEntry(realDDgeo5Hist,"Geometry 5","p");
	legend->AddEntry(daehnick,"Daehnick, PRC 6, 2253 (1980)","L");
	legend->AddEntry(pereyPerey,"Perey-Perey","L");
	legend->Draw();

	TString outDDFileName = "/home/zalewski/Desktop/6He/analysis/dataOut/dd.png";
	canvasDD->Print(outDDFileName.Data());

	delete canvasDD;
	delete ddxSectStack;*/
}

void makeSmallFile()
{
	TString chainName = "smallMC";
	TChain smallChain(chainName.Data());
	TString nameModifier = 5;
	smallChain.Add("/home/zalewski/dataTmp/small/geo1/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo1/mc_out_2H_v" + nameModifier + ".root");

	smallChain.Add("/home/zalewski/dataTmp/small/geo2/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo2/mc_out_2H_v" + nameModifier + ".root");

	smallChain.Add("/home/zalewski/dataTmp/small/geo3/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo3/mc_out_2H_v" + nameModifier + ".root");

	ROOT::RDataFrame smallDF(smallChain);
	TString outFName = "/home/zalewski/dataTmp/small/smallMC_v" + nameModifier + ".root";
	smallDF.Snapshot("small", outFName.Data());
}

void loadGeometryCorrectionParameters(Int_t m_runNo)
{
	parameters.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/experimental2/chosen.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("loadGeometryCorrectionParameters:\tFailed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	float tmpContainer;
    if (verbosity==true) printf("%d.\t", m_runNo);
	for (int iii = 0; iii < 10; iii++)
	{
		outStreamGenerated>>tmpContainer;
		if (iii!=0)
		{
			parameters.push_back(tmpContainer);
			if (verbosity==true) printf("%s = %f\t", parNames[iii-1].c_str(), parameters[iii-1]);
		}
	}
	if (verbosity==true) std::cout<<std::endl;
}

void loadCutsCorrectionParameters(Int_t m_runNo)
{
	vCut.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/dataOut/tarVertexCuts.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("loadCutsCorrectionParameters:\tFailed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	Int_t tmpContainer;
    if (verbosity==true) printf("%d.\t", m_runNo);
	for (int iii = 0; iii < 12; iii++)
	{
		outStreamGenerated>>tmpContainer;
		vCut.push_back(tmpContainer);
		if (verbosity==true) printf("%s = %d\t", corrNames[iii].c_str(), vCut[iii]);
	}
	if (verbosity==true) std::cout<<std::endl;
	
}

std::vector<Double_t> getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y)
{
	std::vector<Double_t> rvecMWPC(6);
	rvecMWPC.at(0) = m_MWPC_1_X;
	rvecMWPC.at(1) = m_MWPC_1_Y;
	rvecMWPC.at(2) = -816.0;
	rvecMWPC.at(3) = m_MWPC_2_X; 
	rvecMWPC.at(4) = m_MWPC_2_Y;
	rvecMWPC.at(5) = -270.0;


	return rvecMWPC;
}

TVector3 getBeamVector(std::vector<Double_t> rvecMWPC, Double_t m_kinE)
{
	TVector3 m_beamVector(rvecMWPC[3] - rvecMWPC[0], rvecMWPC[4] - rvecMWPC[1], rvecMWPC[5] - rvecMWPC[2]);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

TVector3 getTarVertex(std::vector<Double_t> rvecMWPC)
{
	Double_t m_dX = rvecMWPC[3] - rvecMWPC[0];
	Double_t m_dY = rvecMWPC[4] - rvecMWPC[1];
	Double_t m_dZ = rvecMWPC[5] - rvecMWPC[2];
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, tarPos);
	TVector3 m_beamPoint(rvecMWPC[3], rvecMWPC[4], rvecMWPC[5]);
	TVector3 m_tarPerpendicular(sin(tarAngle), 0.0, cos(tarAngle));
	Double_t m_dCoeff = ((m_tarPoint-m_beamPoint).Dot(m_tarPerpendicular))/(m_vBeam.Dot(m_tarPerpendicular));
		
	Double_t m_evX = rvecMWPC[3] + m_dX * m_dCoeff;
	Double_t m_evY = rvecMWPC[4] + m_dY * m_dCoeff;
	Double_t m_evZ = rvecMWPC[5] + m_dZ * m_dCoeff;

	return TVector3(m_evX, m_evY, m_evZ);
}

TVector3 getLeftDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system	
	Double_t X2hDet = cs::widthStripX * m_xStrip * cos(leftAngle);
	Double_t Y2hDet = cs::widthStripY * m_yStrip;
	Double_t Z2hDet = -cs::widthStripX * m_xStrip * sin(leftAngle);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * m_xStrip * cos(rightAngle);
	Double_t Y6HeDet = cs::widthStripY * m_yStrip;
	Double_t Z6HeDet = cs::widthStripX * m_xStrip * sin(rightAngle);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

TVector3 getLeftDetPosition()
{
	Double_t X2Hlab = sqlDist*sin(leftAngle) - (cs::sqlXzero) * cos(leftAngle);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(leftAngle) + (cs::sqlXzero) * sin(leftAngle);
	return TVector3(X2Hlab, Y2Hlab, Z2Hlab);
}

TVector3 getRightDetPosition()
{
	Double_t sqrDistCorrection = 168.0;
	Double_t rightAngleCorrection = rightAngle + 1.65 * TMath::DegToRad();

	Double_t X6HeCorrection = sqrDistCorrection*sin(-rightAngleCorrection);
	Double_t Z6HeCorrection = sqrDistCorrection*cos(-rightAngleCorrection);

	Double_t X6Helab = sqrDist*sin(-rightAngle) + X6HeCorrection - (cs::sqrXzero) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart;
	Double_t Z6Helab = sqrDist*cos(rightAngle) + Z6HeCorrection - (cs::sqrXzero) * sin(rightAngle);
	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

Double_t getMCAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1 + m_gammaSquare2H*m_tanSquare);
	Double_t m_leftAngleCM = (TMath::Pi() - (acos(m_cosLeftAng)-m_ThetaCM2H))*TMath::RadToDeg();
	return m_leftAngleCM;
}

Double_t getMCAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1 + m_gammaSquare2H*m_tanSquare);
	Double_t m_leftAngleCM = (TMath::Pi() - (acos(m_cosLeftAng)-m_ThetaCM2H))*TMath::RadToDeg();
	return m_leftAngleCM;
}

void realDataAnalyzer()
{
	//loadGeometryCorrectionParameters(setNo);
	//loadCutsCorrectionParameters(setNo);

	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small5/small5.root";
	ROOT::RDataFrame smallDF("small", smallFileName.Data());

	auto newDF = smallDF.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
						.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
						.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
						.Define("tarVertex", getTarVertex, {"MWPC"})

						.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip"})
						.Define("leftLabVertex", getLeftDetPosition)
						.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
						.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
						.Define("rightLabVertex", getRightDetPosition)
						.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

						.Define("v2H", "leftGlobVertex-tarVertex")
						.Define("v6He", "rightGlobVertex-tarVertex")
						.Define("leftAngle", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
                        .Define("mc1H", getMCAngle1H, {"leftAngle", "lvBeam"})
					    .Define("mc2H", getMCAngle2H, {"leftAngle", "lvBeam"})
						.Define("rightAngle", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut", "(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9");

	TString ppHistoName5 = "ppReal_geo5_CM";
	auto ppCMHist5 = newDF.Filter("pp && tarCut").Histo1D({ppHistoName5.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName5 = "ddReal_geo5_CM";
	auto ddCMHist5 = newDF.Filter("dd && tarCut").Histo1D({ddHistoName5.Data(), "Real elastic scattering on deuterons", 90,0,180}, {"mc2H"});

	newDF.Snapshot("checkTargetVertex", "/home/zalewski/Desktop/tar.root");
	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/data5Out/realDataCMHisto5.root", "RECREATE");
	ppCMHist5->Write(ppCMHist5->GetName(), 1);
	ddCMHist5->Write(ddCMHist5->GetName(), 1);

	ppCMHist5->SaveAs("/home/zalewski/Desktop/6He/analysis/data5Out/pp_exp_counts.C");
	ddCMHist5->SaveAs("/home/zalewski/Desktop/6He/analysis/data5Out/dd_exp_counts.C");
	histOutputFile.Close();
}

void MCDataAnalyzer()
{
	//loadGeometryCorrectionParameters(setNo);
	//loadCutsCorrectionParameters(setNo);

	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small5/small5_MC.root";
	ROOT::RDataFrame smallDF("small", smallFileName.Data());

	auto newDF = smallDF.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
						.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
						.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
						.Define("tarVertex", getTarVertex, {"MWPC"})

						.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip"})
						.Define("leftLabVertex", getLeftDetPosition)
						.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
						.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
						.Define("rightLabVertex", getRightDetPosition)
						.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

						.Define("v2H", "leftGlobVertex-tarVertex")
						.Define("v6He", "rightGlobVertex-tarVertex")
						.Define("leftAngle", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
                        .Define("mc1H", getMCAngle1H, {"leftAngle", "lvBeam"})
					    .Define("mc2H", getMCAngle2H, {"leftAngle", "lvBeam"})
						.Define("rightAngle", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut", "tarVertex.X()*tarVertex.X()+tarVertex.Y()*tarVertex.Y()<9");


	TString ppHistoName5 = "ppMC_geo5_CM";
	auto ppCMHist5 = newDF.Filter("mcPP && tarCut").Histo1D({ppHistoName5.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName5 = "ddMC_geo5_CM";
	auto ddCMHist5 = newDF.Filter("mcDD && tarCut").Histo1D({ddHistoName5.Data(), "MC elastic scattering on deuterons", 90,0,180}, {"mc2H"});

	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/data5Out/MCDataCMHisto5.root", "UPDATE");
	ppCMHist5->Write(ppCMHist5->GetName(), 1);
	ddCMHist5->Write(ddCMHist5->GetName(), 1);
	histOutputFile.Close();
}

void modifiedAnalysis5()
{
	ROOT::EnableImplicitMT();
	TStopwatch *stopwatch = new TStopwatch();
	verbosity=false;

	
	//makeSmallFile(5);
    //realDataAnalyzer();
	//MCDataAnalyzer();
	makeGraph();

}