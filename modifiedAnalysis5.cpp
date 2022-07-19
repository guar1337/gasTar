#include "modifiedAnalysis5.hh"

std::vector<Double_t> getRelativeError(TH1F *my1DHist)
{
	//prepare relative errors for real data
	Int_t noBins = my1DHist->GetXaxis()->GetNbins();
	std::vector<Double_t> relativeError(noBins);
	for (int iii = 0; iii < noBins; iii++)
	{
		if (my1DHist->GetBinContent(iii)>0.0)
		{
			relativeError.at(iii) = my1DHist->GetBinError(iii)/my1DHist->GetBinContent(iii);
			//printf("Bin: %d\tbin value: %f\tbin error: %f\t%f\n", iii, my1DHist->GetBinContent(iii),my1DHist->GetBinError(iii), relativeError.at(iii)*100.0);
		}
	}
	return relativeError;
}

void setRelativeError(TH1F *myHist, std::vector<Double_t> relativeError, Double_t additionalError = 0.0)
{
	for (int iii = 0; iii < myHist->GetXaxis()->GetNbins(); iii++)
	{
		//realPPgeo1HistError.at(iii) = realPPgeo1Hist->GetBinError(iii)/realPPgeo1Hist->GetBinContent(iii);
		//printf("Error at bin[%d] = %f\n", iii, realPPgeo1HistError.at(iii));
		myHist->SetBinError(iii, myHist->GetBinContent(iii)*(relativeError.at(iii) + additionalError));
	}
}

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
	//open efficiency file, xSect file and file with models
	TFile efficiencyDataFile("/home/zalewski/Desktop/6He/analysis/data5Out/efficiency5.root", "READ");
	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/data5Out/xSection.root", "READ");
	TFile theoHistFile("/home/zalewski/Desktop/6He/analysis/dataOut/source/theoHists.root", "READ");

	//get histograms for protons
	TString efficiencyHistName = "efficiency_geo5_1H";
	TH1F *eff_geo5_1H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());
	TString histoName = "xSect_geo5_1H";	
	TH1F *xSect_geo5_1H = (TH1F*)xSectionFile.Get(histoName.Data());

	//get histograms for deuterons
	efficiencyHistName = "efficiency_geo5_2H";
	TH1F *eff_geo5_2H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());
	histoName = "xSect_geo5_2H";
	TH1F *xSect_geo5_2H = (TH1F*)xSectionFile.Get(histoName.Data());

	//get histograms with models
	TH1F *wolskiExp = (TH1F*)theoHistFile.Get("wolskiExp");
	TH1F *wolskiOM = (TH1F*)theoHistFile.Get("wolskiOM");
	TH1F *expPP = (TH1F*)theoHistFile.Get("expPP");
	TH1F *perey = (TH1F*)theoHistFile.Get("perey");
	TH1F *daehnick = (TH1F*)theoHistFile.Get("daehnick");
	TH1F *Li6d = (TH1F*)theoHistFile.Get("Li6d");
	TH1F *expDD = (TH1F*)theoHistFile.Get("expDD");
	TH1F *expDD1 = (TH1F*)theoHistFile.Get("expDD1");

	TCanvas *canvasPP = new TCanvas("canvasPP","canvasPP",10,10,1200,700);
	TString ppGraphTitle = "Differential cross-section of 6He(p,p)6He reaction";
	THStack *ppxSectStack = new THStack("PP X-section stack", ppGraphTitle.Data());
	canvasPP->SetLogy();


	trimHistogramLeft(eff_geo5_1H, 0);
	trimHistogramRight(eff_geo5_1H, 0);

	ppxSectStack->Add(xSect_geo5_1H);
	expPP->SetMarkerStyle(20);
	expPP->SetMarkerSize(1.2);
	expPP->SetMarkerColor(kRed);
	
	wolskiExp->SetMarkerStyle(3);
	wolskiExp->SetMarkerSize(2);
	wolskiExp->SetMarkerColor(kBlack);


	ppxSectStack->Add(wolskiOM, "L");
	ppxSectStack->Add(wolskiExp, "E1");
	//ppxSectStack->Add(expPP,"E1");
	
	wolskiOM->SetLineWidth(2);
	ppxSectStack->SetMinimum(0.1);
	ppxSectStack->SetMaximum(1000);
	
	ppxSectStack->Draw();
	ppxSectStack->GetXaxis()->SetTitle("#theta_{c.m.}[deg]");
	ppxSectStack->GetXaxis()->SetRangeUser(0,140);
	ppxSectStack->GetXaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->CenterTitle();

	ppxSectStack->GetYaxis()->SetTitle("d#sigma/d#omega [mb/sr]");
	ppxSectStack->Draw("HIST, nostack,p");

	TLegend *legendPP = new TLegend(0.7,0.7,0.9,0.9);
	legendPP->AddEntry(wolskiExp,"Existing data","p");
	legendPP->AddEntry(wolskiOM,"6He+p OM","L");
	//legendPP->AddEntry(expPP,"This measurement","p");
	legendPP->Draw();

	TString outPPFileName = "/home/zalewski/Desktop/6He/analysis/data5Out/pp.png";
	canvasPP->Print(outPPFileName.Data());
	canvasPP->SaveAs("/home/zalewski/Desktop/6He/analysis/data5Out/pp.C");
	delete ppxSectStack;
	delete canvasPP;

	TCanvas *canvasDD = new TCanvas("canvasDD","canvasDD",10,10,1200,700);
	TString ddGraphTitle = "Elastic scattering 6He+d";
	THStack *ddxSectStack = new THStack("DD X-section stack", ddGraphTitle.Data());
	canvasDD->SetLogy();

	trimHistogramLeft(xSect_geo5_2H, 0);
	trimHistogramRight(xSect_geo5_2H, 0);

	ddxSectStack->Add(xSect_geo5_2H);
	ddxSectStack->Add(daehnick, "L");
	ddxSectStack->Add(expDD);

	ddxSectStack->SetMinimum(0.01);
	ddxSectStack->SetMaximum(10000);
	ddxSectStack->Draw();
	ddxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ddxSectStack->GetXaxis()->CenterTitle();
	ddxSectStack->GetXaxis()->SetRangeUser(0,140);
	ddxSectStack->GetYaxis()->CenterTitle();
	ddxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ddxSectStack->Draw("HIST, nostack,p");

	TLegend *legendDD = new TLegend(0.7,0.7,0.9,0.9);
	legendDD->AddEntry(xSect_geo5_2H,"Geometry 5","p");
	legendDD->AddEntry(daehnick,"Daehnick, PRC 6, 2253 (1980)","L");
	legendDD->Draw();

	TString outDDFileName = "/home/zalewski/Desktop/6He/analysis/data5Out/dd.png";
	canvasDD->Print(outDDFileName.Data());
	canvasDD->SaveAs("/home/zalewski/Desktop/6He/analysis/data5Out/dd.C");

	delete canvasDD;
	delete ddxSectStack;
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

Double_t getCMAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi()-acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getCMAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi()-acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
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
                        .Define("mc1H", getCMAngle1H, {"leftAngle", "lvBeam"})
					    .Define("mc2H", getCMAngle2H, {"leftAngle", "lvBeam"})
						.Define("rightAngle", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut", "(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9");

	TString ppHistoName5 = "real_geo5_1H";
	auto ppCMHist5 = newDF.Filter("pp && tarCut").Histo1D({ppHistoName5.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName5 = "real_geo5_2H";
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
                        .Define("mc1H", getCMAngle1H, {"leftAngle", "lvBeam"})
					    .Define("mc2H", getCMAngle2H, {"leftAngle", "lvBeam"})
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

void calculateEfficiency(Int_t mass)
{
	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/data5Out/efficiency5.root", "UPDATE");	
	TString inFName = TString::Format("/home/zalewski/dataTmp/MC/geo5/mc_out5_%dH.root", mass);
	TString histoName = TString::Format("scattered_geo5_%dH", mass);	
	TString histoTitle = TString::Format("Theta angle scattered of {}^{%d}H in geo5 [CM deg]", mass);

	if (mass==1)
	{
		ROOT::RDataFrame inDF("analyzed", inFName.Data());
		auto scattered = inDF.Filter("(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9")
							.Define("thetaCMAngle", getCMAngle1H,{"fsqlang","lvBeam"})
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"thetaCMAngle"});
		scattered->Write(scattered->GetName(), 1);

			//pp in geo1 observed
		histoName.ReplaceAll("scattered", "observed");
		histoTitle.ReplaceAll("scattered", "observed");
		auto observed = inDF.Filter("(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9")
							.Filter("sqlde>0 && sqrde>0 && mcHe6 && mcPP")
							.Define("thetaCMAngle", getCMAngle1H,{"fsqlang","lvBeam"})
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"thetaCMAngle"});
		observed->Write(observed->GetName(), 1);

		TH1D efficiency(*scattered);
		efficiency.Divide(&(observed.GetValue()));
		histoName.ReplaceAll("observed", "efficiency");
		histoTitle = TString::Format("Efficiency of {}^{%d}H detection in geo5 [CM deg]", mass);
		efficiency.SetNameTitle(histoName.Data(), histoTitle.Data());
		efficiency.SetMarkerStyle(3);
		efficiency.SetMarkerSize(2);
		efficiency.SetMarkerColor(kBlack);
		efficiency.Write(efficiency.GetName(), 1);
	}
	
	else if (mass==2)
	{
		ROOT::RDataFrame inDF("analyzed", inFName.Data());
		auto scattered = inDF.Filter("(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9")
							.Define("thetaCMAngle", getCMAngle2H,{"fsqlang","lvBeam"})
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"thetaCMAngle"});
		scattered->Write(scattered->GetName(), 1);

			//pp in geo1 observed
		histoName.ReplaceAll("scattered", "observed");
		histoTitle.ReplaceAll("scattered", "observed");
		auto observed = inDF.Filter("(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+tarVertex.Y()*tarVertex.Y()<9")
							.Filter("sqlde>0 && sqrde>0 && mcHe6 && mcDD")
							.Define("thetaCMAngle", getCMAngle2H,{"fsqlang","lvBeam"})
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"thetaCMAngle"});
		observed->Write(observed->GetName(), 1);

		TH1D efficiency(*scattered);
		efficiency.Divide(&(observed.GetValue()));
		histoName.ReplaceAll("observed", "efficiency");
		histoTitle = TString::Format("Efficiency of {}^{%d}H detection in geo5 [CM deg]", mass);
		efficiency.SetNameTitle(histoName.Data(), histoTitle.Data());
		efficiency.SetMarkerStyle(3);
		efficiency.SetMarkerSize(2);
		efficiency.SetMarkerColor(kBlack);
		efficiency.Write(efficiency.GetName(), 1);
	}


}

void calculateXsection(Int_t mass, Double_t additionalError = 0.0)
{
	//efficiency histogram
	TString efficiencyDataPath = "/home/zalewski/Desktop/6He/analysis/data5Out/efficiency5.root";
	TFile efficiencyDataFile(efficiencyDataPath.Data(), "READ");
	if (!efficiencyDataFile.IsOpen())
	{
		printf("Couldn't open efficiency histogram file in geo5\n");
	}

	TString efficiencyHistName = TString::Format("efficiency_geo5_%dH", mass);
	TH1F *efficiencyHist = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());

	//theta integral histogram
	TH1F thetaIntegral("thetaIntegral", "thetaIntegral", 90,0,180);
	for (int iii = 0; iii < 90; iii++)
	{
		thetaIntegral.SetBinContent(iii, cos(2.0*iii*TMath::DegToRad())-cos(2.0*(iii+1)*TMath::DegToRad()));
	}

	//real data histogram
	TString realDataPath = "/home/zalewski/Desktop/6He/analysis/data5Out/realDataCMHisto5.root";
	TFile realDataFile(realDataPath.Data(), "READ");
	TString realHistName = TString::Format("real_geo5_%dH", mass);
	TH1F *realHist = (TH1F*)realDataFile.Get(realHistName.Data());

	std::vector<Double_t> relativeError = getRelativeError(realHist);

	TString xSectionHistName = TString::Format("xSect_geo5_%dH", mass);
	TString xSectionHistTitle = TString::Format("Elastic scattering cross section for {}^{6}He + {}^{%d}H reaction in geo5 [CM deg]", mass);
	TH1F xSection(*realHist);
	xSection.SetNameTitle(xSectionHistName.Data(), xSectionHistTitle.Data());
	xSection.Divide(&thetaIntegral);
	xSection.Multiply(efficiencyHist);
	if (mass == 1)
	{
		xSection.Scale(proBeamTargetConst);
	}

	else if (mass == 2)
	{
		xSection.Scale(deutBeamTargetConst);
	}
	
	setRelativeError(&xSection, relativeError, additionalError);

	xSection.SetMarkerStyle(20);
	xSection.SetMarkerSize(1);
	xSection.SetMarkerColor(kRed);

	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/data5Out/xSection.root", "UPDATE");
	xSection.Write(xSection.GetName(), 1);

	getRelativeError(&xSection);
}

void modifiedAnalysis5()
{
	ROOT::EnableImplicitMT();
	TStopwatch *stopwatch = new TStopwatch();
	verbosity=false;

	proBeamIntegral = 68.5221322e+07;
	deutBeamIntegral = 252.0620308e+07;
	proTargetThick = 1.75e+20;
	deutTargetThick = 1.71e+20;
	//	1mb = 1e-27cm^2 => 1/(beamIntegral*Tar*2*PI()) [mb]=
	//	protons:	1/(68.5221322 * 1.75 * 2 * PI()) [mb]
	//	deuterons	1/(252.0620308 * 1.71 * 2 * PI()) [mb]
	proBeamTargetConst = 1.0/(68.5221322 * 1.75 * TMath::TwoPi());
	deutBeamTargetConst = 1.0/(252.0620308 * 1.71 * TMath::TwoPi());
	printf("pro: %f\tdeu: %f\n", proBeamTargetConst, deutBeamTargetConst);
	
	//makeSmallFile(5);
    //realDataAnalyzer();
	//calculateEfficiency(protium);
	//calculateEfficiency(deuterium);
	//calculateXsection(protium);
	//calculateXsection(deuterium);
	makeGraph();

}