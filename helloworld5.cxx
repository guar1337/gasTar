// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//	 this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//	 this list of conditions and the following disclaimer in the documentation
//	 and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//	 used to endorse or promote products derived from this software without
//	 specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


#include "helloworld5.hxx"

std::vector<Double_t> getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_1_Z, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t m_MWPC_2_Z)
{
	std::vector<Double_t> rvecMWPC(6);
	rvecMWPC.at(0) = m_MWPC_1_X;
	rvecMWPC.at(1) = m_MWPC_1_Y;
	rvecMWPC.at(2) = m_MWPC_1_Z;
	rvecMWPC.at(3) = m_MWPC_2_X;
	rvecMWPC.at(4) = m_MWPC_2_Y;
	rvecMWPC.at(5) = m_MWPC_2_Z;

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
	Double_t m_tarPos = 0.0;
	Double_t m_tarAngle = 33.0 * TMath::DegToRad();

	Double_t m_dX = rvecMWPC[3] - rvecMWPC[0];
	Double_t m_dY = rvecMWPC[4] - rvecMWPC[1];
	Double_t m_dZ = rvecMWPC[5] - rvecMWPC[2];
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, m_tarPos);
	TVector3 m_beamPoint(rvecMWPC[3], rvecMWPC[4], rvecMWPC[5]);
	TVector3 m_tarPerpendicular(sin(m_tarAngle), 0.0, cos(m_tarAngle));
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
	Double_t X2hDet = cs::widthStripX * xStrip * cos(leftAngle);
	Double_t Y2hDet = cs::widthStripY * yStrip;
	Double_t Z2hDet = -cs::widthStripX * xStrip * sin(leftAngle);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * xStrip * cos(rightAngle);
	Double_t Y6HeDet = cs::widthStripY * yStrip;
	Double_t Z6HeDet = cs::widthStripX * xStrip * sin(rightAngle);
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
	Double_t X6HeCorrection = sqrDistCorrection*sin(-rightAngleCorrection);
	Double_t Z6HeCorrection = sqrDistCorrection*cos(-rightAngleCorrection);

	Double_t X6Helab = sqrDist*sin(-rightAngle) + X6HeCorrection - (cs::sqrXzero) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart;
	Double_t Z6Helab = sqrDist*cos(rightAngle) + Z6HeCorrection - (cs::sqrXzero) * sin(rightAngle);
	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

Double_t myAngAngFit(Double_t m_leftAngle, TLorentzVector m_lvBeam, Double_t tarMass)
{
	TLorentzVector m_lvBeamCopy(m_lvBeam);
	TLorentzVector m_TarCM{0.0,0.0,0.0,tarMass};
	TLorentzVector m_lvCM = m_lvBeam + m_TarCM;
	Double_t m_thetaCM = m_lvCM.Theta();
	TVector3 boostVect = m_lvCM.BoostVector();
	
	Double_t gammaSquare = m_lvCM.Gamma() *  m_lvCM.Gamma();
	Double_t tanSquare = pow(tan(m_leftAngle*TMath::DegToRad()),2);
	Double_t cosLeftAng = (1.0 - gammaSquare*tanSquare)/(1 + gammaSquare*tanSquare);
	Double_t thetaCM = TMath::Pi() - (acos(cosLeftAng)+m_thetaCM);

	m_lvBeam.Boost(-boostVect);
	//printf("TMath::Pi(): %f\tthetaCM: %f\tlAng: %f\n",TMath::Pi()*TMath::RadToDeg(), thetaCM*TMath::RadToDeg(), m_leftAngle);
	m_lvBeam.SetTheta(thetaCM);
	m_lvBeam.Boost(boostVect);
	/*
	m_TarCM.Boost(-boostVect);
	m_TarCM.SetTheta(acos(cosLeftAng)-m_Theta);
	m_TarCM.Boost(boostVect);
	//printf("sqlangIN: %f\tsqlangOUT: %f\tdiff: %f\tm_Theta: %f\n",m_leftAngle, m_TarCM.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg(), m_leftAngle-m_TarCM.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg(),m_Theta*TMath::RadToDeg());
	*/
	return m_lvBeam.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg();
}

struct CostFunctor
{
	bool operator()(const double *par, double *residual) const
	{
		rightAngle = (9.0 + par[sRang]) * TMath::DegToRad();
		rightAngleCorrection =  rightAngle + 1.65 * TMath::DegToRad();

		ROOT::RDataFrame smallDF("small", smallFile);

		auto newDF = smallDF.Define("X1",[&par](Double_t MWPC_1_X){return (MWPC_1_X + par[sMWPC_1_X]);}, {"MWPC_1_X"})
							.Define("Y1",[&par](Double_t MWPC_1_Y){return (MWPC_1_Y + par[sMWPC_1_Y]);}, {"MWPC_1_Y"})
							.Define("Z1",[&par](){return -816.0;})
							.Define("X2",[&par](Double_t MWPC_2_X){return (MWPC_2_X + par[sMWPC_2_X]);}, {"MWPC_2_X"})
							.Define("Y2",[&par](Double_t MWPC_2_Y){return (MWPC_2_Y + par[sMWPC_2_Y]);}, {"MWPC_2_Y"})
							.Define("Z2",[&par](){return -270.0;})
							.Define("MWPC", getMWPC, {"X1", "Y1", "Z1", "X2", "Y2", "Z2"})
							.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
							.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
							.Define("tarVertex", getTarVertex, {"MWPC"})

							.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip"})
							.Define("leftLabVertex", getLeftDetPosition)
							.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
							.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
							.Define("rightLabVertex", getRightDetPosition)
							.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

							.Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
							.Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
							.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
							.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"});

		auto ppDF = newDF.Filter("pp").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
		auto ddDF = newDF.Filter("dd").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});

		// for pp and dd sqlang describes the ang-ang relation. For the pr reaction I might have to divide the data into two ranges
		auto nEventsPP = ppDF.Count().GetValue()/(1000);
		auto sumPP = ppDF.Define("resqrang",[](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass1H);}, {"sqlang", "lvBeam"})
							.Define("difSqrang","pow(resqrang-sqrang,2)")
							.Sum<double>("difSqrang");

		auto nEventsDD = ddDF.Count().GetValue()/(1000);
		auto sumDD = ddDF.Define("resqrang",[](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass2H);}, {"sqlang", "lvBeam"})
							.Define("difSqrang","pow(resqrang-sqrang,2)")
							.Sum<double>("difSqrang");

		auto histTargetX = newDF.Filter("pp || dd").Define("evX", "tarVertex.X()").Histo1D({"tarVertexX", "Fitting X dimension of target", 100,-30,30}, {"evX"});
		histTargetX->Fit("gaus", "Q");
		TF1 *tmpFitX = (TF1*)histTargetX->GetListOfFunctions()->FindObject("gaus");
		Double_t targetXPos = tmpFitX->GetParameter(1);
		delete tmpFitX;

		auto histTargetY = newDF.Filter("pp || dd").Define("evY", "tarVertex.Y()").Histo1D({"tarVertexY", "Fitting Y dimension of target", 100,-30,30}, {"evY"});
		histTargetY->Fit("gaus", "Q");
		TF1 *tmpFitY = (TF1*)histTargetY->GetListOfFunctions()->FindObject("gaus");
		Double_t targetYPos = tmpFitY->GetParameter(1);
		delete tmpFitY;

		Double_t normChi2PP = sumPP.GetValue()/nEventsPP;
		Double_t normChi2DD = sumDD.GetValue()/nEventsDD;

		residual[0] = normChi2PP;
		residual[1] = normChi2DD;
		//residual[2] = targetXPos;
		//residual[3] = targetYPos;

		return true;
	}
};

void drawResults(double *finalPars, int myRunNumber = 0)
{
	rightAngle = (9.0 + finalPars[sRang]) * TMath::DegToRad();
	rightAngleCorrection =  rightAngle + 1.65 * TMath::DegToRad();

	ROOT::EnableImplicitMT();
	ROOT::RDataFrame smallDF("small", smallFile);

	auto newDF = smallDF.Define("X1",[&finalPars](Double_t MWPC_1_X){return (MWPC_1_X + finalPars[sMWPC_1_X]);}, {"MWPC_1_X"})
						.Define("Y1",[&finalPars](Double_t MWPC_1_Y){return (MWPC_1_Y + finalPars[sMWPC_1_Y]);}, {"MWPC_1_Y"})
						.Define("Z1",[&finalPars](){return -816.0;})
						.Define("X2",[&finalPars](Double_t MWPC_2_X){return (MWPC_2_X + finalPars[sMWPC_2_X]);}, {"MWPC_2_X"})
						.Define("Y2",[&finalPars](Double_t MWPC_2_Y){return (MWPC_2_Y + finalPars[sMWPC_2_Y]);}, {"MWPC_2_Y"})
						.Define("Z2",[&finalPars](){return -270.0;})
						.Define("MWPC", getMWPC, {"X1", "Y1", "Z1", "X2", "Y2", "Z2"})
						.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
						.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
						.Define("tarVertex", getTarVertex, {"MWPC"})

						.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip"})
						.Define("leftLabVertex", getLeftDetPosition)
						.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
						.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
						.Define("rightLabVertex", getRightDetPosition)
						.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

						.Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
						.Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
						.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
						.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"});


	auto ppDF = newDF.Filter("pp").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
	auto ddDF = newDF.Filter("dd").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});


	auto nEventsPP = ppDF.Count().GetValue();
	auto sumPP = ppDF.Define("resqrang",[](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass1H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");
						
	auto nEventsDD = ddDF.Count().GetValue();
	auto sumDD = ddDF.Define("resqrang",[](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass2H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");


	auto histTargetX = newDF.Filter("pp || dd").Define("evX", "tarVertex.X()").Histo1D({"tarVertexX", "Fitting X dimension of target", 100,-30,30}, {"evX"});
	histTargetX->Fit("gaus", "Q");
	TF1 *tmpFitX = (TF1*)histTargetX->GetListOfFunctions()->FindObject("gaus");
	Double_t targetXPos = tmpFitX->GetParameter(1);

	auto histTargetY = newDF.Filter("pp || dd").Define("evY", "tarVertex.Y()").Histo1D({"tarVertexY", "Fitting Y dimension of target", 100,-30,30}, {"evY"});
	histTargetY->Fit("gaus", "Q");
	TF1 *tmpFitY = (TF1*)histTargetY->GetListOfFunctions()->FindObject("gaus");
	Double_t targetYPos = tmpFitY->GetParameter(1);

	TGraph ppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	
	TGraph fppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	TGraph fddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();


	TCanvas myCanvas("myCanvas", "Minimize results", 1200, 800);
	myCanvas.Divide(2,2);

	myCanvas.cd(1);
	ppAngAng.GetXaxis()->SetLimits(50.0, 90.0);
	ppAngAng.GetYaxis()->SetRangeUser(0.0, 15.0);
	ppAngAng.GetXaxis()->SetTitle("Angle of {}^{6}He [LAB deg]");
	ppAngAng.GetYaxis()->SetTitle("Angle of {}^{1}H [LAB deg]");
	ppAngAng.SetTitle("Angle-angle relation for elastic scattering");
	ppAngAng.Draw("AP");
	fppAngAng.SetMarkerColor(kRed);
	fppAngAng.Draw("P,same");

/*
	TLatex *tex1 = new TLatex(26.05426,10.62081,"Elastic scattering on protons");
   tex1->SetLineWidth(2);
   tex1->Draw("same");
	TLatex *tex2 = new TLatex(39.47308,20.15142,"Elastic scattering on deuterons");
   tex2->SetLineWidth(2);
   tex2->Draw("same");
*/
	myCanvas.cd(2);
	ddAngAng.GetXaxis()->SetLimits(50.0, 90.0);
	ddAngAng.GetYaxis()->SetRangeUser(0.0, 25.0);
	ddAngAng.GetXaxis()->SetTitle("Angle of {}^{6}He [LAB deg]");
	ddAngAng.GetYaxis()->SetTitle("Angle of {}^{2}H [LAB deg]");
	ddAngAng.SetTitle("Angle-angle relation for elastic scattering");
	ddAngAng.Draw("AP");
	fddAngAng.SetMarkerColor(kRed);
	fddAngAng.Draw("P,same");

	myCanvas.cd(3);
	histTargetX->Draw();

	myCanvas.cd(4);
	histTargetY->Draw();


	myCanvas.Print(TString::Format("/home/zalewski/Desktop/6He/analysis/geo5/MWPC_rAng_pp_dd/pp_dd_%d.png", myRunNumber));
	//myCanvas.SaveAs("/home/zalewski/Desktop/macro.C");


//gPad->Print("/home/zalewski/Desktop/myPP.jpeg");

}

int main(int argc, char** argv)
{
	firstRun = (argc > 1) ? atoi(argv[1]) : 100;
	//Int_t mySeed = (argc > 2) ? atoi(argv[2]) : 0;
	gRandom->SetSeed();
	printf("Passed seed: %d\n", gRandom->GetSeed());
	calculate=true;
	google::InitGoogleLogging(argv[0]);
	// The variable to solve for with its initial value. It will be
	// mutated in place by the solver.
	

	// Build the problem.
	ceres::Problem problem;

	// Set up the only cost function (also known as residual). This uses
	// auto-differentiation to obtain the derivative (jacobian).
	ceres::NumericDiffOptions diffOpts;
	diffOpts.relative_step_size = 0.1;
	double x[numberOfParameters];
	ceres::CostFunction* cost_function =
		new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, 2, numberOfParameters>(new CostFunctor, ceres::TAKE_OWNERSHIP, 2, diffOpts);
//                                                           |        |   |
//                               Finite Differencing Scheme -+        |   |
//                               Dimension of residual ---------------+   |
//                               Dimension of x --------------------------+
	problem.AddResidualBlock(cost_function, nullptr, x);

	lowerParamBound[sMWPC_1_X] = -5.0;
	lowerParamBound[sMWPC_1_Y] = -5.0;
	lowerParamBound[sMWPC_2_X] = -5.0;
	lowerParamBound[sMWPC_2_Y] = -5.0;
	lowerParamBound[sRang] = -2.0;


	upperParamBound[sMWPC_1_X] = 5.0;
	upperParamBound[sMWPC_1_Y] = 5.0;
	upperParamBound[sMWPC_2_X] = 5.0;
	upperParamBound[sMWPC_2_Y] = 5.0;
	upperParamBound[sRang] = 2.0;


	for (int iii = 0; iii < numberOfParameters; iii++)
	{
		problem.SetParameterLowerBound(x, iii, lowerParamBound[iii]);
		problem.SetParameterUpperBound(x, iii, upperParamBound[iii]);
	}

	// Run the solver!
	ceres::Solver::Options options;
	options.trust_region_strategy_type = ceres::DOGLEG;
	options.dogleg_type = ceres::SUBSPACE_DOGLEG;
	options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = true;
	options.num_threads = 8;
	options.max_num_iterations = 50;
	ceres::Solver::Summary summary;

	if (calculate)
	{
		std::ofstream outStreamGenerated(outGeneratedParams, std::ios::out | std::ofstream::app);
		std::ofstream outStreamObtained(outObtainedParams, std::ios::out | std::ofstream::app);
		for (auto &&iii : parNames)
		{
			//outStreamGenerated<<iii<<" ";
			//outStreamObtained<<iii<<" ";
		}

		//outStreamGenerated<<std::endl;
		//outStreamObtained<<std::endl;

		outStreamGenerated<<firstRun<<"\t";
		//randomly set initial value of the parameters
		for (int jjj = 0; jjj < numberOfParameters; jjj++)
		{	
			x[jjj] = gRandom->Uniform(lowerParamBound[jjj], upperParamBound[jjj]);
			//print the generated parameters to the file				
			outStreamGenerated<<x[jjj]<<" ";
		}
		outStreamGenerated<<std::endl;

		//solve for the generated values
		Solve(options, &problem, &summary);

		outStreamObtained<<firstRun<<"\t";
		for (int jjj = 0; jjj < numberOfParameters; jjj++)
		{
			//print the obtained parameters to the file				
			outStreamObtained<<x[jjj]<<" ";
		}
		//save the graph, generated parameters and final parameters
		outStreamObtained<<std::endl;
		outStreamGenerated.close();
		outStreamObtained.close();
	}
	drawResults(x, firstRun);

	//std::cout << summary.FullReport() << "\n";

	return 0;
}