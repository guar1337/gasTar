#include "/home/zalewski/aku/gasTar/ui5.hh"

R__LOAD_LIBRARY(libgsl.so);
R__LOAD_LIBRARY(/home/zalewski/aku/ELC/build/libEloss.so);
//R__LOAD_LIBRARY(/home/zalewski/aku/TELoss/libTELoss.so);

void makeSmallMC()
{
	TChain MCChain("analyzed");
	MCChain.Add("/home/zalewski/dataTmp/MC/geo5/mc_out5*_big*.root");
	ROOT::RDataFrame MCDF(MCChain);
	auto smallDF = MCDF.Filter("he6 && sqlde>0.1 && (mcPP || mcDD)");
	smallDF.Snapshot("small", "/home/zalewski/dataTmp/small5/small5_MC.root", MCcolumnList);
}

void joindE()
{
	printf("\nJoining files in ...dE/geo5\n");
	TChain dEChain("analyzed");
	dEChain.Add("/home/zalewski/dataTmp/dE/geo5/dE_he6_*.root");
	ROOT::RDataFrame dEDF(dEChain);
	dEDF.Snapshot("analyzed", "/home/zalewski/dataTmp/dE/geo5/dE_he6.root");

	//graphical cuts (pp && dd) are based on MWPC shifts = 0.0 and right detector at 9 deg for 132mm and (9+1.65)deg at 168mm
	auto smallDF = dEDF.Filter("he6 && (pp || dd)");
	//smallDF.Snapshot("small", "/home/zalewski/dataTmp/small5/small5_sizeChck.root", columnList);
}

Double_t getMCAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng)-m_ThetaCM2H)*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getMCAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng)-m_ThetaCM2H)*TMath::RadToDeg();
	return m_sqlangCM;
}

//d,t reaction
Double_t getDTLang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass5He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	m_lv3H.SetTheta(m_thetaCM);
	m_lv3H.Boost(m_boostVect);

	return m_lv3H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getDTRang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm5HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm5Heene = m_cm5HekinE + cs::mass5He;
	Double_t m_cm5Hemom = sqrt(m_cm5Heene*m_cm5Heene - cs::mass5He*cs::mass5He);

	TLorentzVector m_lv5He(0.0,0.0,m_cm5Hemom, m_cm5Heene);

	m_lv5He.SetTheta(TMath::Pi() - m_thetaCM);
	m_lv5He.Boost(m_boostVect);
	return m_lv5He.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getDDEne(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0,0.0,0.0,cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_sqlangRad = m_sqlang * TMath::DegToRad();
	Double_t m_gamma2 = m_lvCM.Gamma() * m_lvCM.Gamma();
	Double_t m_tanSqlang2 = tan(m_sqlangRad) * tan(m_sqlangRad);
	Double_t m_thetaCM = acos((1-m_gamma2*m_tanSqlang2)/(1+m_gamma2*m_tanSqlang2));

	m_lvBeam.Boost(-m_vBoost);
	m_lvBeam.SetTheta(TMath::Pi()-m_thetaCM);
	m_lvBeam.Boost(m_vBoost);
	return m_lvBeam.E() - m_lvBeam.M();
}

//p,t reaction
Double_t getPTLang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + m_Q;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);
	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	m_lv3H.SetTheta(m_thetaCM);
	m_lv3H.Boost(m_boostVect);

	return m_lv3H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getPTRang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + m_Q;

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	m_lv4He.SetTheta(TMath::Pi() - m_thetaCM);
	m_lv4He.Boost(m_boostVect);
	return m_lv4He.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t recoPT(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	//generate pt reaction (lv3H actually) with:
	//			theta in range 0:3.14
	//			beam vector
	//			q of the reaction
	Int_t multi = (rnd->Integer(2)==1) ? 1 : -1;
	TLorentzVector m_lv6He(m_lvBeam);
	// /m_lv6He.SetTheta(0.0);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + cs::Qpt;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	Double_t m_betaCM = m_lvCM.Beta();
	Double_t m_gammaCM = m_lvCM.Gamma();
	Double_t m_beta3HCM = m_lv3H.Beta();
	Double_t m_beta4HeCM = m_lv4He.Beta();
	//m_lv3H.SetTheta(m_thetaCM);
	//m_lv3H.Boost(m_boostVect);

	//reconstruct pt reaction (angle) knowing:
	//			beta and gamma of the CM, beta of 3H in CM, beta of 4He in CM
	//			angle between 3H and beam vectors
	Double_t m_betaCM3H = (sqrt(m_cm3HkinE*m_cm3HkinE+2*m_cm3HkinE*cs::mass3H))/(m_cm3HkinE+cs::mass3H);
	Double_t m_betaCM4He = (sqrt(m_cm4HekinE*m_cm4HekinE+2*m_cm4HekinE*cs::mass4He))/(m_cm4HekinE+cs::mass4He);
	Double_t m_beta3HRatio = m_betaCM/m_betaCM3H;
	Double_t m_beta4HeRatio = m_betaCM/m_betaCM4He;
	Double_t m_beta3HRatio2 = m_beta3HRatio*m_beta3HRatio;

	//calculating CM angles
	Double_t m_gammaTan2 = std::pow(m_gammaCM*tan(m_sqlang*TMath::DegToRad()),2);
	Double_t m_cosThetaCM3H = (-m_beta3HRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta3HRatio2)))/(1.0+m_gammaTan2);
	Double_t m_thetaCM3H = acos(m_cosThetaCM3H);
	Double_t m_thetaCM4He = TMath::Pi() - m_thetaCM3H;
	Double_t m_4He_lab = atan(sin(m_thetaCM4He)/(m_gammaCM*(cos(m_thetaCM4He)+m_betaCM/m_betaCM4He)));
	Double_t m_3H_lab = atan(sin(m_thetaCM3H)/(m_gammaCM*(cos(m_thetaCM3H)+m_betaCM/m_betaCM3H)));
	return m_4He_lab*TMath::RadToDeg();
}

//p,p reaction
Double_t getPPLang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0,0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv1H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv1H.SetTheta(m_thetaCM);

	//m_lv6He.Boost(m_boostVect);
	m_lv1H.Boost(m_boostVect);
	Double_t m_sqlangpp = 180.0 * (m_lvBeam.Angle(m_lv1H.Vect()))/double(TMath::Pi());
	//Double_t m_sqrangpp = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqlangpp;
}

Double_t getPPRang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0,0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	//m_lv1H.Boost(-m_boostVect);

	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	//m_lv1H.SetTheta(m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	//m_lv1H.Boost(m_boostVect);
	//Double_t m_sqlangpp = 180.0 * (m_lvBeam.Angle(m_lv1H.Vect()))/double(TMath::Pi());
	Double_t m_sqrangpp = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqrangpp;
}

//d,d reaction
Double_t getDDRang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	//m_lv2H.Boost(-m_boostVect);

	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	//m_lv2H.SetTheta(m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	//m_lv2H.Boost(m_boostVect);
	//Double_t m_sqlangdd = 180.0 * (m_lvBeam.Angle(m_lv2H.Vect()))/double(TMath::Pi());
	Double_t m_sqrangdd = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqrangdd;
}

Double_t getDDLang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv2H.SetTheta(m_thetaCM);

	//m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);
	Double_t m_sqlangdd = 180.0 * (m_lvBeam.Angle(m_lv2H.Vect()))/double(TMath::Pi());
	//Double_t m_sqrangdd = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqlangdd;
}

Int_t getStripNumber(ROOT::VecOps::RVec<Double_t> &inputArray)
{
	return std::distance(inputArray.begin(),std::max_element(inputArray.begin(), inputArray.end()));
}

TVector3 getLeftDetVertex_5(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (rnd->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (rnd->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system	
	Double_t X2hDet = cs::widthStripX * xStrip * cos(sqlAng);
	Double_t Y2hDet = cs::widthStripY * yStrip;
	Double_t Z2hDet = -cs::widthStripX * xStrip * sin(sqlAng);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (rnd->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (rnd->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * xStrip * cos(sqrAng);
	Double_t Y6HeDet = cs::widthStripY * yStrip;
	Double_t Z6HeDet = cs::widthStripX * xStrip * sin(sqrAng);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

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
	Double_t m_tarAngle = 33.0*TMath::DegToRad();

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

Double_t getKineticEnergy(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	return cs::mass6He*(m_gamma-1.0);
}

bool filterLeftDetector(ROOT::VecOps::RVec<Double_t> &leftDetectorArray)
{
	int aboveThresholdEventCounter(0);
	auto leftDetectorArraySize = leftDetectorArray.size();
	for (size_t iii = 0; iii < leftDetectorArraySize; iii++)
	{
		if (leftDetectorArray[iii]>0.1)
		{
			aboveThresholdEventCounter++;
		}
	}

	if (aboveThresholdEventCounter==1)
	{
		return true;
	}

	else
	{
		return false;
	}
	
}

bool filterRightDetector(ROOT::VecOps::RVec<Double_t> &rightDetectorArray)
{
	int aboveThresholdEventCounter(0);
	auto rightDetectorArraySize = rightDetectorArray.size();
	for (size_t iii = 0; iii < rightDetectorArraySize; iii++)
	{
		if (rightDetectorArray[iii]>5.0)
		{
			aboveThresholdEventCounter++;
		}
	}
	if (aboveThresholdEventCounter==1)
	{
		return true;
	}

	else
	{
		return false;
	}
	
}

struct filterMWPC
{
	int MWPC_ID;
	//create object containing ID of the column which is to be calibrated
	filterMWPC(int MWPC_number) 
	{
		MWPC_ID = MWPC_number;
	};

	bool operator()(ROOT::VecOps::RVec<unsigned short> &MWPCHitsArray, unsigned short MWPCHitsCount)
	{		
		if (MWPCHitsCount!=0)
		{
			int sizeof_clust = 1;
			//one wire got signal - simplest solution
			if (MWPCHitsCount==1)
			{
				//printf("simplest solution %i\n", MWPCHitsArray[0]);
				if (MWPC_ID == 0 || MWPC_ID == 2)
				{
					//*MWPC_pos = (15.5 - MWPCHitsArray[0])*1.25;
				}

				else
				{
					//*MWPC_pos = (-15.5 + MWPCHitsArray[0])*1.25;
				}
				return 1;
			}

			//more than one but is it one cluster?
			else
			{	//checking...
				for (int iii = 1; iii < MWPCHitsCount; iii++)
				{
					if ((MWPCHitsArray[iii] - MWPCHitsArray[iii-1])==1)
					{
						sizeof_clust++;
					}
				}
				//
				if (sizeof_clust==MWPCHitsCount)
				{

					if (MWPC_ID == 0 || MWPC_ID == 2)
					{
						//*MWPC_pos = (15.5 - (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25;
					}

					else
					{
						//*MWPC_pos = (-15.5 + (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25;
					}
					return 1;

				}
				else if (sizeof_clust<MWPCHitsCount)
				{
					return 0;
				}
			}
		}
		else
		{
			return 0;
		}
		return 0;
	}

};

struct getMWPCpos
{
	int MWPC_ID;
	//create object containing ID of the column which is to be calibrated
	getMWPCpos(int MWPC_number) 
	{
		MWPC_ID = MWPC_number;
	};

	Double_t operator()(ROOT::VecOps::RVec<unsigned short> &MWPCHitsArray, unsigned short MWPCHitsCount)
	{		
		if (MWPCHitsCount!=0)
		{
			int sizeof_clust = 1;
			//one wire got signal - simplest solution
			if (MWPCHitsCount==1)
			{
				//printf("simplest solution %i\n", MWPCHitsArray[0]);
				if (MWPC_ID == 0 || MWPC_ID == 2)
				{
					return (15.5 - MWPCHitsArray[0])*1.25 + rnd->Uniform(0.0,1.25)-0.625;
				}

				else
				{
					return (-15.5 + MWPCHitsArray[0])*1.25 + rnd->Uniform(0.0,1.25)-0.625;
				}
			}

			//more than one but is it one cluster?
			else
			{	//checking...
				for (int iii = 1; iii < MWPCHitsCount; iii++)
				{
					if ((MWPCHitsArray[iii] - MWPCHitsArray[iii-1])==1)
					{
						sizeof_clust++;
					}
				}
				//
				if (sizeof_clust==MWPCHitsCount)
				{

					if (MWPC_ID == 0 || MWPC_ID == 2)
					{
						return (15.5 - (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25 + rnd->Uniform(0.0,1.25)-0.625;
					}

					else
					{
						return (-15.5 + (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25 + rnd->Uniform(0.0,1.25)-0.625;
					}
				}
				else if (sizeof_clust<MWPCHitsCount)
				{
					printf("Rogue MWPC event got through\tcluster size = %d\tMWPC hits = %d\n",sizeof_clust, MWPCHitsCount);
				}
			}
		}
		else
		{
			std::cout<<"Rogue MWPC event got through"<<std::endl;
		}
		return 0;
	}

};

bool filterToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return (ToF>150 && ToF<185);
}

Double_t calculateToF(ROOT::VecOps::RVec<Double_t> &timeF3, ROOT::VecOps::RVec<Double_t> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1]+timeF5[2]+timeF5[3])/4.0);
	Double_t ToF = (tF5-tF3)+ToFconstant;
	return ToF;
}

bool filterSQL(ROOT::VecOps::RVec<unsigned short> &SQ, Double_t threshold)
{
	Bool_t sqFlag=false;
	Double_t calibratedSQ;

	for (int iii = 0; iii < 32; iii++)
	{
		//vecAcoef.at(0) is array with calibration parameters of SQX_L
		calibratedSQ = vecAcoef.at(0).at(iii) + vecBcoef.at(0).at(iii) * SQ[iii];
		
		if (calibratedSQ>threshold)
		{
			sqFlag=true;
		}
	}
	return sqFlag;
}

bool filterSQR(ROOT::VecOps::RVec<unsigned short> &SQ, Double_t threshold)
{
	Double_t calibratedSQ;
	Bool_t sqFlag=false;

	for (int iii = 0; iii < 32; iii++)
	{
		//vecAcoef.at(3) is array with calibration parameters of SQX_R
		calibratedSQ = vecAcoef.at(3).at(iii) + vecBcoef.at(3).at(iii) * SQ[iii];
				
		if (calibratedSQ>threshold)
		{
			sqFlag=true;
		}
	}
	return sqFlag;
}

void loadCorrectionParameters()
{
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/experimental/good/pics/db";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("Failed to open file: %s\n", fName.c_str());
	}

	int jumpTo = cs::fileToProcess + (int)'a';
	std::getline(outStreamGenerated, line, (char)jumpTo);
	printf("%s\n", line.c_str());

	float tmpContainer;
	for (int iii = 0; iii < 11; iii++)
	{
		outStreamGenerated>>tmpContainer;
		parameters.push_back(tmpContainer);
		printf("%s = %f\n", parNames[iii].c_str(), parameters[iii]);
	}

}

void loadCalibrationParameters()
{
	TString fileName;
	std::string dummy;
	Double_t tmpA, tmpB, tmpC;

	for (auto &&iii : vecCalibratedColumns)
	{
		fileName = "/home/zalewski/dataTmp/calibrationParameters/geo5/" + iii + ".cal";

		//std::cout<<fileName<<std::endl;

		std::ifstream instream(fileName.Data());

		if (instream.is_open()) 
		{
			std::getline(instream, dummy);
			int chNo = std::stoi( dummy );
			Acoef.clear();
			Bcoef.clear();

			for (int jjj = 0; jjj < chNo; jjj++)
			{
				instream>>tmpA>>tmpB>>tmpC;
				Acoef.push_back(tmpA);
				Bcoef.push_back(tmpB);
				//std::cout<<tmpA<< "\t" <<tmpB<<std::endl;
			}

			vecAcoef.push_back(Acoef);
			vecBcoef.push_back(Bcoef);
			instream.close();
		}

		else
		{
			printf ("#Cannot open %s coefficient file\n",
			fileName.Data());
		}
	}
}

void cleaner(TString inFileName)
{
	Double_t MWPCposContainer;
	Double_t sqlThreshold = 0.5;
	Double_t sqrThreshold = 2.0;

	inFileName.ReplaceAll("raw","simp");
	std::cout<<"Cleaning "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("simplified", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("simp","cln");
	//get rid of broken MWPC events
	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
					.Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
					.Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
					.Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
					.Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
					.Filter(filterMWPC(0),{"x1","nx1"})
					.Filter(filterMWPC(1),{"y1","ny1"})
					.Filter(filterMWPC(2),{"x2","nx2"})
					.Filter(filterMWPC(3),{"y2","ny2"})
					.Filter(filterToF,{"tdcF3", "tdcF5"})
					.Filter([sqlThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_L){return filterSQL(SQX_L, sqlThreshold);},{"SQX_L"})
					.Filter([sqrThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_R){return filterSQR(SQX_R, sqrThreshold);},{"SQX_R"});


	outDF.Snapshot("cleaned", outFilename.Data());
}

struct Calibrator
{
	int columnID;
	//create object containing ID of the column which is to be calibrated
	Calibrator(int columnToCalibrate) 
	{
		columnID = columnToCalibrate;
	};

	ROOT::VecOps::RVec<Double_t> operator()(ROOT::VecOps::RVec<unsigned short> &inputArray)
	{	
		auto mySize = inputArray.size();
		ROOT::VecOps::RVec<Double_t> outputArray(mySize);
		for (size_t iii = 0; iii < mySize; iii++)
		{
			//std::cout<<vecAcoef.at(columnID).at(iii)<<" "<<"My ID is: "<<vecCalibratedColumns.at(columnID)<<std::endl;
			outputArray.at(iii) = vecAcoef.at(columnID).at(iii) + vecBcoef.at(columnID).at(iii) * inputArray[iii];
		}

		return outputArray;
	}
};

ROOT::RDF::RNode ApplyDefines(ROOT::RDF::RNode df, Int_t fNo, const std::vector<std::string> &colNames, unsigned int iii = 0)
{
	//printf("passed variable: %d\n", fNo);
	//recursive function updating RDataFrame over scope of std::vector with column names for calibration
	if (iii == colNames.size())
	{

		return df.Define("tof", calculateToF, {"cal_tdcF3", "cal_tdcF5"})
				 .Define("aF5", "(F5[0]+F5[1]+F5[2]+F5[3])/4.0")
				 .Define("MWPC_1_X", getMWPCpos(0),{"x1","nx1"})
				 .Define("MWPC_1_Y", getMWPCpos(1),{"y1","ny1"})
				 .Define("MWPC_2_X", getMWPCpos(2),{"x2","nx2"})
				 .Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"})
				 .Define("geo", "5")
				 .Define("fNo", [fNo](){return fNo;});
	}

	std::string inputColumn = colNames[iii];
	inputColumn.insert(0, "cal_");
	//std::cout<<colNames[iii]<<"\t"<<inputColumn<<std::endl;
	//calling 'Calibrator' so I can pass other info
	return ApplyDefines(df.Define(inputColumn, Calibrator(iii), {colNames[iii]}), fNo, colNames, iii + 1);
}

void calibratorCaller(TString inFileName)
{	
	inFileName.ReplaceAll("raw","cln");
	std::cout<<"Calibrating "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("cleaned", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cln","cal");
	// Print columns' names
	Int_t fileNumber = inFileName.ReplaceAll("/home/zalewski/dataTmp/cal/geo5/cal_he6_","").Remove(2).Atoi();
	auto dfWithDefines = ApplyDefines(inDF, fileNumber, vecCalibratedColumns);
	ROOT::RDF::RSnapshotOptions myOpts;
	dfWithDefines.Snapshot("calibrated", outFilename.Data(), dfWithDefines.GetColumnNames(), myOpts);
}

void translator(TString inFileName)
{
	std::cout<<"Translating "<<inFileName<<std::endl;
	TFile inputFile(inFileName.Data(), "READ");
	TTree *inputTree;
	inputTree = (TTree*)inputFile.Get("AnalysisxTree");
	TString outFilename = inFileName.ReplaceAll("raw","simp");

	inputTree->Process(selector, outFilename.Data());
	
	inputFile.Close();	
}

void analysis5(TString inFileName)
{
	inFileName.ReplaceAll("raw","cal");
	TString treeName = "calibrated";
	ROOT::RDataFrame inDF(treeName.Data(), inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cal","dE").ReplaceAll("mc","mc_out");
	std::cout<<"Analysing "<<inFileName<<std::endl;

	
	tarThickness = 160 + 2*cs::tarThicknessShift;
	tarAngle = 33.0 * TMath::DegToRad();
	sqlAng = (70.0) * TMath::DegToRad();
	sqlDist = 170.0;

	sqrAng = (9.0 + 1.25) * TMath::DegToRad();
	sqrDist = 132.0;
	Double_t sqrDistCorrection = 168.0;	//1.65
	Double_t sqrAngCorrection = sqrAng + (1.65) * TMath::DegToRad();
	

	Double_t X2Hlab = sqlDist*sin(sqlAng) - (cs::sqlXzero) * cos(sqlAng);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(sqlAng) + (cs::sqlXzero) * sin(sqlAng);
	TVector3 leftDetPosition(X2Hlab, Y2Hlab, Z2Hlab);

	Double_t X6HeCorrection = sqrDistCorrection*sin(-sqrAngCorrection);
	Double_t Z6HeCorrection = sqrDistCorrection*cos(-sqrAngCorrection);

	Double_t X6Helab = sqrDist*sin(-sqrAng) + X6HeCorrection - (cs::sqrXzero) * cos(sqrAng);
	Double_t Y6Helab = cs::sqrYstart;
	Double_t Z6Helab = sqrDist*cos(sqrAng) + Z6HeCorrection - (cs::sqrXzero) * sin(sqrAng);
	TVector3 rightDetPosition(X6Helab, Y6Helab, Z6Helab);

	TFile cutgFile("/home/zalewski/aku/gasTar/gcuts5.root","READ");
	GCutdehe6 = (TCutG*)cutgFile.Get("dehe6");
	GCutdAngAng = (TCutG*)cutgFile.Get("dAngAng");
	GCutdAngE = (TCutG*)cutgFile.Get("dAngE");
	GCutpAngAng = (TCutG*)cutgFile.Get("pAngAng");
	GCutpAngE = (TCutG*)cutgFile.Get("pAngE");


	TFile MCcutgFile("/home/zalewski/aku/analysis/mcCuts5.root","READ");
	GCutmcPP = (TCutG*)MCcutgFile.Get("mcPP");
	GCutmcDD = (TCutG*)MCcutgFile.Get("mcDD");
	GCutmcHe6 = (TCutG*)MCcutgFile.Get("mcHe6");

	inDF.Filter(filterLeftDetector, {"cal_SQX_L"})
		.Filter(filterLeftDetector, {"cal_SQY_L"})
		.Filter(filterRightDetector, {"cal_SQX_R"})
		.Filter(filterRightDetector, {"cal_SQY_R"});
	auto outDF = inDF.Define("kinE", getKineticEnergy, {"tof"})
					 .Define("MWPC_1_Z",[](){return -816.0;})
					 .Define("MWPC_2_Z",[](){return -270.0;})
					 .Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_1_Z", "MWPC_2_X", "MWPC_2_Y", "MWPC_2_Z"})
					 .Define("vBeam", getBeamVector, {"MWPC", "kinE"})
					 .Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
					 .Define("tarVertex", getTarVertex, {"MWPC"})
					 .Define("SQX_L_strip", getStripNumber, {"cal_SQX_L"})
					 .Define("SQY_L_strip", getStripNumber, {"cal_SQY_L"})
					 .Define("SQX_R_strip", getStripNumber, {"cal_SQX_R"})
					 .Define("SQY_R_strip", getStripNumber, {"cal_SQY_R"})
					 .Define("CsI_L_strip", getStripNumber, {"cal_CsI_L"})
					 .Define("CsI_R_strip", getStripNumber, {"cal_CsI_R"})
					 .Define("sqlde", [](ROOT::VecOps::RVec<Double_t> &SQX_L, Int_t stripNumber){return SQX_L[stripNumber];}, {"cal_SQX_L", "SQX_L_strip"})
					 .Define("sqrde", [](ROOT::VecOps::RVec<Double_t> &SQX_R, Int_t stripNumber){return SQX_R[stripNumber];}, {"cal_SQX_R", "SQX_R_strip"})
					 .Define("sqletot", [](ROOT::VecOps::RVec<Double_t> &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"cal_CsI_L", "CsI_L_strip"})
					 .Define("sqretot", [](ROOT::VecOps::RVec<Double_t> &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"cal_CsI_R", "CsI_R_strip"})
					 .Define("sqltime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_L, Int_t stripNumber){return tSQX_L[stripNumber];}, {"tSQX_L", "SQX_L_strip"})
					 .Define("sqrtime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_R, Int_t stripNumber){return tSQX_R[stripNumber];}, {"tSQX_R", "SQX_R_strip"})

					 .Define("leftDetVertex", getLeftDetVertex_5, {"SQX_L_strip", "SQY_L_strip"})
					 .Define("leftLabVertex", [leftDetPosition](TVector3 leftDetVertex){return leftDetVertex+leftDetPosition;}, {"leftDetVertex"})
					 .Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
					 .Define("rightLabVertex", [rightDetPosition](TVector3 rightDetVertex){return rightDetVertex+rightDetPosition;}, {"rightDetVertex"})
					 .Define("v2H", [](TVector3 leftLabVertex, TVector3 tarVertex){return TVector3(leftLabVertex-tarVertex);}, {"leftLabVertex", "tarVertex"})
					 .Define("v6He", [](TVector3 rightLabVertex, TVector3 tarVertex){return TVector3(rightLabVertex-tarVertex);}, {"rightLabVertex", "tarVertex"})
					 .Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
					 .Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
					 .Define("mc1H", getMCAngle1H, {"sqlang", "lvBeam"})
					 .Define("mc2H", getMCAngle2H, {"sqlang", "lvBeam"})
					
					 .Define("he6", [](Double_t sqrde, Double_t sqretot){return GCutdehe6->IsInside(sqretot, sqrde);}, {"sqrde", "sqretot"})
					 .Define("pAngAng", [](Double_t sqrang, Double_t sqlang){return GCutpAngAng->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
					 .Define("dAngAng", [](Double_t sqrang, Double_t sqlang){return GCutdAngAng->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
					 .Define("pAngE", [](Double_t sqlang, Double_t sqlde, Double_t sqletot){return GCutpAngE->IsInside(sqlang,sqlde+sqletot);}, {"sqlang","sqlde", "sqletot"})
					 .Define("dAngE", [](Double_t sqlang, Double_t sqlde, Double_t sqletot){return GCutdAngE->IsInside(sqlang,sqlde+sqletot);}, {"sqlang","sqlde", "sqletot"})
					 .Define("pp", "he6 && pAngAng && pAngE && fNo>9")
					 .Define("dd", "he6 && dAngAng && dAngE && fNo<10")
					 .Define("tarCut", "(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+(tarVertex.Y())*(tarVertex.Y())<9")

					 .Define("thetaCM", "rnd->Uniform(0.0,TMath::Pi())")
					 .Define("he6Ene", getDDEne, {"sqlang","lvBeam"})
					 .Define("sqlangpp", getPPLang, {"thetaCM","lvBeam"})
					 .Define("sqrangpp", getPPRang, {"thetaCM","lvBeam"})
					 .Define("sqlangdd", getDDLang, {"thetaCM","lvBeam"})
					 .Define("sqrangdd", getDDRang, {"thetaCM","lvBeam"})
					 .Define("mcPP", [](Double_t sqrang, Double_t sqlang){return GCutmcPP->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
					 .Define("mcDD", [](Double_t sqrang, Double_t sqlang){return GCutmcDD->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
					 .Define("mcHe6", [](Double_t sqretot, Double_t sqrde){return GCutmcHe6->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"});

	outDF.Snapshot("analyzed", outFilename.Data());
}

void ui5()
{
	ROOT::EnableImplicitMT();
	TStopwatch *stopwatch = new TStopwatch();
	rnd = new TRandom3();
	selector = TSelector::GetSelector("/home/zalewski/aku/analysis/simplifier5.C");
	std::cout<<"..."<<std::endl<<".."<<std::endl<<"."<<std::endl;
	ToFconstant = cs::tof_const_5;
	tdcBinning = 0.0625;

	loadCalibrationParameters();
	//loadCorrectionParameters();
	TString str_name, sourceDir;
	sourceDir = "/home/zalewski/dataTmp/raw/geo5/";
	TSystemDirectory *dir_data = new TSystemDirectory("data",sourceDir.Data());
	
	TIter bluster = dir_data->GetListOfFiles();
	while (TObject *obj = bluster())
	{
		str_name = obj->GetName();
		std::cout<<str_name<<std::endl;
		if (str_name.Contains("he6") ||	(str_name.Contains("mc") && !str_name.Contains("out")))
		{
			TString inputFilePath = sourceDir + str_name;
			//printf("fName:\t%s\n",inputFilePath.Data());
			//translator(inputFilePath);
			//cleaner(inputFilePath);
			//calibratorCaller(inputFilePath);
			//analysis5(inputFilePath);
		}
	}

	analysis5("/home/zalewski/dataTmp/MC/geo5/mc5_1H.root");
	//makeSmallMC();
	//joindE();
	stopwatch->Print();
}