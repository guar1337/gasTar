#include "/home/zalewski/aku/analysis/constants.h"
#include "TVector3.h"
#include "TLorentzVector.h"

std::vector<std::string> columnList{"MWPC_1_X",
									"MWPC_1_Y",
									"MWPC_2_X",
									"MWPC_2_Y",
									"lvBeam"};

Double_t ToFconstant = cs::tof_const_5;
Double_t tdcBinning = 0.0625;
TRandom3 *rnd;
Double_t m_MWPC_1_Z = -816.0;
Double_t m_MWPC_2_Z = -270.0;

bool filterToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tdcF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tdcF5 = ((timeF5[0]+timeF5[1]+timeF5[2]+timeF5[3])/4.0);
	Double_t ToF = (tdcF5-tdcF3)*tdcBinning+ToFconstant;
	return (ToF>175 && ToF<185);
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

Double_t getToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tdcF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tdcF5 = ((timeF5[0]+timeF5[1]+timeF5[2]+timeF5[3])/4.0);
	Double_t ToF = (tdcF5-tdcF3)*tdcBinning+ToFconstant;
	return ToF;
}

Double_t getKinE(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	return cs::mass6He*(m_gamma-1.0);
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
				//printdcF("simplest solution %i\n", MWPCHitsArray[0]);
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
				//printdcF("simplest solution %i\n", MWPCHitsArray[0]);
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
					std::cout<<"Rogue MWPC event got through"<<std::endl;
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

TLorentzVector getLVbeam(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t kinE)
{
	TVector3 m_vBeam;
	TLorentzVector m_lvBeam;
	Double_t m_dX = m_MWPC_2_X - m_MWPC_1_X;
	Double_t	m_dY = m_MWPC_2_Y - m_MWPC_1_Y;
	Double_t	m_dZ = m_MWPC_2_Z - m_MWPC_1_Z;

	m_vBeam.SetXYZ(m_dX, m_dY, m_dZ);
	Double_t ene_beam = cs::mass6He + kinE;
	Double_t mom_beam = sqrt(ene_beam*ene_beam - cs::mass6He*cs::mass6He);
						
	m_vBeam.SetMag(mom_beam);
	m_lvBeam.SetVectM(m_vBeam, cs::mass6He);

	return m_lvBeam;
}

int beamCutter5()
{
	rnd = new TRandom3();
	ROOT::EnableImplicitMT();
	Int_t numberOfThreads = 20;
	ROOT::RDataFrame inDF("simplified", "/home/zalewski/dataTmp/simp/d2_1024.root");
	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
					.Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
					.Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
					.Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
					.Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
					.Filter(filterMWPC(0),{"x1","nx1"})
					.Filter(filterMWPC(1),{"y1","ny1"})
					.Filter(filterMWPC(2),{"x2","nx2"})
					.Filter(filterMWPC(3),{"y2","ny2"})
					.Filter(filterToF,{"tdcF3", "tdcF5"});

		 
	auto newDF = outDF.Define("tof", getToF, {"tdcF3", "tdcF5"})
							.Define("kinE", getKinE, {"tof"})
							.Define("aF5", "(F5[0]+F5[1]+F5[2]+F5[3])/4.0")
							.Define("MWPC_1_X", getMWPCpos(0),{"x1","nx1"})
							.Define("MWPC_1_Y", getMWPCpos(1),{"y1","ny1"})
							.Define("MWPC_2_X", getMWPCpos(2),{"x2","nx2"})
							.Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"})
					 		.Define("MWPC_1_Z",[](){return -816.0;})
					 		.Define("MWPC_2_Z",[](){return -270.0;})
					 		.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_1_Z", "MWPC_2_X", "MWPC_2_Y", "MWPC_2_Z"})
							.Define("tarVertex", getTarVertex, {"MWPC"})
							.Define("tarCut", "(tarVertex.X()+4.07)*(tarVertex.X()+4.07)+(tarVertex.Y())*(tarVertex.Y())<9")
							.Define("geo", "5")
							.Define("lvBeam", getLVbeam, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y", "kinE"});

	newDF.Snapshot("beamSource", "/home/zalewski/dataTmp/cln/d2_1024.root"/*, columnList*/);

	return 0;
}