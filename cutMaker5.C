	
int cutMaker5()
{
   TCutG *dehe6 = new TCutG("dehe6",16);
   dehe6->SetVarX("sqretot");
   dehe6->SetVarY("sqrde");
   dehe6->SetTitle("Graph");
   dehe6->SetFillStyle(1000);
   dehe6->SetPoint(0,148.655,17.7452);
   dehe6->SetPoint(1,103.328,22.0552);
   dehe6->SetPoint(2,72.7627,26.7367);
   dehe6->SetPoint(3,48.2758,31.5669);
   dehe6->SetPoint(4,21.1839,41.8217);
   dehe6->SetPoint(5,13.5426,47.2463);
   dehe6->SetPoint(6,5.20665,49.4756);
   dehe6->SetPoint(7,5.72764,47.2463);
   dehe6->SetPoint(8,22.0523,37.4374);
   dehe6->SetPoint(9,44.8025,28.966);
   dehe6->SetPoint(10,68.4211,24.0616);
   dehe6->SetPoint(11,99.3337,19.0085);
   dehe6->SetPoint(12,121.737,16.1104);
   dehe6->SetPoint(13,139.103,14.5499);
   dehe6->SetPoint(14,148.481,15.2187);
   dehe6->SetPoint(15,148.655,17.7452);

   TCutG *pAngAng = new TCutG("pAngAng",13);
   pAngAng->SetVarX("sqrang");
   pAngAng->SetVarY("sqlang");
   pAngAng->SetTitle("Graph");
   pAngAng->SetFillStyle(1000);
   pAngAng->SetPoint(0,6.43748,61.3588);
   pAngAng->SetPoint(1,5.49063,66.6136);
   pAngAng->SetPoint(2,3.70213,73.0361);
   pAngAng->SetPoint(3,2.8154,75.6051);
   pAngAng->SetPoint(4,2.33446,76.7728);
   pAngAng->SetPoint(5,1.92866,76.7728);
   pAngAng->SetPoint(6,1.91363,75.1964);
   pAngAng->SetPoint(7,2.12405,72.7442);
   pAngAng->SetPoint(8,2.45469,69.5329);
   pAngAng->SetPoint(9,4.63395,62.8185);
   pAngAng->SetPoint(10,5.40045,60.0159);
   pAngAng->SetPoint(11,6.01666,59.5488);
   pAngAng->SetPoint(12,6.43748,61.3588);

   TCutG *pAngE = new TCutG("pAngE",13);
   pAngE->SetVarX("sqlang");
   pAngE->SetVarY("sqlde+sqletot");
   pAngE->SetTitle("Graph");
   pAngE->SetFillStyle(1000);
   pAngE->SetPoint(0,60.5036,15.8771);
   pAngE->SetPoint(1,70.5217,8.60679);
   pAngE->SetPoint(2,76.2283,4.8147);
   pAngE->SetPoint(3,79.4303,3.19326);
   pAngE->SetPoint(4,79.3034,2.12102);
   pAngE->SetPoint(5,75.7527,1.65028);
   pAngE->SetPoint(6,71.663,3.0625);
   pAngE->SetPoint(7,68.4293,5.73004);
   pAngE->SetPoint(8,64.5299,8.92062);
   pAngE->SetPoint(9,61.0743,11.1697);
   pAngE->SetPoint(10,59.2038,12.5035);
   pAngE->SetPoint(11,57.3967,14.7264);
   pAngE->SetPoint(12,60.5036,15.8771);

   TCutG *dAngAng = new TCutG("dAngAng",17);
   dAngAng->SetVarX("sqrang");
   dAngAng->SetVarY("sqlang");
   dAngAng->SetTitle("Graph");
   dAngAng->SetFillStyle(1000);
   dAngAng->SetPoint(0,2.74025,83.4139);
   dAngAng->SetPoint(1,4.40851,81.4267);
   dAngAng->SetPoint(2,5.92648,77.5106);
   dAngAng->SetPoint(3,10.0746,67.9835);
   dAngAng->SetPoint(4,12.2388,63.8921);
   dAngAng->SetPoint(5,12.9151,62.0802);
   dAngAng->SetPoint(6,12.6146,59.2163);
   dAngAng->SetPoint(7,11.2769,60.5021);
   dAngAng->SetPoint(8,7.91036,66.9899);
   dAngAng->SetPoint(9,5.77619,71.6073);
   dAngAng->SetPoint(10,4.52875,73.4192);
   dAngAng->SetPoint(11,3.46166,75.4649);
   dAngAng->SetPoint(12,3.0709,76.517);
   dAngAng->SetPoint(13,2.43966,78.212);
   dAngAng->SetPoint(14,2.19919,80.1993);
   dAngAng->SetPoint(15,2.65007,82.3618);
   dAngAng->SetPoint(16,2.74025,83.4139);

   TCutG *dAngE = new TCutG("dAngE",14);
   dAngE->SetVarX("sqlang");
   dAngE->SetVarY("sqlde+sqletot");
   dAngE->SetTitle("Graph");
   dAngE->SetFillStyle(1000);
   dAngE->SetPoint(0,59.365,23.5919);
   dAngE->SetPoint(1,62.25,22.9295);
   dAngE->SetPoint(2,70.683,12.3938);
   dAngE->SetPoint(3,75.8505,6.3373);
   dAngE->SetPoint(4,79.5915,4.98091);
   dAngE->SetPoint(5,84.2835,2.55202);
   dAngE->SetPoint(6,84.6322,0.533203);
   dAngE->SetPoint(7,73.2192,0.470115);
   dAngE->SetPoint(8,73.346,4.82319);
   dAngE->SetPoint(9,69.6051,8.70311);
   dAngE->SetPoint(10,63.0426,15.6428);
   dAngE->SetPoint(11,59.9357,19.365);
   dAngE->SetPoint(12,58.5408,22.9295);
   dAngE->SetPoint(13,59.365,23.5919);


	TFile gcuts("/home/zalewski/aku/gasTar/gcuts5.root", "RECREATE");
	dehe6->Write();
	pAngAng->Write();
	dAngAng->Write();
	pAngE->Write();
	dAngE->Write();

   
   gcuts.Write();
   gcuts.Close();
   return 5;
}
