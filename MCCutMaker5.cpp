int MCCutMaker5()
{
   TFile mcCuts("/home/zalewski/aku/analysis/mcCuts5.root", "RECREATE");

   TCutG *mcPP = new TCutG("mcPP",13);
   mcPP->SetVarX("sqrang");
   mcPP->SetVarY("sqlang");
   mcPP->SetTitle("Graph");
   mcPP->SetFillStyle(1000);
   mcPP->SetPoint(0,9.44293,22.4098);
   mcPP->SetPoint(1,11.2545,31.7516);
   mcPP->SetPoint(2,10.4242,51.1783);
   mcPP->SetPoint(3,7.72569,65.7219);
   mcPP->SetPoint(4,3.30993,81.2208);
   mcPP->SetPoint(5,0.422705,88.4395);
   mcPP->SetPoint(6,0.517059,84.724);
   mcPP->SetPoint(7,4.29121,69.6497);
   mcPP->SetPoint(8,7.53699,54.2569);
   mcPP->SetPoint(9,8.51827,39.1826);
   mcPP->SetPoint(10,7.76344,26.9745);
   mcPP->SetPoint(11,7.49925,22.4098);
   mcPP->SetPoint(12,9.44293,22.4098);

   TCutG *mcDD = new TCutG("mcDD",20);
   mcDD->SetVarX("sqrang");
   mcDD->SetVarY("sqlang");
   mcDD->SetTitle("Graph");
   mcDD->SetFillStyle(1000);
   mcDD->SetPoint(0,15.2962,58.8004);
   mcDD->SetPoint(1,12.6664,65.2229);
   mcDD->SetPoint(2,9.38889,72.9724);
   mcDD->SetPoint(3,6.97494,78.9172);
   mcDD->SetPoint(4,6.32729,80.828);
   mcDD->SetPoint(5,5.24789,83.1635);
   mcDD->SetPoint(6,3.65821,83.482);
   mcDD->SetPoint(7,3.20682,82.155);
   mcDD->SetPoint(8,3.30495,80.6688);
   mcDD->SetPoint(9,3.77597,78.121);
   mcDD->SetPoint(10,4.65912,75.9448);
   mcDD->SetPoint(11,6.09179,73.6624);
   mcDD->SetPoint(12,7.46558,70.3185);
   mcDD->SetPoint(13,9.42814,67.0276);
   mcDD->SetPoint(14,11.0178,63.5775);
   mcDD->SetPoint(15,12.5682,60.4459);
   mcDD->SetPoint(16,13.8635,57.5265);
   mcDD->SetPoint(17,14.256,56.6242);
   mcDD->SetPoint(18,15.0214,57.1019);
   mcDD->SetPoint(19,15.2962,58.8004);

   TCutG *mcHe6 = new TCutG("mcHe6",20);
   mcHe6->SetVarX("sqretot");
   mcHe6->SetVarY("sqrde");
   mcHe6->SetTitle("Graph");
   mcHe6->SetFillStyle(1000);
   mcHe6->SetPoint(0,173.45,14.8639);
   mcHe6->SetPoint(1,152.511,16.6194);
   mcHe6->SetPoint(2,94.4263,22.4713);
   mcHe6->SetPoint(3,78.3635,24.8957);
   mcHe6->SetPoint(4,67.7506,26.9857);
   mcHe6->SetPoint(5,49.3931,31.4164);
   mcHe6->SetPoint(6,30.1751,38.8567);
   mcHe6->SetPoint(7,17.1241,46.297);
   mcHe6->SetPoint(8,3.78623,55.9944);
   mcHe6->SetPoint(9,1.77838,54.824);
   mcHe6->SetPoint(10,9.37953,48.1361);
   mcHe6->SetPoint(11,24.8687,38.7731);
   mcHe6->SetPoint(12,37.4894,33.0048);
   mcHe6->SetPoint(13,52.9786,27.6545);
   mcHe6->SetPoint(14,75.3518,22.6385);
   mcHe6->SetPoint(15,97.4381,18.2914);
   mcHe6->SetPoint(16,120.098,15.6998);
   mcHe6->SetPoint(17,150.79,13.4427);
   mcHe6->SetPoint(18,172.302,12.5231);
   mcHe6->SetPoint(19,173.45,14.8639);


   mcPP->Write();
   mcDD->Write();
   mcHe6->Write();

   mcCuts.Write();
   mcCuts.Close();
   return 5;
}