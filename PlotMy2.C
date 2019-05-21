#include <iostream>
#include <vector>
#include <map>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"


void PlotMy2(


//	  TString Variable = "mu_relIso[0]",
//	  TString Variable = "mu_eta[0]",
	  TString Variable = "puppi_pt",
//	  TString Variable = "mu_dxy[0]",
//	  TString Variable = "met_pt",

	  TString Suffix = "Data2017",
//	  int nBins  =   100,//dxy!!!!!!!!!!!
//	  float xmin =  -0.055,
//	  float xmax =  0.055,
//	  int nBins  =   80,//dz!!!!!!!!!!!
//	  float xmin =  -0.22,
//	  float xmax =  0.22,
	  int nBins  =   80,
	  float xmin =    0,
	  float xmax =  400,

//	  int nBins  =   50,//eta!!!!!!!!!!!!!!!!
//	  float xmin =    -2.5,
//	  float xmax =  2.5,

//	  int nBins  =   100,//telIso!!!!!!!!!!
//	  float xmin =    0,
//	  float xmax =  0.5,

//	  int nBins  =   6,//njets!!!!!
//	  float xmin =    0,
//	  float xmax =  6,

//	  int nBins  =   50,//phi!!!!!
//	  float xmin =    0,
//	  float xmax =  3.14,

//	  int nBins  =   51,//npv!!!!!
//	  float xmin =    -0.5,
//	  float xmax =  50.5,

//	  TString Weight = "gen_weight*LSF_weight*trig_weight*",
	  //TString Cuts = "&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",&& npv < 20 && npv > 15
// fabs(mu_dz[0])<0.01 && fabs(mu_dxy[0])<0.003 

//	  TString Cuts = "MT>60 && MT<120 && mu_relIso[0] < 0.5 && fabs(mu_charge[0])==1", && njets > 0.5 && MT>60
//	  TString Cuts = "MT>60 && MT<120 && fabs(mu_charge[0])==1",
	  TString xtitle = "met_pt",
  	  bool NormalMetForData = false,
	  TString channel = "Wmnu" /////SPISOK FAYLOV!!!!!!!!!!!!!!!!!1




          )
{
  //ModTDRStyle();



      cout << Variable<<endl;
	  TString DataFile1; 
	  TString DataFile2;

   	if (channel=="Wmnu") 
	{
           DataFile1 = "SingleMuon_Run2017F";//-31Mar2018";
	   DataFile2 = "SingleMuon_Run2017F";//-31Mar2018";	
	}

   	if (channel=="Wenu") 
	{
	   DataFile1 = "SingleElectron_RunBCDE";
	   DataFile2 = "SingleElectron__Run2017F-17Nov2017-v1";


	}
//cout <<"AAAAAAAAAAAAAAAAAAAA"<<  DataFile1<<DataFile2 <<endl;
	  //TString Weight = "gen_weight*PU_weightRunF*LSF_weight*trig_weight*"; //run F
	  	  TString Weight = "gen_weight*PU_weightRunBCDE*LSF_weight*trig_weight*"; //run F
if (Suffix == "RunF")  Weight = "gen_weight*PU_weightRunF*LSF_weight*trig_weight*";

//	  TString Weight = "gen_weight*LSF_weight*trig_weight*";
	TString MT = "MT";
	if (Variable.Contains("puppi_pt")) MT = "MTpuppi"; //////?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	if (Variable.Contains("met") && Variable.Contains("smeared")) MT = "MT_smeared"; //////?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	if (Variable.Contains("MTpuppi")) MT = "MTpuppi"; //////?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	if (Variable.Contains("MT") && Variable.Contains("smeared")) MT = "MT_smeared"; //////?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	//if (Variable.Contains("puppi_pt")) MT = "MT";

	  //TString Cuts = "fabs(mu_charge[0])==1 && mu_relIso[0] < 0.5 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && " + MT + ">50 && nbtag<0.5  ";
	TString Cuts;
	  if (channel=="Wmnu")  Cuts = "fabs(mu_charge[0])==1  && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && nbtag<0.5  ";
	  if (channel=="Wenu")  Cuts = "fabs(el_charge[0])==1  && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && nbtag<0.5  ";
//cout <<"AAAAAAAAAAAAAAAAAAAA  " <<Cuts<<endl;
	  TString ytitle = "Events";
	  //TString suffix = "_withPUOfficial";
	  TString suffix = "_Wmnu";
      	if (channel=="Wenu") suffix = "_WenuNewID";
     //	      	if (channel=="Wenu") suffix = "_WenuNewIDnoJet";

	//TString directory = "/nfs/dust/cms/user/dydukhle/METstudy/CMSSW_9_4_0_patch1/src/DesyTauAnalyses/NTupleMaker/test/Wmnu/backupMoriondOfficialPU/";
	TString directory = "/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/Wmnu/";

    if (channel=="Wenu"){
        directory = "/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/Wenu/";
    }

   ////// check for existance
  cout << "start check "<< channel << endl;
  TString FileToWrite;
   if (channel=="Wmnu") FileToWrite = "/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/RootWithWeightsMay"+Suffix+"/PlotsWmnuWithWeight_"+Variable+".root";
   if (channel=="Wenu") FileToWrite = "/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/WenuPlots/RootWithWeights"+Suffix+"/PlotsWenuWithWeight_"+Variable+".root";
  TFile * fileW = new TFile(FileToWrite);
  TH1F * check = NULL; 
  check = (TH1F * ) fileW->Get("data_obs_"+Variable+suffix);
  cout << "check = " << check << endl;
  //if (check != NULL) exit(1);
  bool logY = false;
  bool blindData = false;
  float QCDscale=1.00;
  bool NewMCmethod = false;
  //Run F
  //weightIsoEW = 0.93139 weightIsoQCD = 1.19217
  //weightIsoEW = 0.930434 weightIsoQCD = 0.70528
  float  EWweight = 0.930434;
  float QCDweight = 0.70528;

if (Suffix == "RunF") {
	 EWweight = 0.93139;
	QCDweight = 1.19217;
}
// RunBCDE


// RunBCDE
  //float QCDweight = 0.710319;
  //float EWweight = 0.944619;
  
// RunF
//float  EWweight = 0.983276; 
//float QCDweight = 0.885075;
  //float QCDweight = 0.823658;
  //float EWweight = 0.688938;

  //bool MCDataDriven = false;
  bool MCDataDriven = true;
  int nbMin = 4;
  int nbMax = 11;
  bool plotLeg = true;
  int position = 0; // 0 - right, 1 - left, 2 - central
  bool showSignal = true;
  int nSamples = 30;
   if (channel=="Wenu") nSamples = 26; 
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  //double lumi = 41290; // for test BCDE isoMu24
  // for RunF 
  //| HLT_IsoMu24_v10                                                                | 34    | 82   | 49069 | 11.531            | 10.979           |
  //| HLT_IsoMu24_v11                                                                | 2     | 13   | 4469  | 0.896
  
  //double lumi = 37836;

  double lumi = 37836-10979-896; //runBCDE 
    if (Suffix == "RunF")  lumi = 10979+896; //run F

  double TTnorm = 1.0;
  double Wnorm  = 1.0;

  //  TString topweight("");
  TString topweight2("");
  TString topweight("topptweight*");
  //  TString topweight2("topptweight*topptweight*");
  //  TString qcdweight("2.02*");
  TString qcdweight("qcdweight*");
  TString zptmassweight("zptmassweight*");
  //  TString qcdweight("(qcdweight*qcdweight/qcdweightup)*");
  //  TString qcdweight("qcdweight_nodzeta*");

    TString sampleNames[35] = {
    DataFile1, // data (0)
    DataFile2, // data (0)
        /*
	"DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" //1
	"DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	"DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	"DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	"DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" //1
	"DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	"DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	"DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
	*/
 	"DYJetsToLL_M-5to50_13TeV-12Apr2018",//DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", //1
	"W1JetsToLNu_TuneCP5_13TeV-12Apr2018",
	"W2JetsToLNu_TuneCP5_13TeV-12Apr2018",
	"W3JetsToLNu_TuneCP5_13TeV-12Apr2018",
	"W4JetsToLNu_TuneCP5_13TeV-12Apr2018",
	"WJetsToLNu_TuneCP5_13TeV-12Apr2018",
	"ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
	"ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
	"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
	"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", //10
	"WW_TuneCP5_13TeV-pythia8",
	"WZ_TuneCP5_13TeV-pythia8",
	"ZZ_TuneCP5_13TeV-pythia8",
//	"GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8",
	"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
	"TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
	"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",	
	
	//"QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8", //18
	//"QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia",
//	"QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8", //23
//	"QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
 //   "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8", 
 //   "QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8",//27
	
	
//	"QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8", //23
//	"QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_MINIAODSIM",
//	"QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8",
//	"QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8",//27
    
  };


  double xsec[35] = {1, // data (0)
			1, // data (0)
		    5765,  // DY(50) (1)
		    1.221*9644.5, //W1 (19)
			1.221*3144.5, //W1 (20)
			1.221*954.8,  //W1 (21)
			1.221*485.6,  //W1 (22)
			61526.7,// WJets (3)
//		     TTnorm*831.76,  // TT  (4)
		    136.02, // ST_t-channel_top (5)
		    80.95,  // ST_t-channel_antitop (6)
		    35.6,           // ST_tW_antitop (7)
		    35.6,           // ST_tW_top_5f (8)

		    63.21, //WW
		    47.13, // WZ
		    16.523, //ZZ

	    	//16792*0.0514, // GJet
		    364.4,
	        380.1,
	        87.31, //TT
//             720648000.0*0.00042,  // QCD        (18)
//		   10.4305*0.15544, //QCD
//		   469797*0.05362,
//		   117989*0.07335,
//		   558528000*0.0053,
//		   7820.25*0.10196,
//		   139803000*0.01182,
//		   645.528*0.12242,
//		   19222500*0.02276,
//		   187.109*0.13412,
//		   32.3486*0.14552,
//		   2758420*0.03844,
//		   0   // dummy
  };     


/* 10.4305 0.15544
469797 0.05362
 1273190000 0.003
117989 0.07335
558528000 0.0053
7820.25 0.10196
139803000 0.01182
645.528 0.12242
19222500 0.02276

187.109 0.13412
32.3486 0.14552
2758420 0.03844*/
/*
   if (channel=="Wenu")
{
   TString sampleNames[30] = {
    DataFile1, // data (0)
    DataFile2, // data (0)
   "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",     // isZTT  (1)
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", // isZTT  (2)
    "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (3)
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",               // (4)
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",// (5)
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",// (6)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (7)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (8)
    "VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8",// (9)
    "WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (10)
    //"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (11)
    "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (12)
    "WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8",// (13)
    //"WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8",// (14)
    "ZZTo4L_13TeV_powheg_pythia8",// (15)
    "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (16)
    //"QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8",// (18)
    "W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (19)
    "W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (20)
    "W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (21)
    "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (22)
    "ZZTo2L2Nu_13TeV_powheg_pythia8",// (17)
    "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
    "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
   // "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    //"QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
   // "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8"
    ""    // (24);
  };


  double xsec[30] = {1, // data (0)
			1, // data (0)
		     5765,  // DY(50) (1)
		     18610,   // DY(10to50)     (2)
		     Wnorm*61526.7,// WJets (3)
		     TTnorm*831.76,  // TT  (4)
		     136.02, // ST_t-channel_top (5)
		     80.95,  // ST_t-channel_antitop (6)
		     35.6,           // ST_tW_antitop (7)
		     35.6,           // ST_tW_top_5f (8)
		     11.95,  // VV          (9)
		     49.997, // WWToLNuQQ   (10)
//		     5.595,  // WZTo2L2Q    (11)
		     10.71,  // WZTo1L1Nu2Q (12)
		     3.05,   // WZTo1L3Nu   (13)
//		     4.42965,   // WZTo3L1Nu   (14)
		     1.212,  // ZZTo4L      (15)
		     3.22,   // ZZTo2L2Q    (16)
//             720648000.0*0.00042,  // QCD        (18)
		1.221*9644.5, //W1 (19)
		1.221*3144.5, //W1 (20)
		1.221*954.8,  //W1 (21)
		1.221*485.6,  //W1 (22)
             405.271,  // WGToLNuG        (17)
		16792*0.0514, // GJet
		477000*0.132,
		//114000*0.165,
		//557600000*0.0096,
		//9000*0.15,
		136000000*0.073,
		19800000*0.146,
		2800000*0.125,
		     0   // dummy 
  };     
//}
*/




  TString cuts[35];
  TString cutsInvIso[35];
  TString cutsInvMT[35];
  TString cutsInvMTandIso[35];
  


  for (int i=0; i<35; ++i) {
    cuts[i] = Weight+"("+Cuts+")";
  }

// for Data
	  //if (channel=="Wmnu")  Cuts = "fabs(mu_charge[0])==1  && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && nbtag<0.5  ";

    cuts[0] = "(met_flag>0.5 && fabs(mu_charge[0])==1  && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";

    cuts[1] = "(met_flag>0.5 && fabs(mu_charge[0])==1  && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
    // for Wjets
    cuts[7] = Weight+"("+Cuts+" &&(npartons==0||npartons>4))";


   	if (channel=="Wenu") 
	{
  	for (int i=0; i<35; ++i) {
    	cuts[i] = Weight+"("+Cuts+")";
  	}

// for Data

  	cuts[0] = "(met_flag>0.5 && fabs(el_charge[0])==1  && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";

  	cuts[1] = "(met_flag>0.5 && fabs(el_charge[0])==1  && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
// for Wjets
    	cuts[7] = Weight+"("+Cuts+" &&(npartons==0||npartons>4))";
	}



//cout <<"AAAAAAAAAAAAAAAAAAAA"<<  cuts[0]<< "   "<<cuts[1]<<"  "<<cuts[4] <<endl;

  TH1D * hist[35];
  TH1D * histInvIso[35];
  TH1D * histInvMT[35];
  TH1D * histInvMTandIso[35];



  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;




 /////////////////////////////////////////////////////////////////////////////////////////////////////// // filling histograms


  for (int i=0; i<nSamples; ++i) {
    //    std::cout << sampleNames[i] << std::endl;

//if (i==10 || i==15) continue;
// 
   if (i>=18) continue; 
    cout<<"Sample Name"<<sampleNames[i]<<endl;
    cout<<"Index: "<<i <<endl;
    try{
    	TFile * file = new TFile(directory+sampleNames[i]+".root");
	}
    catch(std::runtime_error){
	continue;
	}
    TFile * file = new TFile(directory+sampleNames[i]+".root");
    if (i>1) file = new TFile(directory+sampleNames[i]+".root");
    TH1D * histWeightsH;
    TTree * tree; 
	cout<<directory+sampleNames[i]+".root"<<endl;
   	histWeightsH = (TH1D*)file->Get("Wmnu/histWeightsH");
    	tree = (TTree*)file->Get("Wmnu/T");
/*
   	if (channel=="Wenu") 
	{
   	histWeightsH = (TH1D*)file->Get("Wenu/histWeightsH");
    	tree = (TTree*)file->Get("Wenu/T");
	}
   	if (channel=="Wmnu") 
	{
   	histWeightsH = (TH1D*)file->Get("Wmnu/histWeightsH");
    	tree = (TTree*)file->Get("Wmnu/T");
	}
*/
/*

Float_t         pu_weight;
TBranch        *b_pu_weight;
tree->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
Float_t         el_pt[20];
TBranch        *b_el_pt;
tree->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
    TH1D * aaa = new TH1D ("aaa","",50,0,500);

Float_t sum=0;
for (int iB=0; iB<=tree->GetEntriesFast(); ++iB) {tree->GetEntry(iB); cout<<pu_weight<<endl; aaa->Fill(el_pt[0],pu_weight); 

//cout<<"pu_weight  "<<pu_weight<<endl;
//sum = sum + pu_weight;
//sum+=pu_weight;
//cout<<"sum  "<< sum<<endl;
if (el_pt[0]<26)cout<<"el_pt[0]!!!!!!!!!!!!!!!!!!!!!!!!!!  "<<el_pt[0]<<endl;

}
cout<<"aaa!!!!!!!!!!   "<<aaa->GetSumOfWeights() <<endl;
cout<<"aaa!!????????????   "<<sum <<endl;
*/

cout << "sampleNames[i]  "<<sampleNames[i] <<endl;
cout << "xsec[i]  "<<xsec[i]<< endl;
cout << "lumi  "<<lumi<< endl;
cout << "histWeightsH->GetSumOfWeights()  "<<histWeightsH->GetSumOfWeights()<< endl;
    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
cout << "norm  "<<norm<< endl;
	TString Variable2;
	if (NormalMetForData && i<1.5 && !Variable.Contains("mu") && !Variable.Contains("el")) 
	{
	Variable2 = Variable;
	Variable = "met_pt";
cout << "?????????????????????????????"<< endl;
	if (Variable2.Contains("puppi_pt")) Variable = "puppi_pt";
	if (Variable2.Contains("MT")) Variable = "MT";
	if (Variable2.Contains("MTpuppi")) Variable = "MTpuppi";
	if (Variable2.Contains("Ut")) Variable = "Ut";
	if (Variable2.Contains("Ut_puppi")) Variable = "Ut_puppi";
	if (Variable2.Contains("Utr")) Variable = "Utr";
	if (Variable2.Contains("Ucol")) Variable = "Ucol";
	}
cout << "Variable   "<<Variable <<endl;

    TString histName = sampleNames[i] + Variable + "_";
    TString histNameInvIso = sampleNames[i] + Variable + "_InvIso";
    TString histNameInvMT = sampleNames[i] + Variable + "_InvMT";
    TString histNameInvMTandIso = sampleNames[i] + Variable + "_InvMTandIso";
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    hist[i]->Sumw2();
    histInvIso[i] = new TH1D(histNameInvIso,"",nBins,xmin,xmax);
    histInvIso[i]->Sumw2();
    histInvMT[i] = new TH1D(histNameInvMT,"",nBins,xmin,xmax);
    histInvMT[i]->Sumw2();
    histInvMTandIso[i] = new TH1D(histNameInvMTandIso,"",nBins,xmin,xmax);
    histInvMTandIso[i]->Sumw2();


    tree->Draw(Variable+">>"+histName,cuts[i]);
    

if (!NewMCmethod)
{
if (MCDataDriven)    tree->Draw(Variable+">>"+histNameInvMT,cutsInvMT[i]);
if (MCDataDriven)    tree->Draw(Variable+">>"+histNameInvMTandIso,cutsInvMTandIso[i]);
    tree->Draw(Variable+">>"+histNameInvIso,cutsInvIso[i]);
}	
	cout << "hist[i]->GetSumOfWeights()  "<<hist[i]->GetSumOfWeights()<< endl;
if (!NewMCmethod)
{
	cout << "histInvIso[i]->GetSumOfWeights()  "<<histInvIso[i]->GetSumOfWeights()<< endl;
if (MCDataDriven)	cout << "histInvMT[i]->GetSumOfWeights()  "<<histInvMT[i]->GetSumOfWeights()<< endl;
if (MCDataDriven)	cout << "histInvMTandIso[i]->GetSumOfWeights()  "<<histInvMTandIso[i]->GetSumOfWeights()<< endl;
}
    if (i>1) {
      for (int iB=1; iB<=nBins; ++iB) {
double x,e;
if (!NewMCmethod)
{
	x= histInvIso[i]->GetBinContent(iB);
	e = histInvIso[i]->GetBinError(iB);
    	histInvIso[i]->SetBinContent(iB,norm*x);
    	histInvIso[i]->SetBinError(iB,norm*e);
if (MCDataDriven)
{
	 x = histInvMT[i]->GetBinContent(iB);
	 e = histInvMT[i]->GetBinError(iB);
    	histInvMT[i]->SetBinContent(iB,norm*x);
    	histInvMT[i]->SetBinError(iB,norm*e);

	 x = histInvMTandIso[i]->GetBinContent(iB);
	 e = histInvMTandIso[i]->GetBinError(iB);
    	histInvMTandIso[i]->SetBinContent(iB,norm*x);
    	histInvMTandIso[i]->SetBinError(iB,norm*e);
}
}
	 x = hist[i]->GetBinContent(iB);
	 e = hist[i]->GetBinError(iB);
    	hist[i]->SetBinContent(iB,norm*x);
    	hist[i]->SetBinError(iB,norm*e);
      }
    }
/*
	TCanvas* canv1; 
	canv1 = new TCanvas("c1"+sampleNames[i], "c1"+sampleNames[i]);
	canv1->cd();
	hist[i]->Draw();
	canv1->Print("/nfs/dust/cms/user/dydukhle/METstudy/CMSSW_9_4_0_patch1/src/DesyTauAnalyses/NTupleMaker/test/WenuPlots/TEST"+sampleNames[i]+".pdf");
*/


	if (NormalMetForData && i<1.5 && !Variable.Contains("mu") && !Variable.Contains("el")) 
	{
	Variable = Variable2;
	}


  }

  delete dummyCanv;
    

    cout << "!!end!!!!!!!!!!!!!!!!!!!!"<<endl;



cout<<"da"<<endl;
  //TH1D * histData;

  
  
  
  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs_"+Variable+suffix);cout<<sampleNames[0]<<endl;
  //  TH1D * histData;
  if (Suffix == "RunF") histData = (TH1D*)hist[1]->Clone("data_obs_"+Variable+suffix);cout<<sampleNames[1]<<endl;
//  histData->Add(histData,hist[1]);cout<<sampleNames[1]<<endl;
  TH1D * DY = (TH1D*)hist[2]->Clone("DY_"+Variable+suffix);cout<<sampleNames[2]<<endl;
  //DY->Add(DY,hist[3]);cout<<sampleNames[3]<<endl;
  TH1D * WJ = (TH1D*)hist[3]->Clone("WJ_"+Variable+suffix);cout<<sampleNames[3]<<endl;
  WJ->Add(WJ,hist[4]);cout<<sampleNames[4]<<endl;
  WJ->Add(WJ,hist[5]);cout<<sampleNames[5]<<endl;
  WJ->Add(WJ,hist[6]);cout<<sampleNames[6]<<endl;
  WJ->Add(WJ,hist[7]);cout<<sampleNames[7]<<endl;

  cout<<"WJ Finish"<<endl;
  //TH1D * TT  = (TH1D*)hist[5]->Clone("TT_"+Variable+suffix);cout<<sampleNames[5]<<endl;
  TH1D * VV  = (TH1D*)hist[12]->Clone("VV_"+Variable+suffix);cout<<sampleNames[12]<<endl;
  VV->Add(VV,hist[13]);cout<<sampleNames[13]<<endl;
  VV->Add(VV,hist[14]);cout<<sampleNames[14]<<endl;
//  VV->Add(VV,hist[15]);cout<<sampleNames[15]<<endl;

  cout<<"WW Finish"<<endl;
  TH1D * TT  = (TH1D*)hist[15]->Clone("TT_"+Variable+suffix);cout<<sampleNames[15]<<endl;
  cout<<"WW Finish"<<endl;
  cout<<sampleNames[16]<<endl;
  TT->Add(TT,hist[16]);cout<<sampleNames[16]<<endl;
  cout<<"WW Finish"<<endl;
  TT->Add(TT,hist[17]);cout<<sampleNames[17]<<endl;
  cout<<"TT Finish"<<endl;
  TH1D * ST  = (TH1D*)hist[8]->Clone("ST_"+Variable+suffix);cout<<sampleNames[8]<<endl;
  ST->Add(ST,hist[9]);cout<<sampleNames[9]<<endl;
  ST->Add(ST,hist[10]);cout<<sampleNames[10]<<endl;
  ST->Add(ST,hist[11]);cout<<sampleNames[11]<<endl;
   cout<<"ST Finish"<<endl;
  TH1D * QCD;
   cout<<"QCD Finish"<<endl; 
if (!NewMCmethod)
{
  histInvIso[0]->Add(histInvIso[1],1);
if (MCDataDriven)  histInvMT[0]->Add(histInvMT[1],1);
if (MCDataDriven)  histInvMTandIso[0]->Add(histInvMTandIso[1],1);

cout<<"QCD troubles"<<endl;
  for (int i=2; i<(nSamples); ++i) 
	
	{ if (i>=18)continue;
	cout<<"da"<<endl;
	histInvIso[0]->Add(histInvIso[i],-1);
        if (MCDataDriven)	histInvMT[0]->Add(histInvMT[i],-1);
        if (MCDataDriven)	histInvMTandIso[0]->Add(histInvMTandIso[i],-1);
	 
	}
cout<<"QCD troubles"<<endl;
QCD  = (TH1D*)histInvIso[0]->Clone("QCD_"+Variable+suffix);
if (MCDataDriven)QCD->Scale(histInvMT[0]->GetSumOfWeights()/histInvMTandIso[0]->GetSumOfWeights());

if (MCDataDriven)cout << "QCD scale = " << histInvMT[0]->GetSumOfWeights()/histInvMTandIso[0]->GetSumOfWeights() << endl;
if (!MCDataDriven)QCD->Scale(QCDscale); 
}

if (NewMCmethod)
{
QCD  = (TH1D*)hist[19]->Clone("QCD_"+Variable+suffix);
QCD->Add(QCD,hist[20]);
QCD->Add(QCD,hist[21]);
QCD->Add(QCD,hist[22]);
QCD->Add(QCD,hist[23]);
QCD->Add(QCD,hist[24]);
QCD->Add(QCD,hist[25]);
QCD->Add(QCD,hist[26]);
QCD->Add(QCD,hist[27]);
QCD->Add(QCD,hist[28]);
QCD->Add(QCD,hist[29]);
	
QCD->Scale(QCDweight);
ST->Scale(EWweight);
VV->Scale(EWweight);
WJ->Scale(EWweight);
TT->Scale(EWweight);
DY->Scale(EWweight);
}	

  std::cout << "ST   : " << ST->GetSumOfWeights() << " : " << ST->Integral(1,nBins+1) <<" GetSumOfWeights= " << ST->GetSumOfWeights()<< std::endl;
  std::cout << "VV   : " << VV->GetSumOfWeights() << " : " << VV->Integral(1,nBins+1) <<" GetSumOfWeights= " << VV->GetSumOfWeights() <<std::endl;
  std::cout << "QCD   : " << QCD->GetSumOfWeights() << " : " << QCD->Integral(1,nBins+1) <<" GetSumOfWeights= " << QCD->GetSumOfWeights()<< std::endl;
  std::cout << "WJ   : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(1,nBins+1) <<" GetSumOfWeights= " << WJ->GetSumOfWeights()<< std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : " << TT->Integral(1,nBins+1) <<" GetSumOfWeights= " << TT->GetSumOfWeights()<< std::endl;
  std::cout << "DY : " << DY->GetSumOfWeights() << " : " << DY->Integral(1,nBins+1) <<" GetSumOfWeights= " << DY->GetSumOfWeights()<< std::endl;


  QCD->Add(QCD,VV);
  ST->Add(ST,QCD);
  DY->Add(DY,ST);
  TT->Add(TT,DY);
  WJ->Add(WJ,DY);
  std::cout << "BKG : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/WJ->GetSumOfWeights() << "+/-" 
	    << TMath::Sqrt(histData->GetSumOfWeights())/WJ->GetSumOfWeights() << std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////
    //ModTDRStyle();
    TH1D * bkgdErr = (TH1D*)WJ->Clone("bkgdErr");
    ModTDRStyle();

  
    TCanvas* canv1 = new TCanvas("c1", "c1");
    canv1->cd();
    std::vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
    pads[0]->SetLogy(logY);
    
    std::vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
    h[0]->Draw();
    
    std::string units="";
    std::string xtitle_ = (std::string) xtitle;
    size_t pos = xtitle_.find("[");
    if(pos!=std::string::npos) {
        units = xtitle_.substr(pos+1, xtitle_.find("]") - pos -1 );
        xtitle_ = xtitle_.substr(0, pos);
    }
    
    pads[1]->cd();
    h[1]->Draw();
    SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.4, 1.6);
    StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
    h[1]->GetYaxis()->SetNdivisions(4);
    h[1]->GetXaxis()->SetTitleOffset(1.2);
    h[1]->GetYaxis()->SetTitleOffset(2.0);
    pads[0]->cd();
    h[0]->GetYaxis()->SetTitleOffset(2.0);
    pads[1]->SetGrid(0,1);
    //it complains if the minimum is set to 0 and you try to set log y scale
    h[0]->SetMinimum(1);
    pads[0]->cd();
    
    // Setup legend
    TLegend *legend = PositionedLegend(0.40, 0.30, 3, 0.03);
    legend->SetTextFont(42);
    
    histData->SetMarkerColor(1);
    histData->SetLineColor(1);
    histData->SetFillColor(1);
    histData->SetFillStyle(0);
    histData->SetLineWidth(2);
    histData->SetMarkerStyle(20);
    histData->SetMarkerSize(1.1);
    
    legend->AddEntry(histData, "Observed", "ple");
    
    InitHist(QCD,TColor::GetColor("#FFCCFF"));
    InitHist(DY,TColor::GetColor("#DE5A6A"));
    InitHist(TT,TColor::GetColor("#9999CC"));
    InitHist(WJ,TColor::GetColor("#6F2D35"));
    InitHist(VV,TColor::GetColor("#4496C8"));
    InitHist(ST,TColor::GetColor("#FFCC66"));
    
    legend->AddEntry(WJ,"WJETS","f");
    legend->AddEntry(TT,"t#bar{t}","f");
    legend->AddEntry(DY,"DY","f");
    //  leg->AddEntry(VV,"VV+VVV","f");
    legend->AddEntry(ST,"ST","f");
    legend->AddEntry(QCD,"QCD","f");
    legend->AddEntry(VV,"VV","f");

  /*  
   WJ->Draw("sameh");
    QCD->Draw("sameh");
    //ZLL->Draw("sameh");
    DY->Draw("sameh");
    VV->Draw("sameh");
    TT->Draw("sameh");
    ST->Draw("sameh");
*/
    WJ->Draw("sameh");
    TT->Draw("sameh");
    DY->Draw("sameh");
    ST->Draw("sameh"); 
    QCD->Draw("sameh");   
    VV->Draw("sameh");

    canv1->Update();
    
   // InitSignal(vbfH,2);
   // InitSignal(SMH,2);
   /* 
    if (showSignal)
    {
        legend->AddEntry(SMH,"SM Higgs(125) #times 10","f");
    }
    
    if (showSignal)
    {
        SMH->Draw("hsame");
    }
    */
    canv1->Update();
    
    if (blindData)
    {
        for (int iB=nbMin; iB<=nbMax; ++iB)
        {
            histData->SetBinContent(iB,-1);
            histData->SetBinError(iB,0);
        }
    }
    
    
    //float errLumi = 0.03;
    //float errMuon = 0.03;
    //float errElectron = 0.04;
    for (int iB=1; iB<=nBins; ++iB) {
        //QCD->SetBinError(iB,0);
        //VV->SetBinError(iB,0);
        //TT->SetBinError(iB,0);
        WJ->SetBinError(iB,0);
        //DY->SetBinError(iB,0);
        //ZTT->SetBinError(iB,0);
        float eStat =  bkgdErr->GetBinError(iB);
        float X = bkgdErr->GetBinContent(iB);
        //float eLumi = errLumi * X;
        //float eMuon = errMuon * X;
        //float eElectron = errElectron * X;
        //float eBkg = dummy->GetBinError(iB);
        //float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg+eMuon*eMuon+eElectron*eElectron);
	float Err = TMath::Sqrt(eStat*eStat);
        bkgdErr->SetBinError(iB,Err);
    }

    
    bkgdErr->SetMarkerSize(0);
    int new_idx = CreateTransparentColor(13,0.4);
    bkgdErr->SetFillColor(new_idx);
    bkgdErr->SetLineWidth(1);
    bkgdErr->Draw("e2same");
    legend->AddEntry(bkgdErr, "Bkg. uncertainty" , "F" );
    canv1->Update();

    TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
    TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    
    for (int iB=1; iB<=nBins; ++iB) {
        float x1 = histData->GetBinContent(iB);
        float x2 = WJ->GetBinContent(iB);
        ratioErrH->SetBinContent(iB,1.0);
        ratioErrH->SetBinError(iB,0.0);
        float xBkg = bkgdErr->GetBinContent(iB);
        float errBkg = bkgdErr->GetBinError(iB);
        if (xBkg>0) {
            float relErr = errBkg/xBkg;
            ratioErrH->SetBinError(iB,relErr);
            
        }
        if (x1>0&&x2>0) {
            float e1 = histData->GetBinError(iB);
            float ratio = x1/x2;
            float eratio = e1/x2;
            ratioH->SetBinContent(iB,ratio);
            ratioH->SetBinError(iB,eratio);
        }
        else {
            ratioH->SetBinContent(iB,1000);
        }
    }
    
    pads[1]->cd();
    ratioErrH->GetYaxis()->SetRangeUser(0.4,1.6);
    ratioErrH->Draw("e2same");
    ratioH->SetMarkerStyle(2);
    ratioH->Draw("pe0same");

    pads[0]->cd();
    histData->Draw("pesame");
    
    FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
    DrawTitle(pads[0], "37.8 fb^{-1} (13 TeV)", 3);
   	  if (channel=="Wmnu") DrawTitle(pads[0], "W to #mu#nu", 1);
    	  if (channel=="Wenu") DrawTitle(pads[0], "W to e#nu", 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();
       	  if (channel=="Wmnu") canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/"+Variable+Suffix+suffix+"May.pdf");
       	  if (channel=="Wenu") canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_5/src/DesyTauAnalyses/NTupleMaker/test/WenuPlots/"+Variable+Suffix+suffix+".pdf");

	TFile *Target = TFile::Open (FileToWrite, "update");    
	
    	WJ->Write();
    	TT->Write();
    	DY->Write();
    	ST->Write(); 
    	QCD->Write();   
    	VV->Write();
	histData->Write();

	Target->Write();
	delete Target;

}
