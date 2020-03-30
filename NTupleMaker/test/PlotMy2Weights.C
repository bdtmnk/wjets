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



void PlotMy2Weights(TString Variable = "mu_relIso[0]",
                    TString Suffix = "_Moriond17",
                    int nBins  =   100,
                    float xmin =    0,
                    float xmax =  0.5,
                    TString xtitle = "iso",
                    bool NormalMetForData = false,
                    TString channel = "Wmnu",
                    TString Run="ABCD")
    {


	  TString DataFile1; 
	  TString DataFile2;
	  TString DataFile3;
	  TString DataFile4;
	  TString DataFile5;
	  TString DataFile6;
	  TString DataFile7;

   	if (channel=="Wmnu") {

       DataFile1 = "SingleMuon_Run2018B";
       DataFile2 = "SingleMuon_Run2018A";
       DataFile3 = "SingleMuon_Run2018C";
       DataFile4 = "SingleMuon_Run2018D-22Jan2019-v2_0_";
       DataFile5 = "SingleMuon_Run2018D-22Jan2019-v2_1_";
       DataFile6 = "SingleMuon_Run2018D-22Jan2019-v2_2_";
       DataFile7 = "SingleMuon_Run2018D_22Jan2019-v2_-1";
	}

  TString directory = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/Wmnu_new/";
	TString Weight = "gen_weight*PUweight*LSF_weight*trig_weight*";
  TString MT = "MT";
  TString Cuts;

	if (Variable.Contains("puppi_pt")) MT = "MTpuppi";
	if (Variable.Contains("met") && Variable.Contains("smeared")) MT = "MT_smeared";
	if (Variable.Contains("MTpuppi")) MT = "MTpuppi";
	if (Variable.Contains("MT") && Variable.Contains("smeared")) MT = "MT_smeared";

	if (channel=="Wmnu")  Cuts = "mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5";

  TString ytitle = "Events";
  TString suffix = "_Wmnu";
    
  if (channel=="Wenu") suffix = "_WenuNewID";

  TString FileToWrite;
  if (channel=="Wmnu") FileToWrite = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/RootWithWeights/"+Run+"/"+Variable+".root";
  TFile * fileW = new TFile(FileToWrite);
  TH1F * check = NULL;
  check = (TH1F * ) fileW->Get("data_obs_"+Variable+suffix);

  bool logY = false;
  TString isLog;
  if (logY) isLog ="Log";
  if (logY==false) isLog ="NonLog";
  bool blindData = false;
  float QCDscale = 1.00;
  int nbMin = 4;
  int nbMax = 11;
  bool plotLeg = true;
  int position = 0;
  bool showSignal = true;
  int nSamples = 34;
  if (channel=="Wenu") nSamples = 29;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  //Lumi Settings:
  double lumi;
  if (Run == "A") {lumi = 13480;}
  else if (Run=="B"){lumi= 6785;}
  else if (Run=="C"){lumi= 6612;}
  else if (Run=="D"){lumi=58822.126 - (13480+6785+6612);}
  else if(Run=="ABC"){lumi = 13480+6785+6612;}
  else{lumi = 58822.126;}

  double TTnorm = 1.0;
  double Wnorm  = 1.0;

  TString topweight2("");
  TString topweight("topptweight*");
  TString zptmassweight("zptmassweight*");
  TString qcdweight("(qcdweight*qcdweight/qcdweightup)*");

  TString sampleNames[35] = {

         DataFile1, // data (0)
        "DYJetsToLL_M-50_new",//1
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_",//2
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_",//3
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_",//4
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_",//5
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_",//6
        "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",//7
        "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",//8
        "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",//9
        "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", //10
        "WW_TuneCP5_13TeV-pythia8",//11
        "WZ_TuneCP5_13TeV-pythia8",//12
        "ZZ_TuneCP5_13TeV-pythia8",//13
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",//14
        "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",//15
        "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",//16
        "QCD_Pt-1000toInf",//17
        "QCD_Pt-120to170",//18
        "QCD_Pt-170to300",//19
        "QCD_Pt-20to30",//20
        "QCD_Pt-300to470",//21
        "QCD_Pt-30to50",//22
        "QCD_Pt-470to600",//23
        "QCD_Pt-50to80",//24
        "QCD_Pt-600to800",//25
        "QCD_Pt-800to1000",//26
        "QCD_Pt-80to120",//27
        DataFile2,//28
        DataFile3,//29
        DataFile4,//30
        DataFile5,//31
        DataFile6,//32
        DataFile7//33

      };

    double xsec[35] = {

      1, // Single Muon data (1)
      2075.14*3 , ////5765,  // DY(50) (2)

      1.221*9644.5, //W1 (3)
      1.221*3144.5, //W1 (4)
      1.221*954.8,  //W1 (5)
      1.221*485.6,  //W1 (6)
      61526.7,// WJets (7)

      44.33,   // ST_t-channel_top (8)
      26.38,   // ST_t-channel_antitop (9)
      35.85,   // ST_tW_antitop (10)
      35.85,   // ST_tW_top_5f (11)

      63.21,   //WW (12)
      47.13,   // WZ (13)
      16.523,  //ZZ (14)

      364.4, // TT (15)
      380.1, // TT (16)
      87.31, // TT (17)

      10.4305*0.15544, //QCD (18)
      469797*0.05362,  //QCD (19)
      117989*0.07335,  //QCD (20)
      558528000*0.0053, //QCD (21)
      7820.25*0.10196,  //QCD (22)
      139803000*0.01182, //QCD (23)
      645.528*0.12242,   //QCD (24)
      19222500*0.02276,  //QCD (25)
      187.109*0.13412,   //QCD (26)
      32.3486*0.14552,   //QCD (27)
      2758420*0.03844,   //QCD (28)

      1,//29
      1,//30
      1,//31
      1,//32
      1,//33
      1,//34
      1//35
    };

  TString cuts[35];
  TString cutsInvIso[35];
  TString cutsInvMT[35];
  TString cutsInvMTandIso[35];

  for (int i=0; i<35; ++i) {
    
    cuts[i] = Weight+"("+Cuts+")";
    cutsInvIso[i] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
    cutsInvMT[i] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";
    cutsInvMTandIso[i] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";

  }

  // For Data
  for (int j=0; j<35;j++){

  	if (j==0||j>=28){
  		cuts[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
  		cutsInvIso[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
  		cutsInvMT[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";
  		cutsInvMTandIso[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";

  	}
  }
  // For Wjets
  cuts[6] = Weight+"("+Cuts+" &&(npartons==0||npartons>4))";
  cutsInvIso[6] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    &&(npartons==0||npartons>4) )";
  cutsInvMT[6] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    &&(npartons==0||npartons>4) )";
  cutsInvMTandIso[6] = Weight+"(mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] > 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    &&(npartons==0||npartons>4))";

	  if (channel=="Wenu"){

		for (int i=0; i<30; ++i) {
            cuts[i] = Weight+"("+Cuts+")";
            cutsInvIso[i] = Weight+"(fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
            cutsInvMT[i] = Weight+"(fabs(el_charge[0])==1 && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";
            cutsInvMTandIso[i] = Weight+"(fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";
         }

        // For Data
  	for (int j=0; j<35;j++){

    	if (j==0||j>=28){ 
    		cuts[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(el_charge[0])==1 && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
    		cutsInvIso[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    )";
    		cutsInvMT[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(el_charge[0])==1 && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";
    		cutsInvMTandIso[j] = "(mu_pt[0]>29 && met_flag>0.5 && fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    )";

      }
    }
      // for Wjets  	
      cuts[6] = Weight+"("+Cuts+" &&(npartons==0||npartons>4))";
    	cutsInvIso[6] = Weight+"(mu_pt[0]>29 && fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    &&(npartons==0||npartons>4) )";
    	cutsInvMT[6] = Weight+"(mu_pt[0]>29 && fabs(el_charge[0])==1 && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    &&(npartons==0||npartons>4) )";
    	cutsInvMTandIso[6] = Weight+"(mu_pt[0]>29 && fabs(el_charge[0])==1 && el_relIso[0] > 0.05 && el_relIso[0] < 0.06 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 &&  nbtag<0.5    &&(npartons==0||npartons>4))";
    }

  TH1D * hist[35];
  TH1D * histInvIso[35];
  TH1D * histInvMT[35];
  TH1D * histInvMTandIso[35];

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
  for (int i=0; i<nSamples; ++i) {
      if (Run == "A"){if(i>=28)continue;}
      else if(Run == "B"){if(i>=28)continue;}
      else if(Run == "C"){if(i>=28)continue;}
      else if(Run == "D"){if(i>=28)continue;}
      else if(Run=="ABC"){if(i>=30)continue;}
      std::cout << sampleNames[i] << std::endl;
      cout<<"Iteration: "<<i<<endl;
      TFile * file = new TFile(directory+sampleNames[i]+".root");
      TH1D * histWeightsH;
      TTree * tree;

      histWeightsH = (TH1D*)file->Get("Wmnu/histWeightsH");
      tree = (TTree*)file->Get("Wmnu/T");
      cout << "sampleNames[i]  "<<sampleNames[i] <<endl;
      cout << "xsec[i]"<<xsec[i]<< endl;
      cout << "lumi"<<lumi<< endl;
      cout << "histWeightsH->GetSumOfWeights()"<<histWeightsH->GetSumOfWeights()<< endl;

      double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
      TString Variable2;

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

      tree->Draw(Variable+">>"+histNameInvIso,cutsInvIso[i]);
      tree->Draw(Variable+">>"+histNameInvMT,cutsInvMT[i]);
      tree->Draw(Variable+">>"+histNameInvMTandIso,cutsInvMTandIso[i]);
      tree->Draw(Variable+">>"+histName,cuts[i]);
  	
      cout << "hist[i]->GetSumOfWeights()  "<<hist[i]->GetSumOfWeights()<< endl;
      cout << "histInvIso[i]->GetSumOfWeights()  "<<histInvIso[i]->GetSumOfWeights()<< endl;
      cout << "histInvMT[i]->GetSumOfWeights()  "<<histInvMT[i]->GetSumOfWeights()<< endl;
      cout << "histInvMTandIso[i]->GetSumOfWeights()  "<<histInvMTandIso[i]->GetSumOfWeights()<< endl;

    if (i>1 &&i<28) {
      for (int iB=1; iB<=nBins; ++iB) {

    	    double x = histInvIso[i]->GetBinContent(iB);
    	    double e = histInvIso[i]->GetBinError(iB);
        	histInvIso[i]->SetBinContent(iB,norm*x);
        	histInvIso[i]->SetBinError(iB,norm*e);

    	    x = hist[i]->GetBinContent(iB);
    	    e = hist[i]->GetBinError(iB);
        	hist[i]->SetBinContent(iB,norm*x);
        	hist[i]->SetBinError(iB,norm*e);
          }
        }
    }

  delete dummyCanv;

  TH1D * histData;
  
  //TODO Set different Run Trigger;
  histData = (TH1D*)hist[0]->Clone("data_obs_"+Variable+suffix);cout<<sampleNames[0]<<endl;
  if (Run=="ABC"){
	histData->Add(histData,hist[28]);cout<<sampleNames[28]<<endl;
	histData->Add(histData,hist[29]);cout<<sampleNames[29]<<endl;
  }
  if(Run=="ABCD"){ 
histData->Add(histData,hist[28]);cout<<sampleNames[28]<<endl;
histData->Add(histData,hist[29]);cout<<sampleNames[29]<<endl;
histData->Add(histData,hist[30]);cout<<sampleNames[30]<<endl;
histData->Add(histData,hist[31]);cout<<sampleNames[31]<<endl;
histData->Add(histData,hist[32]);cout<<sampleNames[32]<<endl;
histData->Add(histData,hist[33]);cout<<sampleNames[33]<<endl;
 }
  cout<<"Hist Data:"<< histData->GetSumOfWeights()<<endl;

  TH1D * DY = (TH1D*)hist[1]->Clone("DY_"+Variable+suffix);cout<<sampleNames[2]<<endl;

  TH1D * WJ = (TH1D*)hist[2]->Clone("WJ_"+Variable+suffix);cout<<sampleNames[3]<<endl;
  WJ->Add(WJ,hist[3]);cout<<sampleNames[3]<<endl;
  WJ->Add(WJ,hist[4]);cout<<sampleNames[4]<<endl;
  WJ->Add(WJ,hist[5]);cout<<sampleNames[5]<<endl;
  WJ->Add(WJ,hist[6]);cout<<sampleNames[6]<<endl;

  TH1D * VV  = (TH1D*)hist[11]->Clone("VV_"+Variable+suffix);cout<<sampleNames[11]<<endl;
  VV->Add(VV,hist[12]);cout<<sampleNames[12]<<endl;
  VV->Add(VV,hist[13]);cout<<sampleNames[13]<<endl;

  TH1D * TT  = (TH1D*)hist[14]->Clone("TT_"+Variable+suffix);cout<<sampleNames[14]<<endl;
  TT->Add(TT,hist[15]);cout<<sampleNames[15]<<endl;
  TT->Add(TT,hist[16]);cout<<sampleNames[16]<<endl;

  TH1D * ST  = (TH1D*)hist[7]->Clone("ST_"+Variable+suffix);cout<<sampleNames[7]<<endl;
  ST->Add(ST,hist[8]);cout<<sampleNames[8]<<endl;
  ST->Add(ST,hist[9]);cout<<sampleNames[9]<<endl;
  ST->Add(ST,hist[10]);cout<<sampleNames[10]<<endl;

  TH1D * QCD  = (TH1D*)hist[17]->Clone("QCD_"+Variable+suffix);
  TH1D * QCDInvIso  = (TH1D*)histInvIso[17]->Clone("QCDInvIso_"+Variable+suffix);
  TH1D * QCDInvMT  = (TH1D*)histInvMT[17]->Clone("QCDInvMT_"+Variable+suffix);
  TH1D * QCDInvMTandIso  = (TH1D*)histInvMTandIso[17]->Clone("QCDInvMTandIso_"+Variable+suffix);

  TH1D * EW  = (TH1D*)hist[1]->Clone("EW_"+Variable+suffix);
  TH1D * EWInvIso  = (TH1D*)histInvIso[1]->Clone("EWInvIso_"+Variable+suffix);
  TH1D * EWInvMT  = (TH1D*)histInvMT[1]->Clone("EWInvMT_"+Variable+suffix);
  TH1D * EWInvMTandIso  = (TH1D*)histInvMTandIso[1]->Clone("EWInvMTandIso_"+Variable+suffix);


  for (int i=2; i<(17); ++i)
	{
  	cout<<"EWK Iteration: "<<i<<endl;
  	EW->Add(hist[i],1);
  	EWInvIso->Add(histInvIso[i],1);
  	cout<<sampleNames[i]<<endl;
  	cout<<sampleNames[i]+"  :"<<endl;
  	cout<<"Sum Of W: "<<hist[i]->GetSumOfWeights()<<endl;
          cout<<"Inv Iso: "<<histInvIso[i]->GetSumOfWeights()<<endl;
  	EWInvMT->Add(histInvMT[i],1);
  	EWInvMTandIso->Add(histInvMTandIso[i],1);
	}
  for (int i=18; i<28; ++i)
	{
        cout<<"QCD Iteration: "<<i<<endl;
        cout<<sampleNames[i]<<endl;
        QCD->Add(hist[i],1);
        QCDInvIso->Add(histInvIso[i],1);
        QCDInvMT->Add(histInvMT[i],1);
        QCDInvMTandIso->Add(histInvMTandIso[i],1);
        cout<<"Sum Of W: "<<hist[i]->GetSumOfWeights()<<endl;
        cout<<"Inv Iso: "<<histInvIso[i]->GetSumOfWeights()<<endl;
	}

  std::cout << "ST   : " << ST->GetSumOfWeights() << " : " << ST->Integral(1,nBins+1) <<" GetSumOfWeights= " << ST->GetSumOfWeights()<< std::endl;
  std::cout << "VV   : " << VV->GetSumOfWeights() << " : " << VV->Integral(1,nBins+1) <<" GetSumOfWeights= " << VV->GetSumOfWeights() <<std::endl;
  std::cout << "QCD   : " << QCD->GetSumOfWeights() << " : " << QCD->Integral(1,nBins+1) <<" GetSumOfWeights= " << QCD->GetSumOfWeights()<< std::endl;
  std::cout << "WJ   : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(1,nBins+1) <<" GetSumOfWeights= " << WJ->GetSumOfWeights()<< std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : " << TT->Integral(1,nBins+1) <<" GetSumOfWeights= " << TT->GetSumOfWeights()<< std::endl;
  std::cout << "DY : " << DY->GetSumOfWeights() << " : " << DY->Integral(1,nBins+1) <<" GetSumOfWeights= " << DY->GetSumOfWeights()<< std::endl;

  //// weight calculation
  std::cout << "WEIGHTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< endl;

   float fData;
   float fDataInvIso;

   if (Run=="A"){
        fData = hist[0]->GetSumOfWeights();
        fDataInvIso = histInvIso[0]->GetSumOfWeights();
   }
   else if(Run=="B"){
        fData = hist[0]->GetSumOfWeights();
        fDataInvIso = histInvIso[0]->GetSumOfWeights();
   }
   else if(Run=="C"){
        fData = hist[0]->GetSumOfWeights();
        fDataInvIso = histInvIso[0]->GetSumOfWeights();
   }
   else if(Run=="D"){
        fData = hist[30]->GetSumOfWeights()+ hist[31]->GetSumOfWeights() + hist[32]->GetSumOfWeights();
        fDataInvIso = histInvIso[30]->GetSumOfWeights()+ hist[31]->GetSumOfWeights() + hist[32]->GetSumOfWeights();

   }
   else if(Run=="ABC"){
        fData = hist[0]->GetSumOfWeights()
               + hist[28]->GetSumOfWeights()
               + hist[29]->GetSumOfWeights();

        fDataInvIso = histInvIso[0]->GetSumOfWeights()
                     + histInvIso[28]->GetSumOfWeights()
                     + histInvIso[29]->GetSumOfWeights();

        }
   else{
        fData = hist[0]->GetSumOfWeights() + hist[28]->GetSumOfWeights() 
              + hist[29]->GetSumOfWeights() + hist[30]->GetSumOfWeights() 
              + hist[31]->GetSumOfWeights() + hist[32]->GetSumOfWeights() 
              + hist[33]->GetSumOfWeights();

        fDataInvIso = histInvIso[0]->GetSumOfWeights()  + histInvIso[28]->GetSumOfWeights() 
              + histInvIso[29]->GetSumOfWeights() + histInvIso[30]->GetSumOfWeights()
              + histInvIso[31]->GetSumOfWeights() + histInvIso[32]->GetSumOfWeights() 
              + histInvIso[33]->GetSumOfWeights();
        }
   float fDataInvIsoFull = fDataInvIso;
   float fDataFull = fData;

   cout <<" fDataFull = "<<fDataFull <<" fDataInvIsoFull = "<<fDataInvIsoFull<<endl;

   float fEW = EW->GetSumOfWeights();
   float fEWInvIso = EWInvIso->GetSumOfWeights();
   float fEWInvIsoFull = fEWInvIso;
   float fEWFull = fEW;

   cout <<" fEWFull = "<<fEWFull <<" fEWInvIsoFull = "<<fEWInvIsoFull<<endl;

   float fQCD = QCD->GetSumOfWeights();
   float fQCDInvIso = QCDInvIso->GetSumOfWeights();
   float fQCDInvMT = QCDInvMT->GetSumOfWeights();
   float fQCDInvMTandIso = QCDInvMTandIso->GetSumOfWeights();
   float fQCDInvIsoFull = fQCDInvIso;
   float fQCDFull = fQCD;

   cout <<" fQCDFull = "<<fQCDFull <<" fQCDInvIsoFull = "<<fQCDInvIsoFull<<endl;

   float weightIsoEW = (1/(fQCDFull*fEWInvIsoFull - fEWFull*fQCDInvIsoFull))*(fQCDFull*fDataInvIsoFull - fQCDInvIsoFull*fDataFull );
   float weightIsoQCD = (1/(fQCDFull*fEWInvIsoFull - fEWFull*fQCDInvIsoFull))*(fEWInvIsoFull*fDataFull - fEWFull*fDataInvIsoFull );

   cout <<" weightIsoEW = "<<weightIsoEW <<" weightIsoQCD = "<<weightIsoQCD<<endl;


//   Calculate Weight Normalisation:
//   float weightMTEW = (1/(fQCD*fEWInvMT - fEW*fQCDInvMT))*(fQCD*fDataInvMT -fQCDInvMT*fData );
//   float weightMTQCD = (1/(fQCD*fEWInvMT - fEW*fQCDInvMT))*(fQCDInvMT*fEW - fEW*fDataInvMT );
//   cout <<" weightMTEW = "<<weightMTEW <<" weightMTQCD = "<<weightMTQCD<<endl;

     QCD->Add(QCD,VV);
     ST->Add(ST,QCD);
     DY->Add(DY,ST);
     TT->Add(TT,DY);
     WJ->Add(WJ,TT);


     std::cout << "BKG : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(0,nBins+1) << std::endl;
     std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
     std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/WJ->GetSumOfWeights() << "+/-"
            << TMath::Sqrt(histData->GetSumOfWeights())/WJ->GetSumOfWeights() << std::endl;

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
    if(logY) h[0]->SetMinimum(1);
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
    //legend->AddEntry(VV,"VV+VVV","f");
    legend->AddEntry(ST,"ST","f");
    legend->AddEntry(QCD,"QCD","f");
    legend->AddEntry(VV,"VV","f");

    WJ->Draw("sameh");
    TT->Draw("sameh");
    DY->Draw("sameh");
    ST->Draw("sameh"); 
    QCD->Draw("sameh");   
    VV->Draw("sameh");

    canv1->Update();
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
    DrawCMSLogo(pads[0], "CMS", "Preliminary Run "+Run, 11, 0.045, 0.035, 1.2);
    DrawTitle(pads[0], to_string(lumi) + "fb^{-1} (13 TeV)", 3);
    if (channel=="Wmnu") DrawTitle(pads[0], "W to #mu#nu", 1);
    if (channel=="Wenu") DrawTitle(pads[0], "W to e#nu", 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();

	  //TODO declare the global Path on the Top of the script;
    if (channel=="Wmnu"){
        canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/"+Run+isLog+"/" +Variable+Suffix+suffix+".pdf");
    }
    if (channel=="Wenu"){
        canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/CMSSW_10_2_10/src/DesyTauAnalyses/NTupleMaker/test/WenuPlots/"+Run+"/"+Variable+Suffix+suffix+".pdf");
    }

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
