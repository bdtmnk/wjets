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

void PlotMy2( TString Variable = "puppi_pt",
              TString Suffix = "RunA",
              int nBins  =   80,
              float xmin =    0,
              float xmax =  400,
              TString xtitle = "met_pt",
              bool NormalMetForData = true,
              TString Run = "ABCD",
              TString channel = "Wmnu" )
{
         cout << Variable<<endl;
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

   	if (channel=="Wenu") {

	   DataFile1 = "SingleElectron_RunBCDE";
	   DataFile2 = "SingleElectron__Run2017F-17Nov2017-v1";
	}

	TString Weight = "gen_weight*PUweight*LSF_weight*trig_weight*";
	//PUWeight
	TString MT = "MT";
    	TString Cuts;

	if (Variable.Contains("puppi_pt")) MT = "MTpuppi";
	if (Variable.Contains("met") && Variable.Contains("smeared")) MT = "MT_smeared";
	if (Variable.Contains("MTpuppi")) MT = "MTpuppi";
	if (Variable.Contains("MT") && Variable.Contains("smeared")) MT = "MT_smeared";

	if (channel=="Wmnu")  Cuts = "mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5    ";
	if (channel=="Wenu")  Cuts = "fabs(el_charge[0])==1 && el_relIso[0] < 0.05 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && nbtag<0.5    ";

	Cuts += "&& TMath::IsNaN("+Variable+")<0.5";// && "+Variable+"!=0";
	Cuts += "&& TMath::IsNaN(puppi_pt_JetResUp)<0.5";
	Cuts += "&& TMath::IsNaN(puppi_pt)<0.5";  
	Cuts += "&& TMath::IsNaN(puppi_pt_JetResDown)<0.5";


    TString ytitle = "Events";
    TString suffix = "_Wmnu";
    

    TString directory = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/Wmnu/";
    if (channel=="Wenu") directory = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/Wenu/";

    TString FileToWrite;
    if (channel=="Wmnu") FileToWrite = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/RootWithWeights/"+Run+"/"+Variable+".root";
    TFile * fileW = new TFile(FileToWrite);

   bool logY =true;
   bool blindData = false;
   float QCDscale=1.00;
   bool NewMCmethod = true;
   float  EWweight = 0.805464;//;0.732697;
   float QCDweight = 0.513706;//0.448448;//0.401957;
   if(Run=="ABC"){
   EWweight = 0.620475; QCDweight = 0.298157;
   //EWweight = 0.80871;
   //QCDweight = 0.514931;
   } 
   if (Run=="A"){
   EWweight = 0.779873; //weightIsoEW = 0.779873 weightIsoQCD = 0.494389
   QCDweight = 0.49438;
   }
   if (Run=="B"){// weightIsoEW = 0.75994 weightIsoQCD = 0.503694
   EWweight = 0.75994;
   QCDweight = 0.503694;
   }
   if (Run=="C"){
   EWweight =  0.757434;//weightIsoEW = 0.757434 weightIsoQCD = 0.527242
   QCDweight = 0.527242;}
   if (Run=="D"){
   EWweight= 0.943622;
   QCDweight= 0.604933;}
   if (Run=="ABCD") {
    EWweight = 0.889838;
    QCDweight = 0.580721;
   }
   else{
   EWweight =  0.811365;
   QCDweight = 0.489164;
   //0.805464;//0.732697;
   //QCDweight = 0.513706;//0.448448;
   }

    bool MCDataDriven = true;
    int nbMin = 4;
    int nbMax = 11;
    bool plotLeg = true;
    int position = 0;
    bool showSignal = true;
    int nSamples = 34;
    if (channel=="Wenu") nSamples = 26;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //Lumi Settings:
    double lumi;
    if (Run == "A") {
        lumi = 13480;
    }
    else if (Run=="B"){
        lumi= 6785;
    }
    else if (Run=="C"){
        lumi= 6612;
    }
    else if (Run=="D"){
        lumi=0;
    }
    else if(Run=="ABC"){
        lumi= 13480+6785+6612;
    }
    else if (Run=="ABCD"){
         lumi = 58822.126;
   }
    else{
        lumi = 58822.126;
    }

    double TTnorm = 1.0;
    double Wnorm  = 1.0;

    TString topweight2("");
    TString topweight("topptweight*");
    //TString qcdweight("qcdweight*");
    TString zptmassweight("zptmassweight*");
    TString qcdweight("(qcdweight*qcdweight/qcdweightup)*");

    TString sampleNames[34] = {

         DataFile1, // data (0)

//        "DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8",//2
        "DYJetsToLL_M-50_new",//1
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",//2
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",//3
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",//4
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",//5
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",//6
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
	 DataFile2, // data (B) 28
         DataFile3, // data (C) 29
 	DataFile4, // data (D0) 30
	DataFile5, // D1 31
	DataFile6, // D2 32
	DataFile7  // D3 33
      };

  double xsec[34] = {

        1, // Single Muon data (0)
        2075.14*3 , ////5765,  // DY(50) (1)

        1.221*9644.5, //W1 (2)
        1.221*3144.5, //W1 (3)
        1.221*954.8,  //W1 (4)
        1.221*485.6,  //W1 (5)
        61526.7,// WJets (6)

        44.33,   // ST_t-channel_top (7)
        26.38,   // ST_t-channel_antitop (8)
        35.85,   // ST_tW_antitop (9)
        35.85,   // ST_tW_top_5f (10)

        63.21,   //WW (11)
        47.13,   // WZ (12)
        16.523,  //ZZ (13)

        364.4, // TT (14)
        380.1, // TT (15)
        87.31, // TT (16)

        10.4305*0.15544, //QCD (17)
        469797*0.05362,  //QCD (18)
        117989*0.07335,  //QCD (19)
        558528000*0.0053, //QCD (20)
        7820.25*0.10196,  //QCD (21)
        139803000*0.01182, //QCD (22)
        645.528*0.12242,   //QCD (23)
        19222500*0.02276,  //QCD (24)
        187.109*0.13412,   //QCD (25)
        32.3486*0.14552,   //QCD (26)
        2758420*0.03844,   //QCD (27)
	1,//SM 28 
	1,//SM 29
	1,//SM 30
	1,//SM 31
	1,//SM 32
	1//SM 33

    };

  TString cuts[35];
  TString cutsInvIso[35];
  TString cutsInvMT[35];
  TString cutsInvMTandIso[35];

  for (int i=0; i<35; ++i) {

    cuts[i] = Weight+"("+Cuts+")";

  }
   // for Data
   if (channel=="Wmnu"){
	   Cuts = "mu_pt[0]>29 && fabs(mu_charge[0])==1 && mu_relIso[0] < 0.15 &&  event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && nbtag<0.5  ";
   }

   // for Wjets
   if (channel=="Wmnu"){
        // for Data
        cuts[0] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
        cuts[28] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
        cuts[29] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
	cuts[30] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
	cuts[31] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
	cuts[32] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";
	cuts[33] = "mu_pt[0]>29&&(mu_relIso[0] < 0.15 && met_flag>0.5 && fabs(mu_charge[0])==1 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5  && nbtag<0.5  )";

         // for Wjets
        cuts[6] = Weight+"("+Cuts+" &&(npartons==0||npartons>4))";
	// for DY
	cuts[2] = Weight +"*zptmassweight*" + "("+Cuts+")";
	}


  TH1D * hist[35];
  TH1D * histInvIso[35];
  TH1D * histInvMT[35];
  TH1D * histInvMTandIso[35];



  TCanvas * dummyCanv = new TCanvas("dummy", "", 500, 500);

  float QCD_Weight = 0;
  for (int i=0; i<nSamples; ++i) {
  cout<<"Iteration: "<<i<<endl;
  //if (i==0) continue;
  //if(i>=28) continue;
    TFile * file = new TFile(directory+sampleNames[i]+".root");
    TH1D * histWeightsH;
    TTree * tree;
    histWeightsH = (TH1D*)file->Get("Wmnu/histWeightsH");

    tree = (TTree*)file->Get("Wmnu/T");
    cout << "sampleNames[i]  "<<sampleNames[i] <<endl;
    cout << "xsec[i]  "<<xsec[i]<< endl;
    cout << "lumi  "<<lumi<< endl;
    cout << "histWeightsH->GetSumOfWeights()  "<<histWeightsH->GetSumOfWeights()<< endl;

    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
    cout << "norm  "<<norm<< endl;
    TString Variable2;
    cout << "Variable   "<<Variable <<endl;
    Variable2 = Variable;//"MT";

    TString histName = sampleNames[i] + Variable2 + "_";
    TString histNameInvIso = sampleNames[i] + Variable2 + "_InvIso";
    TString histNameInvMT = sampleNames[i] + Variable2 + "_InvMT";
    TString histNameInvMTandIso = sampleNames[i] + Variable2 + "_InvMTandIso";
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
    if (i>1 && i<28) {
      for (int iB=1; iB<=nBins; ++iB) {
        double x,e;
    if (!NewMCmethod){
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



  }

      delete dummyCanv;

      TH1D * histData = (TH1D*)hist[0]->Clone("data_obs_"+Variable+suffix);
      histData->Add(histData,hist[28]);cout<<sampleNames[28]<<endl;
      histData->Add(histData,hist[29]);cout<<sampleNames[29]<<endl;
      histData->Add(histData,hist[30]);cout<<sampleNames[30]<<endl;
      histData->Add(histData,hist[31]);cout<<sampleNames[31]<<endl;
      histData->Add(histData,hist[32]);cout<<sampleNames[32]<<endl;
      histData->Add(histData,hist[33]);cout<<sampleNames[33]<<endl;

//      cout<<sampleNames[0]<<endl;
      cout<<"Hist Data:"<< histData->GetSumOfWeights()<<endl;

      TH1D * DY = (TH1D*)hist[1]->Clone("DY_"+Variable+suffix);cout<<sampleNames[1]<<endl;

      TH1D * WJ = (TH1D*)hist[2]->Clone("WJ_"+Variable+suffix);cout<<sampleNames[2]<<endl;
      WJ->Add(WJ,hist[3]);cout<<sampleNames[3]<<endl;
      WJ->Add(WJ,hist[4]);cout<<sampleNames[4]<<endl;
      WJ->Add(WJ,hist[5]);cout<<sampleNames[5]<<endl;
      WJ->Add(WJ,hist[6]);cout<<sampleNames[6]<<endl;

      TH1D * VV  = (TH1D*)hist[11]->Clone("VV_"+Variable+suffix);cout<<sampleNames[11]<<endl;
      VV->Add(VV,hist[12]);cout<<sampleNames[12]<<endl;
      VV->Add(VV,hist[13]);cout<<sampleNames[13]<<endl;

      TH1D * TT  = (TH1D*)hist[14]->Clone("TT_"+Variable+suffix);cout<<sampleNames[14]<<endl;
      TT->Add(TT,hist[15]);cout<<sampleNames[14]<<endl;
      TT->Add(TT,hist[16]);cout<<sampleNames[16]<<endl;

      TH1D * ST  = (TH1D*)hist[7]->Clone("ST_"+Variable+suffix);cout<<sampleNames[7]<<endl;
      ST->Add(ST,hist[8]);cout<<sampleNames[8]<<endl;
      ST->Add(ST,hist[9]);cout<<sampleNames[9]<<endl;
      ST->Add(ST,hist[10]);cout<<sampleNames[10]<<endl;


      TH1D * QCD;

    if (!NewMCmethod){
        histInvIso[0]->Add(histInvIso[1],1);
        if (MCDataDriven)  histInvMT[0]->Add(histInvMT[1],1);
        if (MCDataDriven)  histInvMTandIso[0]->Add(histInvMTandIso[1],1);

    for (int i=1; i<(28); ++i)
        {
        histInvIso[0]->Add(histInvIso[i],-1);
        if (MCDataDriven)	histInvMT[0]->Add(histInvMT[i],-1);
        if (MCDataDriven)	histInvMTandIso[0]->Add(histInvMTandIso[i],-1);
        }
        QCD  = (TH1D*)histInvIso[18]->Clone("QCD_"+Variable+suffix);
    if (MCDataDriven)QCD->Scale(histInvMT[18]->GetSumOfWeights()/histInvMTandIso[18]->GetSumOfWeights());

    if (MCDataDriven)cout << "QCD scale = " << histInvMT[0]->GetSumOfWeights()/histInvMTandIso[0]->GetSumOfWeights() << endl;
    if (!MCDataDriven)QCD->Scale(QCDscale);
    }

    if (NewMCmethod){
    QCD  = (TH1D*)hist[17]->Clone("QCD_"+Variable+suffix);
    QCD->Add(QCD,hist[18]);
    QCD->Add(QCD,hist[19]);
    QCD->Add(QCD,hist[20]);
    QCD->Add(QCD,hist[21]);
    QCD->Add(QCD,hist[22]);
    QCD->Add(QCD,hist[23]);
    QCD->Add(QCD,hist[24]);
    QCD->Add(QCD,hist[25]);
    QCD->Add(QCD,hist[26]);
    QCD->Add(QCD,hist[27]);


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
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
    TString _lumi;
    if (Run=="ABD"){_lumi = "26.87";}
    if (Run=="D"){_lumi = "31.95";}
    if (Run=="ABCD"){_lumi = "58.82";}

    DrawTitle(pads[0], _lumi + " fb^{-1} (13 TeV)", 3);//13480
   	//if (channel=="Wmnu") DrawTitle(pads[0], "W to #mu#nu", 1.5);
    //	  if (channel=="Wenu") DrawTitle(pads[0], "W to e#nu", 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();

    if (channel=="Wmnu"){
    canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/"+Run+"/" +Variable+Suffix+suffix+".pdf");
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
