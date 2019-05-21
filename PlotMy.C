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
//WithoutNPV/
void PlotMy(TString directory = "/nfs/dust/cms/user/bobovnii/new/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/Wmnu/",
	  TString suffix = "",
	  TString DataFile = "SingleMuon",
	  TString Variable = "mu_relIso[0]",
//	  TString Variable = "npv",
//	  TString Variable = "mu_pt[0]",
//	  TString Variable = "mu_dxy[0]",
//	  TString Variable = "mu_dz[0]",

	  TString Suffix = "_Rereco",
//	  int nBins  =   100,//dxy!!!!!!!!!!!
//	  float xmin =  -0.055,
//	  float xmax =  0.055,
//	  int nBins  =   80,//dz!!!!!!!!!!!
//	  float xmin =  -0.22,
//	  float xmax =  0.22,
//	  int nBins  =   50,
//	  float xmin =    0,
//	  float xmax =  400,

//	  int nBins  =   50,//eta!!!!!!!!!!!!!!!!
//	  float xmin =    -2.5,
//	  float xmax =  2.5,

	  int nBins  =   100,//relIso!!!!!!!!!!
	  float xmin =    0,
	  float xmax =  0.5,

//	  int nBins  =   10,//njets!!!!!
//	  float xmin =    0,
//	  float xmax =  10,

//	  int nBins  =   51,//npv!!!!!
//	  float xmin =    -0.5,
//	  float xmax =  50.5,
	  TString Weight = "gen_weight*pu_weight*LSF_weight*trig_weight*",
//	  TString Weight = "gen_weight*LSF_weight*trig_weight*",
	  //TString Cuts = "&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",&& npv < 20 && npv > 15
// fabs(mu_dz[0])<0.01 && fabs(mu_dxy[0])<0.003 && mu_relIso[0] < 0.15 && event_secondLeptonVeto < 0.5 && event_thirdLeptonVeto < 0.5 && MT<30
	  TString Cuts = "fabs(mu_charge[0])==1 ",
//	  TString Cuts = "MT>60 && MT<120 && mu_relIso[0] < 0.15 && fabs(mu_charge[0])==1", && njets > 0.5 && MT>60
//	  TString Cuts = "MT>60 && MT<120 && fabs(mu_charge[0])==1",
//	  TString Cuts = "&&dzeta_mvamet<-60&&mvamet>80&&iso_1<0.15&&iso_2<0.2&&extraelec_veto<0.5&&extramuon_veto<0.5&&pt_1>13&&pt_2>10&&TMath::Max(pt_1,pt_2)>24",
	  TString xtitle = "mu_relIso",
	  TString ytitle = "Events",
      bool logY = true
          )
{
  //ModTDRStyle();

  bool blindData = false;
  int nbMin = 4;
  int nbMax = 11;
  bool plotLeg = true;
  int position = 0; // 0 - right, 1 - left, 2 - central
  bool showSignal = true;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  //double lumi =23912;// for used json
  //double lumi = 27660;
  //double lumi = 4000;
  double lumi = 15700; // for test BCDE
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

  TString sampleNames[30] = {
    DataFile, // data (0)
    "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",     // isZTT  (1)
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", // isZTT  (2)
    "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",      // (3)
    "TT_TuneCUETP8M1_13TeV-powheg-pythia8_ext4-v1",               // (4)
    "ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (5)
    "ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (6)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (7)
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",// (8)
    "VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8",// (9)
    "WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (10)
    "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (11)
    "WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (12)
    "WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8",// (13)
    "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8",// (14)
    "ZZTo4L_13TeV_powheg_pythia8",// (15)
    "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",// (16)
    "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",// (17)
    "QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8",// (18)
/*    "ST_t-channel_top_4f_leptonDecays_13TeV-powheg", // (7)
    "ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg",    // (8)
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg", // (9)
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg",     // (10)
    "VVTo2L2Nu_13TeV_amcatnloFXFX",    // (11)
    "WWTo1L1Nu2Q_13TeV_amcatnloFXFX",  // (12)
    "WZTo2L2Q_13TeV_amcatnloFXFX",     // (13)
    "WZTo1L1Nu2Q_13TeV_amcatnloFXFX",  // (14)
    "WZTo1L3Nu_13TeV_amcatnloFXFX",    // (15)
    "WZTo3LNu_13TeV-powheg",           // (16)
    "ZZTo4L_13TeV_powheg",             // (17)
    "ZZTo2L2Q_13TeV_amcatnloFXFX",     // (18)
    "WGToLNuG_13TeV-madgraphMLM",             // (19)
    "WGstarToLNuMuMu_012Jets_13TeV-madgraph", // (20)
    "WGstarToLNuEE_012Jets_13TeV-madgraph",   // (21)
    "GluGluHToTauTau_M125_13TeV_powheg",    // (22)
    "VBFHToTauTau_M125_13TeV_powheg", // (23)*/
    ""    // (24);
  };


  double xsec[30] = {1, // data (0)
		     5765,  // DY(50) (1)
		     18610,   // DY(10to50)     (2)
		     Wnorm*61526.7,// WJets (3)
		     TTnorm*831.76,  // TT  (4)
		     44.33, // ST_t-channel_top (5)
		     26.38,  // ST_t-channel_antitop (6)
		     35.6,           // ST_tW_antitop (7)
		     35.6,           // ST_tW_top_5f (8)
		     11.95,  // VV          (9)
		     49.997, // WWToLNuQQ   (10)
		     5.595,  // WZTo2L2Q    (11)
		     10.71,  // WZTo1L1Nu2Q (12)
		     3.05,   // WZTo1L3Nu   (13)
		     4.42965,   // WZTo3L1Nu   (14)
		     1.212,  // ZZTo4L      (15)
		     3.22,   // ZZTo2L2Q    (16)
//             489.0,  // WGToLNuG        (19)
             405.271,  // WGToLNuG        (17)
             720648000.0*0.00042,  // QCD        (18)
             /*2.793,  // WGstarToLNuMuMu (20)
             3.526,  // WGstarToLNuEE   (21)
		     43.92*0.0632,  // gg->H (22)
		     3.748*0.0632,  // VBF H   (23)*/
		     0   // dummy 
  };     


  double xsec1[30] = {1, // data (0)
		     0,  // DY(50) (1)
		     0,   // DY(10to50)     (4)
		     Wnorm*0,// WJets (5)
		     TTnorm*0,  // TT  (6)
		     0, // ST_t-channel_top (7)
		     0,  // ST_t-channel_antitop (8)
		     0,           // ST_tW_antitop (9)
		     0,           // ST_tW_top_5f (10)
		     0,  // VV          (11)
		     0, // WWToLNuQQ   (12)
		     0,  // WZTo2L2Q    (13)
		     0,  // WZTo1L1Nu2Q (14)
		     0,   // WZTo1L3Nu   (15)
		     0,   // WZTo3L1Nu   (16)
		     0,  // ZZTo4L      (17)
		     0,   // ZZTo2L2Q    (18)
             0,  // WGToLNuG        (19)
             720648000.0*0.00042,  // QCD        (19)
		     0   // dummy 
  };     
  
  TString cuts[30];
  TString cutsSS[30];
  


  for (int i=0; i<30; ++i) {
    cuts[i] = Weight+"("+Cuts+")";
    //cutsSS[i] = Weight+qcdweight+"(os<0.5"+Cuts+")";
  }
  cuts[0] = "(met_flag>0.5 && "+Cuts+")";
/*
  cutsSS[0] = qcdweight+"(os<0.5"+Cuts+"&&metFilters>0.5)";
  cutsSS[1] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[2] = Weight+zptmassweight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";
  cutsSS[3] = Weight+qcdweight+"(os<0.5"+Cuts+"&&isZTT)";
  cutsSS[4] = Weight+qcdweight+"(os<0.5"+Cuts+"&&!isZTT)";

  cutsSS[6]  = Weight+topweight+qcdweight+"(os<0.5"+Cuts+")";
*/

  TH1D * hist[30];
  TH1D * histSS[30];

  int nSamples = 19;

  TCanvas * dummyCanv = new TCanvas("dummy","",500,500);
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  // filling histograms
  for (int i=0; i<nSamples; ++i) {
    //    std::cout << sampleNames[i] << std::endl;
    TFile * file = new TFile(directory+sampleNames[i]+".root");
    TH1D * histWeightsH = (TH1D*)file->Get("Wmnu/histWeightsH");
    TTree * tree = (TTree*)file->Get("Wmnu/T"); ////?????????????????????????????????????????????????????????????
cout << "sampleNames[i]  "<<sampleNames[i] <<endl;
cout << "xsec[i]  "<<xsec[i]<< endl;
cout << "lumi  "<<lumi<< endl;
cout << "histWeightsH->GetSumOfWeights()  "<<histWeightsH->GetSumOfWeights()<< endl;
    double norm = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
cout << "norm  "<<norm<< endl;
    TString histName = sampleNames[i] + Variable + "_ss";
    //TString histNameSS = sampleNames[i] + Variable + "_os";
    hist[i] = new TH1D(histName,"",nBins,xmin,xmax);
    //histSS[i] = new TH1D(histNameSS,"",nBins,xmin,xmax);

    hist[i]->Sumw2();
    //histSS[i]->Sumw2();

    tree->Draw(Variable+">>"+histName,cuts[i]);
cout << "hist[i]->GetSumOfWeights()  "<<hist[i]->GetSumOfWeights()<< endl;
    //tree->Draw(Variable+">>"+histNameSS,cutsSS[i]);
    std::cout << sampleNames[i] << " " << hist[i]->GetEntries() << " " << hist[i]->Integral(0,nBins+1) <<" GetSumOfWeights= " << hist[i]->GetSumOfWeights() << std::endl;
    if (i>0) {
      for (int iB=1; iB<=nBins; ++iB) {
	double x = hist[i]->GetBinContent(iB);
	double e = hist[i]->GetBinError(iB);
    	hist[i]->SetBinContent(iB,norm*x);
    	hist[i]->SetBinError(iB,norm*e);
	//double xSS = histSS[i]->GetBinContent(iB);
	//double eSS = histSS[i]->GetBinError(iB);
    	//histSS[i]->SetBinContent(iB,norm*xSS);
    	//histSS[i]->SetBinError(iB,norm*eSS);
      }
    }
  }

  delete dummyCanv;
    
    cout << "!!end!!!!!!!!!!!!!!!!!!!!"<<endl;
   
    



/*
  std::cout << "SS region" << std::endl;
  std::cout << "VV (MC)      : " << histSS[7]->GetSumOfWeights() << " : "<< histSS[7]->Integral(0,nBins+1) << std::endl;
  std::cout << "W  (MC)      : " << histSS[5]->GetSumOfWeights() << " : "<< histSS[5]->Integral(0,nBins+1) << std::endl;
  std::cout << "TT (MC)      : " << histSS[6]->GetSumOfWeights() << " : "<< histSS[6]->Integral(0,nBins+1) <<  std::endl;
  std::cout << "DY (MC)      : " << histSS[1]->GetSumOfWeights() << " : "<< histSS[1]->Integral(0,nBins+1) << std::endl;
  std::cout << "non-QCD (MC) : " << nonQCD << " : " << nonQCDfull << std::endl;
  std::cout << "data         : " << dataSS << " : " << dataSSfull << std::endl;
  std::cout << "W+Jets  fraction : " << Wfraction << " : " << WfractionFull << std::endl;
  std::cout << "non-QCD fraction : " << nonQCDfraction << " : " << nonQCDfractionFull << std::endl; 
  std::cout << std::endl;
*/
  TH1D * histData = (TH1D*)hist[0]->Clone("data_obs");
  TH1D * DY = (TH1D*)hist[1]->Clone("DY");
  DY->Add(DY,hist[2]);
  TH1D * WJ = (TH1D*)hist[3]->Clone("WJ");
  TH1D * TT  = (TH1D*)hist[4]->Clone("TT");
  TH1D * VV  = (TH1D*)hist[9]->Clone("VV");
  VV->Add(VV,hist[10]);
  VV->Add(VV,hist[11]);
  VV->Add(VV,hist[12]);
  VV->Add(VV,hist[13]);
  VV->Add(VV,hist[14]);
  VV->Add(VV,hist[15]);
  VV->Add(VV,hist[16]);
  VV->Add(VV,hist[17]);
  TH1D * ST  = (TH1D*)hist[5]->Clone("ST");
  ST->Add(ST,hist[6]);
  ST->Add(ST,hist[7]);
  ST->Add(ST,hist[8]);
  TH1D * QCD  = (TH1D*)hist[18]->Clone("QCD");


  std::cout << "ST   : " << ST->GetSumOfWeights() << " : " << ST->Integral(1,nBins+1) <<" GetSumOfWeights= " << ST->GetSumOfWeights()<< std::endl;
  std::cout << "VV   : " << VV->GetSumOfWeights() << " : " << VV->Integral(1,nBins+1) <<" GetSumOfWeights= " << VV->GetSumOfWeights() <<std::endl;
  std::cout << "QCD   : " << QCD->GetSumOfWeights() << " : " << QCD->Integral(1,nBins+1) <<" GetSumOfWeights= " << QCD->GetSumOfWeights()<< std::endl;
  std::cout << "WJ   : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(1,nBins+1) <<" GetSumOfWeights= " << WJ->GetSumOfWeights()<< std::endl;
  std::cout << "TT  : " << TT->GetSumOfWeights() << " : " << TT->Integral(1,nBins+1) <<" GetSumOfWeights= " << TT->GetSumOfWeights()<< std::endl;
  std::cout << "DY : " << DY->GetSumOfWeights() << " : " << DY->Integral(1,nBins+1) <<" GetSumOfWeights= " << DY->GetSumOfWeights()<< std::endl;
/*
  float nData = histData->GetSumOfWeights();
  float nTT   = TT->GetSumOfWeights();
  float eData = TMath::Sqrt(nData);
  float nNonTT = 
    VV->GetSumOfWeights() + 
    ZTT->GetSumOfWeights() +
    ZLL->GetSumOfWeights() + 
    QCD->GetSumOfWeights() +
    W->GetSumOfWeights(); 

  float ttScale = (nData-nNonTT)/nTT;
  float ttScaleE = eData/nTT;
  float bkgE = 0.3*nNonTT/nTT;

  //  std::cout << "TT scale factor = " << ttScale << " +/- " << ttScaleE << " +/- " << bkgE << std::endl;

*/
  // ***********************************
  // **** Systematic uncertainties *****
  // ***********************************
 /* TH1D * dummy = (TH1D*)ZTT->Clone("dummy");
  float errQCD = 0.15; // normalization of QCD background
  float errVV = 0.15; // normalization of VV background
  float errW = 0.1; // normalization of W+Jets background
  float errTT = 0.07; // normalization of TT background
  for (int iB=1; iB<=nBins; ++iB) {
    float eQCD = errQCD*QCD->GetBinContent(iB);
    float eVV = errVV*VV->GetBinContent(iB);
    float eW = errW*W->GetBinContent(iB);
    float eTT = errTT*TT->GetBinContent(iB);
    float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eTT*eTT;
    float errTot = TMath::Sqrt(err2);
    dummy->SetBinError(iB,errTot);
    SMH->SetBinError(iB,0);
  }
   */ 
  

 /* VV->Add(VV,QCD);
  W->Add(W,VV);
  ZLL->Add(ZLL,W);
  TT->Add(TT,ZLL);
  ZTT->Add(ZTT,TT);
 */ 
 /* TT->Add(TT,ST);
  VV->Add(VV,TT);
  DY->Add(DY,VV);
  QCD->Add(QCD,DY);
  WJ->Add(WJ,QCD);*/

  QCD->Add(QCD,VV);
  ST->Add(ST,QCD);
  DY->Add(DY,ST);
  TT->Add(TT,DY);
  WJ->Add(WJ,TT);
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
    legend->AddEntry(QCD,"QCD","f");
    legend->AddEntry(DY,"DY","f");
    legend->AddEntry(VV,"VV","f");
    legend->AddEntry(TT,"t#bar{t}","f");

    //  leg->AddEntry(VV,"VV+VVV","f");
    legend->AddEntry(ST,"ST","f");


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
    ratioH->Draw("pe0same");

    pads[0]->cd();
    histData->Draw("pesame");
    
    FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);
    DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.035, 1.2);
    DrawTitle(pads[0], "16 fb^{-1} (13 TeV)", 3);
    DrawTitle(pads[0], "W to #mu#nu", 1);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();
    canv1->Print("/nfs/dust/cms/user/bobovnii/new/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/"+Variable+Suffix+suffix+".pdf");
}
