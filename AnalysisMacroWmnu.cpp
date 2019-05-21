#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>

#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
int main(int argc, char * argv[]) {



  // **** configuration
  Config cfg(argv[1]);
  string Channel="Wmnu";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
  
  
  const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

////////////muons

  const double ptMuonCut   = cfg.get<double>("ptMuonCut");
  const double etaMuonCut     = cfg.get<double>("etaMuonCut");
  const double dxyMuonCut     = cfg.get<double>("dxyMuonCut");
  const double dzMuonCut      = cfg.get<double>("dzMuonCut");
  const double isoMuonLowCut  = cfg.get<double>("isoMuonLowCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

////// tau
  const double ptTauCut = cfg.get<double>("ptTauCut"); 
  const double etaTauCut = cfg.get<double>("etaTauCut"); 
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double  againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
  const double  againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double  vertexz =  cfg.get<double>("vertexz");
  const double  byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");


  // veto muons
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // dilemuon veto 
  const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
  const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
  const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
  const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
  const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");


//veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");



  const string dataBaseDir = cfg.get<string>("DataBaseDir");



  const string MuonidIsoEffFile = cfg.get<string>("MuonidIsoEffFile");
  const string MuontrigEffFile = cfg.get<string>("MuontrigEffFile");


  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");


  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");


  // topSingleMuonTriggerFile
  const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  const string TrigLeg  = cfg.get<string>("SingleMuonFilterName") ;
  const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");



  // vertex distributions filenames and histname


  const string jsonFile = cfg.get<string>("jsonFile");

  string cmsswBase = (getenv ("CMSSW_BASE"));
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
 
  //RecoilCorrector recoilMetCorrector("DesyTauAnalyses/NTupleMaker/data/PFMET_Run2016BCDEFGH_Spring16.root");

  MEtSys metSys("HTT-utilities/RecoilCorrections/data/MEtSys.root");

//  const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");

  // Run-lumi selector
  std::vector<Period> periods;  
  if (isData) { // read the good runs 
    std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
    if (inputFileStream.fail() ) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
    }

    for(std::string s; std::getline(inputFileStream, s); ) {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }
  TString MainTrigger(TrigLeg);


  const double bTag   = cfg.get<double>("bTag");

  CutList.clear();
  CutList.push_back("No cut");
  CutList.push_back("No cut after PU");
  CutList.push_back("mu-tau");
  CutList.push_back("2nd lepV");
  CutList.push_back("3rd lepV");
  CutList.push_back("Trigger");
  CutList.push_back("Lepton SF");
  CutList.push_back("TauFakeRate");
  CutList.push_back("topPtRwgt");
  CutList.push_back("Cleaned jets");
  CutList.push_back("b-Veto");


  int CutNumb = (int)CutList.size();
  xs=1;fact=1;fact2=1;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  ifstream ifs("xsecs");
  string line;


      XSec=1.;
   xsecs=XSec;
  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;

// PU reweighting


  PileUp * PUofficial = new PileUp();
  //TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_Cert_271036-277148_13TeV_PromptReco_Collisions16_xsec69p2mb.root","read");
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileup_data_2017Rereco_80bins.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileup_MC_Fall17_80bins.root", "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);


	cout << "PU init"<< endl;
  TFile * filePUdistributionRunBCDE = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/PUforMETrunBCDE.root", "read");
  TFile * filePUdistributionRunF = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PUforMETrunF.root", "read");

  TH1D * PU_RunBCDE = (TH1D *)filePUdistributionRunBCDE->Get("official");
  TH1D * PU_RunF = (TH1D *)filePUdistributionRunF->Get("official");

  TH1D * PU_mc2 = (TH1D *)filePUdistributionRunBCDE->Get("private");
cout << PU_RunBCDE->GetSumOfWeights()<< endl;
cout << PU_RunF->GetSumOfWeights()<< endl;
cout << PU_mc2->GetSumOfWeights()<< endl;


  /*
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);
  */
bool SUSY = false;
  	string BTag_ = "central";

 string BtagCVS = "CSVv2_94XSF_V2_B_F.csv" ;  
  if (SUSY) BtagCVS = "fastsim_csvv2_ttbar_26_1_2017.csv";

  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/"+BtagCVS);


  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,BTag_);
  if (!SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"comb");
  reader_C.load(calib,BTagEntry::FLAV_C,"comb");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");}
  if (SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"fastsim");
  reader_C.load(calib,BTagEntry::FLAV_C,"fastsim");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"fastsim");}



  float etaBTAG[2] = {0.5,2.1};
  float ptBTAG[5] = {25.,35.,50.,100.,200.};

  std::cout << std::endl;
  for (int iEta=0; iEta<2; ++iEta) {
    for (int iPt=0; iPt<5; ++iPt) {
      float sfB = reader_B.eval_auto_bounds(BTag_,BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC = reader_C.eval_auto_bounds(BTag_,BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight = reader_Light.eval_auto_bounds(BTag_,BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
      printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
    }
  }
  std::cout << std::endl;

  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_ichep2016.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;
	
  float MaxBJetPt = 1000.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;

 // ScaleFactor * SF_HIP; 
//    SF_HIP = new ScaleFactor();
//    SF_HIP->init_ScaleFactor(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker//test/HIP.root");
  // Lepton Scale Factors 

  TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);

  cout<<"  Initializing iD SF files....."<<endl;

//  TFile *fileSF_BCDEF = new TFile("/nfs/dust/cms/user/bobovnii/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/SFforW/EfficienciesAndSF_BCDEF.root");
//  TFile *fileSF_GH = new TFile("/nfs/dust/cms/user/bobovnii/CMSSW_8_0_25/src/DesyTauAnalyses/NTupleMaker/test/SFforW/EfficienciesAndSF_GH.root");
//    TH2F * hSF_BCDEF = (TH2F*)fileSF_BCDEF->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
//    TH2F * hSF_GH = (TH2F*)fileSF_GH->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
//cout <<hSF_BCDEF->GetBinContent(1,1) <<endl;

  ScaleFactor * SF_muonIdIso; 
  if (applyLeptonSF) {
    SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));
    //special for RunG
    //SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFileRunGIso));
  }


  cout<<"  Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuontrigEffFile));

cout<<" ended initialization here "<<endl;
  char ff[100];
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  int nTotalFiles = 0;

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  //std::ifstream fileList(argv[2]);
  std::ifstream fileList(ff);
  //std::ifstream fileList0(argv[2]);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString era=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::cout <<" The filename will be "<<TStrName <<"  "<<datasetName<<std::endl;  
  // output fileName with histograms

std::string st1,st2;
//bool SUSY = false;


  TFile * file;
  if (isData) file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");
  if (!isData) file = new TFile(era+"/"+TStrName+TString(".root"),"update");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;


  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;
  //string treename = rootFileName+"_tree.root";
  
// additional branches
  SetupTree(); 
  Float_t         MT;

  T->Branch("MT", &MT, "MT/F");
  Float_t         MTpuppi;
  T->Branch("MTpuppi", &MTpuppi, "MTpuppi/F");
  Float_t         MT_smeared;
  T->Branch("MT_smeared", &MT_smeared, "MT_smeared/F");

  Float_t         MT_PF;
  T->Branch("MT_PF", &MT_PF, "MT_PF/F");

  Float_t      puppi_ex;
  Float_t      puppi_ey;
  Float_t      puppi_ez;
  Float_t      puppi_pt;
  Float_t      puppi_phi;
  Float_t     puppi_ex_JetEnUp;
  Float_t     puppi_ey_JetEnUp;
  Float_t      puppi_ex_JetEnDown;
  Float_t      puppi_ey_JetEnDown;
  Float_t      puppi_ex_UnclusteredEnUp;
  Float_t      puppi_ey_UnclusteredEnUp;
  Float_t      puppi_ex_UnclusteredEnDown;
  Float_t      puppi_ey_UnclusteredEnDown;

  T->Branch("puppi_ex", &puppi_ex, "puppi_ex/F");
  T->Branch("puppi_ey", &puppi_ey, "puppi_ey/F");
  T->Branch("puppi_ez", &puppi_ez, "puppi_ez/F");
  T->Branch("puppi_pt", &puppi_pt, "puppi_pt/F");
  T->Branch("puppi_phi", &puppi_phi, "puppi_phi/F");
  T->Branch("puppi_ex_JetEnUp", &puppi_ex_JetEnUp, "puppi_ex_JetEnUp/F");
  T->Branch("puppi_ey_JetEnUp", &puppi_ey_JetEnUp, "puppi_ey_JetEnUp/F");
  T->Branch("puppi_ex_JetEnDown", &puppi_ex_JetEnDown, "puppi_ex_JetEnDown/F");
  T->Branch("puppi_ey_JetEnDown", &puppi_ey_JetEnDown, "puppi_ey_JetEnDown/F");
  T->Branch("puppi_ex_UnclusteredEnUp", &puppi_ex_UnclusteredEnUp, "puppi_ex_UnclusteredEnUp/F");
  T->Branch("puppi_ey_UnclusteredEnUp", &puppi_ey_UnclusteredEnUp, "puppi_ey_UnclusteredEnUp/F");
  T->Branch("puppi_ex_UnclusteredEnDown", &puppi_ex_UnclusteredEnDown, "puppi_ex_UnclusteredEnDown/F");
  T->Branch("puppi_ey_UnclusteredEnDown", &puppi_ey_UnclusteredEnDown, "puppi_ey_UnclusteredEnDown/F");


//if (!isData)
//{
  Float_t      PFmet_ex;
  Float_t      PFmet_ey;
  Float_t      PFmet_ez;
  Float_t      PFmet_pt;
  Float_t      PFmet_phi;
  Float_t     PFmet_ex_JetEnUp;
  Float_t     PFmet_ey_JetEnUp;
  Float_t      PFmet_ex_JetEnDown;
  Float_t      PFmet_ey_JetEnDown;
  Float_t      PFmet_ex_UnclusteredEnUp;
  Float_t      PFmet_ey_UnclusteredEnUp;
  Float_t      PFmet_ex_UnclusteredEnDown;
  Float_t      PFmet_ey_UnclusteredEnDown;
   Float_t PFmet_pt_JetEnUp;
   Float_t PFmet_pt_JetEnDown;
   Float_t PFmet_pt_UnclusteredEnUp;
   Float_t PFmet_pt_UnclusteredEnDown;

   Float_t MT_PF_JetEnUp;
   Float_t MT_PF_JetEnDown;
   Float_t MT_PF_UnclusteredEnUp;
   Float_t MT_PF_UnclusteredEnDown;






  Float_t      puppi_ex_JetResDown;
  Float_t      puppi_ey_JetResDown;
  Float_t      puppi_ex_JetResUp;
  Float_t      puppi_ey_JetResUp;

  Float_t      met_ex_JetResDown;
  Float_t      met_ey_JetResDown;
  Float_t      met_ex_JetResUp;
  Float_t      met_ey_JetResUp;

  Float_t      met_ex_smeared;
  Float_t      met_ey_smeared;
  Float_t      met_ez_smeared;
  Float_t      met_pt_smeared;
  Float_t      met_phi_smeared;
  Float_t     met_ex_JetEnUp_smeared;
  Float_t     met_ey_JetEnUp_smeared;
  Float_t      met_ex_JetEnDown_smeared;
  Float_t      met_ey_JetEnDown_smeared;
  Float_t      met_ex_UnclusteredEnUp_smeared;
  Float_t      met_ey_UnclusteredEnUp_smeared;
  Float_t      met_ex_UnclusteredEnDown_smeared;
  Float_t      met_ey_UnclusteredEnDown_smeared;
  Float_t      met_ex_JetResUp_smeared;
  Float_t      met_ey_JetResUp_smeared;
  Float_t      met_ex_JetResDown_smeared;
  Float_t      met_ey_JetResDown_smeared;



   Float_t met_pt_JetEnUp_smeared;
   Float_t met_pt_JetEnDown_smeared;
   Float_t met_pt_UnclusteredEnUp_smeared;
   Float_t met_pt_UnclusteredEnDown_smeared;
   Float_t met_pt_JetResUp_smeared;
   Float_t met_pt_JetResDown_smeared;

   Float_t met_pt_JetEnUp;
   Float_t met_pt_JetEnDown;
   Float_t met_pt_UnclusteredEnUp;
   Float_t met_pt_UnclusteredEnDown;
   Float_t met_pt_JetResUp;
   Float_t met_pt_JetResDown;

   Float_t puppi_pt_JetEnUp;
   Float_t puppi_pt_JetEnDown;
   Float_t puppi_pt_UnclusteredEnUp;
   Float_t puppi_pt_UnclusteredEnDown;
   Float_t puppi_pt_JetResUp;
   Float_t puppi_pt_JetResDown;

   Float_t MT_JetEnUp_smeared;
   Float_t MT_JetEnDown_smeared;
   Float_t MT_UnclusteredEnUp_smeared;
   Float_t MT_UnclusteredEnDown_smeared;
   Float_t MT_JetResUp_smeared;
   Float_t MT_JetResDown_smeared;

   Float_t MT_JetEnUp;
   Float_t MT_JetEnDown;
   Float_t MT_UnclusteredEnUp;
   Float_t MT_UnclusteredEnDown;
   Float_t MT_JetResUp;
   Float_t MT_JetResDown;

   Float_t MTpuppi_JetEnUp;
   Float_t MTpuppi_JetEnDown;
   Float_t MTpuppi_UnclusteredEnUp;
   Float_t MTpuppi_UnclusteredEnDown;
   Float_t MTpuppi_JetResUp;
   Float_t MTpuppi_JetResDown;


  T->Branch("puppi_ex_JetResDown", &puppi_ex_JetResDown, "puppi_ex_JetResDown/F");
  T->Branch("puppi_ey_JetResDown", &puppi_ey_JetResDown, "puppi_ey_JetResDown/F");
  T->Branch("puppi_ex_JetResUp", &puppi_ex_JetResUp, "puppi_ex_JetResUp/F");
  T->Branch("puppi_ey_JetResUp", &puppi_ey_JetResUp, "puppi_ey_JetResUp/F");

  T->Branch("met_ex_JetResDown", &met_ex_JetResDown, "met_ex_JetResDown/F");
  T->Branch("met_ey_JetResDown", &met_ey_JetResDown, "met_ey_JetResDown/F");
  T->Branch("met_ex_JetResUp", &met_ex_JetResUp, "met_ex_JetResUp/F");
  T->Branch("met_ey_JetResUp", &met_ey_JetResUp, "met_ey_JetResUp/F");

  T->Branch("met_ex_smeared", &met_ex_smeared, "met_ex_smeared/F");
  T->Branch("met_ey_smeared", &met_ey_smeared, "met_ey_smeared/F");
  T->Branch("met_ez_smeared", &met_ez_smeared, "met_ez_smeared/F");
  T->Branch("met_pt_smeared", &met_pt_smeared, "met_pt_smeared/F");
  T->Branch("met_phi_smeared", &met_phi_smeared, "met_phi_smeared/F");
  T->Branch("met_ex_JetEnUp_smeared", &met_ex_JetEnUp_smeared, "met_ex_JetEnUp_smeared/F");
  T->Branch("met_ey_JetEnUp_smeared", &met_ey_JetEnUp_smeared, "met_ey_JetEnUp_smeared/F");
  T->Branch("met_ex_JetEnDown_smeared", &met_ex_JetEnDown_smeared, "met_ex_JetEnDown_smeared/F");
  T->Branch("met_ey_JetEnDown_smeared", &met_ey_JetEnDown_smeared, "met_ey_JetEnDown_smeared/F");
  T->Branch("met_ex_UnclusteredEnUp_smeared", &met_ex_UnclusteredEnUp_smeared, "met_ex_UnclusteredEnUp_smeared/F");
  T->Branch("met_ey_UnclusteredEnUp_smeared", &met_ey_UnclusteredEnUp_smeared, "met_ey_UnclusteredEnUp_smeared/F");
  T->Branch("met_ex_UnclusteredEnDown_smeared", &met_ex_UnclusteredEnDown_smeared, "met_ex_UnclusteredEnDown_smeared/F");
  T->Branch("met_ey_UnclusteredEnDown_smeared", &met_ey_UnclusteredEnDown_smeared, "met_ey_UnclusteredEnDown_smeared/F");
  T->Branch("met_ex_JetResUp_smeared", &met_ex_JetResUp_smeared, "met_ex_JetResUp_smeared/F");
  T->Branch("met_ey_JetResUp_smeared", &met_ey_JetResUp_smeared, "met_ey_JetResUp_smeared/F");
  T->Branch("met_ex_JetResDown_smeared", &met_ex_JetResDown_smeared, "met_ex_JetResDown_smeared/F");
  T->Branch("met_ey_JetResDown_smeared", &met_ey_JetResDown_smeared, "met_ey_JetResDown_smeared/F");


  T->Branch("met_pt_JetEnUp_smeared", &met_pt_JetEnUp_smeared, "met_pt_JetEnUp_smeared/F");
  T->Branch("met_pt_JetEnDown_smeared", &met_pt_JetEnDown_smeared, "met_pt_JetEnDown_smeared/F");
  T->Branch("met_pt_UnclusteredEnUp_smeared", &met_pt_UnclusteredEnUp_smeared, "met_pt_UnclusteredEnUp_smeared/F");
  T->Branch("met_pt_UnclusteredEnDown_smeared", &met_pt_UnclusteredEnDown_smeared, "met_pt_UnclusteredEnDown_smeared/F");
  T->Branch("met_pt_JetResUp_smeared", &met_pt_JetResUp_smeared, "met_pt_JetResUp_smeared/F");
  T->Branch("met_pt_JetResDown_smeared", &met_pt_JetResDown_smeared, "met_pt_JetResDown_smeared/F");

  T->Branch("met_pt_JetEnUp", &met_pt_JetEnUp, "met_pt_JetEnUp/F");
  T->Branch("met_pt_JetEnDown", &met_pt_JetEnDown, "met_pt_JetEnDown/F");
  T->Branch("met_pt_UnclusteredEnUp", &met_pt_UnclusteredEnUp, "met_pt_UnclusteredEnUp/F");
  T->Branch("met_pt_UnclusteredEnDown", &met_pt_UnclusteredEnDown, "met_pt_UnclusteredEnDown/F");
  T->Branch("met_pt_JetResUp", &met_pt_JetResUp, "met_pt_JetResUp/F");
  T->Branch("met_pt_JetResDown", &met_pt_JetResDown, "met_pt_JetResDown/F");

  T->Branch("puppi_pt_JetEnUp", &puppi_pt_JetEnUp, "puppi_pt_JetEnUp/F");
  T->Branch("puppi_pt_JetEnDown", &puppi_pt_JetEnDown, "puppi_pt_JetEnDown/F");
  T->Branch("puppi_pt_UnclusteredEnUp", &puppi_pt_UnclusteredEnUp, "puppi_pt_UnclusteredEnUp/F");
  T->Branch("puppi_pt_UnclusteredEnDown", &puppi_pt_UnclusteredEnDown, "puppi_pt_UnclusteredEnDown/F");
  T->Branch("puppi_pt_JetResUp", &puppi_pt_JetResUp, "puppi_pt_JetResUp/F");
  T->Branch("puppi_pt_JetResDown", &puppi_pt_JetResDown, "puppi_pt_JetResDown/F");

  T->Branch("MT_JetEnUp_smeared", &MT_JetEnUp_smeared, "MT_JetEnUp_smeared/F");
  T->Branch("MT_JetEnDown_smeared", &MT_JetEnDown_smeared, "MT_JetEnDown_smeared/F");
  T->Branch("MT_UnclusteredEnUp_smeared", &MT_UnclusteredEnUp_smeared, "MT_UnclusteredEnUp_smeared/F");
  T->Branch("MT_UnclusteredEnDown_smeared", &MT_UnclusteredEnDown_smeared, "MT_UnclusteredEnDown_smeared/F");
  T->Branch("MT_JetResUp_smeared", &MT_JetResUp_smeared, "MT_JetResUp_smeared/F");
  T->Branch("MT_JetResDown_smeared", &MT_JetResDown_smeared, "MT_JetResDown_smeared/F");

  T->Branch("MT_JetEnUp", &MT_JetEnUp, "MT_JetEnUp/F");
  T->Branch("MT_JetEnDown", &MT_JetEnDown, "MT_JetEnDown/F");
  T->Branch("MT_UnclusteredEnUp", &MT_UnclusteredEnUp, "MT_UnclusteredEnUp/F");
  T->Branch("MT_UnclusteredEnDown", &MT_UnclusteredEnDown, "MT_UnclusteredEnDown/F");
  T->Branch("MT_JetResUp", &MT_JetResUp, "MT_JetResUp/F");
  T->Branch("MT_JetResDown", &MT_JetResDown, "MT_JetResDown/F");

  T->Branch("MTpuppi_JetEnUp", &MTpuppi_JetEnUp, "MTpuppi_JetEnUp/F");
  T->Branch("MTpuppi_JetEnDown", &MTpuppi_JetEnDown, "MTpuppi_JetEnDown/F");
  T->Branch("MTpuppi_UnclusteredEnUp", &MTpuppi_UnclusteredEnUp, "MTpuppi_UnclusteredEnUp/F");
  T->Branch("MTpuppi_UnclusteredEnDown", &MTpuppi_UnclusteredEnDown, "MTpuppi_UnclusteredEnDown/F");
  T->Branch("MTpuppi_JetResUp", &MTpuppi_JetResUp, "MTpuppi_JetResUp/F");
  T->Branch("MTpuppi_JetResDown", &MTpuppi_JetResDown, "MTpuppi_JetResDown/F");



  T->Branch("PFmet_ex", &PFmet_ex, "PFmet_ex/F");
  T->Branch("PFmet_ey", &PFmet_ey, "PFmet_ey/F");
  T->Branch("PFmet_ez", &PFmet_ez, "PFmet_ez/F");
  T->Branch("PFmet_pt", &PFmet_pt, "PFmet_pt/F");
  T->Branch("PFmet_phi", &PFmet_phi, "PFmet_phi/F");
  T->Branch("PFmet_ex_JetEnUp", &PFmet_ex_JetEnUp, "PFmet_ex_JetEnUp/F");
  T->Branch("PFmet_ey_JetEnUp", &PFmet_ey_JetEnUp, "PFmet_ey_JetEnUp/F");
  T->Branch("PFmet_ex_JetEnDown", &PFmet_ex_JetEnDown, "PFmet_ex_JetEnDown/F");
  T->Branch("PFmet_ey_JetEnDown", &PFmet_ey_JetEnDown, "PFmet_ey_JetEnDown/F");
  T->Branch("PFmet_ex_UnclusteredEnUp", &PFmet_ex_UnclusteredEnUp, "PFmet_ex_UnclusteredEnUp/F");
  T->Branch("PFmet_ey_UnclusteredEnUp", &PFmet_ey_UnclusteredEnUp, "PFmet_ey_UnclusteredEnUp/F");
  T->Branch("PFmet_ex_UnclusteredEnDown", &PFmet_ex_UnclusteredEnDown, "PFmet_ex_UnclusteredEnDown/F");
  T->Branch("PFmet_ey_UnclusteredEnDown", &PFmet_ey_UnclusteredEnDown, "PFmet_ey_UnclusteredEnDown/F");
  T->Branch("PFmet_pt_JetEnUp", &PFmet_pt_JetEnUp, "PFmet_pt_JetEnUp/F");
  T->Branch("PFmet_pt_JetEnDown", &PFmet_pt_JetEnDown, "PFmet_pt_JetEnDown/F");
  T->Branch("PFmet_pt_UnclusteredEnUp", &PFmet_pt_UnclusteredEnUp, "PFmet_pt_UnclusteredEnUp/F");
  T->Branch("PFmet_pt_UnclusteredEnDown", &PFmet_pt_UnclusteredEnDown, "PFmet_pt_UnclusteredEnDown/F");

  T->Branch("MT_PF_JetEnUp", &MT_PF_JetEnUp, "MT_PF_JetEnUp/F");
  T->Branch("MT_PF_JetEnDown", &MT_PF_JetEnDown, "MT_PF_JetEnDown/F");
  T->Branch("MT_PF_UnclusteredEnUp", &MT_PF_UnclusteredEnUp, "MT_PF_UnclusteredEnUp/F");
  T->Branch("MT_PF_UnclusteredEnDown", &MT_PF_UnclusteredEnDown, "MT_PF_UnclusteredEnDown/F");


/////MET performance 

   Float_t Ut;
   Float_t Utr;
   Float_t Ucol;
  T->Branch("Ut", &Ut, "Ut/F");
  T->Branch("Utr", &Utr, "Utr/F");
  T->Branch("Ucol", &Ucol, "Ucol/F");


   Float_t Ut_JetEnUp;
   Float_t Utr_JetEnUp;
   Float_t Ucol_JetEnUp;
  T->Branch("Ut_JetEnUp", &Ut_JetEnUp, "Ut_JetEnUp/F");
  T->Branch("Utr_JetEnUp", &Utr_JetEnUp, "Utr_JetEnUp/F");
  T->Branch("Ucol_JetEnUp", &Ucol_JetEnUp, "Ucol_JetEnUp/F");

   Float_t Ut_JetEnDown;
   Float_t Utr_JetEnDown;
   Float_t Ucol_JetEnDown;
  T->Branch("Ut_JetEnDown", &Ut_JetEnDown, "Ut_JetEnDown/F");
  T->Branch("Utr_JetEnDown", &Utr_JetEnDown, "Utr_JetEnDown/F");
  T->Branch("Ucol_JetEnDown", &Ucol_JetEnDown, "Ucol_JetEnDown/F");

   Float_t Ut_JetResDown;
   Float_t Utr_JetResDown;
   Float_t Ucol_JetResDown;
  T->Branch("Ut_JetResDown", &Ut_JetResDown, "Ut_JetResDown/F");
  T->Branch("Utr_JetResDown", &Utr_JetResDown, "Utr_JetResDown/F");
  T->Branch("Ucol_JetResDown", &Ucol_JetResDown, "Ucol_JetResDown/F");

   Float_t Ut_JetResUp;
   Float_t Utr_JetResUp;
   Float_t Ucol_JetResUp;
  T->Branch("Ut_JetResUp", &Ut_JetResUp, "Ut_JetResUp/F");
  T->Branch("Utr_JetResUp", &Utr_JetResUp, "Utr_JetResUp/F");
  T->Branch("Ucol_JetResUp", &Ucol_JetResUp, "Ucol_JetResUp/F");

   Float_t Ut_UnclusteredEnUp;
   Float_t Utr_UnclusteredEnUp;
   Float_t Ucol_UnclusteredEnUp;
  T->Branch("Ut_UnclusteredEnUp", &Ut_UnclusteredEnUp, "Ut_UnclusteredEnUp/F");
  T->Branch("Utr_UnclusteredEnUp", &Utr_UnclusteredEnUp, "Utr_UnclusteredEnUp/F");
  T->Branch("Ucol_UnclusteredEnUp", &Ucol_UnclusteredEnUp, "Ucol_UnclusteredEnUp/F");

   Float_t Ut_UnclusteredEnDown;
   Float_t Utr_UnclusteredEnDown;
   Float_t Ucol_UnclusteredEnDown;
  T->Branch("Ut_UnclusteredEnDown", &Ut_UnclusteredEnDown, "Ut_UnclusteredEnDown/F");
  T->Branch("Utr_UnclusteredEnDown", &Utr_UnclusteredEnDown, "Utr_UnclusteredEnDown/F");
  T->Branch("Ucol_UnclusteredEnDown", &Ucol_UnclusteredEnDown, "Ucol_UnclusteredEnDown/F");

   Float_t UtMC;
   Float_t UtrMC;
   Float_t UcolMC;
  T->Branch("UtMC", &UtMC, "UtMC/F");
  T->Branch("UtrMC", &UtrMC, "UtrMC/F");
  T->Branch("UcolMC", &UcolMC, "UcolMC/F");


//// puppi


   Float_t Ut_puppi;
   Float_t Utr_puppi;
   Float_t Ucol_puppi;
  T->Branch("Ut_puppi", &Ut_puppi, "Ut_puppi/F");
  T->Branch("Utr_puppi", &Utr_puppi, "Utr_puppi/F");
  T->Branch("Ucol_puppi", &Ucol_puppi, "Ucol_puppi/F");


   Float_t Ut_puppi_JetEnUp;
   Float_t Utr_puppi_JetEnUp;
   Float_t Ucol_puppi_JetEnUp;
  T->Branch("Ut_puppi_JetEnUp", &Ut_puppi_JetEnUp, "Ut_puppi_JetEnUp/F");
  T->Branch("Utr_puppi_JetEnUp", &Utr_puppi_JetEnUp, "Utr_puppi_JetEnUp/F");
  T->Branch("Ucol_puppi_JetEnUp", &Ucol_puppi_JetEnUp, "Ucol_puppi_JetEnUp/F");

   Float_t Ut_puppi_JetEnDown;
   Float_t Utr_puppi_JetEnDown;
   Float_t Ucol_puppi_JetEnDown;
  T->Branch("Ut_puppi_JetEnDown", &Ut_puppi_JetEnDown, "Ut_puppi_JetEnDown/F");
  T->Branch("Utr_puppi_JetEnDown", &Utr_puppi_JetEnDown, "Utr_puppi_JetEnDown/F");
  T->Branch("Ucol_puppi_JetEnDown", &Ucol_puppi_JetEnDown, "Ucol_puppi_JetEnDown/F");

   Float_t Ut_puppi_JetResDown;
   Float_t Utr_puppi_JetResDown;
   Float_t Ucol_puppi_JetResDown;
  T->Branch("Ut_puppi_JetResDown", &Ut_puppi_JetResDown, "Ut_puppi_JetResDown/F");
  T->Branch("Utr_puppi_JetResDown", &Utr_puppi_JetResDown, "Utr_puppi_JetResDown/F");
  T->Branch("Ucol_puppi_JetResDown", &Ucol_puppi_JetResDown, "Ucol_puppi_JetResDown/F");

   Float_t Ut_puppi_JetResUp;
   Float_t Utr_puppi_JetResUp;
   Float_t Ucol_puppi_JetResUp;
  T->Branch("Ut_puppi_JetResUp", &Ut_puppi_JetResUp, "Ut_puppi_JetResUp/F");
  T->Branch("Utr_puppi_JetResUp", &Utr_puppi_JetResUp, "Utr_puppi_JetResUp/F");
  T->Branch("Ucol_puppi_JetResUp", &Ucol_puppi_JetResUp, "Ucol_puppi_JetResUp/F");

   Float_t Ut_puppi_UnclusteredEnUp;
   Float_t Utr_puppi_UnclusteredEnUp;
   Float_t Ucol_puppi_UnclusteredEnUp;
  T->Branch("Ut_puppi_UnclusteredEnUp", &Ut_puppi_UnclusteredEnUp, "Ut_puppi_UnclusteredEnUp/F");
  T->Branch("Utr_puppi_UnclusteredEnUp", &Utr_puppi_UnclusteredEnUp, "Utr_puppi_UnclusteredEnUp/F");
  T->Branch("Ucol_puppi_UnclusteredEnUp", &Ucol_puppi_UnclusteredEnUp, "Ucol_puppi_UnclusteredEnUp/F");

   Float_t Ut_puppi_UnclusteredEnDown;
   Float_t Utr_puppi_UnclusteredEnDown;
   Float_t Ucol_puppi_UnclusteredEnDown;
  T->Branch("Ut_puppi_UnclusteredEnDown", &Ut_puppi_UnclusteredEnDown, "Ut_puppi_UnclusteredEnDown/F");
  T->Branch("Utr_puppi_UnclusteredEnDown", &Utr_puppi_UnclusteredEnDown, "Utr_puppi_UnclusteredEnDown/F");
  T->Branch("Ucol_puppi_UnclusteredEnDown", &Ucol_puppi_UnclusteredEnDown, "Ucol_puppi_UnclusteredEnDown/F");




/////smeared


   Float_t Ut_smeared;
   Float_t Utr_smeared;
   Float_t Ucol_smeared;
  T->Branch("Ut_smeared", &Ut_smeared, "Ut_smeared/F");
  T->Branch("Utr_smeared", &Utr_smeared, "Utr_smeared/F");
  T->Branch("Ucol_smeared", &Ucol_smeared, "Ucol_smeared/F");


   Float_t Ut_smeared_JetEnUp;
   Float_t Utr_smeared_JetEnUp;
   Float_t Ucol_smeared_JetEnUp;
  T->Branch("Ut_smeared_JetEnUp", &Ut_smeared_JetEnUp, "Ut_smeared_JetEnUp/F");
  T->Branch("Utr_smeared_JetEnUp", &Utr_smeared_JetEnUp, "Utr_smeared_JetEnUp/F");
  T->Branch("Ucol_smeared_JetEnUp", &Ucol_smeared_JetEnUp, "Ucol_smeared_JetEnUp/F");

   Float_t Ut_smeared_JetEnDown;
   Float_t Utr_smeared_JetEnDown;
   Float_t Ucol_smeared_JetEnDown;
  T->Branch("Ut_smeared_JetEnDown", &Ut_smeared_JetEnDown, "Ut_smeared_JetEnDown/F");
  T->Branch("Utr_smeared_JetEnDown", &Utr_smeared_JetEnDown, "Utr_smeared_JetEnDown/F");
  T->Branch("Ucol_smeared_JetEnDown", &Ucol_smeared_JetEnDown, "Ucol_smeared_JetEnDown/F");

   Float_t Ut_smeared_JetResDown;
   Float_t Utr_smeared_JetResDown;
   Float_t Ucol_smeared_JetResDown;
  T->Branch("Ut_smeared_JetResDown", &Ut_smeared_JetResDown, "Ut_smeared_JetResDown/F");
  T->Branch("Utr_smeared_JetResDown", &Utr_smeared_JetResDown, "Utr_smeared_JetResDown/F");
  T->Branch("Ucol_smeared_JetResDown", &Ucol_smeared_JetResDown, "Ucol_smeared_JetResDown/F");

   Float_t Ut_smeared_JetResUp;
   Float_t Utr_smeared_JetResUp;
   Float_t Ucol_smeared_JetResUp;
  T->Branch("Ut_smeared_JetResUp", &Ut_smeared_JetResUp, "Ut_smeared_JetResUp/F");
  T->Branch("Utr_smeared_JetResUp", &Utr_smeared_JetResUp, "Utr_smeared_JetResUp/F");
  T->Branch("Ucol_smeared_JetResUp", &Ucol_smeared_JetResUp, "Ucol_smeared_JetResUp/F");

   Float_t Ut_smeared_UnclusteredEnUp;
   Float_t Utr_smeared_UnclusteredEnUp;
   Float_t Ucol_smeared_UnclusteredEnUp;
  T->Branch("Ut_smeared_UnclusteredEnUp", &Ut_smeared_UnclusteredEnUp, "Ut_smeared_UnclusteredEnUp/F");
  T->Branch("Utr_smeared_UnclusteredEnUp", &Utr_smeared_UnclusteredEnUp, "Utr_smeared_UnclusteredEnUp/F");
  T->Branch("Ucol_smeared_UnclusteredEnUp", &Ucol_smeared_UnclusteredEnUp, "Ucol_smeared_UnclusteredEnUp/F");

   Float_t Ut_smeared_UnclusteredEnDown;
   Float_t Utr_smeared_UnclusteredEnDown;
   Float_t Ucol_smeared_UnclusteredEnDown;
  T->Branch("Ut_smeared_UnclusteredEnDown", &Ut_smeared_UnclusteredEnDown, "Ut_smeared_UnclusteredEnDown/F");
  T->Branch("Utr_smeared_UnclusteredEnDown", &Utr_smeared_UnclusteredEnDown, "Utr_smeared_UnclusteredEnDown/F");
  T->Branch("Ucol_smeared_UnclusteredEnDown", &Ucol_smeared_UnclusteredEnDown, "Ucol_smeared_UnclusteredEnDown/F");

   Float_t PU_weightRunF;
   Float_t PU_weightRunBCDE;
  T->Branch("PU_weightRunF", &PU_weightRunF, "PU_weightRunF/F");
  T->Branch("PU_weightRunBCDE", &PU_weightRunBCDE, "PU_weightRunBCDE/F");


////////////////////////////////////// MET performance end declaration
//}

  SetupHists(CutNumb); 
  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
  //if (nTotalFiles>50) nTotalFiles=50;
  //nTotalFiles = 10;
 
for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
//////////// for SUSY!!!
//if (SUSY){
//if (iF+1 != nTotalFiles) continue;}
    TFile * file_ = TFile::Open(TString(filen));

/*    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
    int NE = (int)histoInputEvents->GetEntries();
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    std::cout << "      number of input events         = " << NE << std::endl;
*/


bool WithInit = true;
if (SUSY) WithInit=false;

if (WithInit) cout << "With initroottree"<<endl;
if (!WithInit) cout << "Without initroottree"<<endl;


    TTree * _inittree = NULL;
if (!WithInit)  _inittree = (TTree*)file_->Get(TString(ntupleName));
if (WithInit)  _inittree = (TTree*)file_->Get(TString(initNtupleName));

    if (_inittree==NULL) continue;
    Float_t genweight;
    if (!isData)
      _inittree->SetBranchAddress("genweight",&genweight);
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
    std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; ++iEntry) {
      _inittree->GetEntry(iEntry);
      if (isData)
	histWeightsH->Fill(0.,1.);
      //else
      //histWeightsH->Fill(0.,genweight);
    }


    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      number  of entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);




	if (!isData && !WithInit)
	//if (!isData)
		{    
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
			{
			analysisTree.GetEntry(iEntry);
		/*	if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
			histWeightsH->Fill(0.,analysisTree.genweight);
			}
		}
  	float genweights=1.;

	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	genweightsTree->SetBranchAddress("genweight",&genweights);

    if(!isData && WithInit) 
      {

	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) { 
	  genweightsTree->GetEntry(iEntryInit);
	/*		if (SUSY)
			{
			if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
			&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
			}*/
	  histWeightsH->Fill(0.,genweights);
	}
    
      }




    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 
      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
/*	if (SUSY)
		{
		if (!(SusyMotherMassF < (analysisTree.SusyMotherMass+1) && SusyMotherMassF > (analysisTree.SusyMotherMass - 1) 
		&& SusyLSPMassF <(analysisTree.SusyLSPMass + 1) && SusyLSPMassF > (analysisTree.SusyLSPMass - 1))) continue;
		}*/
      nEvents++;

      //std::cout << "      number of entries in Tree = " << numberOfEntries <<" starting weight "<<weight<< std::endl;

      if (nEvents%50000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;  


      //isData= false;
      bool lumi=false;
      bool CutBasedTauId = false;

      float topPt = 0;
      float antitopPt = 0;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
            PU_weightRunF = 1.;
      PU_weightRunBCDE = 1.;

      
      gen_weight = 1.;
      trig_weight = 1.;

      bool isW = false;
      bool isDY = false;
      bool isZTT = false;
      bool isZMM = false;
      bool isZEE = false;
      bool isTOP = false;
      if (!isData &&  string::npos != filen.find("JetsToLNu") ) isW=true;
      if (!isData &&  string::npos != filen.find("JetsToLL_M") )  isDY=true;
      if (!isData &&  string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8") ) isTOP=true;

      float nuPx = 0;
      float nuPy = 0;
      float nuPz = 0;
      float nuPt = 0;
      float nuPhi = 0;
      
      float nuPx_msv = 0;
      float nuPy_msv = 0;
      float nuPz_msv = 0;
      float nuPt_msv = 0;
      float nuPhi_msv = 0;
      
      float lepPx = 0;
      float lepPy = 0;
      float lepPz = 0;
      float bosonPx = 0;
      float bosonPy = 0;
      float bosonPz = 0;
      float bosonPt = 0;
      float bosonEta = 0;
      float bosonMass = -1;
	  
    /*  bool isZfound = false;
      bool isWfound = false;
      bool isHfound = false;
      bool isGSfound = false;*/
      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      std::vector<TLorentzVector> tauNeutrinos; tauNeutrinos.clear();

      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector tauNeutrinosLV;  tauNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);


	  if (!isData ) {
		genweightsTree->GetEntry(iEntry);
	    weight *= genweights;
	    gen_weight *=genweights;
	   // std::cout <<"analysisTree.genweight "<< float(analysisTree.genweight) << std::endl;
	  lumi=true;
	  }
	//cout << "Data!?!?!?!"  << isData <<endl;
      if (isData)  {
	XSec = 1.;
	histRuns->Fill(analysisTree.event_run);
	int n=analysisTree.event_run;
	int lum = analysisTree.event_luminosityblock;

	std::string num = std::to_string(n);
	std::string lnum = std::to_string(lum);
	for(const auto& a : periods)
	  {

	    if ( num.c_str() ==  a.name ) {
	      //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	      //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {

		//	cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	       //   std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }

	  }
	if (!lumi) continue;
	//if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
      }


      if (analysisTree.event_run<RunMin)
	RunMin = analysisTree.event_run;

      if (analysisTree.event_run>RunMax)
	RunMax = analysisTree.event_run;

      //std::cout << " Run : " << analysisTree.event_run << std::endl;

      bool isNewRun = true;
      if (allRuns.size()>0) {
	for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	  if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }

      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);

      if (!lumi) continue;
//cout << "lumi"<< endl;
	std::vector<TString> metFlags; metFlags.clear();
     //////////////MET filters flag

	 	 metFlags.push_back("Flag_goodVertices");
	 	 metFlags.push_back("Flag_globalTightHalo2016Filter");
	 	 metFlags.push_back("Flag_HBHENoiseFilter");
	 	 metFlags.push_back("Flag_HBHENoiseIsoFilter");
	 	 metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	 	 metFlags.push_back("Flag_BadPFMuonFilter");
	 	 metFlags.push_back("Flag_BadChargedCandidateFilter");
	 	 metFlags.push_back("Flag_eeBadScFilter");
	 	 metFlags.push_back("Flag_ecalBadCalibFilter");



	bool METflag = metFiltersPasses2(analysisTree, metFlags);
	met_flag = METflag;
	//if (!METflag && isData) continue;


	//cout << "nachali"<< endl;
      if (!isData) 
	{
	  if (applyPUreweighting)	 {
	    puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
		//puweight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
	    //puweight = float((PU_data->GetBinContent(analysisTree.primvertex_count)/PU_data->GetSumOfWeights())/(PU_mc->GetBinContent(analysisTree.primvertex_count)/PU_mc->GetSumOfWeights()));
	    if (puweight !=puweight) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}
	    if (puweight >10) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}

	    weight *=puweight; 

	    pu_weight = puweight;
	    
		puweight = float((PU_RunBCDE->GetBinContent(analysisTree.primvertex_count)/PU_RunBCDE->GetSumOfWeights())/(PU_mc2->GetBinContent(analysisTree.primvertex_count)/PU_mc2->GetSumOfWeights()));
	    if (puweight !=puweight) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}
	    if (puweight >10) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}
		PU_weightRunBCDE = puweight;
	    
	    
	    puweight = float((PU_RunF->GetBinContent(analysisTree.primvertex_count)/PU_RunF->GetSumOfWeights())/(PU_mc2->GetBinContent(analysisTree.primvertex_count)/PU_mc2->GetSumOfWeights()));
	    if (puweight !=puweight) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}
	    if (puweight >10) { puweight=1; cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;}
		PU_weightRunF = puweight;
	    
	//cout << PU_data->GetBinLowEdge(analysisTree.primvertex_count+1)<< "   vs    "<< analysisTree.primvertex_count<< endl;
	  }
	}


      bool trigAccept = false;

      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;



      if (1){
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      //  std::cout << "nfiltres = " << nfilters << std::endl;
      for (unsigned int i=0; i<nfilters; ++i) {
	//	std::cout << "HLT Filter : " << i << " = " << analysisTree.run_hltfilters->at(i) << std::endl;
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==MainTrigger) {
	  nMainTrigger = i;
	  isMainTrigger = true;
	}

      }
	}//if isData check for filters


      if (!isMainTrigger) {
	std::cout << "HLT filter for Mu Trigger " << MainTrigger << " not found" << std::endl;
	return(-1);
      }
      /////now clear the Mu.El.Jets again to fill them again after cleaning
	
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	//if (applyMuonId && !analysisTree.muon_isICHEP[im]) continue;
        if ( fabs(analysisTree.muon_charge[im]) != 1) continue;
	muons.push_back((int)im);

      }
      if (muons.size()==0) continue;


      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float isoMuMin  = 1e+10;
      float isoTauMin = 1.; 
      float isoTau = 1.; 
      float ptMu = 0;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;
      
	bool isLegMatch = false;


	for (unsigned int im=0; im<muons.size(); ++im) {
	isLegMatch = false;
	//	bool isMuonTauMuonLegMatch = false;
	//	bool isMuonTauOverlapMuonMatch = false;
	unsigned int mIndex  = muons.at(im);
	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[mIndex];
	float photonIsoMu = analysisTree.muon_photonIso[mIndex];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[mIndex];
	float puIsoMu = analysisTree.muon_puIso[mIndex];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[mIndex];
	  photonIsoMu = analysisTree.muon_r04_sumPhotonEt[mIndex];
	  chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[mIndex];
	  puIsoMu = analysisTree.muon_r04_sumPUPt[mIndex];
	}
	double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	double neutralIsoMu = max(double(0),neutralIsoMuN); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[mIndex];

	//if (relIsoMu < isoMuMin) {isoMuMin = relIsoMu; mu_index = mIndex;}

	if (1)
	{	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {

	  	if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      	&& analysisTree.muon_pt[mIndex]>ptMuonCut&&
	      	analysisTree.trigobject_pt[iT]>SingleMuonTriggerPtCut) { // IsoMu Leg
	    	float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    	if (dRtrig<deltaRTrigMatch) 
	      	isLegMatch = true;
		}
	  

	  	  }
		}
	

	if (!isLegMatch) continue;


          if ((int)mIndex!=(int)mu_index) {
            if (relIsoMu==isoMuMin) {
              if (analysisTree.muon_pt[mIndex]>ptMu) {
                isoMuMin  = relIsoMu;
                ptMu = analysisTree.muon_pt[mIndex];
                mu_index =(int)mIndex;
              }
            }

            else if (relIsoMu<isoMuMin) {
              isoMuMin  = relIsoMu;
              ptMu = analysisTree.muon_pt[mIndex];
              mu_index =(int)mIndex;
            }
          }


    }



      if ((int)mu_index<0) continue;

	
      mu_relIso[0]=isoMuMin;
      



  bool          dilepton_veto=false;
  bool          extraelec_veto=false;
  bool          extramuon_veto=false;

  event_secondLeptonVeto = false;
  event_thirdLeptonVeto = false;

      // looking for extra electron
      bool foundExtraElectron = false;
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie];
	if (!electronMvaId&&applyVetoElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyVetoElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyVetoElectronId) continue;
	float neutralHadIsoEle = analysisTree.electron_neutralHadIso[ie];
	float photonIsoEle = analysisTree.electron_photonIso[ie];
	float chargedHadIsoEle = analysisTree.electron_chargedHadIso[ie];
	float puIsoEle = analysisTree.electron_puIso[ie];
	if (isIsoR03) {
	  neutralHadIsoEle = analysisTree.electron_r03_sumNeutralHadronEt[ie];
	  photonIsoEle = analysisTree.electron_r03_sumPhotonEt[ie];
	  chargedHadIsoEle = analysisTree.electron_r03_sumChargedHadronPt[ie];
	  puIsoEle = analysisTree.electron_r03_sumPUPt[ie];
	}
	double neutralIsoEleN = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	double neutralIsoEle = max(double(0),neutralIsoEleN); 
	float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];	
	if (relIsoEle>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      // looking for extra muon's (dimuon veto)
      bool foundExtraMuon = false;
      vector<int> mu_dimuons; mu_dimuons.clear(); 
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if ((int)im==(int)mu_index) continue;

	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	float photonIsoMu = analysisTree.muon_photonIso[im];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	float puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r04_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r04_sumPUPt[im];
	}
	double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	double neutralIsoMu = max(double(0),neutralIsoMuN); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];

	if (analysisTree.muon_pt[im]>ptDilepMuonCut&&
	    fabs(analysisTree.muon_eta[im])<etaDilepMuonCut&&
	    analysisTree.muon_isGlobal[im]&&
	    analysisTree.muon_isTracker[im]&&
	    analysisTree.muon_isPF[im]&&
	    fabs(analysisTree.muon_dxy[im])<dxyDilepMuonCut&&
	    fabs(analysisTree.muon_dz[im])<dzDilepMuonCut&&
	    relIsoMu<isoDilepMuonCut &&
	    fabs(analysisTree.muon_charge[im]) ==1)  
	{
	    float dRmuons = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				   analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);
	    if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[im]*analysisTree.muon_charge[mu_index]<0)) dilepton_veto = true;
 	  }
	 // mu_dimuons.push_back(im);

	if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	//if (applyVetoMuonId && !analysisTree.muon_isICHEP[im]) continue;
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;
    
/*
      if (mu_dimuons.size()>1) {
	for (unsigned int i1=0; i1<mu_dimuons.size()-1; ++i1) {
	 1unsigned int indx1 = mu_dimuons[i1];
	  for (unsigned int i2=i1+1; i2<mu_dimuons.size(); ++i2 ) {
	    unsigned int indx2 = mu_dimuons[i2];
	    float dRmuons = deltaR(analysisTree.muon_eta[indx1],analysisTree.muon_phi[indx1],
				   analysisTree.muon_eta[indx2],analysisTree.muon_phi[indx2]);
	    if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[indx1]*analysisTree.muon_charge[indx2]<0)) dilepton_veto = true;
 	  }
	}
      }
     */
      //      cout << analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] << endl;
      //      cout << analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index] << endl;



      // applying inclusive selection

   	event_secondLeptonVeto = dilepton_veto;
//	if (dilepton_veto)  continue;


	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
//for now
//extraelec_veto = false;///////////////????????????

      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;




///////////////Trigger weight 
      double ptMu1 = (double)analysisTree.muon_pt[mu_index];
      double etaMu1 = (double)analysisTree.muon_eta[mu_index];
///cout <<"ptMu   " << ptMu1<<endl;
///cout <<"etaMu   " << etaMu1 <<endl;
      float trigweight=1.;

      float EffFromData = (float)SF_muonTrigger->get_EfficiencyData(double(ptMu1),double(etaMu1));
      float Mu17EffMC   = (float)SF_muonTrigger->get_EfficiencyMC(double(ptMu1),double(etaMu1));
	
      bool Signal = true;

	//if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) Signal=true;
      if (!isData) {
	if (Mu17EffMC>1e-6)
	  trigweight = EffFromData / Mu17EffMC;
	//if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) trigweight = EffFromData;
	weight *= trigweight;
	trig_weight = trigweight;
	//	cout<<" Trigger weight "<<trigweight<<endl;
      }
	//if (!isData) trigweight = EffFromData;
	//weight *= trigweight;
	//trig_weight = trigweight;



	///LSF 

      if (!isData && applyLeptonSF) {

	//leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)	
	double IdIsoSF_mu1 = SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
	//double IdIsoSF_mu1 =0;
	//IdIsoSF_mu1=hSF_BCDEF->GetBinContent(hSF_BCDEF->GetXaxis()->FindBin(abs(etaMu1)),hSF_BCDEF->GetYaxis()->FindBin(ptMu1));
	//if (ptMu1>120) IdIsoSF_mu1=hSF_BCDEF->GetBinContent(hSF_BCDEF->GetXaxis()->FindBin(abs(etaMu1)),hSF_BCDEF->GetYaxis()->FindBin(111));
//cout <<"SF!!!!!!!!!!!!!!!!!!!!!!!!!!"<< IdIsoSF_mu1<< endl;
//cout <<"IdIsoSF_mu1   " << IdIsoSF_mu1<<endl;
	MuSF_IdIso_Mu1H->Fill(IdIsoSF_mu1);
	weight *= IdIsoSF_mu1;
	LSF_weight = IdIsoSF_mu1;
      }



	//double HIP_SF1 = SF_HIP->get_ScaleFactor(ptMu1, etaMu1);
	//cout<<"  "<<ptMu1<<"  "<<etaMu1<<"  "<<HIP_SF1<<endl;

      bool isTauMatched = false;
      bool isGenLeptonMatched = false;


      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      mu_count= (int)analysisTree.muon_count;
      //cout<<" here ============================> "<<iEntry<<"  "<<mu_count<<"  "<<(int)analysisTree.muon_count<<"  "<<analysisTree.muon_count<<endl;
      //for (unsigned int im=0;im<analysisTree.muon_count; ++im){
	int im = mu_index;
	mu_px[0]=analysisTree.muon_px[im];
	mu_py[0]=analysisTree.muon_py[im];
	mu_pz[0]=analysisTree.muon_pz[im];
	mu_eta[0]=analysisTree.muon_eta[im];
	mu_pt[0]=analysisTree.muon_pt[im];
	mu_phi[0]=analysisTree.muon_phi[im];
	mu_charge[0]=analysisTree.muon_charge[im];
	mu_dxy[0]=analysisTree.muon_dxy[im];
	mu_dz[0]=analysisTree.muon_dz[im];
	mu_dxyerr[0]=analysisTree.muon_dxyerr[im];
	mu_dzerr[0]=analysisTree.muon_dzerr[im];

        mu_neutralHadIso[0] = analysisTree.muon_r04_sumNeutralHadronEt[im];
        mu_photonIso[0] = analysisTree.muon_r04_sumPhotonEt[im];
        mu_chargedHadIso[0] = analysisTree.muon_r04_sumChargedHadronPt[im];
        mu_puIso[0] = analysisTree.muon_r04_sumPUPt[im];
 
        double neutralIso = mu_neutralHadIso[im] + mu_photonIso[im] - 0.5*mu_puIso[im];
        neutralIso = max(double(0),neutralIso);
	mu_neutralIso[0] = neutralIso;
        mu_absIsoMu[0] = mu_chargedHadIso[im] + neutralIso;
	mu_relIsoMu[0]  = mu_absIsoMu[im]/mu_pt[im] ;
   
     	//}


/*
      el_count=(int)analysisTree.electron_count;
      for (unsigned int ie=0;ie<analysisTree.electron_count; ++ie){
	el_px[ie]=analysisTree.electron_px[ie];
	el_py[ie]=analysisTree.electron_py[ie];
	el_pz[ie]=analysisTree.electron_pz[ie];
	el_eta[ie]=analysisTree.electron_eta[ie];
	el_pt[ie]=analysisTree.electron_pt[ie];
	el_phi[ie]=analysisTree.electron_phi[ie];
	el_charge[ie]=analysisTree.electron_charge[ie];
	el_dxy[ie]=analysisTree.electron_dxy[ie];
	el_dz[ie]=analysisTree.electron_dz[ie];
	el_dxyerr[ie]=analysisTree.electron_dxyerr[ie];
	el_dzerr[ie]=analysisTree.electron_dzerr[ie];

        el_neutralHadIso[ie] = analysisTree.electron_r03_sumNeutralHadronEt[ie];
        el_photonIso[ie] = analysisTree.electron_r03_sumPhotonEt[ie];
        el_chargedHadIso[ie] = analysisTree.electron_r03_sumChargedHadronPt[ie];
        el_puIso[ie] = analysisTree.electron_r03_sumPUPt[ie];
 
        double neutralIso = el_neutralHadIso[ie] + el_photonIso[ie] - 0.5*el_puIso[ie];
        neutralIso = max(double(0),neutralIso);
        el_neutralIso[ie] = neutralIso ;
        el_absIsoEl[ie] = el_chargedHadIso[ie] + mu_neutralIso[ie];
	el_relIsoEl[ie]  = el_absIsoEl[ie]/el_pt[ie] ;

      }

				
      ta_count=(int)analysisTree.tau_count;
      for (unsigned int it=0;it<analysisTree.tau_count; ++it){
	ta_mass[it]=analysisTree.tau_mass[it];
	ta_px[it]=analysisTree.tau_px[it];
	ta_py[it]=analysisTree.tau_py[it];
	ta_pz[it]=analysisTree.tau_pz[it];
	ta_eta[it]=analysisTree.tau_eta[it];
	ta_pt[it]=analysisTree.tau_pt[it];
	ta_phi[it]=analysisTree.tau_phi[it];
	ta_charge[it]=analysisTree.tau_charge[it];
	ta_dxy[it]=analysisTree.tau_dxy[it];
	ta_dz[it]=analysisTree.tau_dz[it];
     	ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[it];
     	ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[it];
     	ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[it];



      }*/
      jet_count=(int)analysisTree.pfjet_count;
      for (unsigned int jj=0;jj<analysisTree.pfjet_count; ++jj){

	jet_e[jj] = analysisTree.pfjet_e[jj];
	jet_px[jj] = analysisTree.pfjet_px[jj];
	jet_py[jj] = analysisTree.pfjet_py[jj];
	jet_pz[jj] = analysisTree.pfjet_pz[jj];
	jet_pt[jj] = analysisTree.pfjet_pt[jj];
	jet_eta[jj] = analysisTree.pfjet_eta[jj];
	jet_phi[jj] = analysisTree.pfjet_phi[jj];
	jet_flavour[jj] = analysisTree.pfjet_flavour[jj];
	jet_btag[jj] = analysisTree.pfjet_btag[jj][0];
      }





  ////////jets cleaning 
      TLorentzVector leptonsV, muonJ, jetsLV;


      //JetsV.SetPxPyPzE(analysisTree.pfjet_px[ij], analysisTree.pfjet_py[ij], analysisTree.pfjet_pz[ij], analysisTree.pfjet_e[ij]);
      float jetEta = 2.4;
      float DRmax = 0.5;
      bool dRmuJet = false;
      bool dRtauJet = false;
      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;

	int counter_cleaned_jets = 0;


      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	if (fabs(analysisTree.pfjet_pt[jet])<30.) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];


	bool isPFJetId = false ; 
      	bool btagged= false;
	isPFJetId =looseJetiD(analysisTree,jet);
	//isPFJetId =tightLepVetoJetiD(analysisTree,jet);

	if (!isPFJetId) continue;
	bool cleanedJet = true;

	double Dr=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Dr  < DRmax)  cleanedJet=false;

	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

	if (analysisTree.pfjet_btag[jet][0]  > bTag) btagged = true;
	

	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);

	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;

	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
	      tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
	    }
	    
	    if (tageff<1e-5)      tageff = 1e-5;
	    if (tageff>0.99999)   tageff = 0.99999;
	    rand.SetSeed((int)((jetEta+5)*100000));
	    double rannum = rand.Rndm();
	    
	    if (jet_scalefactor<1 && btagged) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		btagged = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !btagged) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) { 
		btagged = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }
	  } //is Data

	  if (btagged && cleanedJet) bjets.push_back(jet);
	}


	if (cleanedJet){
		
	//	cout<<"  will push to save now cleaned jet  "<<(int)jet<<"  for counter_cleaned_jet "<<(int)counter_cleaned_jets<<" event "<<iEntry<<endl;

	jets.push_back((int)jet);
	jets_cleaned[counter_cleaned_jets]=(int)jet;
	jet_jecUn[counter_cleaned_jets] = analysisTree.pfjet_jecUncertainty[jet];
	counter_cleaned_jets++;
	}


      }///loop in all jets

      njets = jets.size();
      jet_count = jets.size();
      //njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      //nbtag_nocleaned = bjets_nocleaned.size();

	if (nbtag > 0.5) continue;
	//if (njets < 0.5) continue;



     

/////////////////// Recoil corrections

      int njetsforrecoil = njets;

      float pfmet_corr_x = analysisTree.pfmetcorr_ex;
      float pfmet_corr_y = analysisTree.pfmetcorr_ey;
      float met_x = analysisTree.pfmetcorr_ex;
      float met_y = analysisTree.pfmetcorr_ey;
/*
      if ((isW||isDY) && !isData) {

	  recoilMetCorrector.CorrectByMeanResolution(analysisTree.pfmetcorr_ex,analysisTree.pfmetcorr_ey,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
 
        met_x = pfmet_corr_x;
        met_y = pfmet_corr_y;
 
      // MEt related systematic uncertainties
      int bkgdType = 0;
      if (isDY||isW)
	bkgdType = MEtSys::ProcessType::BOSON;
      else if (isTOP)
	bkgdType = MEtSys::ProcessType::TOP;
      else 
	bkgdType = MEtSys::ProcessType::EWK; 

      float met_scaleUp_x   = met_x;
      float met_scaleUp_y   = met_y;
      float met_scaleDown_x = met_x;
      float met_scaleDown_y = met_y;
      float met_resoUp_x    = met_x;
      float met_resoUp_y    = met_y;
      float met_resoDown_x  = met_x;
      float met_resoDown_y  = met_y;

	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Up,
			   met_scaleUp_x,met_scaleUp_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Response,MEtSys::SysShift::Down,
			   met_scaleDown_x,met_scaleDown_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Up,
			   met_resoUp_x,met_resoUp_y);
	metSys.ApplyMEtSys(met_x,met_y,
			   bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,bkgdType,
			   MEtSys::SysType::Resolution,MEtSys::SysShift::Down,
			   met_resoDown_x,met_resoDown_y);

      
      met_scaleUp = TMath::Sqrt(met_scaleUp_x*met_scaleUp_x+
				   met_scaleUp_y*met_scaleUp_y);
      metphi_scaleUp = TMath::ATan2(met_scaleUp_y,met_scaleUp_x);
      
      met_scaleDown = TMath::Sqrt(met_scaleDown_x*met_scaleDown_x+
				     met_scaleDown_y*met_scaleDown_y);
      metphi_scaleDown = TMath::ATan2(met_scaleDown_y,met_scaleDown_x);
      
      met_resoUp = TMath::Sqrt(met_resoUp_x*met_resoUp_x+
				  met_resoUp_y*met_resoUp_y);
      metphi_resoUp = TMath::ATan2(met_resoUp_y,met_resoUp_x);
      
      met_resoDown = TMath::Sqrt(met_resoDown_x*met_resoDown_x+
				    met_resoDown_y*met_resoDown_y);
      metphi_resoDown = TMath::ATan2(met_resoDown_y,met_resoDown_x);
 
 
 
	

      met_ex_recoil = pfmet_corr_x;
      met_ey_recoil = pfmet_corr_y;

      }//if isW, isDY !isData/
*/
      met_ex = analysisTree.pfmetcorr_ex;
      met_ey = analysisTree.pfmetcorr_ey;


      met_ez = 0;//analysisTree.pfmet_ez;
      //met_pt = analysisTree.pfmet_pt;
      met_pt = TMath::Sqrt(pfmet_corr_x*pfmet_corr_x+pfmet_corr_y*pfmet_corr_y);
      //met_phi = analysisTree.pfmet_phi;
      met_phi = TMath::ATan2(pfmet_corr_y,pfmet_corr_x);
      
      
      genmet = TMath::Sqrt(analysisTree.genmet_ex*analysisTree.genmet_ex+analysisTree.genmet_ey*analysisTree.genmet_ey);

     met_ex_JetEnUp = analysisTree.pfmetcorr_ex_JetEnUp;
     met_ey_JetEnUp = analysisTree.pfmetcorr_ey_JetEnUp;

     met_ex_JetEnDown = analysisTree.pfmetcorr_ex_JetEnDown;
     met_ey_JetEnDown = analysisTree.pfmetcorr_ey_JetEnDown;

     met_ex_UnclusteredEnUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
     met_ey_UnclusteredEnUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
   
     met_ex_UnclusteredEnDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
     met_ey_UnclusteredEnDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;

      met_pt_JetEnUp = TMath::Sqrt(met_ex_JetEnUp*met_ex_JetEnUp+met_ey_JetEnUp*met_ey_JetEnUp);
      met_pt_JetEnDown = TMath::Sqrt(met_ex_JetEnDown*met_ex_JetEnDown+met_ey_JetEnDown*met_ey_JetEnDown);
      met_pt_UnclusteredEnUp = TMath::Sqrt(met_ex_UnclusteredEnUp*met_ex_UnclusteredEnUp+met_ey_UnclusteredEnUp*met_ey_UnclusteredEnUp);
      met_pt_UnclusteredEnDown = TMath::Sqrt(met_ex_UnclusteredEnDown*met_ex_UnclusteredEnDown+met_ey_UnclusteredEnDown*met_ey_UnclusteredEnDown);




      puppi_ex = analysisTree.puppimet_ex;
      puppi_ey = analysisTree.puppimet_ey;
      puppi_ez = 0;//analysisTree.pfmet_ez;
      puppi_pt = TMath::Sqrt(puppi_ex*puppi_ex+puppi_ey*puppi_ey);
      puppi_phi = TMath::ATan2(puppi_ey,puppi_ex);

     puppi_ex_JetEnUp = analysisTree.puppimet_ex_JetEnUp;
     puppi_ey_JetEnUp = analysisTree.puppimet_ey_JetEnUp;
     puppi_ex_JetEnDown = analysisTree.puppimet_ex_JetEnDown;
     puppi_ey_JetEnDown = analysisTree.puppimet_ey_JetEnDown;
     puppi_ex_UnclusteredEnUp = analysisTree.puppimet_ex_UnclusteredEnUp;
     puppi_ey_UnclusteredEnUp = analysisTree.puppimet_ey_UnclusteredEnUp;
     puppi_ex_UnclusteredEnDown = analysisTree.puppimet_ex_UnclusteredEnDown;
     puppi_ey_UnclusteredEnDown = analysisTree.puppimet_ey_UnclusteredEnDown;


      puppi_pt_JetEnUp = TMath::Sqrt(puppi_ex_JetEnUp*puppi_ex_JetEnUp+puppi_ey_JetEnUp*puppi_ey_JetEnUp);
      puppi_pt_JetEnDown = TMath::Sqrt(puppi_ex_JetEnDown*puppi_ex_JetEnDown+puppi_ey_JetEnDown*puppi_ey_JetEnDown);
      puppi_pt_UnclusteredEnUp = TMath::Sqrt(puppi_ex_UnclusteredEnUp*puppi_ex_UnclusteredEnUp+puppi_ey_UnclusteredEnUp*puppi_ey_UnclusteredEnUp);
      puppi_pt_UnclusteredEnDown = TMath::Sqrt(puppi_ex_UnclusteredEnDown*puppi_ex_UnclusteredEnDown+puppi_ey_UnclusteredEnDown*puppi_ey_UnclusteredEnDown);


//if (!isData)
//	{
	

     met_ex_JetResUp = analysisTree.pfmetcorr_ex_JetResUp;
     met_ey_JetResUp = analysisTree.pfmetcorr_ey_JetResUp;

     met_ex_JetResDown = analysisTree.pfmetcorr_ex_JetResDown;
     met_ey_JetResDown = analysisTree.pfmetcorr_ey_JetResDown;
     puppi_ex_JetResUp = analysisTree.puppimet_ex_JetResUp;
     puppi_ey_JetResUp = analysisTree.puppimet_ey_JetResUp;
     puppi_ex_JetResDown = analysisTree.puppimet_ex_JetResDown;
     puppi_ey_JetResDown = analysisTree.puppimet_ey_JetResDown;



      met_ex_smeared = analysisTree.pfmetcorr_ex_smeared;
      met_ey_smeared = analysisTree.pfmetcorr_ey_smeared;


      met_ez_smeared = 0;//analysisTree.pfmet_ez;
      met_pt_smeared = TMath::Sqrt(met_ex_smeared*met_ex_smeared+met_ey_smeared*met_ey_smeared);
      met_phi_smeared = TMath::ATan2(met_ey_smeared,met_ex_smeared);

     met_ex_JetEnUp_smeared = analysisTree.pfmetcorr_ex_JetEnUp_smeared;
     met_ey_JetEnUp_smeared = analysisTree.pfmetcorr_ey_JetEnUp_smeared;
     met_ex_JetEnDown_smeared = analysisTree.pfmetcorr_ex_JetEnDown_smeared;
     met_ey_JetEnDown_smeared = analysisTree.pfmetcorr_ey_JetEnDown_smeared;
     met_ex_UnclusteredEnUp_smeared = analysisTree.pfmetcorr_ex_UnclusteredEnUp_smeared;
     met_ey_UnclusteredEnUp_smeared = analysisTree.pfmetcorr_ey_UnclusteredEnUp_smeared;
     met_ex_UnclusteredEnDown_smeared = analysisTree.pfmetcorr_ex_UnclusteredEnDown_smeared;
     met_ey_UnclusteredEnDown_smeared = analysisTree.pfmetcorr_ey_UnclusteredEnDown_smeared;
     met_ex_JetResUp_smeared = analysisTree.pfmetcorr_ex_JetResUp_smeared;
     met_ey_JetResUp_smeared = analysisTree.pfmetcorr_ey_JetResUp_smeared;
     met_ex_JetResDown_smeared = analysisTree.pfmetcorr_ex_JetResDown_smeared;
     met_ey_JetResDown_smeared = analysisTree.pfmetcorr_ey_JetResDown_smeared;


      met_pt_JetResUp = TMath::Sqrt(met_ex_JetResUp*met_ex_JetResUp+met_ey_JetResUp*met_ey_JetResUp);
      met_pt_JetResDown = TMath::Sqrt(met_ex_JetResDown*met_ex_JetResDown+met_ey_JetResDown*met_ey_JetResDown);
      puppi_pt_JetResUp = TMath::Sqrt(puppi_ex_JetResUp*puppi_ex_JetResUp+puppi_ey_JetResUp*puppi_ey_JetResUp);
      puppi_pt_JetResDown = TMath::Sqrt(puppi_ex_JetResDown*puppi_ex_JetResDown+puppi_ey_JetResDown*puppi_ey_JetResDown);


      met_pt_JetEnUp_smeared = TMath::Sqrt(met_ex_JetEnUp_smeared*met_ex_JetEnUp_smeared+met_ey_JetEnUp_smeared*met_ey_JetEnUp_smeared);
      met_pt_JetEnDown_smeared = TMath::Sqrt(met_ex_JetEnDown_smeared*met_ex_JetEnDown_smeared+met_ey_JetEnDown_smeared*met_ey_JetEnDown_smeared);
      met_pt_UnclusteredEnUp_smeared = TMath::Sqrt(met_ex_UnclusteredEnUp_smeared*met_ex_UnclusteredEnUp_smeared+met_ey_UnclusteredEnUp_smeared*met_ey_UnclusteredEnUp_smeared);
      met_pt_UnclusteredEnDown_smeared = TMath::Sqrt(met_ex_UnclusteredEnDown_smeared*met_ex_UnclusteredEnDown_smeared+met_ey_UnclusteredEnDown_smeared*met_ey_UnclusteredEnDown_smeared);
      met_pt_JetResUp_smeared = TMath::Sqrt(met_ex_JetResUp_smeared*met_ex_JetResUp_smeared+met_ey_JetResUp_smeared*met_ey_JetResUp_smeared);
      met_pt_JetResDown_smeared = TMath::Sqrt(met_ex_JetResDown_smeared*met_ex_JetResDown_smeared+met_ey_JetResDown_smeared*met_ey_JetResDown_smeared);

	



	
//	}





	float dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex,  met_ey );
	MT = TMath::Sqrt(2*mu_pt[0]*met_pt*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetEnUp,  met_ey_JetEnUp );
	MT_JetEnUp = TMath::Sqrt(2*mu_pt[0]*met_pt_JetEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetEnDown,  met_ey_JetEnDown );
	MT_JetEnDown = TMath::Sqrt(2*mu_pt[0]*met_pt_JetEnDown*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_UnclusteredEnUp,  met_ey_UnclusteredEnUp );
	MT_UnclusteredEnUp = TMath::Sqrt(2*mu_pt[0]*met_pt_UnclusteredEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_UnclusteredEnDown,  met_ey_UnclusteredEnDown );
	MT_UnclusteredEnDown = TMath::Sqrt(2*mu_pt[0]*met_pt_UnclusteredEnDown*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetResUp,  met_ey_JetResUp );
	MT_JetResUp = TMath::Sqrt(2*mu_pt[0]*met_pt_JetResUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetResDown,  met_ey_JetResDown );
	MT_JetResDown = TMath::Sqrt(2*mu_pt[0]*met_pt_JetResDown*(1-TMath::Cos(dPhi)));




	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex,  puppi_ey );
	MTpuppi = TMath::Sqrt(2*mu_pt[0]*puppi_pt*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_JetEnUp,  puppi_ey_JetEnUp );
	MTpuppi_JetEnUp = TMath::Sqrt(2*mu_pt[0]*puppi_pt_JetEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_JetEnDown,  puppi_ey_JetEnDown );
	MTpuppi_JetEnDown = TMath::Sqrt(2*mu_pt[0]*puppi_pt_JetEnDown*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_UnclusteredEnUp,  puppi_ey_UnclusteredEnUp );
	MTpuppi_UnclusteredEnUp = TMath::Sqrt(2*mu_pt[0]*puppi_pt_UnclusteredEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_UnclusteredEnDown,  puppi_ey_UnclusteredEnDown );
	MTpuppi_UnclusteredEnDown = TMath::Sqrt(2*mu_pt[0]*puppi_pt_UnclusteredEnDown*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_JetResUp,  puppi_ey_JetResUp );
	MTpuppi_JetResUp = TMath::Sqrt(2*mu_pt[0]*puppi_pt_JetResUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex_JetResDown,  puppi_ey_JetResDown );
	MTpuppi_JetResDown = TMath::Sqrt(2*mu_pt[0]*puppi_pt_JetResDown*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_smeared,  met_ey_smeared );
	MT_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetEnUp_smeared,  met_ey_JetEnUp_smeared );
	MT_JetEnUp_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_JetEnUp_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetEnDown_smeared,  met_ey_JetEnDown_smeared );
	MT_JetEnDown_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_JetEnDown_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_UnclusteredEnUp_smeared,  met_ey_UnclusteredEnUp_smeared );
	MT_UnclusteredEnUp_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_UnclusteredEnUp_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_UnclusteredEnDown_smeared,  met_ey_UnclusteredEnDown_smeared );
	MT_UnclusteredEnDown_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_UnclusteredEnDown_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetResUp_smeared,  met_ey_JetResUp_smeared );
	MT_JetResUp_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_JetResUp_smeared*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], met_ex_JetResDown_smeared,  met_ey_JetResDown_smeared );
	MT_JetResDown_smeared = TMath::Sqrt(2*mu_pt[0]*met_pt_JetResDown_smeared*(1-TMath::Cos(dPhi)));


///////// PF met

      PFmet_ex = analysisTree.pfmet_ex;
      PFmet_ey = analysisTree.pfmet_ey;


      PFmet_ez = 0;
      PFmet_pt = TMath::Sqrt(PFmet_ex*PFmet_ex+PFmet_ey*PFmet_ey);
      PFmet_phi = TMath::ATan2(PFmet_ey,PFmet_ex);

     PFmet_ex_JetEnUp = analysisTree.pfmet_ex_JetEnUp;
     PFmet_ey_JetEnUp = analysisTree.pfmet_ey_JetEnUp;
     PFmet_ex_JetEnDown = analysisTree.pfmet_ex_JetEnDown;
     PFmet_ey_JetEnDown = analysisTree.pfmet_ey_JetEnDown;
     PFmet_ex_UnclusteredEnUp = analysisTree.pfmet_ex_UnclusteredEnUp;
     PFmet_ey_UnclusteredEnUp = analysisTree.pfmet_ey_UnclusteredEnUp;
     PFmet_ex_UnclusteredEnDown = analysisTree.pfmet_ex_UnclusteredEnDown;
     PFmet_ey_UnclusteredEnDown = analysisTree.pfmet_ey_UnclusteredEnDown;

      PFmet_pt_JetEnUp = TMath::Sqrt(PFmet_ex_JetEnUp*PFmet_ex_JetEnUp+PFmet_ey_JetEnUp*PFmet_ey_JetEnUp);
      PFmet_pt_JetEnDown = TMath::Sqrt(PFmet_ex_JetEnDown*PFmet_ex_JetEnDown+PFmet_ey_JetEnDown*PFmet_ey_JetEnDown);
      PFmet_pt_UnclusteredEnUp = TMath::Sqrt(PFmet_ex_UnclusteredEnUp*PFmet_ex_UnclusteredEnUp+PFmet_ey_UnclusteredEnUp*PFmet_ey_UnclusteredEnUp);
      PFmet_pt_UnclusteredEnDown = TMath::Sqrt(PFmet_ex_UnclusteredEnDown*PFmet_ex_UnclusteredEnDown+PFmet_ey_UnclusteredEnDown*PFmet_ey_UnclusteredEnDown);

	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], PFmet_ex,  PFmet_ey );
	MT_PF = TMath::Sqrt(2*mu_pt[0]*PFmet_pt*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], PFmet_ex_JetEnUp,  PFmet_ey_JetEnUp );
	MT_PF_JetEnUp = TMath::Sqrt(2*mu_pt[0]*PFmet_pt_JetEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], PFmet_ex_JetEnDown,  PFmet_ey_JetEnDown );
	MT_PF_JetEnDown = TMath::Sqrt(2*mu_pt[0]*PFmet_pt_JetEnDown*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], PFmet_ex_UnclusteredEnUp,  PFmet_ey_UnclusteredEnUp );
	MT_PF_UnclusteredEnUp = TMath::Sqrt(2*mu_pt[0]*PFmet_pt_UnclusteredEnUp*(1-TMath::Cos(dPhi)));
	dPhi=dPhiFrom2P( mu_px[0], mu_py[0], PFmet_ex_UnclusteredEnDown,  PFmet_ey_UnclusteredEnDown );
	MT_PF_UnclusteredEnDown = TMath::Sqrt(2*mu_pt[0]*PFmet_pt_UnclusteredEnDown*(1-TMath::Cos(dPhi)));
////////////////// PFmet end


////////////////////////////////////// MET performance 
	Float_t Utx,Uty;
	Utx = - mu_px[0] - met_ex;
	Uty = - mu_py[0] - met_ey;

	Ut = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol = TMath::Cos(dPhi)* Ut;
//cout<<"!!!!" <<endl;
//cout << Ut<< endl;
//cout << Utr<< endl;
//cout << Ucol<< endl;


	Utx = - mu_px[0] - met_ex_JetEnUp;
	Uty = - mu_py[0] - met_ey_JetEnUp;
	Ut_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_JetEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_JetEnUp = TMath::Cos(dPhi)* Ut_JetEnUp ;

	Utx = - mu_px[0] - met_ex_JetEnDown;
	Uty = - mu_py[0] - met_ey_JetEnDown;
	Ut_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_JetEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_JetEnDown = TMath::Cos(dPhi)* Ut_JetEnDown ;

	Utx = - mu_px[0] - met_ex_JetResDown;
	Uty = - mu_py[0] - met_ey_JetResDown;
	Ut_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_JetResDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_JetResDown = TMath::Cos(dPhi)* Ut_JetResDown ;

	Utx = - mu_px[0] - met_ex_JetResUp;
	Uty = - mu_py[0] - met_ey_JetResUp;
	Ut_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_JetResUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_JetResUp = TMath::Cos(dPhi)* Ut_JetResUp ;

	Utx = - mu_px[0] - met_ex_UnclusteredEnUp;
	Uty = - mu_py[0] - met_ey_UnclusteredEnUp;
	Ut_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_UnclusteredEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_UnclusteredEnUp ;

	Utx = - mu_px[0] - met_ex_UnclusteredEnDown;
	Uty = - mu_py[0] - met_ey_UnclusteredEnDown;
	Ut_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_UnclusteredEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_UnclusteredEnDown ;



////////// Puppi
	Utx = - mu_px[0] - puppi_ex;
	Uty = - mu_py[0] - puppi_ey;
	Ut_puppi = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi = TMath::Cos(dPhi)* Ut;

	Utx = - mu_px[0] - puppi_ex_JetEnUp;
	Uty = - mu_py[0] - puppi_ey_JetEnUp;
	Ut_puppi_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_JetEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_JetEnUp = TMath::Cos(dPhi)* Ut_puppi_JetEnUp ;

	Utx = - mu_px[0] - puppi_ex_JetEnDown;
	Uty = - mu_py[0] - puppi_ey_JetEnDown;
	Ut_puppi_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_JetEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_JetEnDown = TMath::Cos(dPhi)* Ut_puppi_JetEnDown ;

	Utx = - mu_px[0] - puppi_ex_JetResDown;
	Uty = - mu_py[0] - puppi_ey_JetResDown;
	Ut_puppi_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_JetResDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_JetResDown = TMath::Cos(dPhi)* Ut_puppi_JetResDown ;

	Utx = - mu_px[0] - puppi_ex_JetResUp;
	Uty = - mu_py[0] - puppi_ey_JetResUp;
	Ut_puppi_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_JetResUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_JetResUp = TMath::Cos(dPhi)* Ut_puppi_JetResUp ;

	Utx = - mu_px[0] - puppi_ex_UnclusteredEnUp;
	Uty = - mu_py[0] - puppi_ey_UnclusteredEnUp;
	Ut_puppi_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_UnclusteredEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_puppi_UnclusteredEnUp ;

	Utx = - mu_px[0] - puppi_ex_UnclusteredEnDown;
	Uty = - mu_py[0] - puppi_ey_UnclusteredEnDown;
	Ut_puppi_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_puppi_UnclusteredEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_puppi_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_puppi_UnclusteredEnDown ;



////// smeared met

	Utx = - mu_px[0] - met_ex_smeared;
	Uty = - mu_py[0] - met_ey_smeared;
	Ut_smeared = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared = TMath::Cos(dPhi)* Ut;

	Utx = - mu_px[0] - met_ex_JetEnUp_smeared;
	Uty = - mu_py[0] - met_ey_JetEnUp_smeared;
	Ut_smeared_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_JetEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_JetEnUp = TMath::Cos(dPhi)* Ut_smeared_JetEnUp ;

	Utx = - mu_px[0] - met_ex_JetEnDown_smeared;
	Uty = - mu_py[0] - met_ey_JetEnDown_smeared;
	Ut_smeared_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_JetEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_JetEnDown = TMath::Cos(dPhi)* Ut_smeared_JetEnDown ;

	Utx = - mu_px[0] - met_ex_JetResDown_smeared;
	Uty = - mu_py[0] - met_ey_JetResDown_smeared;
	Ut_smeared_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_JetResDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_JetResDown = TMath::Cos(dPhi)* Ut_smeared_JetResDown ;

	Utx = - mu_px[0] - met_ex_JetResUp_smeared;
	Uty = - mu_py[0] - met_ey_JetResUp_smeared;
	Ut_smeared_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_JetResUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_JetResUp = TMath::Cos(dPhi)* Ut_smeared_JetResUp ;

	Utx = - mu_px[0] - met_ex_UnclusteredEnUp_smeared;
	Uty = - mu_py[0] - met_ey_UnclusteredEnUp_smeared;
	Ut_smeared_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_UnclusteredEnUp = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_smeared_UnclusteredEnUp ;

	Utx = - mu_px[0] - met_ex_UnclusteredEnDown_smeared;
	Uty = - mu_py[0] - met_ey_UnclusteredEnDown_smeared;
	Ut_smeared_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	Utr_smeared_UnclusteredEnDown = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	Ucol_smeared_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_smeared_UnclusteredEnDown ;

	


////////////////////////////////////// MET performance end

////////////////////////////////////// real MC recoil 
if (string::npos != filen.find("JetsToLNu"))

	{
	float mu_pyMC=0;
	float mu_pxMC=0;
	float nu_pxMC=0;
	float nu_pyMC=0;
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) 
		{

	  	TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);
		//cout << fabs(analysisTree.genparticles_mother[igen]) <<endl;
		if (/*fabs(analysisTree.genparticles_mother[igen])==24 && */fabs(analysisTree.genparticles_pdgid[igen])==14 && (genLV.Pt()*genLV.Pt())>(nu_pxMC*nu_pxMC+nu_pyMC*nu_pyMC)) {nu_pyMC=genLV.Py();nu_pxMC=genLV.Px();}
		if (/*fabs(analysisTree.genparticles_mother[igen])==24 && */fabs(analysisTree.genparticles_pdgid[igen])==13 && (genLV.Pt()*genLV.Pt())>(mu_pxMC*mu_pxMC+mu_pyMC*mu_pyMC)) {mu_pyMC=genLV.Py();mu_pxMC=genLV.Px();}

		}



		
	float UtxMC = - mu_pxMC - nu_pxMC;
	float UtyMC = - mu_pyMC - nu_pyMC;
	UtMC = TMath::Sqrt((UtxMC*UtxMC)+(UtyMC*UtyMC));
	dPhi=dPhiFrom2P( UtxMC, UtyMC, mu_pxMC,  mu_pyMC );
	UtrMC = -UtxMC*((mu_pyMC/mu_pxMC)/TMath::Sqrt(1+(mu_pyMC/mu_pxMC)*(mu_pyMC/mu_pxMC))) + UtyMC*(1/TMath::Sqrt(1+(mu_pyMC/mu_pxMC)*(mu_pyMC/mu_pxMC)));
	UcolMC = TMath::Cos(dPhi)* UtMC;
/*
if (UtMC < 1)
{
cout << "mu_pyMC  " << mu_pyMC<< endl; 
cout << "mu_pxMC  " << mu_pxMC<< endl; 
cout << "nu_pyMC  " << nu_pyMC<< endl; 
cout << "nu_pxMC  " << nu_pxMC<< endl; 
cout <<"UtMC  "<< UtMC<<endl;
}
*/

	}


//////////////////////////////////////


      if (!isData) npartons = analysisTree.genparticles_noutgoing;
      npv =  analysisTree.primvertex_count;
      all_weight = weight;
      T->Fill();
      selEvents++;
      continue;
      /////////////////////////////////////////////////



    } // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }


  cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;


  std::cout << std::endl;
  int allEvents = (int)inputEventsH->GetEntries();
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;


  file->cd(Channel.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histTopPt->Write();
  histRuns->Write();
  CutFlowUnW->Write();
  MuSF_IdIso_Mu1H->Write();
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}
