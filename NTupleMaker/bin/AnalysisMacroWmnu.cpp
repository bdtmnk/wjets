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
#include "XYMETCorrection.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"


int main(int argc, char * argv[]) {

    //Configuration=====================================================
    Config cfg(argv[1]);
    string Channel="Wmnu";
    bool usePuppiMET = false;
    bool applySimpleRecoilCorrections = 1;
    
    //TODO declare it!!!
    //Global variables:=================================================
	//const bool isRunMETFilters = cfg.get<bool>("RunMETFilters");
	//const bool isRunZpTreWeighting = cfg.get<bool>("RunZpTreWeighting");
	//const bool isRunTTreWeighting = cfg.get<bool>("RunTTreWeighting");

    //Kinematic cuts on electrons=======================================
    const bool isData = cfg.get<bool>("IsData");
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
    const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

    //Muons=============================================================

    const double ptMuonCut     = cfg.get<double>("ptMuonCut");
    const double etaMuonCut    = cfg.get<double>("etaMuonCut");
    const double dxyMuonCut    = cfg.get<double>("dxyMuonCut");
    const double dzMuonCut     = cfg.get<double>("dzMuonCut");
    const double isoMuonLowCut = cfg.get<double>("isoMuonLowCut");
    const bool   applyMuonId   = cfg.get<bool>("ApplyMuonId");

    //Tau===============================================================
    const double ptTauCut = cfg.get<double>("ptTauCut");
    const double etaTauCut = cfg.get<double>("etaTauCut");
    const double decayModeFinding    = cfg.get<double>("decayModeFinding");
    const double againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
    const double againstMuonTight3  = cfg.get<double>("againstMuonTight3");
    const double vertexz =  cfg.get<double>("vertexz");
    const double byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");


    //Veto muon=========================================================
    const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
    const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
    const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
    const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
    const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
    const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

    //Dilepton Veto=====================================================
    const float ptDilepMuonCut = cfg.get<float>("ptDilepMuonCut");
    const float etaDilepMuonCut = cfg.get<float>("etaDilepMuonCut");
    const float dxyDilepMuonCut = cfg.get<float>("dxyDilepMuonCut");
    const float dzDilepMuonCut = cfg.get<float>("dzDilepMuonCut");
    const float isoDilepMuonCut = cfg.get<float>("isoDilepMuonCut");
    const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCut");

    //Veto Electrons====================================================
    const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
    const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
    const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
    const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
    const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
    const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

    //Vertex Cuts=======================================================
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


    //Kinematic cuts on Jets============================================
    const double etaJetCut   = cfg.get<double>("etaJetCut");
    const double ptJetCut   = cfg.get<double>("ptJetCut");

    //Top Single Muon Trigger File======================================
    const double dRleptonsCutmutau   = cfg.get<double>("dRleptonsCutmutau");
    const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
    const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
    const string TrigLeg  = cfg.get<string>("SingleMuonFilterName") ;
    const double SingleMuonTriggerPtCut = cfg.get<double>("SingleMuonTriggerPtCut");

    //Vertex distributions filenames and histname
    //Golden Json:
    const string jsonFile = cfg.get<string>("jsonFile");

    string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
 
    //const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");
    //Run-lumi selector

    std::vector<Period> periods;



  if (isData) {
    // read only the good runs:)))======================================
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

  xs=1; fact=1; fact2=1;
  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  ifstream ifs("xsecs");
  string line;

  XSec = 1.;
  xsecs = XSec;
  std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;

  //PU Reweighting======================================================
  PileUp * PUofficial = new PileUp();

  //TODO 1 attach PU file to the config=================================
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase) +
		   "/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_data_Autumn18.root","read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase) +
			"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/pileUp_MC_Autumn18.root", "read");

  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");

  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);
  bool SUSY = false;

  // B-tag==============================================================
  // TODO attach B-tag file to the config
  string BTag_ = "central";
  string BtagCVS = "DeepCSV_2017data.csv" ;

  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/"+BtagCVS);
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,BTag_);

  if (!SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"comb");
      reader_C.load(calib,BTagEntry::FLAV_C,"comb");
      reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");
  }

  if (SUSY){

    reader_B.load(calib,BTagEntry::FLAV_B,"fastsim");
    reader_C.load(calib,BTagEntry::FLAV_C,"fastsim");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"fastsim");

  }

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


  float MaxBJetPt = 1000.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;

  //TODO attach the tagging efficiencies to the config file:============
  TFile * fileTagging = new TFile(TString(cmsswBase)+
                        TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_march2018_btageff-all_samp-inc-DeepCSV_medium.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  // TODO add Muon Id Isolation to the config:==========================
  TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);
  cout<<"  Initializing iD SF files....."<<endl;

  ScaleFactor * SF_muonIdIso;

  if (applyLeptonSF) {

    SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonidIsoEffFile));

  }

  cout<<"Reading of Trigger SF files:"<<endl;

  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuontrigEffFile));


  //Recoil path:========================================================
  cout<<"Reading Filenames"<<endl;
  const string recoilFileName   = cfg.get<string>("RecoilFileName");
  TString RecoilFileName(recoilFileName);

  const string recoilFileNamePuppi   = cfg.get<string>("RecoilFileNamePuppi");
  TString RecoilFileNamePuppi(recoilFileNamePuppi);

  const string metSysFileName   = cfg.get<string>("MetSysFileName");
  TString MetSysFileName(metSysFileName);

  const string metSysFileNamePuppi   = cfg.get<string>("MetSysFileNamePuppi");
  TString MetSysFileNamePuppi(metSysFileNamePuppi);

  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
   
  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

/*
  const string zMassPtWeightsFileName   = cfg.get<string>("TopPtWeightsFileName");
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
   
  const string zMassPtWeightsHistName   = cfg.get<string>("TopPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
*/

  //Definition of the Corrector:========================================
  kit::MEtSys metSys(MetSysFileName);
  kit::RecoilCorrector recoilMetCorrector(RecoilFileName);
  kit::MEtSys metSysPuppi(MetSysFileNamePuppi);
  kit::RecoilCorrector recoilMetCorrectorPuppi(RecoilFileNamePuppi);

  cout<<"Start"<<endl;

  char ff[100];
  sprintf(ff,"%s/%s",argv[3],argv[2]);
  int nTotalFiles = 0;

  //File name and tree name=============================================
  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  std::ifstream fileList(ff);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString era=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();
  std::string st1,st2;

  TFile * file;
  file = new TFile(era+"/"+TStrName+TString("_DataDriven.root"),"update");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;
  std::string dummy;
  while (fileList0 >> dummy) nTotalFiles++;

  //additional branches
  SetupTree();

  //Definition of the MT:===============================================
  Float_t         puppiMT;
  Float_t         pfMT;
  Float_t         pfMT_corr;
  Float_t         pfMT_XY;
  Float_t         pfMT_corr_XY;
  Float_t         MT_smeared;

  Float_t         raw_met_pt;
  Float_t         raw_met_phi;


  T->Branch("puppiMT", &puppiMT, "MTpuppi/F");
  T->Branch("pfMT", &pfMT, "pfMT/F");
  T->Branch("pfMT_corr", &pfMT_corr, "pfMT_corr/F");
  T->Branch("pfMT_XY", &pfMT_XY, "pfMT_XY/F");
  T->Branch("pfMT_corr_XY", &pfMT_corr_XY, "pfMT_corr_XY/F");
  T->Branch("MT_smeared", &MT_smeared, "MT_smeared/F");


  Float_t      puppi_ex;
  Float_t      puppi_ey;
  Float_t      puppi_ez;
  Float_t      puppi_pt;
  Float_t      puppi_phi;

  T->Branch("puppi_ex", &puppi_ex, "puppi_ex/F");
  T->Branch("puppi_ey", &puppi_ey, "puppi_ey/F");
  T->Branch("puppi_ez", &puppi_ez, "puppi_ez/F");
  T->Branch("puppi_pt", &puppi_pt, "puppi_pt/F");
  T->Branch("puppi_phi", &puppi_phi, "puppi_phi/F");

  //PF Collections:=====================================================
  Float_t      pfmet_ex;
  Float_t      pfmet_ey;
  Float_t      pfmet_ez;
  Float_t      pfmet_pt;
  Float_t      pfmet_phi;

  T->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
  T->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
  T->Branch("pfmet_ez", &pfmet_ez, "pfmet_ez/F");
  T->Branch("pfmet_pt", &pfmet_pt, "pfmet_pt/F");
  T->Branch("pfmet_phi", &pfmet_phi, "pfmet_phi/F");

   //Corrected MET======================================================
   Float_t      pfmet_pt_corr;
   Float_t      pfmet_phi_corr;
   Float_t      pfmet_ex_corr;
   Float_t      pfmet_ey_corr;

   T->Branch("pfmet_pt_corr", &pfmet_pt_corr, "pfmet_pt_corr/F");
   T->Branch("pfmet_phi_corr", &pfmet_phi_corr, "pfmet_phi_corr/F");
   T->Branch("pfmet_ex_corr", &pfmet_ex_corr, "pfmet_ex_corr/F");
   T->Branch("pfmet_ey_corr", &pfmet_ey_corr, "pfmet_ey_corr/F");

   //XY corrected MET===================================================
   Float_t      pfmet_pt_XY;
   Float_t      pfmet_phi_XY;
   Float_t      pfmet_ex_XY;
   Float_t      pfmet_ey_XY;

   T->Branch("pfmet_pt_XY", &pfmet_pt_XY, "pfmet_pt_XY/F");
   T->Branch("pfmet_phi_XY", &pfmet_phi_XY, "pfmet_phi_XY/F");
   T->Branch("pfmet_ex_XY", &pfmet_ex_XY, "pfmet_ex_XY/F");
   T->Branch("pfmet_ey_XY", &pfmet_ey_XY, "pfmet_ey_XY/F");

   Float_t      pfmet_pt_corr_XY;
   Float_t      pfmet_phi_corr_XY;
   Float_t      pfmet_ex_corr_XY;
   Float_t      pfmet_ey_corr_XY;

   //Add Corrected MET:=================================================
   T->Branch("pfmet_pt_corr_XY", &pfmet_pt_corr_XY, "pfmet_pt_corr_XY/F");
   T->Branch("pfmet_phi_corr_XY", &pfmet_phi_corr_XY, "pfmet_phi_corr_XY/F");
   T->Branch("pfmet_ex_corr_XY", &pfmet_ex_corr_XY, "pfmet_ex_corr_XY/F");
   T->Branch("pfmet_ey_corr_XY", &pfmet_ey_corr_XY, "pfmet_ey_corr_XY/F");

   Float_t      pfmet_ex_corr_smeared;
   Float_t      pfmet_ey_corr_smeared;
   Float_t      pfmet_pt_corr_smeared;
   Float_t      pfmet_phi_corr_smeared;


    //MEt Systematics===================================================
    Float_t pfmetcorr_ex_JetEnUp;
    Float_t pfmetcorr_ey_JetEnUp;

    Float_t pfmetcorr_ex_JetEnDown;
    Float_t pfmetcorr_ey_JetEnDown;

    Float_t pfmetcorr_ex_UnclusteredEnUp;
    Float_t pfmetcorr_ey_UnclusteredEnUp;

    Float_t pfmetcorr_ex_UnclusteredEnDown;
    Float_t pfmetcorr_ey_UnclusteredEnDown;

    Float_t pfmetcorr_ex_JetResUp;
    Float_t pfmetcorr_ey_JetResUp;

    Float_t pfmetcorr_ex_JetResDown;
    Float_t pfmetcorr_ey_JetResDown;

    Float_t puppimet_ex_JetEnUp;
    Float_t puppimet_ey_JetEnUp;

    Float_t puppimet_ex_JetEnDown;
    Float_t puppimet_ey_JetEnDown;

    Float_t puppimet_ex_UnclusteredEnUp;
    Float_t puppimet_ey_UnclusteredEnUp;

    Float_t puppimet_ex_UnclusteredEnDown;
    Float_t puppimet_ey_UnclusteredEnDown;

    Float_t puppimet_ex_JetResUp;
    Float_t puppimet_ey_JetResUp;

    Float_t puppimet_ex_JetResDown;
    Float_t puppimet_ey_JetResDown;

    T->Branch("pfmetcorr_ex_JetEnUp", &pfmetcorr_ex_JetEnUp, "pfmetcorr_ex_JetEnUp/F");
    T->Branch("pfmetcorr_ey_JetEnUp", &pfmetcorr_ey_JetEnUp, "pfmetcorr_ey_JetEnUp/F");

    T->Branch("pfmetcorr_ex_JetEnDown", &pfmetcorr_ex_JetEnDown, "pfmetcorr_ex_JetEnDown/F");
    T->Branch("pfmetcorr_ey_JetEnDown", &pfmetcorr_ey_JetEnDown, "pfmetcorr_ey_JetEnDown/F");

    T->Branch("pfmetcorr_ex_UnclusteredEnUp", &pfmetcorr_ex_UnclusteredEnUp, "pfmetcorr_ex_UnclusteredEnUp/F");
    T->Branch("pfmetcorr_ey_UnclusteredEnUp", &pfmetcorr_ey_UnclusteredEnUp, "pfmetcorr_ey_UnclusteredEnUp/F");

    T->Branch("pfmetcorr_ex_UnclusteredEnDown", &pfmetcorr_ex_UnclusteredEnDown, "pfmetcorr_ex_UnclusteredEnDown/F");
    T->Branch("pfmetcorr_ey_UnclusteredEnDown", &pfmetcorr_ey_UnclusteredEnDown, "pfmetcorr_ey_UnclusteredEnDown/F");

    T->Branch("pfmetcorr_ex_JetResUp", &pfmetcorr_ex_JetResUp, "pfmetcorr_ex_JetResUp/F");
    T->Branch("pfmetcorr_ey_JetResUp", &pfmetcorr_ey_JetResUp, "pfmetcorr_ey_JetResUp/F");

    T->Branch("pfmetcorr_ex_JetResDown", &pfmetcorr_ex_JetResDown, "pfmetcorr_ex_JetResDown/F");
    T->Branch("pfmetcorr_ey_JetResDown", &pfmetcorr_ey_JetResDown, "pfmetcorr_ey_JetResDown/F");

    T->Branch("puppimet_ex_JetEnUp", &puppimet_ex_JetEnUp, "puppimet_ex_JetEnUp/F");
    T->Branch("puppimet_ey_JetEnUp", &puppimet_ey_JetEnUp, "puppimet_ey_JetEnUp/F");

    T->Branch("puppimet_ex_JetEnDown", &puppimet_ex_JetEnDown, "puppimet_ex_JetEnDown/F");
    T->Branch("puppimet_ey_JetEnDown", &puppimet_ey_JetEnDown, "puppimet_ey_JetEnDown/F");

    T->Branch("puppimet_ex_UnclusteredEnUp", &puppimet_ex_UnclusteredEnUp, "puppimet_ex_UnclusteredEnUp/F");
    T->Branch("puppimet_ey_UnclusteredEnUp", &puppimet_ey_UnclusteredEnUp, "puppimet_ey_UnclusteredEnUp/F");

    T->Branch("puppimet_ex_UnclusteredEnDown", &puppimet_ex_UnclusteredEnDown, "puppimet_ex_UnclusteredEnDown/F");
    T->Branch("puppimet_ey_UnclusteredEnDown", &puppimet_ey_UnclusteredEnDown, "puppimet_ey_UnclusteredEnDown/F");

    T->Branch("puppimet_ex_JetResUp", &puppimet_ex_JetResUp, "puppimet_ex_JetResUp/F");
    T->Branch("puppimet_ey_JetResUp", &puppimet_ey_JetResUp, "puppimet_ey_JetResUp/F");

    T->Branch("puppimet_ex_JetResDown", &puppimet_ex_JetResDown, "puppimet_ex_JetResDown/F");
    T->Branch("puppimet_ey_JetResDown", &puppimet_ey_JetResDown, "puppimet_ey_JetResDown/F");


    //Add Corrected MET:================================================
    T->Branch("pfmet_ex_corr_smeared", &pfmet_ex_corr_smeared, "pfmet_ex_corr_smeared/F");
    T->Branch("pfmet_ey_corr_smeared", &pfmet_ey_corr_smeared, "pfmet_ey_corr_smeared/F");
    T->Branch("pfmet_pt_corr_smeared", &pfmet_pt_corr_smeared, "pfmet_pt_corr_smeared/F");
    T->Branch("pfmet_phi_corr_smeared", &pfmet_phi_corr_smeared, "pfmet_phi_corr_smeared/F");

    Float_t         sum_pfjet_energycorr;
    Float_t         sum_pfjet_energycorr_l1fastjet;
    Float_t         sum_pfjet_energycorr_l2relative;
    Float_t         sum_pfjet_energycorr_l3absolute;
    Float_t         sum_pfjet_energycorr_l2l3residual;

    T->Branch("sum_pfjet_energycorr",  &sum_pfjet_energycorr, "sum_pfjet_energycorr/F");
    T->Branch("sum_pfjet_energycorr_l1fastjet", &sum_pfjet_energycorr_l1fastjet, "sum_pfjet_energycorr_l1fastjet/F");
    T->Branch("sum_pfjet_energycorr_l2relative", &sum_pfjet_energycorr_l2relative, "sum_pfjet_energycorr_l2relative//F");
    T->Branch("sum_pfjet_energycorr_l3absolute", &sum_pfjet_energycorr_l3absolute, "sum_pfjet_energycorr_l3absolute/F");
    T->Branch("sum_pfjet_energycorr_l2l3residual", &sum_pfjet_energycorr_l2l3residual, "sum_pfjet_energycorr_l2l3residual/F");

    Float_t         pfjet_energycorr[200];
    Float_t         pfjet_energycorr_l1fastjet[200];
    Float_t         pfjet_energycorr_l2relative[200];
    Float_t         pfjet_energycorr_l3absolute[200];
    Float_t         pfjet_energycorr_l2l3residual[200];

    T->Branch("pfjet_energycorr",  &pfjet_energycorr, "pfjet_energycorr[200]/F");
    T->Branch("pfjet_energycorr_l1fastjet", &pfjet_energycorr_l1fastjet, "pfjet_energycorr_l1fastjet[200]/F");
    T->Branch("pfjet_energycorr_l2relative", &pfjet_energycorr_l2relative, "pfjet_energycorr_l2relative[200]/F");
    T->Branch("pfjet_energycorr_l3absolute", &pfjet_energycorr_l3absolute, "pfjet_energycorr_l3absolute[200]/F");
    T->Branch("pfjet_energycorr_l2l3residual", &pfjet_energycorr_l2l3residual, "pfjet_energycorr_l2l3residual[200]/F");


     //Ut performance===================================================
     Float_t Ut;
     Float_t Utr;
     Float_t Ucol;

     T->Branch("Ut", &Ut, "Ut/F");
     T->Branch("Utr", &Utr, "Utr/F");
     T->Branch("Ucol", &Ucol, "Ucol/F");

     Float_t UtMC;
     Float_t UtrMC;
     Float_t UcolMC;

     T->Branch("UtMC", &UtMC, "UtMC/F");
     T->Branch("UtrMC", &UtrMC, "UtrMC/F");
     T->Branch("UcolMC", &UcolMC, "UcolMC/F");

     Float_t PUweight;
     T->Branch("PUweight", &PUweight, "PUweight/F");

     Float_t raw_dPhi;
     Float_t raw_pfMT;
     Float_t raw_pfmet_ex;
     Float_t raw_pfmet_ey;

     T->Branch("raw_dPhi",  &raw_dPhi, "raw_dPhi/F");
     T->Branch("raw_pfMT", &raw_pfMT, "raw_pfMT/F");
     T->Branch("raw_pfmet_ex", &raw_pfmet_ex, "raw_pfmet_ex//F");
     T->Branch("raw_pfmet_ey", &raw_pfmet_ey, "raw_pfmet_ey/F");

     Float_t raw_dPhi_XY;
     Float_t raw_pfMT_XY;
     Float_t raw_pfmet_ex_XY;
     Float_t raw_pfmet_ey_XY;

     T->Branch("raw_dPhi_XY",  &raw_dPhi_XY, "raw_dPhi_XY/F");
     T->Branch("raw_pfMT_XY", &raw_pfMT_XY, "raw_pfMT_XY/F");
     T->Branch("raw_pfmet_ex_XY", &raw_pfmet_ex_XY, "raw_pfmet_ex_XY//F");
     T->Branch("raw_pfmet_ey_XY", &raw_pfmet_ey_XY, "raw_pfmet_ey_XY/F");

      float met_x_recoilscaleUp;
      float met_x_recoilscaleDown;
      float met_y_recoilscaleUp;
      float met_y_recoilscaleDown;
      float met_x_recoilresoUp;
      float met_x_recoilresoDown;
      float met_y_recoilresoUp;
      float met_y_recoilresoDown;

      float met_unclMetUp_x;
      float met_unclMetUp_y;
      float met_unclMetDown_x;
      float met_unclMetDown_y;

     T->Branch("met_x_recoilscaleUp",  &met_x_recoilscaleUp, "met_x_recoilscaleUp/F");
     T->Branch("met_x_recoilscaleDown",  &met_x_recoilscaleDown, "met_x_recoilscaleDown/F");
     T->Branch("met_y_recoilscaleUp",  &met_y_recoilscaleUp, "met_y_recoilscaleUp/F");
     T->Branch("met_y_recoilscaleDown",  &met_y_recoilscaleDown, "met_y_recoilscaleDown/F");

     T->Branch("met_x_recoilresoUp",  &met_x_recoilresoUp, "met_x_recoilresoUp/F");
     T->Branch("met_x_recoilscaleDown",  &met_x_recoilscaleDown, "met_x_recoilscaleDown/F");
     T->Branch("met_y_recoilscaleUp",  &met_y_recoilscaleUp, "met_y_recoilscaleUp/F");
     T->Branch("met_y_recoilscaleDown",  &met_y_recoilscaleDown, "met_y_recoilscaleDown/F");

     T->Branch("met_unclMetUp_x",  &met_unclMetUp_x, "met_unclMetUp_x/F");
     T->Branch("met_unclMetUp_y",  &met_unclMetUp_y, "met_unclMetUp_y/F");
     T->Branch("met_unclMetDown_x",  &met_unclMetDown_x, "met_unclMetDown_x/F");
     T->Branch("met_unclMetDown_y",  &met_unclMetDown_y, "met_unclMetDown_y/F");


    //////     MET performance end declaration:=========================
   if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);
   for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "# " << iF+1 << " Out of " << nTotalFiles << " Filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));

    bool WithInit = true;
    if (WithInit) cout << "With Init Root Tree"<<endl;
    if (!WithInit) cout << "Without Init Root Tree"<<endl;

    TTree * _inittree = NULL;
    if (!WithInit)  _inittree = (TTree*)file_->Get(TString(ntupleName));
    if (WithInit)  _inittree = (TTree*)file_->Get(TString(initNtupleName));

    if (_inittree==NULL) continue;
    Float_t genweight;
    if (!isData){
      _inittree->SetBranchAddress("genweight",&genweight);
    }
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
    std::cout << "      Number of Entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; ++iEntry) {
            _inittree->GetEntry(iEntry);
            if (isData)histWeightsH->Fill(0.,1.);
    }

    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;

    Long64_t numberOfEntries = _tree->GetEntries();
    std::cout << "      Number  of Entries in Tree      = " << numberOfEntries << std::endl;
    AC1B analysisTree(_tree);

	if (!isData && !WithInit){
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) 
			{
			analysisTree.GetEntry(iEntry);
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
          histWeightsH->Fill(0.,genweights);
        }
      }
      
    //Teresa's Way:=====================================================
	TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/"+ZMassPtWeightsFileName);
	if (fileZMassPtWeights->IsZombie()) {
		std::cout << "File " << TString(cmsswBase) << "/src/" << ZMassPtWeightsFileName << "  does not exist!" << std::endl;
		exit(-1);
	}
	TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get(ZMassPtWeightsHistName);
	if (histZMassPtWeights==NULL) {
		std::cout << "histogram " << ZMassPtWeightsHistName << " is not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << ZMassPtWeightsFileName
		<< std::endl;
		exit(-1);
	}
      
	//Start iteration over the events:==================================
    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) { 

      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (nEvents%50000==0) cout << "      Processed " << nEvents << " events" << endl;
      if (fabs(analysisTree.primvertex_z) > zVertexCut) continue;
      if (analysisTree.primvertex_ndof < ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x
                        + analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex > dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;

      bool lumi=false;
      float topPt = 0;
      float antitopPt = 0;

      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      PUweight = 1.0;
      gen_weight = 1.;
      trig_weight = 1.;


      bool isW = false;
      bool isDY = false;//TODO remove
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
	  
	  // Declaration of variables used for the MET Correction:==========
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



		float met_x;
		float met_y;
		float met;
		TLorentzVector metLV;
		if(!usePuppiMET){
			met_x = analysisTree.pfmetcorr_ex;
			met_y = analysisTree.pfmetcorr_ey;
		}
		else{
			met_x = analysisTree.puppimet_ex;
			met_y = analysisTree.puppimet_ey;
		}
		
		met = TMath::Sqrt(met_x*met_x + met_y*met_y);\
		metLV.SetXYZT(met_x, met_y, 0., met);


	 // METs ===========================================================
/*
	 float met_x_recoilscaleUp;
	 float met_x_recoilscaleDown;
	 float met_y_recoilscaleUp;
	 float met_y_recoilscaleDown;
	 float met_x_recoilresoUp;
	 float met_x_recoilresoDown;
	 float met_y_recoilresoUp;
	 float met_y_recoilresoDown;
	 
	 float met_unclMetUp_x;
	 float met_unclMetUp_y;
	 float met_unclMetDown_x;
	 float met_unclMetDown_y;
 */      
         if (!usePuppiMET){
            met_x_recoilscaleUp = analysisTree.pfmetcorr_ex;
            met_x_recoilscaleDown = analysisTree.pfmetcorr_ex;
            met_y_recoilscaleUp = analysisTree.pfmetcorr_ey;
            met_y_recoilscaleDown = analysisTree.pfmetcorr_ey;
            met_x_recoilresoUp = analysisTree.pfmetcorr_ex;
            met_x_recoilresoDown = analysisTree.pfmetcorr_ex;
            met_y_recoilresoUp = analysisTree.pfmetcorr_ey;
            met_y_recoilresoDown = analysisTree.pfmetcorr_ey;
         
            met_unclMetUp_x    = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            met_unclMetUp_y    = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            met_unclMetDown_x  = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            met_unclMetDown_y  = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
         }
         else{
            met_x_recoilscaleUp = analysisTree.puppimet_ex;
            met_x_recoilscaleDown = analysisTree.puppimet_ex;
            met_y_recoilscaleUp = analysisTree.puppimet_ey;
            met_y_recoilscaleDown = analysisTree.puppimet_ey;
            met_x_recoilresoUp = analysisTree.puppimet_ex;
            met_x_recoilresoDown = analysisTree.puppimet_ex;
            met_y_recoilresoUp = analysisTree.puppimet_ey;
            met_y_recoilresoDown = analysisTree.puppimet_ey;
         
            met_unclMetUp_x    = analysisTree.puppimet_ex_UnclusteredEnUp;
            met_unclMetUp_y    = analysisTree.puppimet_ey_UnclusteredEnUp;
            met_unclMetDown_x  = analysisTree.puppimet_ex_UnclusteredEnDown;
            met_unclMetDown_y  = analysisTree.puppimet_ey_UnclusteredEnDown;
         }
         met = TMath::Sqrt(met_x*met_x + met_y*met_y);
         float metphi = TMath::ATan2(met_y,met_x);
         
		 //MET Significance:============================================
         if (!usePuppiMET){
            float metcov00 = analysisTree.pfmetcorr_sigxx;
            float metcov01 = analysisTree.pfmetcorr_sigxy;
            float metcov10 = analysisTree.pfmetcorr_sigyx;
            float metcov11 = analysisTree.pfmetcorr_sigyy;
         }
         else{
            float metcov00 = analysisTree.puppimet_sigxx;
            float metcov01 = analysisTree.puppimet_sigxy;
            float metcov10 = analysisTree.puppimet_sigyx;
            float metcov11 = analysisTree.puppimet_sigyy;
         }

	  // Apply Simple Recoil Corrections================================
	  int njetsforrecoil = njets;
	  if(isW) njetsforrecoil = njets + 1;
	  
	  float met_uncorr = met;
	  float metphi_uncorr = metphi;
	  float pfmet_corr_x = met_x;
	  float pfmet_corr_y = met_y;
		  
	  
      if ((isW||isDY) && !isData) {

        if (applySimpleRecoilCorrections) {
               
              if (!usePuppiMET){
		 //cout<<met_x<<mer_y<<bosonPx<<bosonPy<<lepPx<<lepPy<<njetsforrecoil<<pfmet_corr_x<<pfmet_corr_y<<endl;
                  //cout<<"met_x<<met_y<<bosonPx<<bosonPy<<lepPx<<lepPy<<njetsforrecoil<<pfmet_corr_x<<pfmet_corr_y"<<endl;
		  //cout<<met_x<<met_y<<bosonPx<<bosonPy<<lepPx<<lepPy<<njetsforrecoil<<pfmet_corr_x<<pfmet_corr_y<<endl;
		  recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                  cout<<met_x<<met_y<<bosonPx<<bosonPy<<lepPx<<lepPy<<njetsforrecoil<<pfmet_corr_x<<pfmet_corr_y<<endl;
		  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                  metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
               }
               else{
                  recoilMetCorrectorPuppi.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                  metSysPuppi.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
               }
            }
            else {
               if (!usePuppiMET)recoilMetCorrector.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
               else recoilMetCorrectorPuppi.Correct(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
            }
			//todo check changes of the correction;
	}

      if (!isData ) {
            genweightsTree->GetEntry(iEntry);
            weight *= genweights;
            gen_weight *=genweights;
            std::cout <<"analysisTree.genweight "<< float(analysisTree.genweight) << std::endl;
            lumi=true;
      }
	/*
	//Zpt Weights:
	// Zpt reweighting for LO DY samples 
	TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
	TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");
	
	ZPtWeightSys* zPtWeightSys = 0;
	TopPtWeightSys* topPtWeightSys = 0;
	
	////////////////////////////////////////////////////////////
	// Z pt weight
	////////////////////////////////////////////////////////////

	TLorentzVector genV( 0., 0., 0., 0.);
	TLorentzVector genL( 0., 0., 0., 0.);

	otree->zptweight = 1.;
	if (!isData && ((isDY && isMG ) || isEWKZ)){
        genV = genTools::genV(analysisTree); // gen Z boson ?
      	float bosonMass = genV.M();
      	float bosonPt = genV.Pt();

        //Merijn determine here some min and max values:
        double massxmin = h_zptweight->GetXaxis()->GetXmin();
        double massxmax = h_zptweight->GetXaxis()->GetXmax();

        double ptxmin = h_zptweight->GetYaxis()->GetXmin();
        double ptxmax = h_zptweight->GetYaxis()->GetXmax();
		//Here is a secret message
      	//Merijn 2019 6 13: adjust to T/M functions, to get boundaries right. Otherwise, for 2017 data we get few outliers that screw up the weight histogram dramatically.
      	Float_t zptmassweight = 1;
      	if (bosonMass > 50.0) {
          float bosonMassX = bosonMass;
          float bosonPtX = bosonPt;
          if (bosonMassX > massxmax) bosonMassX = massxmax - h_zptweight->GetXaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;//Merijn: if doesn't work, lower by 1/2 binwidth..
          if (bosonPtX < ptxmin)     bosonPtX = ptxmin + h_zptweight->GetYaxis()->GetBinWidth(1)*0.5;
          if (bosonPtX > ptxmax)     bosonPtX = ptxmax - h_zptweight->GetYaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;
          zptmassweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(bosonMassX), h_zptweight->GetYaxis()->FindBin(bosonPtX));
          }
          cout<<"zptmassweight"<<zptmassweight<<endl;	
          //otree->zptweight = zptmassweight;
      }
      */
      
	

	TLorentzVector genBosonLV; genBosonLV.SetXYZT(0,0,0,0);
	TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);

	//Iterate over the gen particles and get the gen Boson:=============
	cout<<"genparticles_count:    "<<analysisTree.genparticles_count<<endl;
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
		
		  cout<<analysisTree.genparticles_px[igen]<<endl;
		  cout<<analysisTree.genparticles_py[igen]<<endl;
		  cout<<analysisTree.genparticles_pz[igen]<<endl;
		  cout<<analysisTree.genparticles_e[igen]<<endl;
		  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
      		  				      analysisTree.genparticles_py[igen],
      						      analysisTree.genparticles_pz[igen],
      						      analysisTree.genparticles_e[igen]);
      						      
		 bool fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1;
	     bool isMuon = false;
	     bool isElectron = false;
	     bool isNeutrino = false;
	     bool isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

		
		  if (abs(analysisTree.genparticles_pdgid[igen])==11) {
                  isElectron = true;
                  //if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                  //   promptElectrons.push_back(genLV);
                  //   promptElectronsLV += genLV;
                  //   wDecayProductsLV += genLV;
                  //}
               }
               
               if (abs(analysisTree.genparticles_pdgid[igen])==13) {
                  isMuon = true;
                  //if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
                  //   promptMuons.push_back(genLV);
                  //   promptMuonsLV += genLV;
                  //   wDecayProductsLV += genLV;
                  //}
               }
		
		bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
		bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);

		if (isBoson)genBosonLV += genLV;
		if (isVisibleBoson)genVisBosonLV += genLV;
		
	
	}

	//Later:=^^^
	bosonPx = genBosonLV.Px();
	bosonPy = genBosonLV.Py();
	bosonPz = genBosonLV.Pz();
	bosonPt = genBosonLV.Pt();
	bosonMass = genBosonLV.M();

	lepPx = genVisBosonLV.Px();
	lepPy = genVisBosonLV.Py();
	lepPz = genVisBosonLV.Pz();

	zptmassweight = 1;
	cout<<"zptmassweight"<<zptmassweight<<endl;

	if (isDY) { 
		cout<<"applying Z pt mass weights: "<<endl;
		cout<<"bosonMass: "<<bosonMass<<endl;
		if (bosonMass>50.0) {
		float bosonMassX = bosonMass;
		float bosonPtX = bosonPt;
		if (bosonMassX>1000.) bosonMassX = 1000.;
		if (bosonPtX<1.)      bosonPtX = 1.;
		if (bosonPtX>1000.)   bosonPtX = 1000.;
		zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
		histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
		cout<<"zptmassweight"<<zptmassweight<<endl;
		}
     	}
	/*
	////////////////////////////////////////////////////////////
	// Top pt weight
	////////////////////////////////////////////////////////////

		otree->topptweight = 1.;
		int a_topPtWeight = cfg.get<int>("a_topPtWeight");
		int b_topPtWeight = cfg.get<int>("b_topPtWeight");
		if(!isData)
		// otree->topptweight = genTools::return_topPtWeight(analysisTree, a_topPtWeight, b_topPtWeight);
		otree->topptweight = genTools::topPtWeight(analysisTree, 1); // 1 is for Run1 - use this reweighting as recommended by HTT 17
		counter[11]++;

	*/
	if (isData)  {

		XSec = 1.;
		histRuns->Fill(analysisTree.event_run);
		int n=analysisTree.event_run;
		int lum = analysisTree.event_luminosityblock;

		cout<<"Lumi: "<<lum<<endl;

		std::string num = std::to_string(n);
		std::string lnum = std::to_string(lum);
	    
		for(const auto& a : periods){
	    
			if ( num.c_str() ==  a.name ) {
	           std::cout<< " Run: "<<num<<"  "<<a.name<<" ";
	          for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
			    cout<<b->lower<<"  "<<b->bigger<<endl;
		        if (lum  >= b->lower && lum <= b->bigger ) lumi = true;}
	            auto last = std::prev(a.ranges.end());
	            if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	    }

	  }
	    if (!lumi) {
	    cout<<"No Lumi"<<endl;
	    continue;
	    }
	    if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
      }

      if (analysisTree.event_run<RunMin)RunMin = analysisTree.event_run;
      if (analysisTree.event_run>RunMax)RunMax = analysisTree.event_run;

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

	  std::vector<TString> metFlags; metFlags.clear();

      //TODO Check MET Flags:
      metFlags.push_back("Flag_goodVertices");
      metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
      metFlags.push_back("Flag_HBHENoiseFilter");
      metFlags.push_back("Flag_HBHENoiseIsoFilter");
      metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
      metFlags.push_back("Flag_BadPFMuonFilter");
      if (isData){
      metFlags.push_back("Flag_eeBadScFilter");
      metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
      }
	   
     bool METflag = metFiltersPasses2(analysisTree, metFlags);
	 met_flag = METflag;
	 cout<<"Met flag: "<<met_flag<<endl;

    if (!isData)
	{
	  if (applyPUreweighting)	 {
	    //puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
		  pu_weight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
	    //puweight = float((PU_data->GetBinContent(analysisTree.primvertex_count)/PU_data->GetSumOfWeights())/(PU_mc->GetBinContent(analysisTree.primvertex_count)/PU_mc->GetSumOfWeights()));
	    if (pu_weight !=pu_weight) {
		    pu_weight=1;
		    //cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;
		    }
	    if (pu_weight >10) {
		    pu_weight=1;
		    //cout<<"bad pu for "<< analysisTree.primvertex_count<<"  " <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;
		}

	    weight *=pu_weight;
	    PUweight = float((PU_data->GetBinContent(analysisTree.primvertex_count)/PU_data->GetSumOfWeights())/(PU_mc->GetBinContent(analysisTree.primvertex_count)/PU_mc->GetSumOfWeights()));

	    if (pu_weight !=pu_weight) {
	                               pu_weight=1;
	                               //cout<<"bad pu for "<< analysisTree.primvertex_count<<endl;
	                               //cout<<"" <<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;
	    }
	    if (pu_weight >10) {
	                                pu_weight=1;
	                                //cout<<"bad pu for "<< analysisTree.primvertex_count<<"  "
	                                //<<PU_data->GetBinContent(analysisTree.primvertex_count)<<"  "
	                                //<<PU_mc->GetBinContent(analysisTree.primvertex_count)<<"  " << puweight<<endl;
	    }
	    PUweight = pu_weight;
	 }
	}


      bool trigAccept = false;
      unsigned int nMainTrigger = 0;
      bool isMainTrigger = false;

      if (1){unsigned int nfilters = analysisTree.run_hltfilters->size();
            
            for (unsigned int i=0; i<nfilters; ++i) {
      	            TString HLTFilter(analysisTree.run_hltfilters->at(i));
	                if (HLTFilter==MainTrigger) {
                        nMainTrigger = i;
	                    isMainTrigger = true;
	                }
            }
	  }

      if (!isMainTrigger) {

            return(-1);

      }

    //now clear the Mu.El.Jets again to fill them again after cleaning
    //cout<<"HLT filter for Mu Trigger"<<endl;

    vector<int> muons; muons.clear();
    for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

        if (analysisTree.muon_pt[im]<ptMuonCut) continue;
        if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
        if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
        if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
        if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
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

    //std::cout << "#Muons: " << muons.size() << "#Taus = " << taus.size() << std::endl;

  	bool isLegMatch = false;
  	for (unsigned int im=0; im<muons.size(); ++im) {
  	isLegMatch = false;
  	bool isMuonTauMuonLegMatch = false;
  	bool isMuonTauOverlapMuonMatch = false;
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

	if (1){
	    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {

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

    bool  dilepton_veto = false;
    bool  extraelec_veto = false;
    bool  extramuon_veto = false;
    event_secondLeptonVeto = false;
    event_thirdLeptonVeto = false;


    //Looking for extra electron
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
        break;
        }
    
        extraelec_veto = foundExtraElectron;


    //Looking for extra muon's (dimuon veto)
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
            float dRmuons = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index]);
            if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[im]*analysisTree.muon_charge[mu_index]<0))dilepton_veto = true;
            }
            if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
            if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
            if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
            if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
            if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
            //if (applyVetoMuonId && !analysisTree.muon_isICHEP[im]) continue;
            if (relIsoMu>isoVetoMuonCut) continue;
            foundExtraMuon = true;
            break;
        }
    
        extramuon_veto = foundExtraMuon;


    //Looking for di_muons:
    if (mu_dimuons.size()>1) {

        for (unsigned int i1=0; i1<mu_dimuons.size()-1; ++i1) {
          unsigned int indx1 = mu_dimuons[i1];

          for (unsigned int i2=i1+1; i2<mu_dimuons.size(); ++i2 ) {

            unsigned int indx2 = mu_dimuons[i2];

            float dRmuons = deltaR(analysisTree.muon_eta[indx1],
                                   analysisTree.muon_phi[indx1],
                                   analysisTree.muon_eta[indx2],
                                   analysisTree.muon_phi[indx2]);

            if (dRmuons>dRDilepVetoCut && (analysisTree.muon_charge[indx1]*analysisTree.muon_charge[indx2]<0))dilepton_veto = true;
            break;
          }
          }
    }
    event_secondLeptonVeto = dilepton_veto;


	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
    if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;

    // Trigger weight
    double ptMu1 = (double)analysisTree.muon_pt[mu_index];
    double etaMu1 = (double)analysisTree.muon_eta[mu_index];
    float trigweight=1.;
	//TODO Check the trigger weights;
    float EffFromData = (float)SF_muonTrigger->get_EfficiencyData(double(ptMu1), double(etaMu1));
    float Mu17EffMC   = (float)SF_muonTrigger->get_EfficiencyMC(double(ptMu1), double(etaMu1));
	

      if (!isData) {
	  if (Mu17EffMC>1e-6)
	  trigweight = EffFromData / Mu17EffMC;
	    //if (!isData && (   string::npos != filen.find("stau") || string::npos != filen.find("C1")) ) trigweight = EffFromData;
        weight *= trigweight;
        trig_weight = trigweight;
        cout<<" Trigger weight "<<trigweight<<endl;
     }

	  // LSF
    if (!isData && applyLeptonSF) {

	    //leptonSFweight = SF_yourScaleFactor->get_ScaleFactor(pt, eta)
	    double IdIsoSF_mu1 = SF_muonIdIso->get_ScaleFactor(ptMu1, etaMu1);
        MuSF_IdIso_Mu1H->Fill(IdIsoSF_mu1);
        weight *= IdIsoSF_mu1;
        LSF_weight = IdIsoSF_mu1;
    }

	  //double HIP_SF1 = SF_HIP->get_ScaleFactor(ptMu1, etaMu1);
	  //cout<<"  "<<ptMu1<<"  "<<etaMu1<<"  "<<HIP_SF1<<endl;
    bool isTauMatched = false;
    bool isGenLeptonMatched = false;

    muon_index = (int)mu_index;
    electron_index = (int)el_index;
    taus_index = (int)tau_index;

    mu_count= (int)analysisTree.muon_count;

    //Fill muon collection:
    for (unsigned int im=0;im<analysisTree.muon_count; ++im){
        //int im = mu_index;
        mu_px[im]=analysisTree.muon_px[im];
        mu_py[im]=analysisTree.muon_py[im];
        mu_pz[im]=analysisTree.muon_pz[im];
        mu_eta[im]=analysisTree.muon_eta[im];
        mu_pt[im]=analysisTree.muon_pt[im];
        mu_phi[im]=analysisTree.muon_phi[im];
        mu_charge[im]=analysisTree.muon_charge[im];
        mu_dxy[im]=analysisTree.muon_dxy[im];
        mu_dz[im]=analysisTree.muon_dz[im];
        mu_dxyerr[im]=analysisTree.muon_dxyerr[im];
        mu_dzerr[im]=analysisTree.muon_dzerr[im];

        mu_neutralHadIso[im] = analysisTree.muon_r04_sumNeutralHadronEt[im];
        mu_photonIso[im] = analysisTree.muon_r04_sumPhotonEt[im];
        mu_chargedHadIso[im] = analysisTree.muon_r04_sumChargedHadronPt[im];
        mu_puIso[im] = analysisTree.muon_r04_sumPUPt[im];
 
        double neutralIso = mu_neutralHadIso[im] + mu_photonIso[im] - 0.5*mu_puIso[im];
        neutralIso = max(double(0),neutralIso);
	      mu_neutralIso[im] = neutralIso;
        mu_absIsoMu[im] = mu_chargedHadIso[im] + neutralIso;
	      mu_relIsoMu[im]  = mu_absIsoMu[im]/mu_pt[im] ;
   
      }

    //Fill Electron collection:
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

    //Fill Tau collection:
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
    }

    //Fill jet collection:
    jet_count=(int)analysisTree.pfjet_count;
    sum_pfjet_energycorr = 0;
    sum_pfjet_energycorr_l1fastjet = 0;
    sum_pfjet_energycorr_l2relative = 0;
    sum_pfjet_energycorr_l3absolute = 0;
    sum_pfjet_energycorr_l2l3residual = 0;

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


        pfjet_energycorr[jj] = analysisTree.pfjet_energycorr[jj];
        pfjet_energycorr_l1fastjet[jj] = analysisTree.pfjet_energycorr_l1fastjet[jj];
        pfjet_energycorr_l2relative[jj] = analysisTree.pfjet_energycorr_l2relative[jj];
        pfjet_energycorr_l3absolute[jj] = analysisTree.pfjet_energycorr_l3absolute[jj];
        pfjet_energycorr_l2l3residual[jj] = analysisTree.pfjet_energycorr_l2l3residual[jj];


        //Corrections:
        sum_pfjet_energycorr += analysisTree.pfjet_energycorr[jj];
        sum_pfjet_energycorr_l1fastjet += analysisTree.pfjet_energycorr_l1fastjet[jj];
        sum_pfjet_energycorr_l2relative += analysisTree.pfjet_energycorr_l2relative[jj];
        sum_pfjet_energycorr_l3absolute += analysisTree.pfjet_energycorr_l3absolute[jj];
        sum_pfjet_energycorr_l2l3residual += analysisTree.pfjet_energycorr_l2l3residual[jj];

     }


    //jets cleaning    TLorentzVector leptonsV, muonJ, jetsLV;
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

     //Loop in all jets:
     for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
        if (fabs(analysisTree.pfjet_pt[jet])<30.) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
        if (absJetEta > etaJetCut) continue;
	    float jetPt = analysisTree.pfjet_pt[jet];
	    bool ispfJetId = false ;
      	bool btagged= false;
	    ispfJetId =looseJetiD(analysisTree,jet);
        if (!ispfJetId) continue;
        bool cleanedJet = true;

	    double Dr=deltaR(analysisTree.muon_eta[mu_index],analysisTree.muon_phi[mu_index], analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
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
                  }
                }
                if (jet_scalefactor>1 && !btagged) { // upgrade
                  double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                  if (rannum<fraction) {
                    btagged = true;
                        std::cout << "upgrading " << std::endl;
                  }
                }
	  } //is Data
	  if (btagged && cleanedJet) bjets.push_back(jet);
	}

	if (cleanedJet){
        jets.push_back((int)jet);
        jets_cleaned[counter_cleaned_jets]=(int)jet;
        jet_jecUn[counter_cleaned_jets] = analysisTree.pfjet_jecUncertainty[jet];
        counter_cleaned_jets++;
	    }
    }

      //Loop in all jets:
      njets = jets.size();
      jet_count = jets.size();
      nbtag = bjets.size();
	  if (nbtag > 0.5) continue;

      // Recoil corrections:
        pfmet_ex_corr = pfmet_corr_x;//analysisTree.pfmetcorr_ex;
        pfmet_ey_corr = pfmet_corr_y;//analysisTree.pfmetcorr_ey;

	//TODO recalulate pt and phi
        pfmet_pt_corr = TMath::Sqrt(pfmet_ex_corr*pfmet_ex_corr + pfmet_ey_corr*pfmet_ey_corr);
	if(pfmet_ex_corr==0 && pfmet_ey_corr>0) pfmet_phi_corr = TMath::Pi();

        else if(pfmet_ex_corr==0 && pfmet_ey_XY<0 )pfmet_phi_corr= -TMath::Pi();

        else if(pfmet_ex_corr >0) pfmet_phi_corr = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr);

        else if(pfmet_ex_corr <0 && pfmet_ey_corr>0) pfmet_phi_corr = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr) + TMath::Pi();

        else if(pfmet_ex_corr <0 && pfmet_ey_corr<0) pfmet_phi_corr = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr) - TMath::Pi();

        else pfmet_phi_corr =0;
	
        //Systematic HERE:
      	pfmetcorr_ex_JetEnUp = analysisTree.pfmetcorr_ex_JetEnUp;
      	pfmetcorr_ey_JetEnUp = analysisTree.pfmetcorr_ey_JetEnUp;

      	pfmetcorr_ex_JetEnDown = analysisTree.pfmetcorr_ex_JetEnDown;
      	pfmetcorr_ey_JetEnDown = analysisTree.pfmetcorr_ey_JetEnDown;

      	pfmetcorr_ex_UnclusteredEnUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
      	pfmetcorr_ey_UnclusteredEnUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;

      	pfmetcorr_ex_UnclusteredEnDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
      	pfmetcorr_ey_UnclusteredEnDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;

      	pfmetcorr_ex_JetResUp = analysisTree.pfmetcorr_ex_JetResUp;
      	pfmetcorr_ey_JetResUp = analysisTree.pfmetcorr_ey_JetResUp;

      	pfmetcorr_ex_JetResDown = analysisTree.pfmetcorr_ex_JetResDown;
      	pfmetcorr_ey_JetResDown = analysisTree.pfmetcorr_ey_JetResDown;


        puppimet_ex_JetEnUp = analysisTree.puppimet_ex_JetEnUp;
        puppimet_ey_JetEnUp = analysisTree.puppimet_ey_JetEnUp;

        puppimet_ex_JetEnDown = analysisTree.puppimet_ex_JetEnDown;
        puppimet_ey_JetEnDown = analysisTree.puppimet_ey_JetEnDown;

        puppimet_ex_UnclusteredEnUp = analysisTree.puppimet_ex_UnclusteredEnUp;
        puppimet_ey_UnclusteredEnUp = analysisTree.puppimet_ey_UnclusteredEnUp;

        puppimet_ex_UnclusteredEnDown = analysisTree.puppimet_ex_UnclusteredEnDown;
        puppimet_ey_UnclusteredEnDown = analysisTree.puppimet_ey_UnclusteredEnDown;

        puppimet_ex_JetResUp = analysisTree.puppimet_ex_JetResUp;
        puppimet_ey_JetResUp = analysisTree.puppimet_ey_JetResUp;

        puppimet_ex_JetResDown = analysisTree.puppimet_ex_JetResDown;
        puppimet_ey_JetResDown = analysisTree.puppimet_ey_JetResDown;

        //Raw met:
        raw_pfmet_ex = pfmet_ex_corr - sum_pfjet_energycorr_l1fastjet;
        raw_pfmet_ey = pfmet_ey_corr - sum_pfjet_energycorr_l1fastjet;

      //Systematic for pfMET:
      
      //Generator level based MET:
      genmet = TMath::Sqrt(analysisTree.genmet_ex*analysisTree.genmet_ex+analysisTree.genmet_ey*analysisTree.genmet_ey);

      //PF MET:
      pfmet_ex = analysisTree.pfmet_ex;
      pfmet_ey = analysisTree.pfmet_ey;
      pfmet_ez = 0;
      pfmet_pt = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      pfmet_phi = TMath::ATan2(pfmet_ey,pfmet_ex);

      //PuppiMET:
      puppi_ex = analysisTree.puppimet_ex;
      puppi_ey = analysisTree.puppimet_ey;
      puppi_ez = 0;//analysisTree.pfmet_ez;
      puppi_pt = TMath::Sqrt(puppi_ex*puppi_ex+puppi_ey*puppi_ey);
      puppi_phi = TMath::ATan2(puppi_ey,puppi_ex);

      //Smeared MET:
      pfmet_ex_corr_smeared = analysisTree.pfmetcorr_ex_smeared;
      pfmet_ey_corr_smeared = analysisTree.pfmetcorr_ey_smeared;
      pfmet_pt_corr_smeared = analysisTree.pfmetcorr_pt_smeared;
      pfmet_phi_corr_smeared = analysisTree.pfmetcorr_phi_smeared;

      //Systematic for Puppi MET
      Float_t dPhi;

      dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex,  pfmet_ey );
	    pfMT = TMath::Sqrt(2*mu_pt[0]*pfmet_pt*(1-TMath::Cos(dPhi)));

      dPhi=dPhiFrom2P( mu_px[0], mu_py[0], puppi_ex,  puppi_ey );
	    puppiMT = TMath::Sqrt(2*mu_pt[0]*puppi_pt*(1-TMath::Cos(dPhi)));

  	  dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex_corr_smeared,  pfmet_ey_corr_smeared );
	    MT_smeared = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_corr_smeared*(1-TMath::Cos(dPhi)));

      dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex_corr,  pfmet_ey_corr);
      pfMT_corr = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_corr*(1-TMath::Cos(dPhi)));

      //raw_dPhi=dPhiFrom2P( mu_px[0], mu_py[0], raw_pfmet_ex,  raw_pfmet_ey);
      //raw_pfMT = TMath::Sqrt(2*mu_pt[0]*raw_pfmet_pt*(1-TMath::Cos(raw_dPhi)));

     //MET XY Corrections:
     double uncormet;
     double uncormet_phi;
     int year;
     bool isMC;
     isMC = false;
     year = 2018;

     npv  =  analysisTree.primvertex_count;
     uncormet = analysisTree.pfmet_pt;
     uncormet_phi = analysisTree.pfmet_phi;
     int runnb = 0;
     if (isData) runnb = analysisTree.event_run;
     if (!isData) isMC=true;

     //XY Corrections:
     std::pair<double,double> METXY = METXYCorr_Met_MetPhi(uncormet,
                                                           uncormet_phi,
                                                           runnb,
                                                           year,
                                                           isMC,
                                                           npv);
     pfmet_ex_XY = METXY.first;
     pfmet_ey_XY = METXY.second;
     pfmet_pt_XY = sqrt(pfmet_ex_XY*pfmet_ex_XY + pfmet_ey_XY*pfmet_ey_XY);

     if(pfmet_ex_XY==0 && pfmet_ey_XY>0) pfmet_phi_XY = TMath::Pi();

     else if(pfmet_ex_XY==0 && pfmet_ey_XY<0 )pfmet_phi_XY = -TMath::Pi();

     else if(pfmet_ex_XY >0) pfmet_phi_XY = TMath::ATan(pfmet_ey_XY/pfmet_ex_XY);

     else if(pfmet_ex_XY <0&& pfmet_ey_XY>0) pfmet_phi_XY = TMath::ATan(pfmet_ey_XY/pfmet_ex_XY) + TMath::Pi();

     else if(pfmet_ex_XY <0&& pfmet_ey_XY<0) pfmet_phi_XY = TMath::ATan(pfmet_ey_XY/pfmet_ex_XY) - TMath::Pi();

     else pfmet_phi_XY =0;

     //MT for XY Corrections:
     dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex_XY,  pfmet_ey_XY );
	   pfMT_XY = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_XY*(1-TMath::Cos(dPhi)));

     //The same XY correction for JEC corrected met:
     uncormet = analysisTree.pfmetcorr_pt;
     uncormet_phi = analysisTree.pfmetcorr_phi;

     METXY = METXYCorr_Met_MetPhi( uncormet,  uncormet_phi,  runnb,  year,  isMC,  npv);
     pfmet_ex_corr_XY = METXY.first;
     pfmet_ey_corr_XY = METXY.second;

     pfmet_pt_corr_XY = sqrt(pfmet_ex_corr_XY*pfmet_ex_corr_XY + pfmet_ey_corr_XY*pfmet_ey_corr_XY);

     if(pfmet_ex_corr_XY==0 && pfmet_ey_corr_XY>0) pfmet_phi_corr_XY = TMath::Pi();

     else if(pfmet_ex_corr_XY==0 && pfmet_ey_corr_XY<0 )pfmet_phi_corr_XY = -TMath::Pi();

     else if(pfmet_ex_corr_XY >0) pfmet_phi_corr_XY = TMath::ATan(pfmet_ey_corr_XY/pfmet_ex_corr_XY);

     else if(pfmet_ex_corr_XY <0 && pfmet_ey_corr_XY>0) pfmet_phi_corr_XY = TMath::ATan(pfmet_ey_corr_XY/pfmet_ex_corr_XY) + TMath::Pi();

     else if(pfmet_ex_corr_XY <0 && pfmet_ey_corr_XY<0) pfmet_phi_corr_XY = TMath::ATan(pfmet_ey_corr_XY/pfmet_ex_corr_XY) - TMath::Pi();

     else pfmet_phi_corr_XY =0;

     //MT for XY Corrections:
     dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex_corr_XY,  pfmet_ey_corr_XY );
     pfMT_corr_XY = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_corr_XY*(1-TMath::Cos(dPhi)));

     //MET performance:
	   Float_t Utx,Uty;
	   Utx = - mu_px[0] - met_ex;
	   Uty = - mu_py[0] - met_ey;

	   Ut = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	   dPhi=dPhiFrom2P( Utx, Uty, mu_px[0],  mu_py[0] );
	   Utr = -Utx*((mu_py[0]/mu_px[0])/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0]))) + Uty*(1/TMath::Sqrt(1+(mu_py[0]/mu_px[0])*(mu_py[0]/mu_px[0])));
	   Ucol = TMath::Cos(dPhi)* Ut;


    // MET performance ends but not at this place!

    if (!isData) npartons = analysisTree.genparticles_noutgoing;
    all_weight = weight;
    T->Fill();
    selEvents++;
    continue;

    }
    // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }

  //cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;
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

  delete file;

}
