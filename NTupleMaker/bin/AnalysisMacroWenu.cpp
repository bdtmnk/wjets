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

float calculate_dphi(float pfmet_ex_corr, float pfmet_ey_corr){

      float dphi = 0;

      if(pfmet_ex_corr==0 && pfmet_ey_corr>0) dphi = TMath::Pi();

      else if(pfmet_ex_corr==0 && pfmet_ey_corr<0 )dphi= -TMath::Pi();

      else if(pfmet_ex_corr >0) dphi = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr);

      else if(pfmet_ex_corr <0 && pfmet_ey_corr>0) dphi = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr) + TMath::Pi();

      else if(pfmet_ex_corr <0 && pfmet_ey_corr<0) dphi = TMath::ATan(pfmet_ey_corr/pfmet_ex_corr) - TMath::Pi();

      else dphi =0;

      return dphi;
}


int main(int argc, char * argv[]) {

    //Configuration=====================================================
    Config cfg(argv[1]);
    string _year = argv[2];
    string Channel="Wmnu";
    bool usePuppiMET = true;
    bool applySimpleRecoilCorrections = true;
    bool isSampleForRecoilCorrection = false;
    bool ApplyZptReweighting = false;

    //string era = "2018";
    //TODO declare it!!!
    //Global variables:=================================================
    //const bool isRunMETFilters = cfg.get<bool>("RunMETFilters");
    //const bool isRunZpTreWeighting = cfg.get<bool>("RunZpTreWeighting");
    //const bool isRunTTreWeighting = cfg.get<bool>("RunTTreWeighting");

    //Kinematic cuts on electrons=======================================
    const bool isData = cfg.get<bool>("IsData");
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
    const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");

    // electron

    const float  ptElectronCut       = cfg.get<float>("ptElectronCuteltau");
    const double etaElectronCut     = cfg.get<double>("etaElectronCuteltau");
    const double dxyElectronCut     = cfg.get<double>("dxyElectronCuteltau");
    const double dzElectronCut      = cfg.get<double>("dzElectronCuteltau");
    const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCuteltau");
    const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

    // tau
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

  //  veto electrons
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  // dilepton veto
  const float ptDilepElectronCut = cfg.get<float>("ptDilepElectronCuteltau");
  const float etaDilepElectronCut = cfg.get<float>("etaDilepElectronCuteltau");
  const float dxyDilepElectronCut = cfg.get<float>("dxyDilepElectronCuteltau");
  const float dzDilepElectronCut = cfg.get<float>("dzDilepElectronCuteltau");
  const float isoDilepElectronCut = cfg.get<float>("isoDilepElectronCuteltau");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCuteltau");

  // vertex cuts
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");   
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");


  const string dataBaseDir = cfg.get<string>("DataBaseDir");
  const string ElectronIdIsoFile = cfg.get<string>("ElectronIdIsoEffFile");
  const string SingleElectronTriggerFile = cfg.get<string>("ElectrontrigEffFile");

  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");


  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");

  // kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");


  // topSingleElonTriggerFile
  const double dRleptonsCuteltau   = cfg.get<double>("dRleptonsCuteltau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  string TrigLeg = cfg.get<string>("SingleElectronFilterName");
  const float SingleElectronTriggerPtCut = cfg.get<float>("SingleElectronTriggerPtCut");

    // jets
   const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
   const string bTagFile = cfg.get<string>("BTagFile");
   const string bTagEffFile = cfg.get<string>("BTagEffFile");
   const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
   const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
   const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");
   const float jetEtaCut = cfg.get<float>("JetEtaCut");
   const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
   const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
   const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
   const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
   const float btagCut = cfg.get<float>("btagCut");
   const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
   const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");

   TString BTagDiscriminator1(bTagDiscriminator1);
   TString BTagDiscriminator2(bTagDiscriminator2);
   TString BTagDiscriminator3(bTagDiscriminator3);

   const string jec_UncertaintySources = cfg.get<string>("JEC_UncertaintySources"+_year);
   const string jer_resolution = cfg.get<string>("JER_Resolution"+_year);
   const string jer_scalefactor = cfg.get<string>("JER_ScaleFactor"+_year);

    //Vertex distributions filenames and histname
    //Golden Json:
    const string jsonFile = cfg.get<string>("jsonFile"+_year);

    string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

    //const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");
    //Run-lumi selector
    std::vector<Period> periods;


 // JEC uncertainties (NEW) ==========================================================================================================================================
   for (int isrc = 0; isrc < nsrc_Eta0To5; isrc++) {
      const char *name = srcnames_Eta0To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To5[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_Eta0To3; isrc++) {
      const char *name = srcnames_Eta0To3[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta0To3[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_Eta3To5; isrc++) {
      const char *name = srcnames_Eta3To5[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_Eta3To5[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_RelativeBal; isrc++) {
      const char *name = srcnames_RelativeBal[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_RelativeBal[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_RelativeSampleYear; isrc++) {
      const char *name = srcnames_RelativeSampleYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_RelativeSampleYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_EC2; isrc++) {
      const char *name = srcnames_EC2[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_EC2[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_FlavorQCD; isrc++) {
      const char *name = srcnames_FlavorQCD[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_FlavorQCD[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_AbsoluteYear; isrc++) {
      const char *name = srcnames_AbsoluteYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_AbsoluteYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_HFYear; isrc++) {
      const char *name = srcnames_HFYear[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_HFYear[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_EC2Year; isrc++) {
      const char *name = srcnames_EC2Year[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_EC2Year[isrc] = unc;
   }
   for (int isrc = 0; isrc < nsrc_BBEC1Year; isrc++) {
      const char *name = srcnames_BBEC1Year[isrc];
      JetCorrectorParameters const *p = new JetCorrectorParameters(cmsswBase+"/src/"+jec_UncertaintySources, name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      vsrc_BBEC1Year[isrc] = unc;
   }

   map< TString , vector<JetCorrectionUncertainty*> > jec_unc_map = {
      { "jecUncEta0To5"     , vsrc_Eta0To5 },
      { "jecUncEta0To3"     , vsrc_Eta0To3 },
      { "jecUncEta3To5"     , vsrc_Eta3To5 },
      { "jecUncRelativeBal" , vsrc_RelativeBal },
      { "jecUncRelativeSampleYear", vsrc_RelativeSampleYear },
      { "jecUncEC2"         , vsrc_EC2 },
      { "jecUncFlavorQCD"   , vsrc_FlavorQCD },
      { "jecUncAbsoluteYear"   , vsrc_AbsoluteYear },
      { "jecUncHFYear"   , vsrc_HFYear },
      { "jecUncEC2Year"   , vsrc_EC2Year },
      { "jecUncBBEC1Year"   , vsrc_BBEC1Year },
   };

   // JER corrections and uncertainties (NEW) ==========================================================================================================================

   m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/"+jer_resolution));
   m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/"+jer_scalefactor));

   resolution = *m_resolution_from_file;
   resolution_sf = *m_scale_factor_from_file;


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

  //PU Reweighting======================================================
  PileUp * PUofficial = new PileUp();

  //TODO 1 attach PU file to the config=================================

  TString Pileup_file_name_data = "";
  TString Pileup_file_name_mc = "";

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
  TString tagging_efficiencies_file_name = "";

  TFile * fileTagging = new TFile(TString(cmsswBase)+
                        TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_march2018_btageff-all_samp-inc-DeepCSV_medium.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  // Z pt mass weights
/*
  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016.root"); 
  if (fileZMassPtWeights->IsZombie()) {
    std::cout << "File " << TString(cmsswBase) << "src/DesyTauAnalyses/NTupleMaker/data/zpt_weights_2016.root" << "  does not exist!" << std::endl;
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get("zptmass_histo"); 
  if (histZMassPtWeights==NULL) {
    std::cout << " ZMassPT Weights histogram cannot found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << std::endl;
    exit(-1);
  }*/


 // ScaleFactor * SF_HIP; 
//    SF_HIP = new ScaleFactor();
//    SF_HIP->init_ScaleFactor(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker//test/HIP.root");
  // Lepton Scale Factors 


  // Lepton Scale Factors 


cout<<" Initializing iD SF files....."<<endl;

  TFile *fileSF = new TFile("/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/SFforW/egammaEffi.txt_EGM2D.root");
  TH2F * hSF = (TH2F*)fileSF->Get("EGamma_SF2D");
//    ScaleFactor * SF_elIdIso = new ScaleFactor();
//    SF_elIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoFile));
  cout<<" Initializing Trigger SF files....."<<endl;
  ScaleFactor * SF_electronTrigger = new ScaleFactor();
  SF_electronTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(SingleElectronTriggerFile));

  //Recoil path:========================================================
  const string recoilFileName   = cfg.get<string>("RecoilFileName"+year);
  TString RecoilFileName(recoilFileName);

  const string recoilFileNamePuppi   = cfg.get<string>("RecoilFileNamePuppi"+year);
  TString RecoilFileNamePuppi(recoilFileNamePuppi);

  const string metSysFileName   = cfg.get<string>("MetSysFileName"+year);
  TString MetSysFileName(metSysFileName);

  const string metSysFileNamePuppi   = cfg.get<string>("MetSysFileNamePuppi"+year);
  TString MetSysFileNamePuppi(metSysFileNamePuppi);

  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName"+year);
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);

  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);


  const string TopPtWeightsFileName   = cfg.get<string>("TopPtWeightsFileName");
  TString TopPtWeightsFileName(zMassPtWeightsFileName);

  const string TopPtWeightsHistName   = cfg.get<string>("TopPtWeightsHistName");
  TString TopPtWeightsHistName(zMassPtWeightsHistName);


  //Definition of the Corrector:========================================
  kit::MEtSys metSys(MetSysFileName);
  kit::RecoilCorrector recoilMetCorrector(RecoilFileName);
  kit::MEtSys metSysPuppi(MetSysFileNamePuppi);
  kit::RecoilCorrector recoilMetCorrectorPuppi(RecoilFileNamePuppi);

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
  std::string st1,st2;


  TFile * file;
  cout<<"Era"<< era<<endl;
  file = new TFile("./"+era+"/"+TStrName+TString("_DataDriven.root"),"update");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  bool lumi=false;
  std::string dummy;
  while (fileList0 >> dummy) nTotalFiles++;

  //additional branches:=======================================================================================
  SetupTree();

  //Definition of the MT:======================================================================================
  Float_t MT_JetEnUp;
  Float_t MT_JetEnDown;
  Float_t MT_UnclusteredEnUp;
  Float_t MT_UnclusteredEnDown;
  Float_t MT_JetResUp;
  Float_t MT_JetResDown;

  Float_t MT_XY_JetEnUp;
  Float_t MT_XY_JetEnDown;
  Float_t MT_XY_UnclusteredEnUp;
  Float_t MT_XY_UnclusteredEnDown;
  Float_t MT_XY_JetResUp;
  Float_t MT_XY_JetResDown;

  Float_t MTpuppi;
  Float_t MTpuppi_JetEnUp;
  Float_t MTpuppi_JetEnDown;
  Float_t MTpuppi_UnclusteredEnUp;
  Float_t MTpuppi_UnclusteredEnDown;
  Float_t MTpuppi_JetResUp;
  Float_t MTpuppi_JetResDown;

  Float_t MT_JetEnUp_smeared;
  Float_t MT_JetEnDown_smeared;
  Float_t MT_UnclusteredEnUp_smeared;
  Float_t MT_UnclusteredEnDown_smeared;
  Float_t MT_JetResUp_smeared;
  Float_t MT_JetResDown_smeared;

  T->Branch("MT_JetEnUp", &MT_JetEnUp, "MT_JetEnUp/F");
  T->Branch("MT_JetEnDown", &MT_JetEnDown, "MT_JetEnDown/F");
  T->Branch("MT_UnclusteredEnUp", &MT_UnclusteredEnUp, "MT_UnclusteredEnUp/F");
  T->Branch("MT_UnclusteredEnDown", &MT_UnclusteredEnDown, "MT_UnclusteredEnDown/F");
  T->Branch("MT_JetResUp", &MT_JetResUp, "MT_JetResUp/F");
  T->Branch("MT_JetResDown", &MT_JetResDown, "MT_JetResDown/F");

  T->Branch("MT_JetEnUp_smeared", &MT_JetEnUp_smeared, "MT_JetEnUp_smeared/F");
  T->Branch("MT_JetEnDown_smeared", &MT_JetEnDown_smeared, "MT_JetEnDown_smeared/F");
  T->Branch("MT_UnclusteredEnUp_smeared", &MT_UnclusteredEnUp_smeared, "MT_UnclusteredEnUp_smeared/F");
  T->Branch("MT_UnclusteredEnDown_smeared", &MT_UnclusteredEnDown_smeared, "MT_UnclusteredEnDown_smeared/F");
  T->Branch("MT_JetResUp_smeared", &MT_JetResUp_smeared, "MT_JetResUp_smeared/F");
  T->Branch("MT_JetResDown_smeared", &MT_JetResDown_smeared, "MT_JetResDown_smeared/F");

  T->Branch("MT_XY_JetEnUp", &MT_XY_JetEnUp, "MT_XY_JetEnUp/F");
  T->Branch("MT_XY_JetEnDown", &MT_XY_JetEnDown, "MT_XY_JetEnDown/F");
  T->Branch("MT_XY_UnclusteredEnUp", &MT_XY_UnclusteredEnUp, "MT_XY_UnclusteredEnUp/F");
  T->Branch("MT_XY_UnclusteredEnDown", &MT_XY_UnclusteredEnDown, "MT_XY_UnclusteredEnDown/F");
  T->Branch("MT_XY_JetResUp", &MT_XY_JetResUp, "MT_XY_JetResUp/F");
  T->Branch("MT_XY_JetResDown", &MT_XY_JetResDown, "MT_XY_JetResDown/F");

  T->Branch("MTpuppi", &MTpuppi, "MTpuppi/F");
  T->Branch("MTpuppi_JetEnUp", &MTpuppi_JetEnUp, "MTpuppi_JetEnUp/F");
  T->Branch("MTpuppi_JetEnDown", &MTpuppi_JetEnDown, "MTpuppi_JetEnDown/F");
  T->Branch("MTpuppi_UnclusteredEnUp", &MTpuppi_UnclusteredEnUp, "MTpuppi_UnclusteredEnUp/F");
  T->Branch("MTpuppi_UnclusteredEnDown", &MTpuppi_UnclusteredEnDown, "MTpuppi_UnclusteredEnDown/F");
  T->Branch("MTpuppi_JetResUp", &MTpuppi_JetResUp, "MTpuppi_JetResUp/F");
  T->Branch("MTpuppi_JetResDown", &MTpuppi_JetResDown, "MTpuppi_JetResDown/F");

  Float_t         raw_met_pt;
  Float_t         raw_met_phi;

  Float_t         puppiMT;
  Float_t         pfMT;
  Float_t         pfMT_corr;
  Float_t         pfMT_XY;
  Float_t         pfMT_corr_XY;
  Float_t         MT_smeared;

  T->Branch("puppiMT", &puppiMT, "MTpuppi/F");
  T->Branch("pfMT", &pfMT, "pfMT/F");
  T->Branch("pfMT_corr", &pfMT_corr, "pfMT_corr/F");
  T->Branch("pfMT_XY", &pfMT_XY, "pfMT_XY/F");
  T->Branch("pfMT_corr_XY", &pfMT_corr_XY, "pfMT_corr_XY/F");
  T->Branch("MT_smeared", &MT_smeared, "MT_smeared/F");


  //Puppi Met declaration:========================================================================================

  Float_t      puppi_ex;
  Float_t      puppi_ey;
  Float_t      puppi_ez;
  Float_t      puppi_pt;
  Float_t      puppi_phi;

  Float_t      puppi_ex_JetEnUp;
  Float_t      puppi_ey_JetEnUp;
  Float_t      puppi_ex_JetEnDown;
  Float_t      puppi_ey_JetEnDown;
  Float_t      puppi_ex_UnclusteredEnUp;
  Float_t      puppi_ey_UnclusteredEnUp;
  Float_t      puppi_ex_UnclusteredEnDown;
  Float_t      puppi_ey_UnclusteredEnDown;
  Float_t      puppi_ex_JetResUp;
  Float_t      puppi_ey_JetResUp
  Float_t      puppi_ex_JetResDown;
  Float_t      puppi_ey_JetResDown;

  //Systematic for the Puppi Pt:
  Float_t puppi_pt_JetEnUp;
  Float_t puppi_pt_JetEnDown;
  Float_t puppi_pt_UnclusteredEnUp;
  Float_t puppi_pt_UnclusteredEnDown;
  Float_t puppi_pt_JetResUp;
  Float_t puppi_pt_JetResDown;

  //Systematic for the Puppi Phi
  Float_t puppi_phi_JetEnUp;
  Float_t puppi_phi_JetEnDown;
  Float_t puppi_phi_UnclusteredEnUp;
  Float_t puppi_phi_UnclusteredEnDown;
  Float_t puppi_phi_JetResUp;
  Float_t puppi_phi_JetResDown;

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

  T->Branch("puppi_pt_JetEnUp", &puppi_pt_JetEnUp, "puppi_pt_JetEnUp/F");
  T->Branch("puppi_pt_JetEnDown", &puppi_pt_JetEnDown, "puppi_pt_JetEnDown/F");
  T->Branch("puppi_pt_UnclusteredEnUp", &puppi_pt_UnclusteredEnUp, "puppi_pt_UnclusteredEnUp/F");
  T->Branch("puppi_pt_UnclusteredEnDown", &puppi_pt_UnclusteredEnDown, "puppi_pt_UnclusteredEnDown/F");
  T->Branch("puppi_pt_JetResUp", &puppi_pt_JetResUp, "puppi_pt_JetResUp/F");
  T->Branch("puppi_pt_JetResDown", &puppi_pt_JetResDown, "puppi_pt_JetResDown/F");

  T->Branch("puppi_phi_JetEnUp", &puppi_phi_JetEnUp, "puppi_phi_JetEnUp/F");
  T->Branch("puppi_phi_JetEnDown", &puppi_phi_JetEnDown, "puppi_phi_JetEnDown/F");
  T->Branch("puppi_phi_UnclusteredEnUp", &puppi_phi_UnclusteredEnUp, "puppi_phi_UnclusteredEnUp/F");
  T->Branch("puppi_phi_UnclusteredEnDown", &puppi_phi_UnclusteredEnDown, "puppi_phi_UnclusteredEnDown/F");
  T->Branch("puppi_phi_JetResUp", &puppi_phi_JetResUp, "puppi_phi_JetResUp/F");
  T->Branch("puppi_phi_JetResDown", &puppi_phi_JetResDown, "puppi_phi_JetResDown/F");

  //PF Collections:=============================================================================================
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

   //Corrected MET==============================================================================================
   Float_t      pfmet_pt_corr;
   Float_t      pfmet_phi_corr;
   Float_t      pfmet_ex_corr;
   Float_t      pfmet_ey_corr;

   Float_t      pfmet_ex_JetEnUp;
   Float_t      pfmet_ey_JetEnUp;

   Float_t      pfmet_ex_JetEnDown;
   Float_t      pfmet_ey_JetEnDown;

   Float_t      pfmet_ex_UnclusteredEnUp;
   Float_t      pfmet_ey_UnclusteredEnUp;
   Float_t      pfmet_ex_UnclusteredEnDown;
   Float_t      pfmet_ey_UnclusteredEnDown;
   Float_t      pfmet_ex_JetResUp;
   Float_t      pfmet_ey_JetResUp;
   Float_t      pfmet_ex_JetResDown;
   Float_t      pfmet_ey_JetResDown;
   Float_t      pfmet_pt_JetEnUp;
   Float_t      pfmet_pt_JetEnDown;
   Float_t      pfmet_pt_UnclusteredEnUp;
   Float_t      pfmet_pt_UnclusteredEnDown;
   Float_t      pfmet_pt_JetResUp;
   Float_t      pfmet_pt_JetResDown;
   Float_t      pfmet_phi_JetEnUp;
   Float_t      pfmet_phi_JetEnDown;
   Float_t      pfmet_phi_UnclusteredEnUp;
   Float_t      pfmet_phi_UnclusteredEnDown;
   Float_t      pfmet_phi_JetResUp;
   Float_t      pfmet_phi_JetResDown;


    T->Branch("pfmet_pt_corr", &pfmet_pt_corr, "pfmet_pt_corr/F");
    T->Branch("pfmet_phi_corr", &pfmet_phi_corr, "pfmet_phi_corr/F");
    T->Branch("pfmet_ex_corr", &pfmet_ex_corr, "pfmet_ex_corr/F");
    T->Branch("pfmet_ey_corr", &pfmet_ey_corr, "pfmet_ey_corr/F");
    //JetEnUp
    T->Branch("pfmet_ex_JetEnUp", &pfmet_ex_JetEnUp, "pfmet_ex_JetEnUp/F");
    T->Branch("pfmet_ey_JetEnUp", &pfmet_ey_JetEnUp, "pfmet_ey_JetEnUp/F");
    T->Branch("pfmet_ex_JetEnDown", &pfmet_ex_JetEnDown, "pfmet_ex_JetEnDown/F");
    T->Branch("pfmet_ey_JetEnDown", &pfmet_ey_JetEnDown, "pfmet_ey_JetEnDown/F");
    //JetRes
    T->Branch("pfmet_ex_JetResUp", &pfmet_ex_JetResUp, "pfmet_ex_JetResUp/F");
    T->Branch("pfmet_ey_JetResUp", &pfmet_ey_JetResUp, "pfmet_ey_JetResUp/F");
    T->Branch("pfmet_ex_JetResDown", &pfmet_ex_JetResDown, "pfmet_ex_JetResDown/F");
    T->Branch("pfmet_ey_JetResDown", &pfmet_ey_JetResDown, "pfmet_ey_JetResDown/F");
    //Unclustered
    T->Branch("pfmet_ex_UnclusteredEnUp", &pfmet_ex_UnclusteredEnUp, "pfmet_ex_UnclusteredEnUp/F");
    T->Branch("pfmet_ey_UnclusteredEnUp", &pfmet_ey_UnclusteredEnUp, "pfmet_ey_UnclusteredEnUp/F");
    T->Branch("pfmet_ex_UnclusteredEnDown", &pfmet_ex_UnclusteredEnDown, "pfmet_ex_UnclusteredEnDown/F");
    T->Branch("pfmet_ey_UnclusteredEnDown", &pfmet_ey_UnclusteredEnDown, "pfmet_ey_UnclusteredEnDown/F");
    //Unc for the pfmet pt:
    T->Branch("pfmet_pt_JetEnUp", &pfmet_pt_JetEnUp, "pfmet_pt_JetEnUp/F");
    T->Branch("pfmet_pt_JetEnDown", &pfmet_pt_JetEnDown, "pfmet_pt_JetEnDown/F");
    T->Branch("pfmet_pt_UnclusteredEnUp", &pfmet_pt_UnclusteredEnUp, "pfmet_pt_UnclusteredEnUp/F");
    T->Branch("pfmet_pt_UnclusteredEnDown", &pfmet_pt_UnclusteredEnDown, "pfmet_pt_UnclusteredEnDown/F");
    T->Branch("pfmet_pt_JetResUp", &pfmet_pt_JetResUp, "pfmet_pt_JetResUp/F");
    T->Branch("pfmet_pt_JetResDown", &pfmet_pt_JetResDown, "pfmet_pt_JetResDown/F");
    //Unc for the pfmet phi:
    T->Branch("pfmet_phi_JetEnUp", &pfmet_phi_JetEnUp, "pfmet_phi_JetEnUp/F");
    T->Branch("pfmet_phi_JetEnDown", &pfmet_phi_JetEnDown, "pfmet_phi_JetEnDown/F");
    T->Branch("pfmet_phi_UnclusteredEnUp", &pfmet_phi_UnclusteredEnUp, "pfmet_phi_UnclusteredEnUp/F");
    T->Branch("pfmet_phi_UnclusteredEnDown", &pfmet_phi_UnclusteredEnDown, "pfmet_phi_UnclusteredEnDown/F");
    T->Branch("pfmet_phi_JetResUp", &pfmet_phi_JetResUp, "pfmet_phi_JetResUp/F");
    T->Branch("pfmet_phi_JetResDown", &pfmet_phi_JetResDown, "pfmet_phi_JetResDown/F");

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

    //MEt Systematics===================================================
    Float_t pfmetcorr_ex_JetEnUp;
    Float_t pfmetcorr_ey_JetEnUp;
    Float_t pfmetcorr_pt_JetEnUp;
    Float_t pfmetcorr_phi_JetEnUp;

    Float_t pfmetcorr_ex_JetEnDown;
    Float_t pfmetcorr_ey_JetEnDown;
    Float_t pfmetcorr_pt_JetEnDown;
    Float_t pfmetcorr_phi_JetEnDown;

    Float_t pfmetcorr_ex_UnclusteredEnUp;
    Float_t pfmetcorr_ey_UnclusteredEnUp;
    Float_t pfmetcorr_pt_UnclusteredEnUp;
    Float_t pfmetcorr_phi_UnclusteredEnUp;

    Float_t pfmetcorr_ex_UnclusteredEnDown;
    Float_t pfmetcorr_ey_UnclusteredEnDown;
    Float_t pfmetcorr_pt_UnclusteredEnDown;
    Float_t pfmetcorr_phi_UnclusteredEnDown;

    Float_t pfmetcorr_ex_JetResUp;
    Float_t pfmetcorr_ey_JetResUp;
    Float_t pfmetcorr_pt_JetResUp;
    Float_t pfmetcorr_phi_JetResUp;

    Float_t pfmetcorr_ex_JetResDown;
    Float_t pfmetcorr_ey_JetResDown;
    Float_t pfmetcorr_pt_JetResDown;
    Float_t pfmetcorr_phi_JetResDown;

    T->Branch("pfmetcorr_ex_JetEnUp", &pfmetcorr_ex_JetEnUp, "pfmetcorr_ex_JetEnUp/F");
    T->Branch("pfmetcorr_ey_JetEnUp", &pfmetcorr_ey_JetEnUp, "pfmetcorr_ey_JetEnUp/F");
    T->Branch("pfmetcorr_pt_JetEnUp", &pfmetcorr_pt_JetEnUp, "pfmetcorr_pt_JetEnUp/F");
    T->Branch("pfmetcorr_phi_JetEnUp", &pfmetcorr_phi_JetEnUp, "pfmetcorr_phi_JetEnUp/F");

    T->Branch("pfmetcorr_ex_JetEnDown", &pfmetcorr_ex_JetEnDown, "pfmetcorr_ex_JetEnDown/F");
    T->Branch("pfmetcorr_ey_JetEnDown", &pfmetcorr_ey_JetEnDown, "pfmetcorr_ey_JetEnDown/F");
    T->Branch("pfmetcorr_pt_JetEnDown", &pfmetcorr_pt_JetEnDown, "pfmetcorr_pt_JetEnDown/F");
    T->Branch("pfmetcorr_phi_JetEnDown", &pfmetcorr_phi_JetEnDown, "pfmetcorr_phi_JetEnDown/F");

    T->Branch("pfmetcorr_ex_UnclusteredEnUp", &pfmetcorr_ex_UnclusteredEnUp, "pfmetcorr_ex_UnclusteredEnUp/F");
    T->Branch("pfmetcorr_ey_UnclusteredEnUp", &pfmetcorr_ey_UnclusteredEnUp, "pfmetcorr_ey_UnclusteredEnUp/F");
    T->Branch("pfmetcorr_pt_UnclusteredEnUp", &pfmetcorr_pt_UnclusteredEnUp, "pfmetcorr_pt_UnclusteredEnUp/F");
    T->Branch("pfmetcorr_phi_UnclusteredEnUp", &pfmetcorr_phi_UnclusteredEnUp, "pfmetcorr_phi_UnclusteredEnUp/F");

    T->Branch("pfmetcorr_ex_UnclusteredEnDown", &pfmetcorr_ex_UnclusteredEnDown, "pfmetcorr_ex_UnclusteredEnDown/F");
    T->Branch("pfmetcorr_ey_UnclusteredEnDown", &pfmetcorr_ey_UnclusteredEnDown, "pfmetcorr_ey_UnclusteredEnDown/F");
    T->Branch("pfmetcorr_pt_UnclusteredEnDown", &pfmetcorr_pt_UnclusteredEnDown, "pfmetcorr_pt_UnclusteredEnDown/F");
    T->Branch("pfmetcorr_phi_UnclusteredEnDown", &pfmetcorr_phi_UnclusteredEnDown, "pfmetcorr_phi_UnclusteredEnDown/F");

    T->Branch("pfmetcorr_ex_JetResUp", &pfmetcorr_ex_JetResUp, "pfmetcorr_ex_JetResUp/F");
    T->Branch("pfmetcorr_ey_JetResUp", &pfmetcorr_ey_JetResUp, "pfmetcorr_ey_JetResUp/F");
    T->Branch("pfmetcorr_pt_JetResUp", &pfmetcorr_pt_JetResUp, "pfmetcorr_pt_JetResUp/F");
    T->Branch("pfmetcorr_phi_JetResUp", &pfmetcorr_phi_JetResUp, "pfmetcorr_phi_JetResUp/F");

    T->Branch("pfmetcorr_ex_JetResDown", &pfmetcorr_ex_JetResDown, "pfmetcorr_ex_JetResDown/F");
    T->Branch("pfmetcorr_ey_JetResDown", &pfmetcorr_ey_JetResDown, "pfmetcorr_ey_JetResDown/F");
    T->Branch("pfmetcorr_pt_JetResDown", &pfmetcorr_pt_JetResDown, "pfmetcorr_pt_JetResDown/F");
    T->Branch("pfmetcorr_phi_JetResDown", &pfmetcorr_phi_JetResDown, "pfmetcorr_phi_JetResDown/F");

    Float_t pfmet_ex_XY_JetEnUp;
    Float_t pfmet_ey_XY_JetEnUp;
    Float_t pfmet_ex_XY_JetEnDown;
    Float_t pfmet_ey_XY_JetEnDown;
    Float_t pfmet_ex_XY_UnclusteredEnUp;
    Float_t pfmet_ey_XY_UnclusteredEnUp;
    Float_t pfmet_ex_XY_UnclusteredEnDown;
    Float_t pfmet_ey_XY_UnclusteredEnDown;
    Float_t pfmet_ex_XY_JetResUp;
    Float_t pfmet_ey_XY_JetResUp;
    Float_t pfmet_ex_XY_JetResDown;
    Float_t pfmet_ey_XY_JetResDown;

    Float_t pfmet_pt_XY_JetEnUp;
    Float_t pfmet_phi_XY_JetEnUp;
    Float_t pfmet_pt_XY_JetEnDown;
    Float_t pfmet_phi_XY_JetEnDown;
    Float_t pfmet_pt_XY_UnclusteredEnUp;
    Float_t pfmet_phi_XY_UnclusteredEnUp;
    Float_t pfmet_pt_XY_UnclusteredEnDown;
    Float_t pfmet_phi_XY_UnclusteredEnDown;
    Float_t pfmet_pt_XY_JetResUp;
    Float_t pfmet_phi_XY_JetResUp;
    Float_t pfmet_pt_XY_JetResDown;
    Float_t pfmet_phi_XY_JetResDown;


    T->Branch("pfmet_ex_XY_JetEnUp", &pfmet_ex_XY_JetEnUp, "pfmet_ex_XY_JetEnUp/F");
    T->Branch("pfmet_ey_XY_JetEnUp", &pfmet_ey_XY_JetEnUp, "pfmet_ey_XY_JetEnUp/F");
    T->Branch("pfmet_ex_XY_JetEnDown", &pfmet_ex_XY_JetEnDown, "pfmet_ex_XY_JetEnDown/F");
    T->Branch("pfmet_ey_XY_JetEnDown", &pfmet_ey_XY_JetEnDown, "pfmet_ey_XY_JetEnDown/F");
    T->Branch("pfmet_ex_XY_UnclusteredEnUp", &pfmet_ex_XY_UnclusteredEnUp, "pfmet_ex_XY_UnclusteredEnUp/F");
    T->Branch("pfmet_ey_XY_UnclusteredEnUp", &pfmet_ey_XY_UnclusteredEnUp, "pfmet_ey_XY_UnclusteredEnUp/F");
    T->Branch("pfmet_ex_XY_UnclusteredEnDown", &pfmet_ex_XY_UnclusteredEnDown, "pfmet_ex_XY_UnclusteredEnDown/F");
    T->Branch("pfmet_ey_XY_UnclusteredEnDown", &pfmet_ey_XY_UnclusteredEnDown, "pfmet_ey_XY_UnclusteredEnDown/F");
    T->Branch("pfmet_ex_XY_JetResUp", &pfmet_ex_XY_JetResUp, "pfmet_ex_XY_JetResUp/F");
    T->Branch("pfmet_ey_XY_JetResUp", &pfmet_ey_XY_JetResUp, "pfmet_ey_XY_JetResUp/F");
    T->Branch("pfmet_ex_XY_JetResDown", &pfmet_ex_XY_JetResDown, "pfmet_ex_XY_JetResDown/F");
    T->Branch("pfmet_ey_XY_JetResDown", &pfmet_ey_XY_JetResDown, "pfmet_ey_XYs_JetResDown/F");

    T->Branch("pfmet_pt_XY_JetEnUp", &pfmet_pt_XY_JetEnUp, "pfmet_pt_XY_JetEnUp/F");
    T->Branch("pfmet_phi_XY_JetEnUp", &pfmet_phi_XY_JetEnUp, "pfmet_phi_XY_JetEnUp/F");
    T->Branch("pfmet_pt_XY_JetEnDown", &pfmet_pt_XY_JetEnDown, "pfmet_pt_XY_JetEnDown/F");
    T->Branch("pfmet_phi_XY_JetEnDown", &pfmet_phi_XY_JetEnDown, "pfmet_phi_XY_JetEnDown/F");
    T->Branch("pfmet_pt_XY_UnclusteredEnUp", &pfmet_pt_XY_UnclusteredEnUp, "pfmet_pt_XY_UnclusteredEnUp/F");
    T->Branch("pfmet_phi_XY_UnclusteredEnUp", &pfmet_phi_XY_UnclusteredEnUp, "pfmet_phi_XY_UnclusteredEnUp/F");
    T->Branch("pfmet_pt_XY_UnclusteredEnDown", &pfmet_pt_XY_UnclusteredEnDown, "pfmet_pt_XY_UnclusteredEnDown/F");
    T->Branch("pfmet_phi_XY_UnclusteredEnDown", &pfmet_phi_XY_UnclusteredEnDown, "pfmet_phi_XY_UnclusteredEnDown/F");
    T->Branch("pfmet_pt_XY_JetResUp", &pfmet_pt_XY_JetResUp, "pfmet_pt_XY_JetResUp/F");
    T->Branch("pfmet_phi_XY_JetResUp", &pfmet_phi_XY_JetResUp, "pfmet_phi_XY_JetResUp/F");
    T->Branch("pfmet_pt_XY_JetResDown", &pfmet_pt_XY_JetResDown, "pfmet_pt_XY_JetResDown/F");
    T->Branch("pfmet_phi_XY_JetResDown", &pfmet_phi_XY_JetResDown, "pfmet_phi_XYs_JetResDown/F");

    //Smeared:===============================================================================================
    Float_t      pfmet_ex_smeared;
    Float_t      pfmet_ey_smeared;
    Float_t      pfmet_ez_smeared;

    Float_t      pfmet_pt_smeared;
    Float_t      pfmet_phi_smeared;

    Float_t      pfmet_ex_JetEnUp_smeared;
    Float_t      pfmet_ey_JetEnUp_smeared;
    Float_t      pfmet_ex_JetEnDown_smeared;
    Float_t      pfmet_ey_JetEnDown_smeared;
    Float_t      pfmet_ex_UnclusteredEnUp_smeared;
    Float_t      pfmet_ey_UnclusteredEnUp_smeared;
    Float_t      pfmet_ex_UnclusteredEnDown_smeared;
    Float_t      pfmet_ey_UnclusteredEnDown_smeared;
    Float_t      pfmet_ex_JetResUp_smeared;
    Float_t      pfmet_ey_JetResUp_smeared;
    Float_t      pfmet_ex_JetResDown_smeared;
    Float_t      pfmet_ey_JetResDown_smeared;
    Float_t      pfmet_pt_JetEnUp_smeared;
    Float_t      pfmet_pt_JetEnDown_smeared;
    Float_t      pfmet_pt_UnclusteredEnUp_smeared;
    Float_t      pfmet_pt_UnclusteredEnDown_smeared;
    Float_t      pfmet_pt_JetResUp_smeared;
    Float_t      pfmet_pt_JetResDown_smeared;

    //Add Corrected MET:======================================================================================
   Float_t      pfmet_ex_corr_smeared;
   Float_t      pfmet_ey_corr_smeared;
   Float_t      pfmet_pt_corr_smeared;
   Float_t      pfmet_phi_corr_smeared;

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
//    T->Branch("sum_pfjet_energycorr_l2relative", &sum_pfjet_energycorr_l2relative, "sum_pfjet_energycorr_l2relative//F");
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


    Float_t PUweight;
    T->Branch("PUweight", &PUweight, "PUweight/F");

    Float_t raw_dPhi;
    Float_t raw_pfMT;
    Float_t raw_pfmet_ex;
    Float_t raw_pfmet_ey;
    Float_t raw_dPhi_XY;
    Float_t raw_pfMT_XY;
    Float_t raw_pfmet_ex_XY;
    Float_t raw_pfmet_ey_XY;

    T->Branch("raw_dPhi",  &raw_dPhi, "raw_dPhi/F");
    T->Branch("raw_pfMT", &raw_pfMT, "raw_pfMT/F");
    T->Branch("raw_pfmet_ex", &raw_pfmet_ex, "raw_pfmet_ex/F");
    T->Branch("raw_pfmet_ey", &raw_pfmet_ey, "raw_pfmet_ey/F");
    T->Branch("raw_dPhi_XY",  &raw_dPhi_XY, "raw_dPhi_XY/F");
    T->Branch("raw_pfMT_XY", &raw_pfMT_XY, "raw_pfMT_XY/F");
    T->Branch("raw_pfmet_ex_XY", &raw_pfmet_ex_XY, "raw_pfmet_ex_XY/F");
    T->Branch("raw_pfmet_ey_XY", &raw_pfmet_ey_XY, "raw_pfmet_ey_XY/F");


    float zptmassweight;
    T->Branch("zptmassweight",  &zptmassweight, "zptmassweight/F");

        /////MET performance

    Float_t Ut;
    Float_t Utr;
    Float_t Ucol;
    Float_t Ut_JetEnUp;
    Float_t Utr_JetEnUp;
    Float_t Ucol_JetEnUp;
    Float_t Ut_JetEnDown;
    Float_t Utr_JetEnDown;
    Float_t Ucol_JetEnDown;
    Float_t Ut_JetResDown;
    Float_t Utr_JetResDown;
    Float_t Ucol_JetResDown;
    Float_t Ut_JetResUp;
    Float_t Utr_JetResUp;
    Float_t Ucol_JetResUp;
    Float_t Ut_UnclusteredEnUp;
    Float_t Utr_UnclusteredEnUp;
    Float_t Ucol_UnclusteredEnUp;
    Float_t Ut_UnclusteredEnDown;
    Float_t Utr_UnclusteredEnDown;
    Float_t Ucol_UnclusteredEnDown;
    Float_t UtMC;
    Float_t UtrMC;
    Float_t UcolMC;


    T->Branch("Ut", &Ut, "Ut/F");
    T->Branch("Utr", &Utr, "Utr/F");
    T->Branch("Ucol", &Ucol, "Ucol/F");
    T->Branch("Ut_JetEnUp", &Ut_JetEnUp, "Ut_JetEnUp/F");
    T->Branch("Utr_JetEnUp", &Utr_JetEnUp, "Utr_JetEnUp/F");
    T->Branch("Ucol_JetEnUp", &Ucol_JetEnUp, "Ucol_JetEnUp/F");
    T->Branch("Ut_JetEnDown", &Ut_JetEnDown, "Ut_JetEnDown/F");
    T->Branch("Utr_JetEnDown", &Utr_JetEnDown, "Utr_JetEnDown/F");
    T->Branch("Ucol_JetEnDown", &Ucol_JetEnDown, "Ucol_JetEnDown/F");
    T->Branch("Ut_JetResDown", &Ut_JetResDown, "Ut_JetResDown/F");
    T->Branch("Utr_JetResDown", &Utr_JetResDown, "Utr_JetResDown/F");
    T->Branch("Ucol_JetResDown", &Ucol_JetResDown, "Ucol_JetResDown/F");
    T->Branch("Ut_JetResUp", &Ut_JetResUp, "Ut_JetResUp/F");
    T->Branch("Utr_JetResUp", &Utr_JetResUp, "Utr_JetResUp/F");
    T->Branch("Ucol_JetResUp", &Ucol_JetResUp, "Ucol_JetResUp/F");
    T->Branch("Ut_UnclusteredEnUp", &Ut_UnclusteredEnUp, "Ut_UnclusteredEnUp/F");
    T->Branch("Utr_UnclusteredEnUp", &Utr_UnclusteredEnUp, "Utr_UnclusteredEnUp/F");
    T->Branch("Ucol_UnclusteredEnUp", &Ucol_UnclusteredEnUp, "Ucol_UnclusteredEnUp/F");
    T->Branch("Ut_UnclusteredEnDown", &Ut_UnclusteredEnDown, "Ut_UnclusteredEnDown/F");
    T->Branch("Utr_UnclusteredEnDown", &Utr_UnclusteredEnDown, "Utr_UnclusteredEnDown/F");
    T->Branch("Ucol_UnclusteredEnDown", &Ucol_UnclusteredEnDown, "Ucol_UnclusteredEnDown/F");
    T->Branch("UtMC", &UtMC, "UtMC/F");
    T->Branch("UtrMC", &UtrMC, "UtrMC/F");
    T->Branch("UcolMC", &UcolMC, "UcolMC/F");


    Float_t Ut_puppi;
    Float_t Utr_puppi;
    Float_t Ucol_puppi;
    Float_t Ut_puppi_JetEnUp;
    Float_t Utr_puppi_JetEnUp;
    Float_t Ucol_puppi_JetEnUp;
    Float_t Ut_puppi_JetEnDown;
    Float_t Utr_puppi_JetEnDown;
    Float_t Ucol_puppi_JetEnDown;
    Float_t Ut_puppi_JetResDown;
    Float_t Utr_puppi_JetResDown;
    Float_t Ucol_puppi_JetResDown;
    Float_t Ut_puppi_JetResUp;
    Float_t Utr_puppi_JetResUp;
    Float_t Ucol_puppi_JetResUp;
    Float_t Ut_puppi_UnclusteredEnUp;
    Float_t Utr_puppi_UnclusteredEnUp;
    Float_t Ucol_puppi_UnclusteredEnUp;
    Float_t Ut_puppi_UnclusteredEnDown;
    Float_t Utr_puppi_UnclusteredEnDown;
    Float_t Ucol_puppi_UnclusteredEnDown;

    T->Branch("Ut_puppi", &Ut_puppi, "Ut_puppi/F");
    T->Branch("Utr_puppi", &Utr_puppi, "Utr_puppi/F");
    T->Branch("Ucol_puppi", &Ucol_puppi, "Ucol_puppi/F");
    T->Branch("Ut_puppi_JetEnUp", &Ut_puppi_JetEnUp, "Ut_puppi_JetEnUp/F");
    T->Branch("Utr_puppi_JetEnUp", &Utr_puppi_JetEnUp, "Utr_puppi_JetEnUp/F");
    T->Branch("Ucol_puppi_JetEnUp", &Ucol_puppi_JetEnUp, "Ucol_puppi_JetEnUp/F");
    T->Branch("Ut_puppi_JetEnDown", &Ut_puppi_JetEnDown, "Ut_puppi_JetEnDown/F");
    T->Branch("Utr_puppi_JetEnDown", &Utr_puppi_JetEnDown, "Utr_puppi_JetEnDown/F");
    T->Branch("Ucol_puppi_JetEnDown", &Ucol_puppi_JetEnDown, "Ucol_puppi_JetEnDown/F");
    T->Branch("Ut_puppi_JetResDown", &Ut_puppi_JetResDown, "Ut_puppi_JetResDown/F");
    T->Branch("Utr_puppi_JetResDown", &Utr_puppi_JetResDown, "Utr_puppi_JetResDown/F");
    T->Branch("Ucol_puppi_JetResDown", &Ucol_puppi_JetResDown, "Ucol_puppi_JetResDown/F");
    T->Branch("Ut_puppi_JetResUp", &Ut_puppi_JetResUp, "Ut_puppi_JetResUp/F");
    T->Branch("Utr_puppi_JetResUp", &Utr_puppi_JetResUp, "Utr_puppi_JetResUp/F");
    T->Branch("Ucol_puppi_JetResUp", &Ucol_puppi_JetResUp, "Ucol_puppi_JetResUp/F");
    T->Branch("Ut_puppi_UnclusteredEnUp", &Ut_puppi_UnclusteredEnUp, "Ut_puppi_UnclusteredEnUp/F");
    T->Branch("Utr_puppi_UnclusteredEnUp", &Utr_puppi_UnclusteredEnUp, "Utr_puppi_UnclusteredEnUp/F");
    T->Branch("Ucol_puppi_UnclusteredEnUp", &Ucol_puppi_UnclusteredEnUp, "Ucol_puppi_UnclusteredEnUp/F");
    T->Branch("Ut_puppi_UnclusteredEnDown", &Ut_puppi_UnclusteredEnDown, "Ut_puppi_UnclusteredEnDown/F");
    T->Branch("Utr_puppi_UnclusteredEnDown", &Utr_puppi_UnclusteredEnDown, "Utr_puppi_UnclusteredEnDown/F");
    T->Branch("Ucol_puppi_UnclusteredEnDown", &Ucol_puppi_UnclusteredEnDown, "Ucol_puppi_UnclusteredEnDown/F");


    //Smeared
    Float_t Ut_smeared;
    Float_t Utr_smeared;
    Float_t Ucol_smeared;

    Float_t Ut_smeared_JetEnUp;
    Float_t Utr_smeared_JetEnUp;
    Float_t Ucol_smeared_JetEnUp;
    Float_t Ut_smeared_JetEnDown;
    Float_t Utr_smeared_JetEnDown;
    Float_t Ucol_smeared_JetEnDown;
    Float_t Ut_smeared_JetResDown;
    Float_t Utr_smeared_JetResDown;
    Float_t Ucol_smeared_JetResDown;
    Float_t Ut_smeared_JetResUp;
    Float_t Utr_smeared_JetResUp;
    Float_t Ucol_smeared_JetResUp;

    Float_t Ut_smeared_UnclusteredEnUp;
    Float_t Utr_smeared_UnclusteredEnUp;
    Float_t Ucol_smeared_UnclusteredEnUp;
    Float_t Ut_smeared_UnclusteredEnDown;
    Float_t Utr_smeared_UnclusteredEnDown;
    Float_t Ucol_smeared_UnclusteredEnDown;

    T->Branch("Ut_smeared", &Ut_smeared, "Ut_smeared/F");
    T->Branch("Utr_smeared", &Utr_smeared, "Utr_smeared/F");
    T->Branch("Ucol_smeared", &Ucol_smeared, "Ucol_smeared/F");
    T->Branch("Ut_smeared_JetEnUp", &Ut_smeared_JetEnUp, "Ut_smeared_JetEnUp/F");
    T->Branch("Utr_smeared_JetEnUp", &Utr_smeared_JetEnUp, "Utr_smeared_JetEnUp/F");
    T->Branch("Ucol_smeared_JetEnUp", &Ucol_smeared_JetEnUp, "Ucol_smeared_JetEnUp/F");
    T->Branch("Ut_smeared_JetEnDown", &Ut_smeared_JetEnDown, "Ut_smeared_JetEnDown/F");
    T->Branch("Utr_smeared_JetEnDown", &Utr_smeared_JetEnDown, "Utr_smeared_JetEnDown/F");
    T->Branch("Ucol_smeared_JetEnDown", &Ucol_smeared_JetEnDown, "Ucol_smeared_JetEnDown/F");
    T->Branch("Ut_smeared_JetResDown", &Ut_smeared_JetResDown, "Ut_smeared_JetResDown/F");
    T->Branch("Utr_smeared_JetResDown", &Utr_smeared_JetResDown, "Utr_smeared_JetResDown/F");
    T->Branch("Ucol_smeared_JetResDown", &Ucol_smeared_JetResDown, "Ucol_smeared_JetResDown/F");
    T->Branch("Ut_smeared_JetResUp", &Ut_smeared_JetResUp, "Ut_smeared_JetResUp/F");
    T->Branch("Utr_smeared_JetResUp", &Utr_smeared_JetResUp, "Utr_smeared_JetResUp/F");
    T->Branch("Ucol_smeared_JetResUp", &Ucol_smeared_JetResUp, "Ucol_smeared_JetResUp/F");
    T->Branch("Ut_smeared_UnclusteredEnUp", &Ut_smeared_UnclusteredEnUp, "Ut_smeared_UnclusteredEnUp/F");
    T->Branch("Utr_smeared_UnclusteredEnUp", &Utr_smeared_UnclusteredEnUp, "Utr_smeared_UnclusteredEnUp/F");
    T->Branch("Ucol_smeared_UnclusteredEnUp", &Ucol_smeared_UnclusteredEnUp, "Ucol_smeared_UnclusteredEnUp/F");
    T->Branch("Ut_smeared_UnclusteredEnDown", &Ut_smeared_UnclusteredEnDown, "Ut_smeared_UnclusteredEnDown/F");
    T->Branch("Utr_smeared_UnclusteredEnDown", &Utr_smeared_UnclusteredEnDown, "Utr_smeared_UnclusteredEnDown/F");
    T->Branch("Ucol_smeared_UnclusteredEnDown", &Ucol_smeared_UnclusteredEnDown, "Ucol_smeared_UnclusteredEnDown/F");

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


    Float_t rho = analysisTree.rho;
    // searching for btagging discriminant ========================================================================================================================
    unsigned int nBTagDiscriminant1 = 0;
    unsigned int nBTagDiscriminant2 = 0;
    unsigned int nBTagDiscriminant3 = 0;
    SearchForBtagDiscriminant(analysisTree, BTagDiscriminator1, BTagDiscriminator2, BTagDiscriminator3, nBTagDiscriminant1, nBTagDiscriminant2, nBTagDiscriminant3, era);

    //JET Iteration and collection filling:
    TLorentzVector metLV;

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
            met_x = raw_pfmet_ex;
            met_y = raw_pfmet_ey;
        }
        else{
            met_x = analysisTree.puppimet_ex;
            met_y = analysisTree.puppimet_ey;
        }

        met = TMath::Sqrt(met_x*met_x + met_y*met_y);
        metLV.SetXYZT(met_x, met_y, 0., met);


	 // METs ===========================================================

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


    if (!isData ) {
            genweightsTree->GetEntry(iEntry);
            weight *= genweights;
            gen_weight *=genweights;
            lumi=true;
      }

    if (ApplyZptReweighting){

	 //Zpt Weights:
	 // Zpt reweighting for LO DY samples
        TString ZptweightFile = "";
	TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
	TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

	//ZPtWeightSys* zPtWeightSys = 0;
	//TopPtWeightSys* topPtWeightSys = 0;

	////////////////////////////////////////////////////////////
	// Z pt weight
	////////////////////////////////////////////////////////////

	TLorentzVector genV( 0., 0., 0., 0.);
	TLorentzVector genL( 0., 0., 0., 0.);

	//otree->zptweight = 1.;
	if (!isData && ((isDY  ))){
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
          //cout<<"zptmassweight"<<zptmassweight<<endl;
          //otree->zptweight = zptmassweight;
      }

    }

	TLorentzVector genBosonLV; genBosonLV.SetXYZT(0,0,0,0);
	TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);

	//Iterate over the gen particles and get the gen Boson:=============
	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

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
		if (bosonMass>50.0) {
		float bosonMassX = bosonMass;
		float bosonPtX = bosonPt;
		if (bosonMassX>1000.) bosonMassX = 1000.;
		if (bosonPtX<1.)      bosonPtX = 1.;
		if (bosonPtX>1000.)   bosonPtX = 1000.;
		zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
		histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
		}
     	}
	/*TODO ADD!!
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
      //TODO Check MET Flags for different years:
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

    if (!isData){
	  if (applyPUreweighting)	 {
        pu_weight = float(PUofficial->get_PUweight(double(analysisTree.primvertex_count)));
        puweight = float((PU_data->GetBinContent(analysisTree.primvertex_count)/PU_data->GetSumOfWeights())/(PU_mc->GetBinContent(analysisTree.primvertex_count)/PU_mc->GetSumOfWeights()));
	    if (pu_weight !=pu_weight) {
		    pu_weight=1;
		    }
	    if (pu_weight >10) {
		    pu_weight=1;
		}

	    weight *=pu_weight;
	    PUweight = float((PU_data->GetBinContent(analysisTree.primvertex_count)/PU_data->GetSumOfWeights())/(PU_mc->GetBinContent(analysisTree.primvertex_count)/PU_mc->GetSumOfWeights()));

	    if (pu_weight !=pu_weight) {
	                               pu_weight=1;
	    }
	    if (pu_weight >10) {
	                                pu_weight=1;
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
      /////now clear the Mu.El.Jets again to fill them again after cleaning
	
    vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<SingleElectronTriggerPtCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	//bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[ie];
	//if (applyElectronId && !electronMvaId) continue;
	bool electronCutBasedId = analysisTree.electron_cutId_tight_Summer16[ie];
	if (applyElectronId && !electronCutBasedId) continue;
	if (applyElectronId && !analysisTree.electron_pass_conversion[ie]) continue;
	if (applyElectronId && analysisTree.electron_nmissinginnerhits[ie]>1) continue;
	if (fabs(analysisTree.electron_charge[ie]) !=1) continue;
	 electrons.push_back((int)ie);

      }
      if (electrons.size()==0 ) continue;


      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float isoElMin  = 1e+10;
      float isoTauMin = 1.; 
      float isoTau = 1.; 
      float ptEl = 0;
      //      if (muons.size()>1||electrons.size()>1)
      //      std::cout << "muons = " << muons.size() << "  taus = " << taus.size() << std::endl;
      
	bool isLegMatch = false;


      for (unsigned int im=0; im<electrons.size(); ++im) {
	isLegMatch = false;
	//	bool isMuonTauMuonLegMatch = false;
	//	bool isMuonTauOverlapMuonMatch = false;
	unsigned int eIndex  = electrons.at(im);
	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[eIndex];
	float photonIsoElec = analysisTree.electron_photonIso[eIndex];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[eIndex];
	float puIsoElec = analysisTree.electron_puIso[eIndex];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[eIndex];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[eIndex];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/float(analysisTree.electron_pt[eIndex]);


	if (1)
	{	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {

	  	if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      	&& analysisTree.electron_pt[eIndex]>ptElectronCut&&
	      	analysisTree.trigobject_pt[iT]>SingleElectronTriggerPtCut) { 
	    	float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    	if (dRtrig<deltaRTrigMatch) 
	      	isLegMatch = true;
		}
	  

	  	  }
		}
	

	if (!isLegMatch) continue;


          if ((int)eIndex!=(int)el_index) {
            if (relIsoElec==isoElMin) {
              if (analysisTree.electron_pt[eIndex]>ptEl) {
                isoElMin  = relIsoElec;
                ptEl = analysisTree.electron_pt[eIndex];
                el_index =(int)eIndex;
              }
            }

            else if (relIsoElec<isoElMin) {
              isoElMin  = relIsoElec;
              ptEl = analysisTree.electron_pt[eIndex];
              el_index =(int)eIndex;
            }
          }


    }



      if ((int)el_index<0) continue;

	
      el_relIso[0]=isoElMin;
      

    bool   dilepton_veto=false;
    bool   extraelec_veto=false;
    bool   extramuon_veto=false;
    event_secondLeptonVeto = false;
    event_thirdLeptonVeto = false;



      // looking for extra muon
      bool foundExtraMuon = false;
      for (unsigned int ie = 0; ie<analysisTree.muon_count; ++ie) {
	if (isData && analysisTree.muon_isDuplicate[ie]) continue;
	if (isData && analysisTree.muon_isBad[ie]) continue;
	if (analysisTree.muon_pt[ie]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[ie])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[ie])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[ie])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[ie]) continue;

	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[ie];
        float photonIsoMu = analysisTree.muon_photonIso[ie];
        float chargedHadIsoMu = analysisTree.muon_chargedHadIso[ie];
        float puIsoMu = analysisTree.muon_puIso[ie];
        if (isIsoR03) {
          neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[ie];
          photonIsoMu = analysisTree.muon_r04_sumPhotonEt[ie];
          chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[ie];
          puIsoMu = analysisTree.muon_r04_sumPUPt[ie];
        }
        double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
        double neutralIsoMu = max(double(0),neutralIsoMuN);
        float absIsoMu = chargedHadIsoMu + neutralIsoMu;
        float relIsoMu = absIsoMu/analysisTree.muon_pt[ie];
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      // looking for extra electron's (dielectron veto)
      bool foundExtraElectron = false;
      vector<unsigned int> e_dielectrons; e_dielectrons.clear(); 
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {

      if ((int)im==(int)el_index) continue;

	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[im];
	float photonIsoElec = analysisTree.electron_photonIso[im];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[im];
	float puIsoElec = analysisTree.electron_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[im];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[im];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[im];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[im];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/analysisTree.electron_pt[im];

	if (analysisTree.electron_pt[im]>ptDilepElectronCut&&
	    fabs(analysisTree.electron_eta[im])<etaDilepElectronCut&&
	    fabs(analysisTree.electron_dxy[im])<dxyDilepElectronCut&&
	    fabs(analysisTree.electron_dz[im])<dzDilepElectronCut&&
	    analysisTree.electron_cutId_veto_Spring15[im]&&
	    relIsoElec<isoDilepElectronCut && 
	    fabs(analysisTree.electron_charge[im]) ==1)
	{
	    
	float dRelectrons = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			   analysisTree.electron_eta[im],analysisTree.electron_phi[im]);

	    if (dRelectrons>dRDilepVetoCut && (analysisTree.electron_charge[el_index]*analysisTree.electron_charge[im]<0.)) 
	      dilepton_veto = true;

	}

	 // e_dielectrons.push_back(im);

	if (analysisTree.electron_pt[im]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_general_Spring16_v1[im];
	if (applyVetoElectronId && !electronMvaId) continue;
	if (applyVetoElectronId && !analysisTree.electron_pass_conversion[im]) continue;
	if (applyVetoElectronId && analysisTree.electron_nmissinginnerhits[im]>1) continue;
	if (relIsoElec>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

   	event_secondLeptonVeto = dilepton_veto;
//	if (dilepton_veto)  continue;


	if (extraelec_veto) continue;
	if (extramuon_veto) continue;
//for now
//extraelec_veto = false;///////////////????????????

      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;



    /////////////////////////////////
      // Apply trigger SF
      // ////////////////////////////
      double ptEl1 = analysisTree.electron_pt[el_index];
      double etaEl1 = analysisTree.electron_eta[el_index];
      float trigweight=1.;
      float EffFromData = 1.;
      float EffFromMC = 1.;


      if (!isData) {
      EffFromData = (float)SF_electronTrigger->get_EfficiencyData(double(ptEl1),double(etaEl1));
      EffFromMC = (float)SF_electronTrigger->get_EfficiencyMC(double(ptEl1),double(etaEl1));

	if (EffFromMC>1e-6)
        trigweight = EffFromData/EffFromMC;

	weight *= trigweight;
	trig_weight = trigweight;


      //double IdIsoSF_el = SF_elIdIso->get_ScaleFactor(ptEl1, etaEl1);
	double IdIsoSF_el =0;
	IdIsoSF_el=hSF->GetBinContent(hSF->GetXaxis()->FindBin(etaEl1),hSF->GetYaxis()->FindBin(ptEl1));
//cout <<"SF!!!!!!!!!!!!!!!!!!!!!!!!!!"<< IdIsoSF_el<< endl;

	LSF_weight = IdIsoSF_el;
	weight *= LSF_weight;
      }



      bool isTauMatched = false;
      bool isGenLeptonMatched = false;


      //////////////////////////////////////////////
      muon_index = (int)mu_index;
      electron_index = (int)el_index;
      taus_index = (int)tau_index;

    el_count= (int)analysisTree.electron_count;
    for (unsigned int im=0;im<analysisTree.electron_count; ++im){
    int im = el_index;
    el_px[im]=analysisTree.electron_px[im];
    el_py[im]=analysisTree.electron_py[im];
    el_pz[im]=analysisTree.electron_pz[im];
    el_eta[im]=analysisTree.electron_eta[im];
    el_pt[im]=analysisTree.electron_pt[im];
    el_phi[im]=analysisTree.electron_phi[im];
    el_charge[im]=analysisTree.electron_charge[im];
    el_dxy[im]=analysisTree.electron_dxy[im];
    el_dz[im]=analysisTree.electron_dz[im];
    el_dxyerr[im]=analysisTree.electron_dxyerr[im];
    el_dzerr[im]=analysisTree.electron_dzerr[im];

    el_neutralHadIso[im] = analysisTree.electron_r03_sumNeutralHadronEt[im];
    el_photonIso[im] = analysisTree.electron_r03_sumPhotonEt[im];
    el_chargedHadIso[im] = analysisTree.electron_r03_sumChargedHadronPt[im];
    el_puIso[im] = analysisTree.electron_r03_sumPUPt[im];
 
    double neutralIso = el_neutralHadIso[im] + el_photonIso[im] - 0.5*el_puIso[im];
    neutralIso = max(double(0),neutralIso);
    el_neutralIso[im] = neutralIso;
    el_absIsoEl[im] = el_chargedHadIso[im] + neutralIso;
    el_relIsoEl[im]  = el_absIsoEl[im]/el_pt[im] ;
   
     }


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


  for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

     double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, analysisTree.pfjet_pt[jet]}, {JME::Binning::JetEta, analysisTree.pfjet_eta[jet]}, {JME::Binning::Rho, rho}});
     double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, analysisTree.pfjet_pt[jet]}, {JME::Binning::JetEta, analysisTree.pfjet_eta[jet]}}, Variation::NOMINAL);
     double jer_sf_up = resolution_sf.getScaleFactor({{JME::Binning::JetPt, analysisTree.pfjet_pt[jet]}, {JME::Binning::JetEta, analysisTree.pfjet_eta[jet]}}, Variation::UP);
     double jer_sf_down = resolution_sf.getScaleFactor({{JME::Binning::JetPt, analysisTree.pfjet_pt[jet]}, {JME::Binning::JetEta, analysisTree.pfjet_eta[jet]}}, Variation::DOWN);

     randm.SetSeed(static_cast<int>((analysisTree.pfjet_eta[jet] + 5) * 1000) * 1000 + static_cast<int>((analysisTree.pfjet_phi[jet] + 4) * 1000) + 10000);
     bool isEmbedded = false;
     if (!isData && !isEmbedded){
        shift = randm.Gaus(0, jet_resolution) * std::sqrt(std::max(jer_sf * jer_sf - 1, 0.0));
        shift_up = randm.Gaus(0, jet_resolution) * std::sqrt(std::max(jer_sf_up * jer_sf_up - 1, 0.0));
        shift_down = randm.Gaus(0, jet_resolution) * std::sqrt(std::max(jer_sf_down * jer_sf_down - 1, 0.0));
        if (shift < -1.0) shift = -1.0;
        if (shift_up < -1.0) shift_up = -1.0;
        if (shift_down < -1.0) shift_down = -1.0;
           }
       jetLV.SetPxPyPzE(analysisTree.pfjet_px[jet] * (1.0 + shift),
                        analysisTree.pfjet_py[jet] * (1.0 + shift),
                        analysisTree.pfjet_pz[jet] * (1.0 + shift),
                        analysisTree.pfjet_e[jet] * (1.0 + shift) );
       jetLVJERUp.SetPxPyPzE(analysisTree.pfjet_px[jet] * (1.0 + shift_up),
                        analysisTree.pfjet_py[jet] * (1.0 + shift_up),
                        analysisTree.pfjet_pz[jet] * (1.0 + shift_up),
                        analysisTree.pfjet_e[jet] * (1.0 + shift_up) );
       jetLVJERDown.SetPxPyPzE(analysisTree.pfjet_px[jet] * (1.0 + shift_down),
                        analysisTree.pfjet_py[jet] * (1.0 + shift_down),
                        analysisTree.pfjet_pz[jet] * (1.0 + shift_down),
                        analysisTree.pfjet_e[jet] * (1.0 + shift_down) );


        double absJetEta = fabs(analysisTree.pfjet_eta[jet]);
        double jetEta = analysisTree.pfjet_eta[jet];
        if (absJetEta>jetEtaCut) continue;


        float jetPt = analysisTree.pfjet_pt[jet] * (1.0 + shift);
        map<TString,TLorentzVector> jetLV_jecUnc;

        // Include variations for jec uncertainties
        for (auto uncer_split : jec_unc_map) {
           float sum_unc   = 0;
           for (auto single_jec_unc : uncer_split.second){
              JetCorrectionUncertainty *unc = single_jec_unc;
              unc->setJetPt(jetPt);
              unc->setJetEta(jetEta);
              double unc_ = unc->getUncertainty(true);
              sum_unc  += pow(unc_,2);
           }
           float unc_total  = TMath::Sqrt(sum_unc);
           jetLV_jecUnc[uncer_split.first + "Up"]   = jetLV * ( 1 + unc_total);
           jetLV_jecUnc[uncer_split.first + "Down"] = jetLV * ( 1 - unc_total);
           // Propagate jec uncertainties to met
           if( metLV_jecUnc.find(uncer_split.first+"Up") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Up"] = metLV;
           if( metLV_jecUnc.find(uncer_split.first+"Down") == metLV_jecUnc.end()) metLV_jecUnc[uncer_split.first+"Down"] = metLV;
           if (!isSampleForRecoilCorrection) {
              metLV_jecUnc[uncer_split.first + "Up"]   -= jetLV* unc_total;
              metLV_jecUnc[uncer_split.first + "Down"] += jetLV* unc_total;
           }
        }
        bool sync=true;
        float jetPt_tocheck = jetPt;
        if (sync) jetPt_tocheck = jetPt;
        else{
           float jetPtmax = jetPt;
           if(jetLV_jecUnc.at("jecUncEta0To5Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To5Down").Pt();
           if(jetLV_jecUnc.at("jecUncEta0To5Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To5Up").Pt();
           if(jetLV_jecUnc.at("jecUncEta0To3Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To3Down").Pt();
           if(jetLV_jecUnc.at("jecUncEta0To3Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta0To3Up").Pt();
           if(jetLV_jecUnc.at("jecUncEta3To5Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta3To5Down").Pt();
           if(jetLV_jecUnc.at("jecUncEta3To5Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEta3To5Up").Pt();
           if(jetLV_jecUnc.at("jecUncRelativeBalDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeBalDown").Pt();
           if(jetLV_jecUnc.at("jecUncRelativeBalUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeBalUp").Pt();
           if(jetLV_jecUnc.at("jecUncEC2Down").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2Down").Pt();
           if(jetLV_jecUnc.at("jecUncEC2Up").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2Up").Pt();
           if(jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt();
           if(jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt();
           if(jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt();
           if(jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt();
           if(jetLV_jecUnc.at("jecUncHFYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncHFYearDown").Pt();
           if(jetLV_jecUnc.at("jecUncHFYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncHFYearUp").Pt();
           if(jetLV_jecUnc.at("jecUncEC2YearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2YearDown").Pt();
           if(jetLV_jecUnc.at("jecUncEC2YearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncEC2YearUp").Pt();
           if(jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt();
           if(jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt();
           if(jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt();
           if(jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt() > jetPtmax) jetPtmax = jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt();
           if (jetLVJERUp.Pt() > jetPtmax) jetPtmax = jetLVJERUp.Pt();
           if (jetLVJERDown.Pt() > jetPtmax) jetPtmax = jetLVJERDown.Pt();
           jetPt_tocheck = jetPtmax;
        }
        if (sync) {if (jetPt<jetPtLowCut) continue; }
        else if (jetPt_tocheck<jetPtLowCut) continue;

        bool isPFJetId =false;
        if (era=="2016") isPFJetId= looseJetiD_2016(analysisTree,int(jet));
        else if (era=="2017") isPFJetId = tightJetiD_2017(analysisTree,int(jet));
        else if (era=="2018") isPFJetId = tightJetiD_2018(analysisTree,int(jet));
        if (!isPFJetId) continue;

        if (era=="2017" && jetPt < 50 && absJetEta > 2.65 && absJetEta < 3.139) continue;

        bool cleanedJet = true;
        float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           eta_1,phi_1);
        if (dR1<dRJetLeptonCut) cleanedJet = false;
        float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           eta_2,phi_2);
        if (dR2<dRJetLeptonCut) cleanedJet = false;
        if (!cleanedJet) continue;

        if (jetPt>jetPtLowCut) jetspt20.push_back(jet);
        if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

         bool tagged = false;
         bool tagged_mistagUp = false;
         bool tagged_mistagDown = false;
         bool tagged_btagUp = false;
         bool tagged_btagDown = false;
         tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3]>btagCut;
         tagged_mistagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
         tagged_mistagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
         tagged_btagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
         tagged_btagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] >btagCut;
         bool taggedRaw = tagged;

         if (!isData) {
            int flavor = abs(analysisTree.pfjet_flavour[jet]);

            double jet_scalefactor      = 1;
            double jet_scalefactor_up   = 1;
            double jet_scalefactor_down = 1;
            double JetPtForBTag = jetPt;
            double tageff = 1;

            if (flavor==5) {
               if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
               if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
               jet_scalefactor      = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
               jet_scalefactor_up   = reader_B.eval_auto_bounds("up"     ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
               jet_scalefactor_down = reader_B.eval_auto_bounds("down"   ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
               tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
            }
            else if (flavor==4) {
               if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
               if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
               jet_scalefactor      = reader_B.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
               jet_scalefactor_up   = reader_B.eval_auto_bounds("up"     , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
               jet_scalefactor_down = reader_B.eval_auto_bounds("down"   , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
               tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
            }
            else {
               if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
               if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
               jet_scalefactor      = reader_B.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
               jet_scalefactor_up   = reader_B.eval_auto_bounds("up"     , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
               jet_scalefactor_down = reader_B.eval_auto_bounds("down"   , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
               tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
            }
            if (tageff<1e-5)      tageff = 1e-5;
            if (tageff>0.99999)   tageff = 0.99999;
            r.SetSeed((int)((jetEta+5)*100000));
            double rannum = r.Rndm();

            if (tagged) { // demote
               if(jet_scalefactor<1){
                  double fraction = 1-jet_scalefactor;
                  if (rannum<fraction) tagged = false;
               }
               if(jet_scalefactor_up<1){
                  double fraction_up = 1-jet_scalefactor_up;
                  if (rannum<fraction_up) tagged_mistagUp = false;
               }
               if(jet_scalefactor_down<1){
                  double fraction_down = 1-jet_scalefactor_down;
                  if (rannum<fraction_down) tagged_mistagDown = false;
               }
               tagged_btagUp   = tagged;
               tagged_btagDown = tagged;
            }
               else if (!tagged) { // promote
                  if(jet_scalefactor>1){
                     double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                     if (rannum<fraction) tagged = true;
                  }
                  if(jet_scalefactor_up>1){
                     double fraction_up = (jet_scalefactor_up-1.0)/(1.0/tageff-1.0);
                     if (rannum<fraction_up) tagged_btagUp = true;
                  }
                  if(jet_scalefactor_down>1){
                     double fraction_down = (jet_scalefactor_down-1.0)/(1.0/tageff-1.0);
                     if (rannum<fraction_down) tagged_btagDown = true;
                  }
                  tagged_mistagUp   = tagged;
                  tagged_mistagDown = tagged;
               }
         }

         if (taggedRaw) bjetsRaw.push_back(jet);

         if (tagged) {
            bjets.push_back(jet);
               if (jetPt>ptLeadingBJet) {
                  ptLeadingBJet = jetPt;
                  indexLeadingBJet = jet;
                  shift_bjet = shift;
               }
         }

            if(tagged_mistagUp)   bjets_mistagUp.push_back(jet);
            if(tagged_mistagDown) bjets_mistagDown.push_back(jet);
            if(tagged_btagUp)     bjets_btagUp.push_back(jet);
            if(tagged_btagDown)   bjets_btagDown.push_back(jet);
      }
      if (jetLV_jecUnc.at("jecUncEta0To5Up").Pt()>jetPtHighCut) njets_jecUncEta0To5Up += 1;
      if (jetLV_jecUnc.at("jecUncEta0To5Down").Pt()>jetPtHighCut) njets_jecUncEta0To5Down += 1;
      if (jetLV_jecUnc.at("jecUncEta0To3Up").Pt()>jetPtHighCut) njets_jecUncEta0To3Up += 1;
      if (jetLV_jecUnc.at("jecUncEta0To3Down").Pt()>jetPtHighCut) njets_jecUncEta0To3Down += 1;
      if (jetLV_jecUnc.at("jecUncEta3To5Up").Pt()>jetPtHighCut) njets_jecUncEta3To5Up += 1;
      if (jetLV_jecUnc.at("jecUncEta3To5Down").Pt()>jetPtHighCut) njets_jecUncEta3To5Down += 1;
      if (jetLV_jecUnc.at("jecUncRelativeBalUp").Pt()>jetPtHighCut) njets_jecUncRelativeBalUp += 1;
      if (jetLV_jecUnc.at("jecUncRelativeBalDown").Pt()>jetPtHighCut) njets_jecUncRelativeBalDown += 1;
      if (jetLV_jecUnc.at("jecUncEC2Up").Pt()>jetPtHighCut) njets_jecUncEC2Up += 1;
      if (jetLV_jecUnc.at("jecUncEC2Down").Pt()>jetPtHighCut) njets_jecUncEC2Down += 1;
      if (jetLV_jecUnc.at("jecUncFlavorQCDUp").Pt()>jetPtHighCut) njets_jecUncFlavorQCDUp += 1;
      if (jetLV_jecUnc.at("jecUncFlavorQCDDown").Pt()>jetPtHighCut) njets_jecUncFlavorQCDDown += 1;
      if (jetLV_jecUnc.at("jecUncRelativeSampleYearUp").Pt()>jetPtHighCut) njets_jecUncRelativeSampleYearUp += 1;
      if (jetLV_jecUnc.at("jecUncRelativeSampleYearDown").Pt()>jetPtHighCut) njets_jecUncRelativeSampleYearDown += 1;
      if (jetLV_jecUnc.at("jecUncEC2YearUp").Pt()>jetPtHighCut) njets_jecUncEC2YearUp += 1;
      if (jetLV_jecUnc.at("jecUncEC2YearDown").Pt()>jetPtHighCut) njets_jecUncEC2YearDown += 1;
      if (jetLV_jecUnc.at("jecUncHFYearUp").Pt()>jetPtHighCut) njets_jecUncHFYearUp += 1;
      if (jetLV_jecUnc.at("jecUncHFYearDown").Pt()>jetPtHighCut) njets_jecUncHFYearDown += 1;
      if (jetLV_jecUnc.at("jecUncAbsoluteYearUp").Pt()>jetPtHighCut) njets_jecUncAbsoluteYearUp += 1;
      if (jetLV_jecUnc.at("jecUncAbsoluteYearDown").Pt()>jetPtHighCut) njets_jecUncAbsoluteYearDown += 1;
      if (jetLV_jecUnc.at("jecUncBBEC1YearUp").Pt()>jetPtHighCut) njets_jecUncBBEC1YearUp += 1;
      if (jetLV_jecUnc.at("jecUncBBEC1YearDown").Pt()>jetPtHighCut) njets_jecUncBBEC1YearDown += 1;
      if (jetLVJERUp.Pt() > jetPtHighCut) njets_jerUp += 1;
      if (jetLVJERDown.Pt() > jetPtHighCut) njets_jerDown += 1;

      if (jetPt>jetPtHighCut) jets.push_back(jet);

      //if (jetPt_tocheck<jetPtHighCut) continue; cut is done later to make sure that the uncertainties are treated correcly
      if (sync) if (jetPt <jetPtHighCut) continue;

            if (indexLeadingJet>=0) {
               if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
                  indexSubLeadingJet = jet;
                  ptSubLeadingJet = jetPt;
                  shift_jet2 =  shift;
                  shift_jet2_up = shift_up;
                  shift_jet2_down = shift_down;
               }
            }
            if (jetPt>ptLeadingJet) {
               indexSubLeadingJet = indexLeadingJet;
               ptSubLeadingJet = ptLeadingJet;
               indexLeadingJet = jet;
               ptLeadingJet = jetPt;
               shift_jet2 =  shift_jet1;
               shift_jet2_up = shift_jet1_up;
               shift_jet2_down = shift_jet1_down;
               shift_jet1 =  shift;
               shift_jet1_up = shift_up;
               shift_jet1_down = shift_down;
            }
         }
         njets = jets.size();
         int njetsMax = njets;

         njetspt20 = jetspt20.size();
         nbtag = bjets.size();
         nbtag_mistagUp   = bjets_mistagUp.size();
         nbtag_mistagDown = bjets_mistagDown.size();
         nbtag_btagUp   = bjets_btagUp.size();
         nbtag_btagDown = bjets_btagDown.size();
         nbtag_noSF = bjetsRaw.size();

         if (indexLeadingBJet>=0) {
            bpt = analysisTree.pfjet_pt[indexLeadingBJet]* (1.0 + shift_bjet);
            beta_1 = analysisTree.pfjet_eta[indexLeadingBJet];
            bphi = analysisTree.pfjet_phi[indexLeadingBJet];
         }

         if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
            cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;

         if (indexLeadingJet>=0) {
            jpt_1 = analysisTree.pfjet_pt[indexLeadingJet] * (1.0 + shift_jet1);
            jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
            jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
            jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
            jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet] * (1.0 + shift_jet1),
                            analysisTree.pfjet_py[indexLeadingJet] * (1.0 + shift_jet1),
                            analysisTree.pfjet_pz[indexLeadingJet] * (1.0 + shift_jet1),
                            analysisTree.pfjet_e[indexLeadingJet] * (1.0 + shift_jet1));

            jet1LV_jerUp.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet] * (1.0 + shift_jet1_up),
                            analysisTree.pfjet_py[indexLeadingJet] * (1.0 + shift_jet1_up),
                            analysisTree.pfjet_pz[indexLeadingJet] * (1.0 + shift_jet1_up),
                            analysisTree.pfjet_e[indexLeadingJet] * (1.0 + shift_jet1_up));
            jet1LV_jerDown.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet] * (1.0 + shift_jet1_down),
                            analysisTree.pfjet_py[indexLeadingJet] * (1.0 + shift_jet1_down),
                            analysisTree.pfjet_pz[indexLeadingJet] * (1.0 + shift_jet1_down),
                            analysisTree.pfjet_e[indexLeadingJet] * (1.0 + shift_jet1_down));

            for (auto uncer_split : jec_unc_map) {
               float sum_unc   = 0;
               for (auto single_jec_unc : uncer_split.second){
                  JetCorrectionUncertainty *unc = single_jec_unc;
                  unc->setJetPt(jet1.Pt());
                  unc->setJetEta(jet1.Eta());
                  double unc_ = unc->getUncertainty(true);
                  sum_unc  += pow(unc_,2);
               }
               float unc_total = TMath::Sqrt(sum_unc);
               jet1LV_jecUnc[uncer_split.first+"Up"]   = jet1*(1+unc_total);
               jet1LV_jecUnc[uncer_split.first+"Down"] = jet1*(1-unc_total);
            }
         }

         if (indexSubLeadingJet>=0) {
            jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet]  * (1.0 + shift_jet2);
            jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
            jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
            jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];

            jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet] * (1.0 + shift_jet2),
                            analysisTree.pfjet_py[indexSubLeadingJet] * (1.0 + shift_jet2),
                            analysisTree.pfjet_pz[indexSubLeadingJet] * (1.0 + shift_jet2),
                            analysisTree.pfjet_e[indexSubLeadingJet] * (1.0 + shift_jet2));
            jet2LV_jerUp.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet] * (1.0 + shift_jet2_up),
                            analysisTree.pfjet_py[indexSubLeadingJet] * (1.0 + shift_jet2_up),
                            analysisTree.pfjet_pz[indexSubLeadingJet] * (1.0 + shift_jet2_up),
                            analysisTree.pfjet_e[indexSubLeadingJet] * (1.0 + shift_jet2_up));
            jet2LV_jerDown.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet] * (1.0 + shift_jet2_down),
                            analysisTree.pfjet_py[indexSubLeadingJet] * (1.0 + shift_jet2_down),
                            analysisTree.pfjet_pz[indexSubLeadingJet] * (1.0 + shift_jet2_down),
                            analysisTree.pfjet_e[indexSubLeadingJet] * (1.0 + shift_jet2_down));

            for (auto uncer_split : jec_unc_map) {
               float sum_unc   = 0;
               for (auto single_jec_unc : uncer_split.second){
                  JetCorrectionUncertainty *unc = single_jec_unc;
                  unc->setJetPt(jet2.Pt());
                  unc->setJetEta(jet2.Eta());
                  double unc_ = unc->getUncertainty(true);
                  sum_unc  += pow(unc_,2);
               }
               float unc_total = TMath::Sqrt(sum_unc);
               jet2LV_jecUnc[uncer_split.first+"Up"]   = jet2*(1+unc_total);
               jet2LV_jecUnc[uncer_split.first+"Down"] = jet2*(1-unc_total);
            }
         }

         if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {

            mjj      = (jet1+jet2).M();
            dijetpt  = (jet1+jet2).Pt();
            dijetphi = (jet1+jet2).Phi();
            jdeta = fabs(analysisTree.pfjet_eta[indexLeadingJet]-analysisTree.pfjet_eta[indexSubLeadingJet]);
         }

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

	    double Dr=deltaR(analysisTree.electron_eta[el_index],
	                    analysisTree.electron_phi[el_index],
	                    analysisTree.pfjet_eta[jet],
	                    analysisTree.pfjet_phi[jet]);
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
        pfmet_ex_corr = analysisTree.pfmetcorr_ex;
        pfmet_ey_corr = analysisTree.pfmetcorr_ey;

        pfmet_phi_corr = calculate_dphi(pfmet_ex_corr, pfmet_ey_corr);
        pfmet_pt_corr = TMath::Sqrt(pfmet_ex_corr*pfmet_ex_corr + pfmet_ey_corr*pfmet_ey_corr);



       //Raw met:
       raw_pfmet_ex = pfmet_ex_corr - sum_pfjet_energycorr_l1fastjet;
       raw_pfmet_ey = pfmet_ey_corr - sum_pfjet_energycorr_l1fastjet;

      //Systematic for pfMET:

      //Generator level based MET:###############################################################################
      genmet = TMath::Sqrt(analysisTree.genmet_ex*analysisTree.genmet_ex+analysisTree.genmet_ey*analysisTree.genmet_ey);

      //PF MET:
      pfmet_ex = raw_pfmet_ex;
      pfmet_ey = raw_pfmet_ey;
      pfmet_ez = 0;
      pfmet_pt = TMath::Sqrt(pfmet_ex*pfmet_ex+pfmet_ey*pfmet_ey);
      pfmet_phi = calculate_dphi(pfmet_ex,pfmet_ey);

      pfmetcorr_ex_JetResUp = analysisTree.pfmetcorr_ex_JetResUp;
      pfmetcorr_ey_JetResUp = analysisTree.pfmetcorr_ey_JetResUp;

      pfmetcorr_ex_JetResDown = analysisTree.pfmetcorr_ex_JetResDown;
      pfmetcorr_ey_JetResDown = analysisTree.pfmetcorr_ey_JetResDown;

      pfmetcorr_ex_JetEnUp = analysisTree.pfmetcorr_ex_JetEnUp;
      pfmetcorr_ey_JetEnUp = analysisTree.pfmetcorr_ey_JetEnUp;

      pfmetcorr_ex_JetEnDown = analysisTree.pfmetcorr_ex_JetEnDown;
      pfmetcorr_ey_JetEnDown = analysisTree.pfmetcorr_ey_JetEnDown;

      pfmetcorr_ex_UnclusteredEnUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
      pfmetcorr_ey_UnclusteredEnUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;

      pfmetcorr_ex_UnclusteredEnDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
      pfmetcorr_ey_UnclusteredEnDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;

    pfmet_phi_JetEnUp = calculate_dphi(pfmet_ex_JetResUp,   pfmet_ey_JetResUp);
    pfmet_phi_JetEnDown = calculate_dphi(pfmet_ex_JetResDown,   pfmet_ey_JetResDown);
    pfmet_phi_UnclusteredEnUp = calculate_dphi(pfmet_ex_UnclusteredEnUp,pfmet_ey_UnclusteredEnUp);
    pfmet_phi_UnclusteredEnDown = calculate_dphi(pfmet_ex_UnclusteredEnDown,pfmet_ey_UnclusteredEnDown);
    pfmet_phi_JetResUp = calculate_dphi(pfmet_ex_JetResUp,pfmet_ey_JetResUp);
    pfmet_phi_JetResDown = calculate_dphi(pfmet_ex_JetResDown,pfmet_ey_JetResDown);

    //Smeared MET:============================================================================================
    pfmet_ex_corr_smeared = analysisTree.pfmetcorr_ex_smeared;
    pfmet_ey_corr_smeared = analysisTree.pfmetcorr_ey_smeared;
    pfmet_pt_corr_smeared = analysisTree.pfmetcorr_pt_smeared;
    pfmet_phi_corr_smeared = analysisTree.pfmetcorr_phi_smeared;


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

    //MET Significance:============================================================================================:
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

    // Apply Simple Recoil Corrections===========================================:
    int njetsforrecoil = njets;
    if(isW) njetsforrecoil = njets + 1;

    float met_uncorr = met;
    float metphi_uncorr = metphi;
    float pfmet_corr_x = analysisTree.pfmetcorr_ex;
    float pfmet_corr_y = analysisTree.pfmetcorr_ey;

    float pfmet_ex = analysisTree.pfmetcorr_ex;
    float pfmet_ey = analysisTree.pfmetcorr_ey;

    float puppimet_ex = analysisTree.puppimet_ex;
    float puppimet_ey = analysisTree.puppimet_ey;


        if ((isW||isDY) && !isData) {

            if (applySimpleRecoilCorrections) {

                  if (!usePuppiMET){
                    recoilMetCorrector.CorrectByMeanResolution(pfmet_corr_x,pfmet_corr_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
                    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                    metSys.ApplyMEtSys(pfmet_corr_x, pfmet_corr_y, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);

	             pfmet_ex_corr = pfmet_corr_x;
                pfmet_ey_corr = pfmet_corr_y;
          }
                   else{
                      recoilMetCorrectorPuppi.CorrectByMeanResolution(puppimet_ex,puppimet_ey,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,puppimet_ex, puppimet_ey);
                      metSysPuppi.ApplyMEtSys(puppimet_ex, puppimet_ey, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Up, met_x_recoilscaleUp, met_y_recoilscaleUp);
                      metSysPuppi.ApplyMEtSys(puppimet_ex, puppimet_ey, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Response, kit::MEtSys::SysShift::Down, met_x_recoilscaleDown, met_y_recoilscaleDown);
                      metSysPuppi.ApplyMEtSys(puppimet_ex, puppimet_ey, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Up, met_x_recoilresoUp, met_y_recoilresoUp);
                      metSysPuppi.ApplyMEtSys(puppimet_ex, puppimet_ey, bosonPx, bosonPy, lepPx, lepPy, njetsforrecoil, kit::MEtSys::SysType::Resolution, kit::MEtSys::SysShift::Down, met_x_recoilresoDown, met_y_recoilresoDown);
                   }
                }
                //todo check changes of the correction;
        }

       //TODO recalulate pt and phi
       pfmet_pt_corr = TMath::Sqrt(pfmet_ex_corr*pfmet_ex_corr + pfmet_ey_corr*pfmet_ey_corr);
       pfmet_phi_corr = calculate_dphi(pfmet_ex_corr, pfmet_ey_corr);

      //PuppiMET:============================================================================:
      puppi_ex = analysisTree.puppimet_ex;
      puppi_ey = analysisTree.puppimet_ey;
      puppi_ez = 0;//analysisTree.pfmet_ez;

      puppi_pt = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);
      puppi_phi = calculate_dphi(puppimet_ey,puppimet_ex);//TODO Recalculate via function;

      puppi_ex_JetEnUp = analysisTree.puppimet_ex_JetEnUp;
      puppi_ey_JetEnUp = analysisTree.puppimet_ey_JetEnUp;

      puppi_ex_JetEnDown = analysisTree.puppimet_ex_JetEnDown;
      puppi_ey_JetEnDown = analysisTree.puppimet_ey_JetEnDown;

      puppi_ex_UnclusteredEnUp = analysisTree.puppimet_ex_UnclusteredEnUp;
      puppi_ey_UnclusteredEnUp = analysisTree.puppimet_ey_UnclusteredEnUp;

      puppi_ex_UnclusteredEnDown = analysisTree.puppimet_ex_UnclusteredEnDown;
      puppi_ey_UnclusteredEnDown = analysisTree.puppimet_ey_UnclusteredEnDown;

      puppi_ex_JetResUp = analysisTree.puppimet_ex_JetResUp;
      puppi_ey_JetResUp = analysisTree.puppimet_ey_JetResUp;

      puppi_ex_JetResDown = analysisTree.puppimet_ex_JetResDown;
      puppi_ey_JetResDown = analysisTree.puppimet_ey_JetResDown;


      //==========================================================================================;

      puppi_pt_JetEnUp = TMath::Sqrt(puppi_ex_JetEnUp*puppi_ex_JetEnUp+
                                     puppi_ey_JetEnUp*puppi_ey_JetEnUp);

      puppi_pt_JetEnDown = TMath::Sqrt(puppi_ex_JetEnDown*puppi_ex_JetEnDown+
                                       puppi_ey_JetEnDown*puppi_ey_JetEnDown);

      puppi_pt_JetResUp = TMath::Sqrt(puppi_ex_JetResUp*puppi_ex_JetResUp+
                                      puppi_ey_JetResUp*puppi_ey_JetResUp);

      puppi_pt_JetResDown = TMath::Sqrt(puppi_ex_JetResDown*puppi_ex_JetResDown+
                                        puppi_ey_JetResDown*puppi_ey_JetResDown);

      puppi_pt_UnclusteredEnUp = TMath::Sqrt(puppi_ex_UnclusteredEnUp*puppi_ex_UnclusteredEnUp+
                                             puppi_ey_UnclusteredEnUp*puppi_ey_UnclusteredEnUp);

      puppi_pt_UnclusteredEnDown = TMath::Sqrt(puppi_ex_UnclusteredEnDown*puppi_ex_UnclusteredEnDown+
                                               puppi_ey_UnclusteredEnDown*puppi_ey_UnclusteredEnDown);

     //Unc for the puppi met phi;

    puppi_phi_JetEnUp = calculate_dphi(puppi_ex_JetEnUp,puppi_ey_JetEnUp);

    puppi_phi_JetEnDown = calculate_dphi(puppi_ex_JetEnDown,puppi_ey_JetEnDown);

    puppi_phi_JetResUp = calculate_dphi(puppi_ex_JetResUp,puppi_ey_JetResUp);

    puppi_phi_JetResDown = calculate_dphi(puppi_ex_JetResDown,puppi_ey_JetResDown);

    puppi_phi_UnclusteredEnUp = calculate_dphi(puppi_ex_UnclusteredEnUp, puppi_ey_UnclusteredEnUp);

    puppi_phi_UnclusteredEnDown = calculate_dphi(puppi_ex_UnclusteredEnDown,puppi_ey_UnclusteredEnDown);

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

    //=============================================================================//
    pfmetcorr_pt_JetEnUp = TMath::Sqrt(pfmetcorr_ex_JetEnUp*pfmetcorr_ex_JetEnUp +
                                       pfmetcorr_ey_JetEnUp*pfmetcorr_ey_JetEnUp);

    pfmetcorr_pt_JetEnDown = TMath::Sqrt(pfmetcorr_ex_JetEnDown*pfmetcorr_ex_JetEnDown +
                                         pfmetcorr_ey_JetEnDown*pfmetcorr_ey_JetEnDown);

    pfmetcorr_pt_UnclusteredEnUp = TMath::Sqrt(pfmetcorr_ex_UnclusteredEnUp*pfmetcorr_ex_UnclusteredEnUp +
                                               pfmetcorr_ey_UnclusteredEnUp*pfmetcorr_ey_UnclusteredEnUp);

    pfmetcorr_pt_UnclusteredEnDown = TMath::Sqrt(pfmetcorr_ex_UnclusteredEnDown*pfmetcorr_ex_UnclusteredEnDown +
                                                 pfmetcorr_ey_UnclusteredEnDown*pfmetcorr_ey_UnclusteredEnDown);

    pfmetcorr_pt_JetResUp = TMath::Sqrt(pfmetcorr_ex_JetResUp*pfmetcorr_ex_JetResUp +
                                        pfmetcorr_ey_JetResUp*pfmetcorr_ey_JetResUp);

    pfmetcorr_pt_JetResDown = TMath::Sqrt(pfmetcorr_ex_JetResDown*pfmetcorr_ex_JetResDown +
                                          pfmetcorr_ey_JetResDown*pfmetcorr_ey_JetResDown);
    //Uncertainties for the MET Phi Corr:
    pfmetcorr_phi_JetEnUp = calculate_dphi(pfmetcorr_ex_JetEnUp,pfmetcorr_ey_JetEnUp);

    pfmetcorr_phi_JetEnDown = calculate_dphi(pfmetcorr_ex_JetEnDown, pfmetcorr_ey_JetEnDown);

    pfmetcorr_phi_UnclusteredEnUp = calculate_dphi(pfmetcorr_ex_UnclusteredEnUp, pfmetcorr_ey_UnclusteredEnUp);

    pfmetcorr_phi_UnclusteredEnDown = calculate_dphi(pfmetcorr_ex_UnclusteredEnDown,pfmetcorr_ey_UnclusteredEnDown);

    pfmetcorr_phi_JetResUp = calculate_dphi(pfmetcorr_ex_JetResUp, pfmetcorr_ey_JetResUp);

    pfmetcorr_phi_JetResDown = calculate_dphi(pfmetcorr_ex_JetResDown, pfmetcorr_ey_JetResDown);

        //::+++:://
        // Uncertainties of the smeared MET:=====================================================================;
        /*
        float pfmetcorr_ex_JetEnUp_smeared = analysisTree.pfmetcorr_ex_JetEnUp_smeared;
        float pfmetcorr_ey_JetEnUp_smeared = analysisTree.pfmetcorr_ey_JetEnUp_smeared;

        pfmetcorr_ex_JetEnDown_smeared = analysisTree.pfmetcorr_ex_JetEnDown_smeared;
        pfmetcorr_ey_JetEnDown_smeared = analysisTree.pfmetcorr_ey_JetEnDown_smeared;

        pfmetcorr_ex_UnclusteredEnUp_smeared = analysisTree.pfmetcorr_ex_UnclusteredEnUp_smeared;
        pfmetcorr_ey_UnclusteredEnUp_smeared = analysisTree.pfmetcorr_ey_UnclusteredEnUp_smeared;

        pfmetcorr_ex_UnclusteredEnDown_smeared = analysisTree.pfmetcorr_ex_UnclusteredEnDown_smeared;
        pfmetcorr_ey_UnclusteredEnDown_smeared = analysisTree.pfmetcorr_ey_UnclusteredEnDown_smeared;

        pfmetcorr_ex_JetResUp_smeared = analysisTree.pfmetcorr_ex_JetResUp_smeared;
        pfmetcorr_ey_JetResUp_smeared = analysisTree.pfmetcorr_ey_JetResUp_smeared;

        pfmetcorr_ex_JetResDown_smeared = analysisTree.pfmetcorr_ex_JetResDown_smeared;
        pfmetcorr_ey_JetResDown_smeared = analysisTree.pfmetcorr_ey_JetResDown_smeared;


        pfmetcorr_pt_JetEnUp_smeared = TMath::Sqrt(pfmetcorr_ex_JetResUp_smeared*pfmetcorr_ex_JetResUp_smeared+
                                                    pfmetcorr_ey_JetResUp_smeared*pfmetcorr_ey_JetResUp_smeared);

        pfmetcorr_pt_JetEnDown_smeared = TMath::Sqrt(pfmetcorr_ex_JetResDown_smeared*pfmetcorr_ex_JetResDown_smeared
                                                    + pfmetcorr_ey_JetResDown_smeared*pfmetcorr_ey_JetResDown_smeared);

        pfmetcorr_pt_UnclusteredEnUp_smeared = TMath::Sqrt(pfmetcorr_ex_UnclusteredEnUp_smeared*pfmetcorr_ex_UnclusteredEnUp_smeared
                                                        + pfmetcorr_ey_UnclusteredEnUp_smeared*pfmetcorr_ey_UnclusteredEnUp_smeared);

        pfmetcorr_pt_UnclusteredEnDown_smeared = TMath::Sqrt(pfmetcorr_ex_UnclusteredEnDown_smeared*pfmetcorr_ex_UnclusteredEnDown_smeared
                                                            +pfmetcorr_ey_UnclusteredEnDown_smeared*pfmetcorr_ey_UnclusteredEnDown_smeared);

        pfmetcorr_pt_JetResUp_smeared = TMath::Sqrt(pfmetcorr_ex_JetResUp_smeared*pfmetcorr_ex_JetResUp_smeared
                                                  + pfmetcorr_ey_JetResUp_smeared*pfmetcorr_ey_JetResUp_smeared);

        pfmetcorr_pt_JetResDown_smeared = TMath::Sqrt(pfmetcorr_ex_JetResDown_smeared*pfmetcorr_ex_JetResDown_smeared +
                                                      pfmetcorr_ey_JetResDown_smeared*pfmetcorr_ey_JetResDown_smeared);

        */

     //MET XY Corrections:=============================================================================================:
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

     //XY Corrections: #########################################################################################
     std::pair<double,double> METXY = METXYCorr_Met_MetPhi(uncormet,
                                                           uncormet_phi,
                                                           runnb,
                                                           year,
                                                           isMC,
                                                           npv);
     pfmet_ex_XY = METXY.first;
     pfmet_ey_XY = METXY.second;
     pfmet_pt_XY = sqrt(pfmet_ex_XY*pfmet_ex_XY + pfmet_ey_XY*pfmet_ey_XY);
     pfmet_phi_XY = calculate_dphi(pfmet_ex_XY, pfmet_ey_XY);

     //MT for XY Corrections:
     double dPhi;
     dPhi = dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY,  pfmet_ey_XY );
	   pfMT_XY = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_XY*(1-TMath::Cos(dPhi)));

    //Get the XY Uncertainties:
    //1

    //1 JetEnUp
    METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_JetEnUp,
                                 pfmetcorr_phi_JetEnUp,
                                 runnb,
                                 year,
                                 isMC,
                                 npv);

     pfmet_ex_XY_JetEnUp = METXY.first;
     pfmet_ey_XY_JetEnUp = METXY.second;
     pfmet_pt_XY_JetEnUp = sqrt(pfmet_ex_XY_JetEnUp*pfmet_ex_XY_JetEnUp +
                                pfmet_ey_XY_JetEnUp*pfmet_ey_XY_JetEnUp);
     pfmet_phi_XY_JetEnUp = calculate_dphi(pfmet_ex_XY_JetEnUp, pfmet_ey_XY_JetEnUp);
     //2 JetEnDow
     METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_JetEnDown,
                                  pfmetcorr_phi_JetEnDown,
                                   runnb,
                                   year,
                                   isMC,
                                   npv);

     pfmet_ex_XY_JetEnDown = METXY.first;
     pfmet_ey_XY_JetEnDown = METXY.second;
     pfmet_pt_XY_JetEnDown = sqrt(pfmet_ex_XY_JetEnDown*pfmet_ex_XY_JetEnDown
                                + pfmet_ey_XY_JetEnDown*pfmet_ey_XY_JetEnDown);
     pfmet_phi_XY_JetEnDown = calculate_dphi(pfmet_ex_XY_JetEnDown, pfmet_ey_XY_JetEnDown);

     //3 UnclusteredEnUp
     METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_UnclusteredEnUp,
                                  pfmetcorr_phi_UnclusteredEnUp,
                                                           runnb,
                                                           year,
                                                           isMC,
                                                           npv);
     pfmet_ex_XY_UnclusteredEnUp = METXY.first;
     pfmet_ey_XY_UnclusteredEnUp = METXY.second;
     pfmet_pt_XY_UnclusteredEnUp = sqrt(pfmet_ex_XY_UnclusteredEnUp*pfmet_ex_XY_UnclusteredEnUp +
                                        pfmet_ey_XY_UnclusteredEnUp*pfmet_ey_XY_UnclusteredEnUp);
     pfmet_phi_XY_UnclusteredEnUp = calculate_dphi(pfmet_ex_XY_UnclusteredEnUp, pfmet_ey_XY_UnclusteredEnUp);

     //4 UnclusteredEnDown
     METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_UnclusteredEnDown,
                                   pfmetcorr_phi_UnclusteredEnDown,
                                   runnb,
                                   year,
                                   isMC,
                                   npv);

     pfmet_ex_XY_UnclusteredEnDown = METXY.first;
     pfmet_ey_XY_UnclusteredEnDown = METXY.second;
     pfmet_pt_XY_UnclusteredEnDown = sqrt(pfmet_ex_XY_UnclusteredEnDown*pfmet_ex_XY_UnclusteredEnDown
                                        + pfmet_ey_XY_UnclusteredEnDown*pfmet_ey_XY_UnclusteredEnDown);
     pfmet_phi_XY_UnclusteredEnDown = calculate_dphi(pfmet_ex_XY_UnclusteredEnDown, pfmet_ey_XY_UnclusteredEnDown);


     //5 JetResUp
     METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_JetResUp,
                                   pfmetcorr_phi_JetResUp,
                                   runnb,
                                   year,
                                   isMC,
                                   npv);

     pfmet_ex_XY_JetResUp = METXY.first;
     pfmet_ey_XY_JetResUp = METXY.second;
     pfmet_pt_XY_JetResUp = sqrt(pfmet_ex_XY_JetResUp*pfmet_ex_XY_JetResUp
                               + pfmet_ey_XY_JetResUp*pfmet_ey_XY_JetResUp);
     pfmet_phi_XY_JetResUp = calculate_dphi(pfmet_ex_XY_JetResUp, pfmet_ey_XY_JetResUp);


     //6 JetResDown
     METXY = METXYCorr_Met_MetPhi(pfmetcorr_pt_JetResDown,
                                  pfmetcorr_phi_JetResDown,
                                   runnb,
                                   year,
                                   isMC,
                                   npv);
     pfmet_ex_XY_JetResDown = METXY.first;
     pfmet_ey_XY_JetResDown = METXY.second;

     pfmet_pt_XY_JetResDown = sqrt(pfmet_ex_XY_JetResDown*pfmet_ex_XY_JetResDown
                                 + pfmet_ey_XY_JetResDown*pfmet_ey_XY_JetResDown);
     pfmet_phi_XY_JetResDown = calculate_dphi(pfmet_ex_XY_JetResDown, pfmet_ey_XY_JetResDown);


     //The same XY correction for JECed MET:===========================================================================:

     uncormet = analysisTree.pfmetcorr_pt;
     uncormet_phi = analysisTree.pfmetcorr_phi;

     METXY = METXYCorr_Met_MetPhi( uncormet,  uncormet_phi,  runnb,  year,  isMC,  npv);

     pfmet_ex_corr_XY = METXY.first;
     pfmet_ey_corr_XY = METXY.second;
     pfmet_pt_corr_XY = sqrt(pfmet_ex_corr_XY*pfmet_ex_corr_XY + pfmet_ey_corr_XY*pfmet_ey_corr_XY);
     pfmet_phi_corr_XY = calculate_dphi(pfmet_ex_corr_XY, pfmet_ey_corr_XY);


    //MT Definition and calculation:===================================================================================:

     dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex,  pfmet_ey );
	 pfMT = TMath::Sqrt(2*el_pt[0]*pfmet_pt*(1-TMath::Cos(dPhi)));

     dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex,  puppi_ey );
	 puppiMT = TMath::Sqrt(2*el_pt[0]*puppi_pt*(1-TMath::Cos(dPhi)));

  	 dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_corr_smeared,  pfmet_ey_corr_smeared );
	 MT_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_corr_smeared*(1-TMath::Cos(dPhi)));

     dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_corr,  pfmet_ey_corr);
     pfMT_corr = TMath::Sqrt(2*el_pt[0]*pfmet_pt_corr*(1-TMath::Cos(dPhi)));


    //========================================================
    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY_JetEnUp,  pfmet_ey_XY_JetEnUp );
    MT_XY_JetEnUp = TMath::Sqrt(2*el_pt[0]*pfmet_pt_XY_JetEnUp*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY_JetEnDown,  pfmet_ey_XY_JetEnDown );
    MT_XY_JetEnDown = TMath::Sqrt(2*el_pt[0]*pfmet_pt_XY_JetEnDown*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY_UnclusteredEnUp,  pfmet_ey_XY_UnclusteredEnUp );
    MT_XY_UnclusteredEnUp = TMath::Sqrt(2*el_pt[0]*pfmet_pt_XY_UnclusteredEnUp*(1-TMath::Cos(dPhi)));

        dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY_UnclusteredEnDown,  pfmet_ey_XY_UnclusteredEnDown );
        MT_XY_UnclusteredEnDown = TMath::Sqrt(2*el_pt[0]*pfmet_pt_XY_UnclusteredEnDown*(1-TMath::Cos(dPhi)));

        dPhi=dPhiFrom2P( mu_px[0], mu_py[0], pfmet_ex_XY_JetResUp,  pfmet_ey_XY_JetResUp );
        MT_XY_JetResUp = TMath::Sqrt(2*mu_pt[0]*pfmet_pt_XY_JetResUp*(1-TMath::Cos(dPhi)));

        dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_XY_JetResDown,  pfmet_ey_XY_JetResDown );
        MT_XY_JetResDown = TMath::Sqrt(2*el_pt[0]*pfmet_pt_XY_JetResDown*(1-TMath::Cos(dPhi)));


       float mu_unc = 0.3;

  	dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_JetEnUp,  pfmetcorr_ey_JetEnUp );
	MT_JetEnUp = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_JetEnUp*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_JetEnDown,  pfmetcorr_ey_JetEnDown );
  	MT_JetEnDown = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_JetEnDown*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_UnclusteredEnUp,  pfmetcorr_ey_UnclusteredEnUp );
  	MT_UnclusteredEnUp = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_UnclusteredEnUp*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_UnclusteredEnDown,  pfmetcorr_ey_UnclusteredEnDown );
  	MT_UnclusteredEnDown = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_UnclusteredEnDown*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_JetResUp,  pfmetcorr_ey_JetResUp );
  	MT_JetResUp = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_JetResUp*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmetcorr_ex_JetResDown,  pfmetcorr_ey_JetResDown );
  	MT_JetResDown = TMath::Sqrt(2*el_pt[0]*pfmetcorr_pt_JetResDown*(1-TMath::Cos(dPhi)));
    //=================================================================================================================:

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex,  puppi_ey );
    MTpuppi = TMath::Sqrt(2*el_pt[0]*puppi_pt*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_JetEnUp,  puppi_ey_JetEnUp );
    MTpuppi_JetEnUp = TMath::Sqrt(2*el_pt[0]*puppi_pt_JetEnUp*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_JetEnDown,  puppi_ey_JetEnDown );
    MTpuppi_JetEnDown = TMath::Sqrt(2*el_pt[0]*puppi_pt_JetEnDown*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_UnclusteredEnUp,  puppi_ey_UnclusteredEnUp );
    MTpuppi_UnclusteredEnUp = TMath::Sqrt(2*el_pt[0]*puppi_pt_UnclusteredEnUp*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_UnclusteredEnDown,  puppi_ey_UnclusteredEnDown );
    MTpuppi_UnclusteredEnDown = TMath::Sqrt(2*el_pt[0]*puppi_pt_UnclusteredEnDown*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_JetResUp,  puppi_ey_JetResUp );
   	MTpuppi_JetResUp = TMath::Sqrt(2*el_pt[0]*puppi_pt_JetResUp*(1-TMath::Cos(dPhi)));

	dPhi=dPhiFrom2P( el_px[0], el_py[0], puppi_ex_JetResDown,  puppi_ey_JetResDown );
    MTpuppi_JetResDown = TMath::Sqrt(2*el_pt[0]*puppi_pt_JetResDown*(1-TMath::Cos(dPhi)));


    //=================================================================================================================:

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_smeared,  pfmet_ey_smeared );
    MT_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_JetEnUp_smeared,  pfmet_ey_JetEnUp_smeared );
    MT_JetEnUp_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_JetEnUp_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_JetEnDown_smeared,  pfmet_ey_JetEnDown_smeared );
    MT_JetEnDown_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_JetEnDown_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_UnclusteredEnUp_smeared,  pfmet_ey_UnclusteredEnUp_smeared );
    MT_UnclusteredEnUp_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_UnclusteredEnUp_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_UnclusteredEnDown_smeared,  pfmet_ey_UnclusteredEnDown_smeared );
    MT_UnclusteredEnDown_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_UnclusteredEnDown_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_JetResUp_smeared,  pfmet_ey_JetResUp_smeared );
    MT_JetResUp_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_JetResUp_smeared*(1-TMath::Cos(dPhi)));

    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_JetResDown_smeared,  pfmet_ey_JetResDown_smeared );
    MT_JetResDown_smeared = TMath::Sqrt(2*el_pt[0]*pfmet_pt_JetResDown_smeared*(1-TMath::Cos(dPhi)));

    //MT for XY Corrections:===========================================================================================:
    dPhi=dPhiFrom2P( el_px[0], el_py[0], pfmet_ex_corr_XY,  pfmet_ey_corr_XY );
    pfMT_corr_XY = TMath::Sqrt(2*el_pt[0]*pfmet_pt_corr_XY*(1-TMath::Cos(dPhi)));

    //MET performance:=================================================================================================:
    Float_t Utx,Uty;
    Utx = - el_px[0] - pfmet_corr_x;
    Uty = - el_py[0] - pfmet_corr_y;

	Ut = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol = TMath::Cos(dPhi)* Ut;



	Utx = - el_px[0] - met_ex_JetEnUp;
	Uty = - el_py[0] - met_ey_JetEnUp;
	Ut_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_JetEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_JetEnUp = TMath::Cos(dPhi)* Ut_JetEnUp ;

	Utx = - el_px[0] - met_ex_JetEnDown;
	Uty = - el_py[0] - met_ey_JetEnDown;
	Ut_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_JetEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_JetEnDown = TMath::Cos(dPhi)* Ut_JetEnDown ;

	Utx = - el_px[0] - pfmet_ex_JetResDown;
	Uty = - el_py[0] - pfmet_ey_JetResDown;
	Ut_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_JetResDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_JetResDown = TMath::Cos(dPhi)* Ut_JetResDown ;

	Utx = - el_px[0] - pfmet_ex_JetResUp;
	Uty = - el_py[0] - pfmet_ey_JetResUp;
	Ut_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_JetResUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_JetResUp = TMath::Cos(dPhi)* Ut_JetResUp ;

	Utx = - el_px[0] - pfmet_ex_UnclusteredEnUp;
	Uty = - el_py[0] - pfmet_ey_UnclusteredEnUp;
	Ut_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_UnclusteredEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_UnclusteredEnUp ;

	Utx = - el_px[0] - pfmet_ex_UnclusteredEnDown;
	Uty = - el_py[0] - pfmet_ey_UnclusteredEnDown;
	Ut_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_UnclusteredEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_UnclusteredEnDown ;



////////// Puppi
	Utx = - el_px[0] - puppi_ex;
	Uty = - el_py[0] - puppi_ey;
	Ut_puppi = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi = TMath::Cos(dPhi)* Ut;

	Utx = - el_px[0] - puppi_ex_JetEnUp;
	Uty = - el_py[0] - puppi_ey_JetEnUp;
	Ut_puppi_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_JetEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_JetEnUp = TMath::Cos(dPhi)* Ut_puppi_JetEnUp ;

	Utx = - el_px[0] - puppi_ex_JetEnDown;
	Uty = - el_py[0] - puppi_ey_JetEnDown;
	Ut_puppi_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_JetEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_JetEnDown = TMath::Cos(dPhi)* Ut_puppi_JetEnDown ;

	Utx = - el_px[0] - puppi_ex_JetResDown;
	Uty = - el_py[0] - puppi_ey_JetResDown;
	Ut_puppi_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_JetResDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_JetResDown = TMath::Cos(dPhi)* Ut_puppi_JetResDown ;

	Utx = - el_px[0] - puppi_ex_JetResUp;
	Uty = - el_py[0] - puppi_ey_JetResUp;
	Ut_puppi_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_JetResUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_JetResUp = TMath::Cos(dPhi)* Ut_puppi_JetResUp ;

	Utx = - el_px[0] - puppi_ex_UnclusteredEnUp;
	Uty = - el_py[0] - puppi_ey_UnclusteredEnUp;
	Ut_puppi_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_UnclusteredEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_puppi_UnclusteredEnUp ;

	Utx = - el_px[0] - puppi_ex_UnclusteredEnDown;
	Uty = - el_py[0] - puppi_ey_UnclusteredEnDown;
	Ut_puppi_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_puppi_UnclusteredEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_puppi_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_puppi_UnclusteredEnDown ;

    //Smeared MET:
	Utx = - el_px[0] - pfmet_ex_smeared;
	Uty = - el_py[0] - pfmet_ey_smeared;
	Ut_smeared = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared = TMath::Cos(dPhi)* Ut;

	Utx = - el_px[0] - pfmet_ex_JetEnUp_smeared;
	Uty = - el_py[0] - pfmet_ey_JetEnUp_smeared;
	Ut_smeared_JetEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_JetEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_JetEnUp = TMath::Cos(dPhi)* Ut_smeared_JetEnUp ;

	Utx = - el_px[0] - pfmet_ex_JetEnDown_smeared;
	Uty = - el_py[0] - pfmet_ey_JetEnDown_smeared;
	Ut_smeared_JetEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_JetEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_JetEnDown = TMath::Cos(dPhi)* Ut_smeared_JetEnDown ;

	Utx = - el_px[0] - pfmet_ex_JetResDown_smeared;
	Uty = - el_py[0] - pfmet_ey_JetResDown_smeared;
	Ut_smeared_JetResDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_JetResDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_JetResDown = TMath::Cos(dPhi)* Ut_smeared_JetResDown ;

	Utx = - el_px[0] - pfmet_ex_JetResUp_smeared;
	Uty = - el_py[0] - pfmet_ey_JetResUp_smeared;
	Ut_smeared_JetResUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_JetResUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_JetResUp = TMath::Cos(dPhi)* Ut_smeared_JetResUp ;

	Utx = - el_px[0] - pfmet_ex_UnclusteredEnUp_smeared;
	Uty = - el_py[0] - pfmet_ey_UnclusteredEnUp_smeared;
	Ut_smeared_UnclusteredEnUp = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_UnclusteredEnUp = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_UnclusteredEnUp = TMath::Cos(dPhi)* Ut_smeared_UnclusteredEnUp ;

	Utx = - el_px[0] - pfmet_ex_UnclusteredEnDown_smeared;
	Uty = - el_py[0] - pfmet_ey_UnclusteredEnDown_smeared;
	Ut_smeared_UnclusteredEnDown = TMath::Sqrt((Utx*Utx)+(Uty*Uty));
	dPhi=dPhiFrom2P( Utx, Uty, el_px[0],  el_py[0] );
	Utr_smeared_UnclusteredEnDown = -Utx*((el_py[0]/el_px[0])/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0]))) + Uty*(1/TMath::Sqrt(1+(el_py[0]/el_px[0])*(el_py[0]/el_px[0])));
	Ucol_smeared_UnclusteredEnDown = TMath::Cos(dPhi)* Ut_smeared_UnclusteredEnDown ;

    //MET performance ends but not at this place!######

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
  file->Write();
  file->Close();

  delete file;

}
