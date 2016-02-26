
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Spring15Tree.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>



#define pi 3.14159265358979312
#define d2r 1.74532925199432955e-02
#define r2d 57.2957795130823229

#define electronMass 0
#define muonMass 0.10565837
#define tauMass 1.77682
#define pionMass 0.1396


typedef std::vector<std::pair<int,int> > lumi_json;

struct compare_lumi {
  bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
    if (left.first < right.first)
      return 1;
    else if (left.first > right.first)
      return 0;
    else
      return left.second < right.second;
  }
};

int read_json(std::string filename, lumi_json& json){

  std::pair <int,int> lumi;

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);

  BOOST_FOREACH(boost::property_tree::ptree::value_type &json_run, pt.get_child("")){
    int irun = atoi(json_run.first.data());
    BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_ranges, json_run.second.get_child("")){
      int ilumi[2] = {};

      int count = 0;
      BOOST_FOREACH(boost::property_tree::ptree::value_type &lumi_boundaries, lumi_ranges.second.get_child("")){
	ilumi[count] = atoi(lumi_boundaries.second.data().data());
	count++;
      }
      
      for (;ilumi[0] <= ilumi[1]; ilumi[0]++){
	lumi = std::make_pair(irun, ilumi[0]);
	json.push_back(lumi);
      }
    }
  }

  sort( json.begin(), json.end(),  compare_lumi());
  json.erase( unique( json.begin(), json.end() ), json.end() );
  
  return 0;
}

bool isGoodLumi(const std::pair<int, int>& lumi, const lumi_json& json){
  static compare_lumi compare;
  static std::pair<int,int> oldlumi = lumi;
  static bool old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  
  if(lumi.first != oldlumi.first || lumi.second != oldlumi.second){
    oldlumi = lumi;
    old = std::binary_search(json.begin(), json.end(), lumi, compare_lumi());
  }

  return old;
}

bool isGoodLumi(int run, int lumi, const lumi_json& json){
  std::pair<int, int> run_lumi = std::make_pair(run, lumi);
  return isGoodLumi(run_lumi, json);
}


int main(int argc, char * argv[]) {

  // first argument - config file for analysis
  // second argument - config file for process
  // third argument - index of first file to run on (optional, ignored if only one file is used)
  // fourth argument - index of last file to run on (optional, ignored if only one file is used)
  
  using namespace std;

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // **** configuration analysis
  Config cfg(argv[1]);

  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];

  float xs = -1.;
  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase)+"/src/"+TString(json_name)).Data(), json);
  }

  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");

  //lep eff
  const string idIsoEffFile = cfg.get<string>("idIsoEffFile");
  const string trigEffFile = cfg.get<string>("trigEffFile");

  //svfit
  const string svFitPtResFile = cfg.get<string>("svFitPtResFile");

  // HLT filters
  string isoLeg;
  if (isData){
    isoLeg = cfg.get<string>("isoLegData");
  }
  else {isoLeg = cfg.get<string>("isoLegMC");}

  const float ptTrigObjCut  = cfg.get<float>("ptTrigObjCut");
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");

  // muon cuts
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  
  // tau cuts
  const float ptTauLowCut    = cfg.get<float>("ptTauLowCut");
  const float ptTauHighCut   = cfg.get<float>("ptTauHighCut");  
  const float etaTauCut      = cfg.get<float>("etaTauCut");
  const float dzTauCut      = cfg.get<float>("dzTauCut");
  const bool applyTauId      = cfg.get<bool>("ApplyTauId");

  // pair selection
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");

  // dimuon veto
  const float ptDiMuonVeto     = cfg.get<float>("ptDiMuonVeto");  
  const float etaDiMuonVeto    = cfg.get<float>("etaDiMuonVeto");
  const float dxyDiMuonVeto    = cfg.get<float>("dxyDiMuonVeto");  
  const float dzDiMuonVeto     = cfg.get<float>("dzDiMuonVeto"); 
  const bool applyDiMuonVetoId = cfg.get<bool>("applyDiMuonVetoId");
  const bool applyDiMuonOS     = cfg.get<bool>("applyDiMuonOS");
  const float isoDiMuonVeto    = cfg.get<float>("isoDiMuonVeto");
  const float drDiMuonVeto     = cfg.get<float>("drDiMuonVeto");  
  
  // extra electron veto
  const float ptVetoElectronCut  = cfg.get<float>("ptVetoElectronCut");  
  const float etaVetoElectronCut = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut = cfg.get<float>("dxyVetoElectronCut");  
  const float dzVetoElectronCut  = cfg.get<float>("dzVetoElectronCut"); 
  const bool applyVetoElectronId = cfg.get<bool>("applyVetoElectronId");
  const float isoVetoElectronCut = cfg.get<float>("isoVetoElectronCut");  
  
  // extra muon veto
  const float ptVetoMuonCut  = cfg.get<float>("ptVetoMuonCut");  
  const float etaVetoMuonCut = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut = cfg.get<float>("dxyVetoMuonCut");  
  const float dzVetoMuonCut  = cfg.get<float>("dzVetoMuonCut"); 
  const bool applyVetoMuonId = cfg.get<bool>("applyVetoMuonId");
  const float isoVetoMuonCut = cfg.get<float>("isoVetoMuonCut");  

    
  // topological cuts
  //const float dZetaCut       = cfg.get<float>("dZetaCut");
  //const bool oppositeSign    = cfg.get<bool>("oppositeSign");

  const float deltaRTrigMatch = cfg.get<float>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  
  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
  const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
  const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
  const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
  const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");

  const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
  const float btagCut = cfg.get<float>("btagCut");

  // check overlap
  const bool checkOverlap = cfg.get<bool>("CheckOverlap");
  const bool debug = cfg.get<bool>("debug");
  
  // **** end of configuration analysis
    
  int ifile = 0;
  int jfile = -1;

  if (argc > 3)
    ifile = atoi(argv[3]);
  if (argc > 4)
    jfile = atoi(argv[4]);
  
  // create input files list
  std::vector<std::string> fileList;  
  if (infiles.find(".root") != std::string::npos){
    ifile = 0;
    jfile = 1;

    fileList.push_back(infiles);
  }
  else{
    ifstream input;
    std::string infile;
    
    input.open(infiles);

    while(true){
      input>>infile;
      if(!input.eof()){
	if (infile.length() > 0)
	  fileList.push_back(infile);
      }
      else
	break;
    }

    if(jfile < 0)
      jfile = fileList.size();   
  }

  for (int iF=ifile; iF<jfile; ++iF) {
    std::cout<<fileList[iF]<<std::endl;
  }
  
  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");

  // PU reweighting - initialization
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(pileUpInDataFile),"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(pileUpInMCFile), "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Lepton Scale Factors
  // Electron Id+Iso scale factor
  ScaleFactor * SF_muIdIso = new ScaleFactor();
  SF_muIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(idIsoEffFile));

  // Electron SingleElectron trigger scale factor
  ScaleFactor * SF_muTrigger = new ScaleFactor();
  SF_muTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(trigEffFile));

  // output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_mt_Sync.root";
  std::cout <<rootFileName <<std::endl;  

  TFile * file = new TFile( rootFileName ,"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5,0.5);
  
  TTree * tree = new TTree("TauCheck","TauCheck");

  Spring15Tree *otree = new Spring15Tree(tree);

  int nTotalFiles = 0;

  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;
  
  vector<int> runList; runList.clear();
  vector<int> eventList; eventList.clear();

  int nonOverlap = 0;

  std::ifstream fileEvents("overlap.txt");
  int Run, Event, Lumi;
  if (checkOverlap) {
    std::cout << "Non-overlapping events ->" << std::endl;
    while (fileEvents >> Run >> Event >> Lumi) {
      runList.push_back(Run);
      eventList.push_back(Event);
      std::cout << Run << ":" << Event << std::endl;
    }
    std::cout << std::endl;
  }
  std::ofstream fileOutput("overlap.out");

  //svFit
  TH1::AddDirectory(false);  
  //TFile* inputFile_visPtResolution = new TFile(svFitPtResFile.data());

  for (int iF=ifile; iF<jfile; ++iF) {

    std::cout << "file " << iF+1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
    TFile * file_ = TFile::Open(fileList[iF].data());
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {       
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (isData)
	nWeightedEventsH->Fill(0., 1.);
      else
	nWeightedEventsH->Fill(0., analysisTree.genweight);

      unsigned int nIsoLeg = 0;
      bool checkIsoLeg = false;

      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
	TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==isoLeg) {
	  nIsoLeg = i;
	  checkIsoLeg = true;
	}
      }
      if (!checkIsoLeg) {
	std::cout << "HLT filter " << isoLeg << " not found" << std::endl;
	exit(-1);
      }
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      otree->run = int(analysisTree.event_run);
      otree->lumi = int(analysisTree.event_luminosityblock);
      otree->evt = int(analysisTree.event_nr);
      
      bool overlapEvent = true;
      for (unsigned int iEvent=0; iEvent<runList.size(); ++iEvent) {
	if (runList.at(iEvent)==otree->run && eventList.at(iEvent)==otree->evt) {
	  overlapEvent = false;	  
	}
      }

      if (overlapEvent&&checkOverlap) continue;
      nonOverlap++;

      if (debug) {
	fileOutput << std::endl;
	fileOutput << "Run = " << otree->run << "   Event = " << otree->evt << std::endl;
      }

      if (isData && !isGoodLumi(otree->run, otree->lumi, json))
	continue;
      
      // weights
      otree->mcweight = 1.;
      otree->pu_weight = 0;


      if(!isData)
	otree->mcweight = analysisTree.genweight;
	otree->pu_weight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));

      otree->xs = xs;
      otree->trigweight_1 = 0;
      otree->trigweight_2 = 1;
      otree->idisoweight_1 = 0;
      otree->idisoweight_2 = 1;

      otree->effweight = 0;
      otree->fakeweight = 0;
      otree->embeddedWeight = 0;
      otree->signalWeight = 0;
      otree->weight = 1;
      otree->lheHt = analysisTree.genparticles_lheHt;
      
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;
	

      // vertex cuts
      //if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      //if (analysisTree.primvertex_ndof<=ndofVertexCut) continue;
      //float dVertex = sqrt(analysisTree.primvertex_x*analysisTree.primvertex_x+
      //		   analysisTree.primvertex_y*analysisTree.primvertex_y);
      //if (dVertex>dVertexCut)
      //continue;
      
      if (debug)
	fileOutput << "Vertex cuts are passed " << std::endl;

      // muon selection
      vector<int> muons; muons.clear();
      if (debug)
	fileOutput << "# muons = " << analysisTree.muon_count << std::endl;
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {

	bool muonMediumId = analysisTree.muon_isMedium[im];
	
	if (debug)
	  fileOutput << "  " << im
		     << " pt = " << analysisTree.muon_pt[im] 
		     << " eta = " << analysisTree.muon_eta[im]
		     << " dxy = " << analysisTree.muon_dxy[im]
		     << " dz  = " << analysisTree.muon_dz[im];
		     //<< " passConv = " << analysisTree.electron_pass_conversion[ie]
		     //<< " nmisshits = " << int(analysisTree.electron_nmissinginnerhits[ie])
		     //<< " mvaTight = " << electronMvaId << std::endl;
	if (analysisTree.muon_pt[im]<=ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>=etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>=dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>=dzMuonCut) continue;
	if (!muonMediumId&&applyMuonId) continue;
	//if (!analysisTree.electron_pass_conversion[ie]&&applyElectronId) continue;
	//if (analysisTree.electron_nmissinginnerhits[ie]>1&&applyElectronId) continue;
	muons.push_back(im);
      }

      // tau selection
      if (debug)
	fileOutput << "# taus = " << analysisTree.tau_count << std::endl;
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {
	if (debug)
	  fileOutput << "  " << it 
		     << " pt = " << analysisTree.tau_pt[it] 
		     << " eta = " << analysisTree.tau_eta[it]
		     << " decayModeFinding = " << analysisTree.tau_decayModeFinding[it]
		     << " decayModeFindingNewDMs = " << analysisTree.tau_decayModeFindingNewDMs[it]
		     << " tau_vertexz = " << analysisTree.tau_vertexz[it] <<std::endl;
	if (analysisTree.tau_pt[it]<=ptTauLowCut) continue;
	if (fabs(analysisTree.tau_eta[it])>=etaTauCut) continue;
	if (fabs(fabs(analysisTree.tau_charge[it])-1)>0.001) continue;
	if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>=dzTauCut) continue;
	if (applyTauId &&
	    //analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) 
	    analysisTree.tau_decayModeFinding[it] < 0.5) continue;
	
	taus.push_back(it);
      }

      if (debug) {
	fileOutput << " # selected muons = " << muons.size() << std::endl;
	fileOutput << " # selected taus = " << taus.size() << std::endl;	
      }
      
      if (muons.size()==0) continue;
      if (taus.size()==0) continue;

      // selecting electron and tau pair (OS or SS);
      int muonIndex = -1;
      int tauIndex = -1;
      
      float isoMuMin = 1e+10;
      float isoTauMin = 1e+10;      
      for (unsigned int im=0; im<muons.size(); ++im) {
	unsigned int mIndex  = muons.at(im);

	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	float photonIsoMu = analysisTree.muon_photonIso[im];
	float chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	float puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[im];
	}
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu =  chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
	
	if (debug)
	  fileOutput << "Muon " << mIndex << " -> relIso = "<<relIsoMu<<" absIso = "<<absIsoMu<<std::endl;

	bool isSingleLepTrig = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);

	  if (dRtrig < deltaRTrigMatch){
	    if (analysisTree.trigobject_filters[iT][nIsoLeg] && ( isData || analysisTree.trigobject_pt[iT] > ptTrigObjCut)) // Ele23 Leg
	      isSingleLepTrig = true;
	  }
	}
	  
	if (debug)
	  fileOutput << "Muon " << mIndex << " -> isTrigMu18 = " << isSingleLepTrig << std::endl;

	if (!isSingleLepTrig) continue;
      
	for (unsigned int it=0; it<taus.size(); ++it) {
	  unsigned int tIndex = taus.at(it);
	  //float absIsoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tIndex];
	  //float absIsoTau = analysisTree.tau_byIsolationMVArun2v1DBnewDMwLTraw[tIndex];
	  float absIsoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
	  float relIsoTau = absIsoTau / analysisTree.tau_pt[tIndex];

	  if (debug)
	    fileOutput << "tau" << tIndex << " -> relIso = "<<relIsoTau<<" absIso = "<<absIsoTau<<std::endl;
	  
	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (debug)
	    fileOutput << "dR(mu,tau) = " << dR << std::endl;

	  if (dR<dRleptonsCut) continue;

	  // kinematic match
	  if (analysisTree.muon_pt[mIndex]<=ptMuonHighCut) continue;

	  // change pair
	  bool changePair =  false;
	  if (relIsoMu<isoMuMin) {
	    if(debug)
	      fileOutput<<"ChangePair mu iso"<<std::endl;
	    changePair = true;
	  }
	  else if (fabs(relIsoMu - isoMuMin) < 1.e-5) {
	    if (analysisTree.muon_pt[mIndex] > analysisTree.muon_pt[muonIndex]) {
	      if(debug)
		fileOutput<<"ChangePair ele pt"<<std::endl; 
	      changePair = true;
	    }	    
	    else if (fabs(analysisTree.muon_pt[mIndex] - analysisTree.muon_pt[muonIndex]) < 1.e-5) {
	      if (absIsoTau > isoTauMin) {    // --------------->>> changed here, before was "<" 
		if(debug)
		  fileOutput<<"ChangePair tau iso"<<std::endl; 
		changePair = true;
	      }
	      else if ((absIsoTau - isoTauMin) < 1.e-5){
		if (analysisTree.tau_pt[tIndex] > analysisTree.tau_pt[tauIndex]) {
		  if(debug)
		    fileOutput<<"ChangePair tau pt"<<std::endl; 
		  changePair = true;
		}
	      }
	    }
	  }
	  
	  if (changePair){
	    isoMuMin  = relIsoMu;
	    muonIndex = mIndex;
	    isoTauMin = absIsoTau;
	    tauIndex = tIndex;
	  }
	}
      }
    
      if (muonIndex<0) continue;
      if (tauIndex<0) continue;
      
      if(debug)
	fileOutput << "Selected Pair (m,tau) = "<<muonIndex<<", "<<tauIndex<<std::endl;

      // filling electron variables
      otree->pt_1 = analysisTree.muon_pt[muonIndex];
      otree->eta_1 = analysisTree.muon_eta[muonIndex];
      otree->phi_1 = analysisTree.muon_phi[muonIndex];
      otree->m_1 = muonMass;
      otree->q_1 = -1;
      if (analysisTree.muon_charge[muonIndex]>0)
        otree->q_1 = 1;
      otree->gen_match_1 = analysisTree.muon_genmatch[muonIndex];
      
      float neutralHadIsoMu = analysisTree.muon_neutralHadIso[muonIndex];
      float photonIsoMu = analysisTree.muon_photonIso[muonIndex];
      float chargedHadIsoMu = analysisTree.muon_chargedHadIso[muonIndex];
      float puIsoMu = analysisTree.muon_puIso[muonIndex];
      if (isIsoR03) {
	neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[muonIndex];
	photonIsoMu = analysisTree.muon_r03_sumPhotonEt[muonIndex];
	chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[muonIndex];
	puIsoMu = analysisTree.muon_r03_sumPUPt[muonIndex];
      }
      float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
      neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
      otree->iso_1 = (chargedHadIsoMu + neutralIsoMu)/analysisTree.muon_pt[muonIndex];
      //otree->mva_1 = analysisTree.electron_mva_value_nontrig_Spring15_v1[electronIndex];

      otree->d0_1 = analysisTree.muon_dxy[muonIndex];
      otree->dZ_1 = analysisTree.muon_dz[muonIndex];

      otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = 0;
      otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
      otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
      otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = 0;
      otree->againstElectronLooseMVA5_1 = 0;
      otree->againstElectronMediumMVA5_1 = 0;
      otree->againstElectronTightMVA5_1 = 0;
      otree->againstElectronVLooseMVA5_1 = 0;
      otree->againstElectronVTightMVA5_1 = 0;
      otree->againstElectronVLooseMVA6_2 = 0;
      otree->againstElectronTightMVA6_2 = 0;
      otree->againstMuonLoose3_1 = 0;
      otree->againstMuonTight3_1 = 0;

      // filling electron weights 
      if (!isData) {
      			// Scale Factor SingleEle trigger SF_eleTrigger
      otree->trigweight_1 = (SF_muTrigger->get_ScaleFactor(double(analysisTree.muon_pt[muonIndex]),double(analysisTree.muon_eta[muonIndex])));

      			// Scale Factor Id+Iso SF_eleIdIso
      otree->idisoweight_1 = (SF_muIdIso->get_ScaleFactor(double(analysisTree.muon_pt[muonIndex]),double(analysisTree.muon_eta[muonIndex])));

      otree->effweight = (otree->trigweight_1)*(otree->idisoweight_1)*(otree->trigweight_2)*(otree->idisoweight_2);
      }

      // filling tau variables
      otree->pt_2 = analysisTree.tau_pt[tauIndex];
      otree->eta_2 = analysisTree.tau_eta[tauIndex];
      otree->phi_2 = analysisTree.tau_phi[tauIndex];
      otree->q_2 = analysisTree.tau_charge[tauIndex];
      otree->gen_match_2 = analysisTree.tau_genmatch[tauIndex];
      //if (analysisTree.tau_charge[tauIndex]>0)
      //otree->q_2 = 1;
      //otree->mva_2 = log(0);
      otree->mva_2 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex]; //before, tau_byIsolationMVArun2v1DBoldDMwLTraw
      otree->d0_2 = analysisTree.tau_leadchargedhadrcand_dxy[tauIndex];
      otree->dZ_2 = analysisTree.tau_leadchargedhadrcand_dz[tauIndex];      
      otree->iso_2 = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
      otree->m_2 = analysisTree.tau_mass[tauIndex];
      otree->tau_decay_mode_2 = analysisTree.tau_decayMode[tauIndex];

      otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
      otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree.tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
      otree->againstElectronLooseMVA5_2 = analysisTree.tau_againstElectronLooseMVA5[tauIndex];
      otree->againstElectronMediumMVA5_2 = analysisTree.tau_againstElectronMediumMVA5[tauIndex];
      otree->againstElectronTightMVA5_2 = analysisTree.tau_againstElectronTightMVA5[tauIndex];
      otree->againstElectronVLooseMVA5_2 = analysisTree.tau_againstElectronVLooseMVA5[tauIndex];
      otree->againstElectronVTightMVA5_2 = analysisTree.tau_againstElectronVTightMVA5[tauIndex];
      otree->againstMuonLoose3_2 = analysisTree.tau_againstMuonLoose3[tauIndex];
      otree->againstMuonTight3_2 = analysisTree.tau_againstMuonTight3[tauIndex];
      otree->againstElectronVLooseMVA6_2 = analysisTree.tau_againstElectronVLooseMVA6[tauIndex];
      otree->againstElectronTightMVA6_2 = analysisTree.tau_againstElectronTightMVA6[tauIndex];

      otree->byTightIsolationMVArun2v1DBoldDMwLT_2 = analysisTree.tau_byTightIsolationMVArun2v1DBoldDMwLT[tauIndex];

      // ditau system
      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector tauLV; tauLV.SetXYZM(analysisTree.tau_px[tauIndex],
					  analysisTree.tau_py[tauIndex],
					  analysisTree.tau_pz[tauIndex],
					  tauMass);

      TLorentzVector metLV; metLV.SetXYZT(analysisTree.pfmet_ex,
					  analysisTree.pfmet_ey,
					  0,
					  TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey));
      
      TLorentzVector dileptonLV = muonLV + tauLV;

      // visible mass
      otree->m_vis = dileptonLV.M();
      // visible ditau pt 
      otree->pt_tt = (dileptonLV+metLV).Pt();

      // opposite charge
      otree->os = (otree->q_1 * otree->q_2) < 0.;

      // dilepton veto
      otree->dilepton_veto = 0;

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<=ptDiMuonVeto) continue;
	if (fabs(analysisTree.muon_eta[im])>=etaDiMuonVeto) continue;	
	
	if (fabs(analysisTree.muon_dxy[im])>=dxyDiMuonVeto) continue;
	if (fabs(analysisTree.muon_dz[im])>=dzDiMuonVeto) continue;

	neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	photonIsoMu = analysisTree.muon_photonIso[im];
	chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[im];
	}
	neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu =  chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];

	if(relIsoMu >= isoDiMuonVeto) continue;
	
	bool passedVetoId =  analysisTree.muon_isMedium[im]; 
	if (!passedVetoId && applyDiMuonVetoId) continue;
	
	for (unsigned int je = im+1; je<analysisTree.muon_count; ++je) {
	  if (analysisTree.muon_pt[je]<=ptDiMuonVeto) continue;
	  if (fabs(analysisTree.muon_eta[je])>=etaDiMuonVeto) continue;	
	  
	  if (fabs(analysisTree.muon_dxy[je])>=dxyDiMuonVeto) continue;
	  if (fabs(analysisTree.muon_dz[je])>=dzDiMuonVeto) continue;
	  
	  if (analysisTree.muon_charge[im] * analysisTree.muon_charge[je] > 0. && applyDiMuonOS) continue;

	  neutralHadIsoMu = analysisTree.muon_neutralHadIso[je];
	  photonIsoMu = analysisTree.muon_photonIso[je];
	  chargedHadIsoMu = analysisTree.muon_chargedHadIso[je];
	  puIsoMu = analysisTree.muon_puIso[je];
	  if (isIsoR03) {
	    neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[je];
	    photonIsoMu = analysisTree.muon_r03_sumPhotonEt[je];
	    chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[je];
	    puIsoMu = analysisTree.muon_r03_sumPUPt[je];
	  }
	  neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	  neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	  absIsoMu =  chargedHadIsoMu + neutralIsoMu;
	  relIsoMu = absIsoMu/analysisTree.muon_pt[je];

	  if(relIsoMu >= isoDiMuonVeto) continue;	

	  passedVetoId =  analysisTree.muon_isMedium[je];
	  if (!passedVetoId && applyDiMuonVetoId) continue;
	  
	  float dr = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
			    analysisTree.muon_eta[je],analysisTree.muon_phi[je]);

	  if(dr<=drDiMuonVeto) continue;

	  otree->dilepton_veto = 1;
	}
      }
      
      // extra electron veto
      otree->extraelec_veto = 0;

      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	
	if (analysisTree.electron_pt[ie]<=ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>=etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>=dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>=dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_nontrig_Spring15_v1[ie];
	if (!electronMvaId && applyVetoElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie] && applyVetoElectronId) continue;
	if (analysisTree.electron_nmissinginnerhits[ie]>1 && applyVetoElectronId) continue;
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
	float neutralIsoEle = neutralHadIsoEle + photonIsoEle - 0.5*puIsoEle;
	neutralIsoEle = TMath::Max(float(0),neutralIsoEle); 
	float absIsoEle =  chargedHadIsoEle + neutralIsoEle;
	float relIsoEle = absIsoEle/analysisTree.electron_pt[ie];	

	if (relIsoEle>=isoVetoElectronCut) continue;
	
	otree->extraelec_veto = 1;
      }
            
      // extra muon veto
      otree->extramuon_veto = 0;

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (int(im)==muonIndex) continue;
	if (analysisTree.muon_pt[im]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[im]) continue;
	 neutralHadIsoMu = analysisTree.muon_neutralHadIso[im];
	 photonIsoMu = analysisTree.muon_photonIso[im];
	 chargedHadIsoMu = analysisTree.muon_chargedHadIso[im];
	 puIsoMu = analysisTree.muon_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoMu = analysisTree.muon_r03_sumNeutralHadronEt[im];
	  photonIsoMu = analysisTree.muon_r03_sumPhotonEt[im];
	  chargedHadIsoMu = analysisTree.muon_r03_sumChargedHadronPt[im];
	  puIsoMu = analysisTree.muon_r03_sumPUPt[im];
	}
	 neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu); 
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/analysisTree.muon_pt[im];
	if (relIsoMu>isoVetoMuonCut) continue;

	otree->extramuon_veto = 1;
      }

      // svfit variables
      otree->m_sv = -9999;
      otree->pt_sv = -9999;
      otree->eta_sv = -9999;
      otree->phi_sv = -9999;
      
      // pfmet variables
      otree->met = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);
      otree->metphi = TMath::ATan2(analysisTree.pfmet_ey,analysisTree.pfmet_ex);
      otree->metcov00 = analysisTree.pfmet_sigxx;
      otree->metcov01 = analysisTree.pfmet_sigxy;
      otree->metcov10 = analysisTree.pfmet_sigyx;
      otree->metcov11 = analysisTree.pfmet_sigyy;

      float met_x = analysisTree.pfmet_ex;
      float met_y = analysisTree.pfmet_ey;
      float met_x2 = met_x * met_x;
      float met_y2 = met_y * met_y;     	

      // puppimet variables
      otree->puppimet = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex +
				    analysisTree.puppimet_ey*analysisTree.puppimet_ey);
      otree->puppimetphi = TMath::ATan2(analysisTree.puppimet_ey,analysisTree.puppimet_ex);      
      
      // choosing mva met
      unsigned int iMet = 0;
      float mvamet_x = 0;
      float mvamet_y = 0;
      otree->mvacov00 = 0.;
      otree->mvacov01 = 0.;
      otree->mvacov10 = 0.;
      otree->mvacov11 = 0.;
      
      for (; iMet<analysisTree.mvamet_count; ++iMet) {
	if ((int)analysisTree.mvamet_channel[iMet]==3) break;
      }
      
      if (iMet>=analysisTree.mvamet_count){
	//if (true) {
	if(debug)
	  fileOutput<<"MVA MET channel ploblems.."<<std::endl;

	otree->mvamet = log(0);
	otree->mvametphi = log(0);
 	otree->mvacov00 = log(0);
	otree->mvacov01 = log(0);
	otree->mvacov10 = log(0);
	otree->mvacov11 = log(0);
      }
      else {
	// choosing mva met
	unsigned int iMet = 0;
	for (; iMet<analysisTree.mvamet_count; ++iMet) {
	  if (analysisTree.mvamet_channel[iMet]==3){
	    if( ((int)analysisTree.mvamet_lep1[iMet])==tauIndex && ((int)analysisTree.mvamet_lep2[iMet])==muonIndex)
	    break;
	  }
	}
/*
	if ( fabs(analysisTree.mvamet_lep1_pt[iMet] - otree->pt_2) > 0.0001 )
	  std::cout<<"tau pt does not match"<<std::endl;
	if ( fabs(analysisTree.mvamet_lep2_pt[iMet] - otree->pt_1) > 0.0001 )
	  std::cout<<"ele pt does not match"<<std::endl;
*/
	float mvamet_x = 0;
	float mvamet_y = 0;
	otree->mvacov00 = 0.;
	otree->mvacov01 = 0.;
	otree->mvacov10 = 0.;
	otree->mvacov11 = 0.;
	
	if(iMet < analysisTree.mvamet_count){
	  mvamet_x = analysisTree.mvamet_ex[iMet];
	  mvamet_y = analysisTree.mvamet_ey[iMet];
	  otree->mvacov00 = analysisTree.mvamet_sigxx[iMet];
	  otree->mvacov01 = analysisTree.mvamet_sigxy[iMet];
	  otree->mvacov10 = analysisTree.mvamet_sigyx[iMet];
	  otree->mvacov11 = analysisTree.mvamet_sigyy[iMet];
	}
	else{
	  std::cout<<"MVA MEt not found! :: "<< analysisTree.mvamet_count << std::endl;
	  iMet = 0;
	  //std::cout<<"tau = "<<tauIndex<<" mu = "<<muonIndex<<std::endl;
	  //for (; iMet<analysisTree.mvamet_count; ++iMet) {
	    //std::cout<<"imet = "<<analysisTree.mvamet_channel[iMet]<<" itau = "<<analysisTree.mvamet_lep1[iMet]<<" imu = "<<analysisTree.mvamet_lep2[iMet]<<std::endl;
	  //}
	  
	}
	
	float mvamet_x2 = mvamet_x * mvamet_x;
	float mvamet_y2 = mvamet_y * mvamet_y;
	otree->mvamet = TMath::Sqrt(mvamet_x2+mvamet_y2);
	otree->mvametphi = TMath::ATan2(mvamet_y,mvamet_x);
      }
      
      
      // define MET covariance
      TMatrixD covMET(2, 2);
      covMET[0][0] = otree->mvacov00;
      covMET[1][0] = otree->mvacov10;
      covMET[0][1] = otree->mvacov01;
      covMET[1][1] = otree->mvacov11;

      // mt calculation
      otree->mt_1 = sqrt(2*otree->pt_1*otree->mvamet*(1.-cos(otree->phi_1-otree->mvametphi)));
      otree->mt_2 = sqrt(2*otree->pt_2*otree->mvamet*(1.-cos(otree->phi_2-otree->mvametphi)));
      
      otree->pfmt_1 = sqrt(2*otree->pt_1*otree->met*(1.-cos(otree->phi_1-otree->metphi)));
      otree->pfmt_2 = sqrt(2*otree->pt_2*otree->met*(1.-cos(otree->phi_2-otree->metphi)));

      otree->puppimt_1 = sqrt(2*otree->pt_1*otree->puppimet*(1.-cos(otree->phi_1-otree->puppimetphi)));
      otree->puppimt_2 = sqrt(2*otree->pt_2*otree->puppimet*(1.-cos(otree->phi_2-otree->puppimetphi)));      

      // bisector of muon and tau transverse momenta
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();
	
      float tauUnitX = tauLV.Px()/tauLV.Pt();
      float tauUnitY = tauLV.Py()/tauLV.Pt();

      float zetaX = muonUnitX + tauUnitX;
      float zetaY = muonUnitY + tauUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorVisX = muonLV.Px() + tauLV.Px();
      float vectorVisY = muonLV.Py() + tauLV.Py();

      otree->pzetavis  = vectorVisX*zetaX + vectorVisY*zetaY;
      otree->pzetamiss = otree->mvamet*TMath::Cos(otree->mvametphi)*zetaX + otree->mvamet*TMath::Sin(otree->mvametphi)*zetaY;
      otree->pfpzetamiss = analysisTree.pfmet_ex*zetaX + analysisTree.pfmet_ey*zetaY;      
      otree->puppipzetamiss = analysisTree.puppimet_ex*zetaX + analysisTree.puppimet_ey*zetaY;

      metLV.SetXYZT(otree->mvamet*TMath::Cos(otree->mvametphi),
		    otree->mvamet*TMath::Sin(otree->mvametphi),
		    0,
		    TMath::Sqrt( otree->mvamet*TMath::Sin(otree->mvametphi)*otree->mvamet*TMath::Sin(otree->mvametphi) +
				 otree->mvamet*TMath::Cos(otree->mvametphi)*otree->mvamet*TMath::Cos(otree->mvametphi)));
      otree->pt_tt = (dileptonLV+metLV).Pt();      
      
      // define lepton four vectors
      //std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
      //measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, otree->pt_1, otree->eta_1,  otree->phi_1, 0.51100e-3)); // tau -> electron decay (Pt, eta, phi, mass)
      //measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, otree->pt_1, otree->eta_1,  otree->phi_1, otree->m_2, 0));//analysisTree.tau_decayMode[tauIndex])); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
      //SVfitStandaloneAlgorithm algo(measuredTauLeptons, analysisTree.mvamet_ex[iMet], analysisTree.mvamet_ey[iMet], covMET, 0);
      //algo.addLogM(false);  
      //algo.shiftVisPt(true, inputFile_visPtResolution);
      //algo.integrateMarkovChain();

      //otree->m_sv = algo.getMass(); // return value is in units of GeV
      
      // counting jets
      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();

      int indexLeadingJet = -1;
      float ptLeadingJet = -1;

      int indexSubLeadingJet = -1;
      float ptSubLeadingJet = -1;
      
      int indexLeadingBJet = -1;
      float ptLeadingBJet = -1;

      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta>=jetEtaCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];
	if (jetPt<=jetPtLowCut) continue;

	float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			   otree->eta_1,otree->phi_1);
	if (dR1<=dRJetLeptonCut) continue;

	float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                           otree->eta_2,otree->phi_2);
        if (dR2<=dRJetLeptonCut) continue;

	// jetId
	float energy = analysisTree.pfjet_e[jet];
        energy *= analysisTree.pfjet_energycorr[jet];
        float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
        float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
        float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
        float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
	float muf = analysisTree.pfjet_muonenergy[jet]/energy;
        float chm = analysisTree.pfjet_chargedmulti[jet];
	float nm = analysisTree.pfjet_neutralmulti[jet];
        float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
	//bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
	bool isPFJetId = false;
	if (absJetEta<=3.0)
	  isPFJetId = (nhf < 0.99 && phf < 0.99 && npr > 1) && (absJetEta>2.4 || (chf>0 && chm > 0 && elf < 0.99));
	else
	  isPFJetId = phf < 0.9 && nm > 10;
	//isPFJetId = (npr>1 && phf<0.99 && nhf<0.99 && muf < 0.8) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
	//isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
	
	if (!isPFJetId) continue;

	jetspt20.push_back(jet);

	if (absJetEta<bJetEtaCut && analysisTree.pfjet_btag[jet][6]>btagCut) { // b-jet
	  bjets.push_back(jet);
	  if (jetPt>ptLeadingBJet) {
	    ptLeadingBJet = jetPt;
	    indexLeadingBJet = jet;
	  }
	} 

	if (indexLeadingJet>=0) {
	  if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
	    indexSubLeadingJet = jet;
	    ptSubLeadingJet = jetPt;
	  }
	}

	if (jetPt>ptLeadingJet) {
	  indexLeadingJet = jet;
	  ptLeadingJet = jetPt;
	}

	if (jetPt<jetPtHighCut) continue;
	jets.push_back(jet);
      }
	
      otree->njets = jets.size();
      otree->njetspt20 = jetspt20.size();
      otree->nbtag = bjets.size();
      
      otree->bpt = -9999;
      otree->beta = -9999;
      otree->bphi = -9999;
      
      if (indexLeadingBJet>=0) {
	otree->bpt = analysisTree.pfjet_pt[indexLeadingBJet];
	otree->beta = analysisTree.pfjet_eta[indexLeadingBJet];
	otree->bphi = analysisTree.pfjet_phi[indexLeadingBJet];
      }

      otree->jpt_1 = -9999;
      otree->jeta_1 = -9999;
      otree->jphi_1 = -9999;
      otree->jptraw_1 = -9999;
      otree->jptunc_1 = -9999;
      otree->jmva_1 = -9999;
      otree->jlrm_1 = -9999;
      otree->jctm_1 = -9999;

      if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
	cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;

      if (indexLeadingJet>=0) {
	otree->jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
	otree->jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
	otree->jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
	otree->jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
	otree->jmva_1 = analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
      }

      otree->jpt_2 = -9999;
      otree->jeta_2 = -9999;
      otree->jphi_2 = -9999;
      otree->jptraw_2 = -9999;
      otree->jptunc_2 = -9999;
      otree->jmva_2 = -9999;
      otree->jlrm_2 = -9999;
      otree->jctm_2 = -9999;

      if (indexSubLeadingJet>=0) {
	otree->jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
	otree->jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
	otree->jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
	otree->jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];
	otree->jmva_2 = analysisTree.pfjet_pu_jet_full_mva[indexSubLeadingJet];
      }

      otree->mjj =  -9999;
      otree->jdeta =  -9999;
      otree->njetingap = -1;

      if (indexLeadingJet>=0 && indexSubLeadingJet>=0) {
	otree->njetingap = 0;
	TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
					     analysisTree.pfjet_py[indexLeadingJet],
					     analysisTree.pfjet_pz[indexLeadingJet],
					     analysisTree.pfjet_e[indexLeadingJet]);

	TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
					     analysisTree.pfjet_py[indexSubLeadingJet],
					     analysisTree.pfjet_pz[indexSubLeadingJet],
					     analysisTree.pfjet_e[indexSubLeadingJet]);

	otree->mjj = (jet1+jet2).M();
	otree->jdeta = abs(analysisTree.pfjet_eta[indexLeadingJet]-
		    analysisTree.pfjet_eta[indexSubLeadingJet]);
 
	float etamax = analysisTree.pfjet_eta[indexLeadingJet];
	float etamin = analysisTree.pfjet_eta[indexSubLeadingJet];
	if (etamax<etamin) {
	  float tmp = etamax;
	  etamax = etamin;
	  etamin = tmp;
	}
	for (unsigned int jet=0; jet<jets.size(); ++jet) {
	  int index = jets.at(jet);
	  float etaX = analysisTree.pfjet_eta[index];
	  if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax) 
	    otree->njetingap++;
	}
      }
      otree->Fill();
      selEvents++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}
