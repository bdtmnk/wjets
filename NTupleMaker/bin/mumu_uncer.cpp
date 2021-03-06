//
//  mumu_uncer.cpp
//  
//
//  Created by Vallary on 2/1/17.
//
//

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <TF1.h>
#include <iomanip>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TCut.h"

#include "TRandom.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/EventWeight.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/rochcor2015.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"

#include "DesyTauAnalyses/NTupleMaker/interface/btagSF.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"

float topPtWeight(float pt1,
                  float pt2) {
    
    float a = 0.0615;    // Run1 a parameter
    float b = -0.0005;  // Run1 b parameter
    
    //    if (pt1>400) pt1 = 400;
    //if (pt2>400) pt2 = 400;
    
    float w1 = TMath::Exp(a+b*pt1);
    float w2 = TMath::Exp(a+b*pt2);
    
    return TMath::Sqrt(w1*w2);
    
}

SVfitStandaloneAlgorithm SVFitMassComputation(svFitStandalone::MeasuredTauLepton svFitMu1,
                                              svFitStandalone::MeasuredTauLepton svFitMu2,
                                              double measuredMETx,
                                              double measuredMETy,
                                              TMatrixD covMVAMET,
					      TFile * inputFile_visPtResolution){
  
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitMu1);
  measuredTauLeptons.push_back(svFitMu2);
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMVAMET, 0);
  algo.addLogM(false);
  algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();
  
  return algo;
  
}

bool metFiltersPasses(AC1B &tree_, std::vector<TString> metFlags) {
    
    bool passed = true;
    unsigned int nFlags = metFlags.size();
    //  std::cout << "MEt filters : " << std::endl;
    for (std::map<string,int>::iterator it=tree_.flags->begin(); it!=tree_.flags->end(); ++it) {
        TString flagName(it->first);
        //    std::cout << it->first << " : " << it->second << std::endl;
        for (unsigned int iFilter=0; iFilter<nFlags; ++iFilter) {
            if (flagName.Contains(metFlags[iFilter])) {
                if (it->second==0) {
                    passed = false;
                    break;
                }
            }
        }
    }
    //  std::cout << "Passed : " << passed << std::endl;
    return passed;
    
}

float nJetsWeight(int nJets) {
    
    float weight = 1;
    if (nJets==0)
        weight = 1.02;
    else if (nJets==1)
        weight = 0.95;
    else
        weight = 0.93;
    
    return weight;
    
}

void computeDzeta(float metX,  float metY,
                  float zetaX, float zetaY,
                  float pzetavis,
                  float & pzetamiss,
                  float & dzeta) {
    
    pzetamiss = metX*zetaX + metY*zetaY;
    dzeta = pzetamiss - 0.85*pzetavis;
    
}


void computeRecoil(float metx, float mety,
                   float unitX,float unitY,
                   float perpUnitX, float perpUnitY,
                   float dimuonPt,
                   float & recoilParal,
                   float & recoilPerp,
                   float & responseHad) {
    
    recoilParal = metx*unitX + mety*unitY;
    recoilPerp = metx*perpUnitX + mety*perpUnitY;
    responseHad = 1 + recoilParal/dimuonPt;
}

bool isICHEPmed(unsigned int Index, AC1B analysisTree) {
    bool goodGlob =
    analysisTree.muon_isGlobal[Index] &&
    analysisTree.muon_normChi2[Index] < 3 &&
    analysisTree.muon_combQ_chi2LocalPosition[Index] < 12 &&
    analysisTree.muon_combQ_trkKink[Index] < 20;
    bool isICHEPmedium  =
    analysisTree.muon_isLoose[Index] &&
    analysisTree.muon_validFraction[Index] >0.49 &&
    analysisTree.muon_segmentComp[Index] > (goodGlob ? 0.303 : 0.451);
    return isICHEPmedium;
}
//------tracking eff--------------------->

double eff(double eta){
    static const double eff[10] = {0.982399, 0.991747, 0.995945, 0.993413, 0.991461, 0.99468, 0.996666, 0.994934, 0.991187,0.976812};
    int region = 0;
    
    if (eta > -2.2309) region++;
    if (eta > -1.827) region++;
    if (eta > -1.34607 ) region++;
    if (eta > -0.843046) region++;
    if (eta > -0.297941) region++;
    if (eta > 0.298253) region++;
    if (eta > 0.843136) region++;
    if (eta > 1.34753) region++;
    if (eta > 1.82701) region++;
    if (eta > 2.2333) region++;
    return eff[region];
}
//-------------main function starts here---------------->

int main(int argc, char * argv[]) {
    //first argument - config file
    // second argument - filelist
    
    using namespace std;
    TH1::SetDefaultSumw2(true);
    TH2::SetDefaultSumw2(true);
    
    //---------configuration----------------------------->
    Config cfg(argv[1]);
    
    const bool isData = cfg.get<bool>("IsData");
    
    //Run-lumi selector
    const string jsonFile = cfg.get<string>("jsonFile");
    
    //---------corrections---------->
    
    //pile up reweighting
    const bool applyPUreweighting_official = cfg.get<bool>("ApplyPUreweighting_official");
    
    const bool applyMEtRecoilCorrections = cfg.get<bool>("ApplyMEtRecoilCorrections");
    const bool applyLeptonSF = cfg.get<bool>("ApplyLeptonSF");
    const bool applyTrackEff = cfg.get<bool>("ApplyTrackEff");
    
    const bool applyTopPtReweighting = cfg.get<bool>("ApplyTopPtReweighting");
    const bool applyBDT = cfg.get<bool>("ApplyBDT");
    const bool applyBDTvariations = cfg.get<bool>("ApplyBDTvariation");
    
    const bool applySVFit = cfg.get<bool>("ApplySVFit");
    const bool applyMSVvariation = cfg.get<bool>("ApplyMSVvariation");
    
    const bool applyZptmassCorr = cfg.get<bool>("ApplyZptmassCorr");
    
    const bool applyCategoryWeights = cfg.get<bool>("ApplyCategoryWeights");
    
    //ztotautautomumu selection
    const bool  applyTauTauSelection = cfg.get<bool>("ApplyTauTauSelection");
    const bool  selectZToTauTauMuMu = cfg.get<bool>("SelectZToTauTauMuMu");
    
    // kinematic cuts on muons
    const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
    const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
    const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
    const float etaMuonLowCut = cfg.get<float>("etaMuonLowCut");
    const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
    const float dzMuonCut      = cfg.get<float>("dzMuonCut");
    const float dxyMuonLooseCut     = cfg.get<float>("dxyMuonLooseCut");
    const float dzMuonLooseCut      = cfg.get<float>("dzMuonLooseCut");
    const float isoMuonCut     = cfg.get<float>("isoMuonCut");
    
    // vertex cuts
    const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
    const float zVertexCut     = cfg.get<float>("ZVertexCut");
    const float dVertexCut     = cfg.get<float>("DVertexCut");
    
    // topological cuts
    const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
    const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
    //const bool oppositeSign    = cfg.get<bool>("OppositeSign");
    const float dimuonMassCut = cfg.get<float>("DimuonMassCut"); // aug 17
    
    // jets
    const string bTagDiscriminator = cfg.get<string>("BTagDiscriminator");
    TString BTagDiscriminator(bTagDiscriminator);
    const float jetEtaCut = cfg.get<float>("JetEtaCut");
    const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
    const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
    const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
    const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
    const float btagCut = cfg.get<float>("btagCut");
    const bool applyJetPfId = cfg.get<bool>("ApplyJetPfId");
    const bool applyJetPuId = cfg.get<bool>("ApplyJetPuId");
    const float jetEtaTrkCut   = cfg.get<float>("JetEtaTrkCut");
    
    // trigger
    const string muonTriggerName = cfg.get<string>("MuonTriggerName");
    const string muonFilterName = cfg.get<string>("MuonFilterName");
    TString MuonTriggerName(muonTriggerName);
    TString MuonFilterName(muonFilterName);
    
    const string MuonTrigFile = cfg.get<string>("MuonTrigEff");
    const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoEff");
    
    const string singleMuonFilterName = cfg.get<string>("SingleMuonFilterName");
    TString SingleMuonFilterName(singleMuonFilterName);
    const float singleMuonTriggerPtCut = cfg.get<float>("SingleMuonTriggerPtCut");
    const float singleMuonTriggerEtaCut = cfg.get<float>("SingleMuonTriggerEtaCut");
    
    const string recoilFileName   = cfg.get<string>("RecoilFileName");
    TString RecoilFileName(recoilFileName);
    
    const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
    const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");
    
    TString PileUpDataFile(pileUpDataFile);
    TString PileUpMCFile(pileUpMCFile);
    
    const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName");
    TString ZMassPtWeightsFileName(zMassPtWeightsFileName);
    
    const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName");
    TString ZMassPtWeightsHistName(zMassPtWeightsHistName);
    
    const string categoryWeightsFileName   = cfg.get<string>("CategoryWeightsFileName");
    TString CategoryWeightsFileName(categoryWeightsFileName);
    
    const string jet0WeightsHist = cfg.get<string>("Jet0WeightsHist");
    TString Jet0WeightsHist(jet0WeightsHist);
    const string boostedWeightsHist  =  cfg.get<string>("BoostedWeightsHist");
    TString BoostedWeightsHist(boostedWeightsHist);
    const string vbfWeightsHist      =  cfg.get<string>("VBFWeightsHist");
    TString VBFWeightsHist(vbfWeightsHist);
    
    
    // Run range
    const unsigned int RunRangeMin = cfg.get<unsigned int>("RunRangeMin");
    const unsigned int RunRangeMax = cfg.get<unsigned int>("RunRangeMax");
    
    const float muonScale = cfg.get<float>("MuonScale");
    
    //--------------end of configuration---------------->
    
    string cmsswBase = (getenv ("CMSSW_BASE"));
    
    TString fullDir = TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/");
    
    //---------file name and tree name----------->
    std::string rootFileName(argv[2]);
    std::ifstream fileList(argv[2]);
    std::ifstream fileList0(argv[2]);
    std::string ntupleName("makeroottree/AC1B");
    
    TString TStrName(rootFileName);
    std::cout <<TStrName <<std::endl;
    
    TFile * file = new TFile(TStrName+TString(".root"),"recreate");
    file->cd("");
    
    //----------defining histograms---------------->
    TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
    TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
    TH1D * histZTTGenWeightsH = new TH1D("histZTTGenWeightsH","",1,-0.5,0.5);
    TH1D * histGenCutsWeightsH = new TH1D("histGenCutsWeightsH","",1,-0.5,0.5);
    TH1D * histGenCutsGenWeightsH = new TH1D("histGenCutsGenWeightsH","",1,-0.5,0.5);
    TH1D * histRecCutsGenWeightsH = new TH1D("histRecCutsGenWeightsH","",1,-0.5,0.5);
    TH1D * histRecCutsWeightsH = new TH1D("histRecCutsWeightsH","",1,-0.5,0.5);
    TH1D * histBDTCutGenWeightsH = new TH1D("histBDTCutGenWeightsH","",1,-0.5,0.5);
    TH1D * histBDTCutWeightsH = new TH1D("histBDTCutWeightsH","",1,-0.5,0.5);
    TH1D * histGenWeightH = new TH1D("histGenWeightH","",1,-0.5,0.5);
    
    // fout
    TH1D * noutDenRecCutsWeightsH = new TH1D("noutDenRecCutsWeightsH","",1,-0.5,0.5);
    TH1D * noutNumRecCutsWeightsH = new TH1D("noutNumRecCutsWeightsH","",1,-0.5,0.5);
    TH1D * noutDenBDTCutWeightsH = new TH1D("noutDenBDTCutWeightsH","",1,-0.5,0.5);
    TH1D * noutNumBDTCutWeightsH = new TH1D("noutNumBDTCutWeightsH","",1,-0.5,0.5);
    TH1D * noutDenRecCutsGenWeightsH = new TH1D("noutDenRecCutsGenWeightsH","",1,-0.5,0.5);
    TH1D * noutNumRecCutsGenWeightsH = new TH1D("noutNumRecCutsGenWeightsH","",1,-0.5,0.5);
    TH1D * noutDenBDTCutGenWeightsH = new TH1D("noutDenBDTCutGenWeightsH","",1,-0.5,0.5);
    TH1D * noutNumBDTCutGenWeightsH = new TH1D("noutNumBDTCutGenWeightsH","",1,-0.5,0.5);
    
    TH1D * MuSF_IdIso_Mu1H = new TH1D("MuIdIsoSF_Mu1H", "MuIdIsoSF_Mu1", 100, 0.5,1.5);
    TH1D * MuSF_IdIso_Mu2H = new TH1D("MuIdIsoSF_Mu2H", "MuIdIsoSF_Mu2", 100, 0.5,1.5);
    
    TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
    TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

    //Declaration of branch types
    Int_t           run;
    Int_t           lumi;
    Int_t           evt;
    Int_t           npv;
    Int_t           npu;
    Float_t         rho;
    
    
    Bool_t          genAccept;
    Float_t         noOfvertices;
    Float_t         genweight;
    Float_t         mcweight;
    Float_t         puweight;
    Float_t         trigweight;
    Float_t         weightTrig_gen;
    Float_t         effweight_gen;
    Float_t         idweight_1;
    Float_t         idweight_2;
    Float_t         isoweight_1;
    Float_t         isoweight_2;
    Float_t         effweight;
    Float_t         topptweight;
    Float_t         trkeffweight;
    Float_t         zptmassweight;
    Float_t         jet0weight;
    Float_t         boostedweight;
    Float_t         vbfweight;
    Float_t         weight;
    Float_t         btag0weight;
    Float_t         btag0weight_Up;
    Float_t         btag0weight_Down;

    Bool_t metFilters_;

    Bool_t badChargedCandidateFilter_;
    Bool_t badPFMuonFilter_;
    Bool_t badGlobalMuonFilter_;
    Bool_t muonBadTrackFilter_;
    Bool_t chargedHadronTrackResolutionFilter_;

    Bool_t badMuonFilter_;
    Bool_t duplicateMuonFilter_;
    
    Float_t         n_genZ_mass;
    Float_t         n_genZ_Pt;
    Float_t         n_genZ_Eta;
    Float_t         n_genZ_Phi;
    Float_t         n_genZTT_mass;
    Float_t         n_genZTT_Pt;
    Float_t         n_genZTT_Eta;
    Float_t         n_genZTT_Phi;
    Float_t         n_genV_mass;
    Float_t         n_genV_Pt;
    Float_t         n_genV_Eta;
    Float_t         n_genV_Phi;
    Float_t         n_gen_mu1_Pt;
    Float_t         n_gen_mu1_Eta;
    Float_t         n_gen_mu1_Phi;
    Float_t         n_gen_mu2_Pt;
    Float_t         n_gen_mu2_Eta;
    Float_t         n_gen_mu2_Phi;
    UInt_t          n_gen_taus;
    UInt_t          n_gen_mutaus;
    
    Float_t         m_vis;            //BDT discriminator (dimuonMass)
    Float_t         m_vis_Up;
    Float_t         m_vis_Down;
    Float_t         m_vis_scaleUp;
    Float_t         m_vis_scaleDown;
    Float_t         m_vis_resoUp;
    Float_t         m_vis_resoDown;
    
    Float_t         m_sv;
    Float_t         pt_sv;            //BDT discriminator (boosted)
    Float_t         eta_sv;
    Float_t         phi_sv;
    
    Float_t         m_sv_muUp;
    Float_t         m_sv_muDown;
    Float_t         m_sv_scaleUp;
    Float_t         m_sv_scaleDown;
    Float_t         m_sv_resoUp;
    Float_t         m_sv_resoDown;

    Float_t         pt_sv_muUp;            
    Float_t         eta_sv_muUp;
    Float_t         phi_sv_muUp;
	
		Float_t         pt_sv_scaleUp;
		Float_t         eta_sv_scaleUp;
		Float_t         phi_sv_scaleUp;
	
		Float_t         pt_sv_resoUp;
		Float_t         eta_sv_resoUp;
		Float_t         phi_sv_resoUp;

    Float_t         pt_sv_muDown;
    Float_t         eta_sv_muDown;
    Float_t         phi_sv_muDown;
	
		Float_t         pt_sv_scaleDown;
		Float_t         eta_sv_scaleDown;
		Float_t         phi_sv_scaleDown;
		
		Float_t         pt_sv_resoDown;
		Float_t         eta_sv_resoDown;
		Float_t         phi_sv_resoDown;
	
    Float_t         pt_1;
    Float_t         pt_1_Up;
    Float_t         pt_1_Down;
    
    Float_t         eta_1;
    Float_t         phi_1;
    Int_t           q_1;
    Float_t         d0_1;
    Float_t         dZ_1;
    Float_t         dcaSigd0_1;
    Float_t         dcaSigdZ_1;
    
    Float_t         pt_2;
    Float_t         pt_2_Up;
    Float_t         pt_2_Down;
    
    Float_t         eta_2;
    Float_t         phi_2;
    Int_t           q_2;
    Float_t         d0_2;
    Float_t         dZ_2;
    Float_t         dcaSigd0_2;
    Float_t         dcaSigdZ_2;
    
    Bool_t          os;
    
    Float_t         met;              //BDT discriminator
    Float_t         metphi;
    Float_t         metcov00;
    Float_t         metcov01;
    Float_t         metcov10;
    Float_t         metcov11;
    
    Float_t         met_uncorr;
    Float_t         metphi_uncorr;

    Float_t         met_Up;
    Float_t         met_Down;

    Float_t         met_JES_Up_uncorr;
    Float_t         met_JES_Down_uncorr;
    
    Float_t         met_JES_Up;
    Float_t         met_JES_Down;
    
    Float_t         met_UnclusteredJES_Up_uncorr;
    Float_t         met_UnclusteredJES_Down_uncorr;
    
    Float_t         met_UnclusteredJES_Up;
    Float_t         met_UnclusteredJES_Down;
    
    Float_t         genmet;
    Float_t         genmetphi;
    
    Float_t         msvmet;
    Float_t         msvmetphi;

    Float_t         met_x;
    Float_t         met_y;
    Float_t         px_mu1;
    Float_t         py_mu1;
    Float_t         px_mu2;
    Float_t         py_mu2;
    Float_t         pt_tot_x;
    Float_t         pt_tot_y;
    Float_t         pt_tot;
    Float_t         pt_tot_Up;
    Float_t         pt_tot_Down;
    Float_t         pt_tot_JES_Up;
    Float_t         pt_tot_JES_Down;
    Float_t         pt_tot_UnclusteredJES_Up;
    Float_t         pt_tot_UnclusteredJES_Down;
    
    //dimuon system
    Float_t         pt_tt;
    Float_t         dr_tt;
    Float_t         dphi_tt;
    Float_t         eta_tt;             //BDT discriminator (dimuonEta)
    Float_t         ptRatio;          //BDT discriminator

    Float_t         pt_tt_Up;
    Float_t         dr_tt_Up;
    Float_t         dphi_tt_Up;
    Float_t         eta_tt_Up;             
    Float_t         ptRatio_Up;

    Float_t         pt_tt_Down;
    Float_t         dr_tt_Down;
    Float_t         dphi_tt_Down;
    Float_t         eta_tt_Down;             
    Float_t         ptRatio_Down;
    
    Float_t         dcasig2Mu2D;
    Float_t         dcasig2Mu3D;
    Float_t         sig2Mu2D;
    Float_t         sig2Mu3D;
    
    Float_t         pzetavis;
    
    Float_t         pzetamiss;
    Float_t         dzeta;
    
    Float_t         pzetamiss_uncorr;
    Float_t         dzeta_uncorr;
    
    Float_t         pzetamiss_Up;
    Float_t         dzeta_Up;
    Float_t         pzetamiss_Down;
    Float_t         dzeta_Down;
    Float_t         pzetamiss_JES_Up;
    Float_t         dzeta_JES_Up;
    Float_t         pzetamiss_JES_Down;
    Float_t         dzeta_JES_Down;
    
    Float_t         pzetamiss_UnclusteredJES_Up;
    Float_t         dzeta_UnclusteredJES_Up;
    
    Float_t         pzetamiss_UnclusteredJES_Down;
    Float_t         dzeta_UnclusteredJES_Down;
    
    Float_t         pzetamiss_genmet;
    Float_t         dzeta_genmet;
    
    Float_t         dphi_mumet_1;
    Float_t         dphi_mumet_uncorr_1;
    Float_t         dphi_mumet_Up_1;
    Float_t         dphi_mumet_Down_1;
    Float_t         dphi_mumet_JES_Up_1;
    Float_t         dphi_mumet_JES_Down_1;
    Float_t         dphi_mumet_UnclusteredJES_Up_1;
    Float_t         dphi_mumet_UnclusteredJES_Down_1;
    
    Float_t         dphi_mumet_2;
    Float_t         dphi_mumet_uncorr_2;
    Float_t         dphi_mumet_Up_2;
    Float_t         dphi_mumet_Down_2;
    Float_t         dphi_mumet_JES_Up_2;
    Float_t         dphi_mumet_JES_Down_2;
    Float_t         dphi_mumet_UnclusteredJES_Up_2;
    Float_t         dphi_mumet_UnclusteredJES_Down_2;
    
    Float_t         dphi_posmu_met;
    Float_t         dphi_posmu_met_uncorr;
    Float_t         dphi_posmu_met_Up;
    Float_t         dphi_posmu_met_Down;
    Float_t         dphi_posmu_met_JES_Up;
    Float_t         dphi_posmu_met_JES_Down;
    Float_t         dphi_posmu_met_UnclusteredJES_Up;
    Float_t         dphi_posmu_met_UnclusteredJES_Down;
    
    Float_t         dphi_twomu;
    
    Float_t         costheta_1;
    Float_t         costheta_2;
    Float_t         costheta;
    
    Float_t         pos_mupt1;
    Float_t         neg_mupt1;
    Float_t         pos_mueta1;
    Float_t         neg_mueta1;
    Float_t         pos_mupt2;
    Float_t         neg_mupt2;
    Float_t         pos_mueta2;
    Float_t         neg_mueta2;
    
    Int_t           njets;
    Int_t           njets_Up;
    Int_t           njets_Down;
    
    Int_t           njetspt20;
    
    Float_t         jpt_1;
    Float_t         jpt_1_Up;
    Float_t         jpt_1_Down;
    
    Float_t         jeta_1;
    Float_t         jphi_1;
    Float_t         jptraw_1;
    Float_t         jptunc_1;
    Float_t         jmva_1;
    Float_t         jlrm_1;
    Int_t           jctm_1;
    Int_t           gen_match_1;
    
    Float_t         jpt_2;
    Float_t         jpt_2_Up;
    Float_t         jpt_2_Down;
    
    Float_t         jeta_2;
    Float_t         jphi_2;
    Float_t         jptraw_2;
    Float_t         jptunc_2;
    Float_t         jmva_2;
    Float_t         jlrm_2;
    Int_t           jctm_2;
    Int_t           gen_match_2;
    
    Float_t         mjj;
    Float_t         mjj_Up;
    Float_t         mjj_Down;
    
    Float_t         jdeta;
    Int_t           njetingap;
    
    Int_t           nbtag;
    Int_t           nbtag_noSF;
    Int_t           nbtag_Up;
    Int_t           nbtag_Down;
    Float_t         bpt;
    Float_t         beta;
    Float_t         bphi;
    
    UInt_t          npartons;

    Float_t         bdt;
    Float_t         bdt_Up;
    Float_t         bdt_Down;
    Float_t         bdt_JES_Up;
    Float_t         bdt_JES_Down;
    Float_t         bdt_UnclusteredJES_Up;
    Float_t         bdt_UnclusteredJES_Down;
    
    Float_t         bdt_0jets;
    Float_t         bdt_0jets_Up;
    Float_t         bdt_0jets_Down;
    Float_t         bdt_0jets_JES_Up;
    Float_t         bdt_0jets_JES_Down;
    Float_t         bdt_0jets_UnclusteredJES_Up;
    Float_t         bdt_0jets_UnclusteredJES_Down;
    
    Float_t         bdt_boosted;
    Float_t         bdt_boosted_Up;
    Float_t         bdt_boosted_Down;
    Float_t         bdt_boosted_JES_Up;
    Float_t         bdt_boosted_JES_Down;
    Float_t         bdt_boosted_UnclusteredJES_Up;
    Float_t         bdt_boosted_UnclusteredJES_Down;
    
    Float_t         bdt_vbf;
    Float_t         bdt_vbf_Up;
    Float_t         bdt_vbf_Down;
    Float_t         bdt_vbf_JES_Up;
    Float_t         bdt_vbf_JES_Down;
    Float_t         bdt_vbf_UnclusteredJES_Up;
    Float_t         bdt_vbf_UnclusteredJES_Down;

    
    TTree * TW = new TTree("TW","Weights");
    TW->Branch("genweight",&genweight,"genweight/F");
    TW->Branch("genZ_mass", &n_genZ_mass, "n_genZ_mass/F");
    TW->Branch("genZ_Pt", &n_genZ_Pt, "n_genZ_Pt/F");
    TW->Branch("genZ_Eta", &n_genZ_Eta, "n_genZ_Eta/F");
    TW->Branch("genZ_Phi", &n_genZ_Phi, "n_genZ_Phi/F");
    TW->Branch("genZTT_mass", &n_genZTT_mass, "n_genZTT_mass/F");
    TW->Branch("genZTT_Pt", &n_genZTT_Pt, "n_genZTT_Pt/F");
    TW->Branch("genZTT_Eta", &n_genZTT_Eta, "n_genZTT_Eta/F");
    TW->Branch("genZTT_Phi", &n_genZTT_Phi, "n_genZTT_Phi/F");
    TW->Branch("genV_mass", &n_genV_mass, "n_genV_mass/F");
    TW->Branch("genV_Pt", &n_genV_Pt, "n_genV_Pt/F");
    TW->Branch("genV_Eta", &n_genV_Eta, "n_genV_Eta/F");
    TW->Branch("genV_Phi", &n_genV_Phi, "n_genV_Phi/F");
    TW->Branch("gen_mu1_Pt", &n_gen_mu1_Pt, "n_gen_mu1_Pt/F");
    TW->Branch("gen_mu1_Eta", &n_gen_mu1_Eta, "n_gen_mu1_Eta/F");
    TW->Branch("gen_mu1_Phi", &n_gen_mu1_Phi, "n_gen_mu1_Phi/F");
    TW->Branch("gen_mu2_Pt", &n_gen_mu2_Pt, "n_gen_mu2_Pt/F");
    TW->Branch("gen_mu2_Eta", &n_gen_mu2_Eta, "n_gen_mu2_Eta/F");
    TW->Branch("gen_mu2_Phi", &n_gen_mu2_Phi, "n_gen_mu2_Phi/F");
    TW->Branch("gen_taus",&n_gen_taus,"n_gen_taus/i");
    TW->Branch("gen_mutaus",&n_gen_mutaus,"n_gen_mutaus/i");

    TTree * T = new TTree("T","mumu channel");
    T->Branch("run", &run, "run/I");
    T->Branch("lumi", &lumi, "lumi/I");
    T->Branch("evt", &evt, "evt/I");
    T->Branch("npv", &npv, "npv/I");
    T->Branch("npu", &npu, "npu/I");
    T->Branch("rho", &rho, "rho/F");
    
    T->Branch("genAccept", &genAccept, "genAccept/O");
    T->Branch("noOfvertices", &noOfvertices, "noOfvertices/F");
    T->Branch("genweight",&genweight,"genweight/F");
    T->Branch("mcweight", &mcweight, "mcweight/F");
    T->Branch("puweight", &puweight, "puweight/F");
    T->Branch("trigweight", &trigweight, "trigweight/F");
    T->Branch("weightTrig_gen", &weightTrig_gen, "weightTrig_gen/F");
    T->Branch("effweight_gen", &effweight_gen, "effweight_gen/F");
    T->Branch("idweight_1", &idweight_1, "idweight_1/F");
    T->Branch("idweight_2", &idweight_2, "idweight_2/F");
    T->Branch("isoweight_1", &isoweight_1, "isoweight_1/F");
    T->Branch("isoweight_2", &isoweight_2, "isoweight_2/F");
    T->Branch("effweight", &effweight, "effweight/F");
    T->Branch("topptweight", &topptweight, "topptweight/F");
    T->Branch("trkeffweight", &trkeffweight, "trkeffweight/F");
    T->Branch("zptmassweight",&zptmassweight,"zptmassweight/F");
    T->Branch("jet0weight", &jet0weight, "jet0weight/F");
    T->Branch("boostedweight", &boostedweight, "boostedweight/F");
    T->Branch("vbfweight", &vbfweight, "vbfweight/F");
    T->Branch("weight", &weight, "weight/F");
    T->Branch("btag0weight",&btag0weight,"btag0weight/F");
    T->Branch("btag0weight_Up",&btag0weight_Up,"btag0weight_Up/F");
    T->Branch("btag0weight_Down",&btag0weight_Down,"btag0weight_Down/F");

    T->Branch("metFilters",&metFilters_,"metFilters/O");
   
    T->Branch("badChargedCandidateFilter",&badChargedCandidateFilter_,"badChargedCandidateFilter/O");
    T->Branch("badPFMuonFilter",&badPFMuonFilter_,"badPFMuonFilter/O");
    T->Branch("badGlobalMuonFilter",&badGlobalMuonFilter_,"badGlobalMuonFilter/O");
    T->Branch("muonBadTrackFilter",&muonBadTrackFilter_,"muonBadTrackFilter/O");
    T->Branch("chargedHadronTrackResolutionFilter",&chargedHadronTrackResolutionFilter_,"chargedHadronTrackResolutionFilter/O");

    T->Branch("badMuonFilter",&badMuonFilter_,"badMuonFilter/O");
    T->Branch("duplicateMuonFilter",&duplicateMuonFilter_,"duplicateMuonFilter/O");
    
    T->Branch("m_vis",&m_vis,"m_vis/F");
    T->Branch("m_vis_Up",&m_vis_Up,"m_vis_Up/F");
    T->Branch("m_vis_Down",&m_vis_Down,"m_vis_Down/F");
    T->Branch("m_vis_scaleUp",&m_vis_scaleUp,"m_vis_scaleUp/F");
    T->Branch("m_vis_scaleDown",&m_vis_scaleDown,"m_vis_scaleDown/F");
    T->Branch("m_vis_resoUp",&m_vis_resoUp,"m_vis_resoUp/F");
    T->Branch("m_vis_resoDown",&m_vis_resoDown,"m_vis_resoDown/F");
    
    T->Branch("m_sv", &m_sv, "m_sv/F");
    T->Branch("pt_sv", &pt_sv, "pt_sv/F");
    T->Branch("eta_sv", &eta_sv, "eta_sv/F");
    T->Branch("phi_sv", &phi_sv, "phi_sv/F");
    
    T->Branch("m_sv_muUp", &m_sv_muUp, "m_sv_muUp/F");
    T->Branch("m_sv_muDown", &m_sv_muDown, "m_sv_muDown/F");
    T->Branch("m_sv_scaleUp", &m_sv_scaleUp, "m_sv_scaleUp/F");
    T->Branch("m_sv_scaleDown", &m_sv_scaleDown, "m_sv_scaleDown/F");
    T->Branch("m_sv_resoUp", &m_sv_resoUp, "m_sv_resoUp/F");
    T->Branch("m_sv_resoDown", &m_sv_resoDown, "m_sv_resoDown/F");

    T->Branch("pt_sv_muUp", &pt_sv_muUp, "pt_sv_muUp/F");
    T->Branch("eta_sv_muUp", &eta_sv_muUp, "eta_sv_muUp/F");
    T->Branch("phi_sv_muUp", &phi_sv_muUp, "phi_sv_muUp/F");
	
		T->Branch("pt_sv_scaleUp", &pt_sv_scaleUp, "pt_sv_scaleUp/F");
		T->Branch("eta_sv_scaleUp", &eta_sv_scaleUp, "eta_sv_scaleUp/F");
		T->Branch("phi_sv_scaleUp", &phi_sv_scaleUp, "phi_sv_scaleUp/F");
	
		T->Branch("pt_sv_resoUp", &pt_sv_resoUp, "pt_sv_resoUp/F");
		T->Branch("eta_sv_resoUp", &eta_sv_resoUp, "eta_sv_resoUp/F");
		T->Branch("phi_sv_resoUp", &phi_sv_resoUp, "phi_sv_resoUp/F");

    T->Branch("pt_sv_muDown", &pt_sv_muDown, "pt_sv_muDown/F");
    T->Branch("eta_sv_muDown", &eta_sv_muDown, "eta_sv_muDown/F");
    T->Branch("phi_sv_muDown", &phi_sv_muDown, "phi_sv_muDown/F");

		T->Branch("pt_sv_scaleDown", &pt_sv_scaleDown, "pt_sv_scaleDown/F");
		T->Branch("eta_sv_scaleDown", &eta_sv_scaleDown, "eta_sv_scaleDown/F");
		T->Branch("phi_sv_scaleDown", &phi_sv_scaleDown, "phi_sv_scaleDown/F");
		
		T->Branch("pt_sv_resoDown", &pt_sv_resoDown, "pt_sv_resoDown/F");
		T->Branch("eta_sv_resoDown", &eta_sv_resoDown, "eta_sv_resoDown/F");
		T->Branch("phi_sv_resoDown", &phi_sv_resoDown, "phi_sv_resoDown/F");

    T->Branch("pt_1", &pt_1, "pt_1/F");
    T->Branch("pt_1_Up", &pt_1_Up, "pt_1_Up/F");
    T->Branch("pt_1_Down", &pt_1_Down, "pt_1_Down/F");
    
    T->Branch("eta_1", &eta_1, "eta_1/F");
    T->Branch("phi_1", &phi_1, "phi_1/F");
    T->Branch("q_1", &q_1, "q_1/I");
    T->Branch("d0_1",&d0_1,"d0_1/F");
    T->Branch("dZ_1",&dZ_1,"dZ_1/F");
    T->Branch("dcaSigd0_1",&dcaSigd0_1,"dcaSigd0_1/F");
    T->Branch("dcaSigdZ_1",&dcaSigdZ_1,"dcaSigdZ_1/F");
    
    T->Branch("pt_2", &pt_2, "pt_2/F");
    T->Branch("pt_2_Up", &pt_2_Up, "pt_2_Up/F");
    T->Branch("pt_2_Down", &pt_2_Down, "pt_2_Down/F");
    
    T->Branch("eta_2", &eta_2, "eta_2/F");
    T->Branch("phi_2", &phi_2, "phi_2/F");
    T->Branch("q_2", &q_2, "q_2/I");
    T->Branch("d0_2",&d0_2,"d0_2/F");
    T->Branch("dZ_2",&dZ_2,"dZ_2/F");
    T->Branch("dcaSigd0_2",&dcaSigd0_2,"dcaSigd0_2/F");
    T->Branch("dcaSigdZ_2",&dcaSigdZ_2,"dcaSigdZ_2/F");
    
    T->Branch("os", &os, "os/O");
    
    T->Branch("met",&met,"met/F");
    T->Branch("metphi",&metphi,"metphi/F");
    T->Branch("metcov00", &metcov00, "metcov00/F");
    T->Branch("metcov01", &metcov01, "metcov01/F");
    T->Branch("metcov10", &metcov10, "metcov10/F");
    T->Branch("metcov11", &metcov11, "metcov11/F");
    
    T->Branch("met_uncorr", &met_uncorr, "met_uncorr/F");
    T->Branch("metphi_uncorr", &metphi_uncorr, "metphi_uncorr/F");

    T->Branch("met_Up", &met_Up, "met_Up/F");
    T->Branch("met_Down", &met_Down, "met_Down/F");
    
    T->Branch("met_JES_Up_uncorr", &met_JES_Up_uncorr, "met_JES_Up_uncorr/F");
    T->Branch("met_JES_Down_uncorr", &met_JES_Down_uncorr, "met_JES_Down_uncorr/F");
    
    T->Branch("met_JES_Up", &met_JES_Up, "met_JES_Up/F");
    T->Branch("met_JES_Down", &met_JES_Down, "met_JES_Down/F");
    
    T->Branch("met_UnclusteredJES_Up_uncorr", &met_UnclusteredJES_Up_uncorr, "met_UnclusteredJES_Up_uncorr/F");
    T->Branch("met_UnclusteredJES_Down_uncorr", &met_UnclusteredJES_Down_uncorr, "met_UnclusteredJES_Down_uncorr/F");
    
    T->Branch("met_UnclusteredJES_Up", &met_UnclusteredJES_Up, "met_UnclusteredJES_Up/F");
    T->Branch("met_UnclusteredJES_Down", &met_UnclusteredJES_Down, "met_UnclusteredJES_Down/F");
    
    T->Branch("genmet", &genmet, "genmet/F");
    T->Branch("genmetphi", &genmetphi, "genmetphi/F");

    T->Branch("msvmet", &msvmet, "msvmet/F");
    T->Branch("msvmetphi", &msvmetphi, "msvmetphi/F");

    T->Branch("met_x", &met_x, "met_x/F");
    T->Branch("met_y", &met_y, "met_y/F");
    T->Branch("px_mu1", &px_mu1, "px_mu1/F");
    T->Branch("py_mu1", &py_mu1, "py_mu1/F");
    T->Branch("px_mu2", &px_mu2, "px_mu2/F");
    T->Branch("py_mu2", &py_mu2, "py_mu2/F");
    T->Branch("pt_tot_x", &pt_tot_x, "pt_tot_x/F");
    T->Branch("pt_tot_y", &pt_tot_y, "pt_tot_y/F");
    T->Branch("pt_tot", &pt_tot, "pt_tot/F");
    T->Branch("pt_tot_Up", &pt_tot_Up, "pt_tot_Up/F");
    T->Branch("pt_tot_Down", &pt_tot_Down, "pt_tot_Down/F");
    T->Branch("pt_tot_JES_Up", &pt_tot_JES_Up, "pt_tot_JES_Up/F");
    T->Branch("pt_tot_JES_Down", &pt_tot_JES_Down, "pt_tot_JES_Down/F");
    T->Branch("pt_tot_UnclusteredJES_Up", &pt_tot_UnclusteredJES_Up, "pt_tot_UnclusteredJES_Up/F");
    T->Branch("pt_tot_UnclusteredJES_Down", &pt_tot_UnclusteredJES_Down, "pt_tot_UnclusteredJES_Down/F");

    T->Branch("pt_tt", &pt_tt, "pt_tt/F");
    T->Branch("dr_tt", &dr_tt, "dr_tt/F");
    T->Branch("dphi_tt", &dphi_tt, "dphi_tt/F");
    T->Branch("eta_tt", &eta_tt, "eta_tt/F");
    T->Branch("ptRatio",&ptRatio,"ptRatio/F");

    T->Branch("pt_tt_Up", &pt_tt_Up, "pt_tt_Up/F");
    T->Branch("dr_tt_Up", &dr_tt_Up, "dr_tt_Up/F");
    T->Branch("dphi_tt_Up", &dphi_tt_Up, "dphi_tt_Up/F");
    T->Branch("eta_tt_Up", &eta_tt_Up, "eta_tt_Up/F");
    T->Branch("ptRatio_Up",&ptRatio_Up,"ptRatio_Up/F");

    T->Branch("pt_tt_Down", &pt_tt_Down, "pt_tt_Down/F");
    T->Branch("dr_tt_Down", &dr_tt_Down, "dr_tt_Down/F");
    T->Branch("dphi_tt_Down", &dphi_tt_Down, "dphi_tt_Down/F");
    T->Branch("eta_tt_Down", &eta_tt_Down, "eta_tt_Down/F");
    T->Branch("ptRatio_Down",&ptRatio_Down,"ptRatio_Down/F");
    
    T->Branch("dcasig2Mu2D", &dcasig2Mu2D, "dcasig2Mu2D/F");
    T->Branch("dcasig2Mu3D", &dcasig2Mu3D, "dcasig2Mu3D/F");
    T->Branch("sig2Mu2D", &sig2Mu2D, "sig2Mu2D/F");
    T->Branch("sig2Mu3D", &sig2Mu3D, "sig2Mu3D/F");
    
    T->Branch("pzetavis", &pzetavis, "pzetavis/F");
    
    T->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
    T->Branch("dzeta", &dzeta, "dzeta/F");
    
    T->Branch("pzetamiss_uncorr", &pzetamiss_uncorr, "pzetamiss_uncorr/F");
    T->Branch("dzeta_uncorr", &dzeta_uncorr, "dzeta_uncorr/F");
    
    T->Branch("pzetamiss_Up", &pzetamiss_Up, "pzetamiss_Up/F");
    T->Branch("dzeta_Up", &dzeta_Up, "dzeta_Up/F");
    T->Branch("pzetamiss_Down", &pzetamiss_Down, "pzetamiss_Down/F");
    T->Branch("dzeta_Down", &dzeta_Down, "dzeta_Down/F");
    
    T->Branch("pzetamiss_JES_Up", &pzetamiss_JES_Up, "pzetamiss_JES_Up/F");
    T->Branch("dzeta_JES_Up", &dzeta_JES_Up, "dzeta_JES_Up/F");
    T->Branch("pzetamiss_JES_Down", &pzetamiss_JES_Down, "pzetamiss_JES_Down/F");
    T->Branch("dzeta_JES_Down", &dzeta_JES_Down, "dzeta_JES_Down/F");
    
    T->Branch("pzetamiss_UnclusteredJES_Up", &pzetamiss_UnclusteredJES_Up, "pzetamiss_UnclusteredJES_Up/F");
    T->Branch("dzeta_UnclusteredJES_Up", &dzeta_UnclusteredJES_Up, "dzeta_UnclusteredJES_Up/F");
    T->Branch("pzetamiss_UnclusteredJES_Down", &pzetamiss_UnclusteredJES_Down, "pzetamiss_UnclusteredJES_Down/F");
    T->Branch("dzeta_UnclusteredJES_Down", &dzeta_UnclusteredJES_Down, "dzeta_UnclusteredJES_Down/F");
    
    T->Branch("pzetamiss_genmet", &genmet, "genmet/F");
    T->Branch("dzeta_genmet", &genmet, "genmet/F");

    T->Branch("dphi_mumet_1", &dphi_mumet_1, "dphi_mumet_1/F");
    T->Branch("dphi_mumet_uncorr_1", &dphi_mumet_uncorr_1, "dphi_mumet_uncorr_1/F");
    T->Branch("dphi_mumet_JES_Up_1", &dphi_mumet_JES_Up_1, "dphi_mumet_JES_Up_1/F");
    T->Branch("dphi_mumet_JES_Down_1", &dphi_mumet_JES_Down_1, "dphi_mumet_JES_Down_1/F");
    T->Branch("dphi_mumet_UnclusteredJES_Up_1", &dphi_mumet_UnclusteredJES_Up_1, "dphi_mumet_UnclusteredJES_Up_1/F");
    T->Branch("dphi_mumet_UnclusteredJES_Down_1", &dphi_mumet_UnclusteredJES_Down_1, "dphi_mumet_UnclusteredJES_Down_1/F");
    
    T->Branch("dphi_mumet_2", &dphi_mumet_2, "dphi_mumet_2/F");
    T->Branch("dphi_mumet_uncorr_2", &dphi_mumet_uncorr_2, "dphi_mumet_uncorr_2/F");
    T->Branch("dphi_mumet_JES_Up_2", &dphi_mumet_JES_Up_2, "dphi_mumet_JES_Up_2/F");
    T->Branch("dphi_mumet_JES_Down_2", &dphi_mumet_JES_Down_2, "dphi_mumet_JES_Down_2/F");
    T->Branch("dphi_mumet_UnclusteredJES_Up_2", &dphi_mumet_UnclusteredJES_Up_2, "dphi_mumet_UnclusteredJES_Up_2/F");
    T->Branch("dphi_mumet_UnclusteredJES_Down_2", &dphi_mumet_UnclusteredJES_Down_2, "dphi_mumet_UnclusteredJES_Down_2/F");
    
    T->Branch("dphi_posmu_met", &dphi_posmu_met, "dphi_posmu_met/F");
    T->Branch("dphi_posmu_met_uncorr", &dphi_posmu_met_uncorr, "dphi_posmu_met_uncorr/F");
    T->Branch("dphi_posmu_met_Up", &dphi_posmu_met_Up, "dphi_posmu_met_Up/F");
    T->Branch("dphi_posmu_met_Down", &dphi_posmu_met_Down, "dphi_posmu_met_Down/F");
    T->Branch("dphi_posmu_met_JES_Up", &dphi_posmu_met_JES_Up, "dphi_posmu_met_JES_Up/F");
    T->Branch("dphi_posmu_met_JES_Down", &dphi_posmu_met_JES_Down, "dphi_posmu_met_JES_Down/F");
    T->Branch("dphi_posmu_met_UnclusteredJES_Up", &dphi_posmu_met_UnclusteredJES_Up, "dphi_posmu_met_UnclusteredJES_Up/F");
    T->Branch("dphi_posmu_met_UnclusteredJES_Down", &dphi_posmu_met_UnclusteredJES_Down, "dphi_posmu_met_UnclusteredJES_Down/F");
    
    T->Branch("dphi_twomu", &dphi_twomu, "dphi_twomu/F");
    
    T->Branch("costheta_1", &costheta_1, "costheta_1/F");
    T->Branch("costheta_2", &costheta_2, "costheta_2/F");
    T->Branch("costheta", &costheta, "costheta/F");
    
    T->Branch("pos_mupt1", &pos_mupt1, "pos_mupt1/F");
    T->Branch("neg_mupt1", &neg_mupt1, "neg_mupt1/F");
    T->Branch("pos_mueta1", &pos_mueta1, "pos_mueta1/F");
    T->Branch("neg_mueta1", &neg_mueta1, "neg_mueta1/F");
    T->Branch("pos_mupt2", &pos_mupt2, "pos_mupt2/F");
    T->Branch("neg_mupt2", &neg_mupt2, "neg_mupt2/F");
    T->Branch("pos_mueta2", &pos_mueta2, "pos_mueta2/F");
    T->Branch("neg_mueta2", &neg_mueta2, "neg_mueta2/F");
    
    T->Branch("genZ_mass", &n_genZ_mass, "n_genZ_mass/F");
    T->Branch("genZ_Pt", &n_genZ_Pt, "n_genZ_Pt/F");
    T->Branch("genZ_Eta", &n_genZ_Eta, "n_genZ_Eta/F");
    T->Branch("genZ_Phi", &n_genZ_Phi, "n_genZ_Phi/F");
    T->Branch("genZTT_mass", &n_genZTT_mass, "n_genZTT_mass/F");
    T->Branch("genZTT_Pt", &n_genZTT_Pt, "n_genZTT_Pt/F");
    T->Branch("genZTT_Eta", &n_genZTT_Eta, "n_genZTT_Eta/F");
    T->Branch("genZTT_Phi", &n_genZTT_Phi, "n_genZTT_Phi/F");
    T->Branch("genV_mass", &n_genV_mass, "n_genV_mass/F");
    T->Branch("genV_Pt", &n_genV_Pt, "n_genV_Pt/F");
    T->Branch("genV_Eta", &n_genV_Eta, "n_genV_Eta/F");
    T->Branch("genV_Phi", &n_genV_Phi, "n_genV_Phi/F");
    T->Branch("gen_mu1_Pt", &n_gen_mu1_Pt, "n_gen_mu1_Pt/F");
    T->Branch("gen_mu1_Eta", &n_gen_mu1_Eta, "n_gen_mu1_Eta/F");
    T->Branch("gen_mu1_Phi", &n_gen_mu1_Phi, "n_gen_mu1_Phi/F");
    T->Branch("gen_mu2_Pt", &n_gen_mu2_Pt, "n_gen_mu2_Pt/F");
    T->Branch("gen_mu2_Eta", &n_gen_mu2_Eta, "n_gen_mu2_Eta/F");
    T->Branch("gen_mu2_Phi", &n_gen_mu2_Phi, "n_gen_mu2_Phi/F");
    T->Branch("gen_taus",&n_gen_taus,"n_gen_taus/i");
    T->Branch("gen_mutaus",&n_gen_mutaus,"n_gen_mutaus/i");
    
    T->Branch("njets", &njets, "njets/I");
    T->Branch("njets_Up", &njets_Up, "njets_Up/I");
    T->Branch("njets_Down", &njets_Down, "njets_Down/I");
    T->Branch("njetspt20", &njetspt20, "njetspt20/I");
    
    T->Branch("jpt_1", &jpt_1, "jpt_1/F");
    T->Branch("jpt_1_Up", &jpt_1_Up, "jpt_1_Up/F");
    T->Branch("jpt_1_Down",&jpt_1_Down, "jpt_1_Down/F");
    
    T->Branch("jeta_1", &jeta_1, "jeta_1/F");
    T->Branch("jphi_1", &jphi_1, "jphi_1/F");
    T->Branch("jptraw_1", &jptraw_1, "jptraw_1/F");
    T->Branch("jptunc_1", &jptunc_1, "jptunc_1/F");
    T->Branch("jmva_1", &jmva_1, "jmva_1/F");
    T->Branch("jlrm_1", &jlrm_1, "jlrm_1/F");
    T->Branch("jctm_1", &jctm_1, "jctm_1/I");
    
    T->Branch("jpt_2", &jpt_2, "jpt_2/F");
    T->Branch("jpt_2_Up", &jpt_2_Up, "jpt_2_Up/F");
    T->Branch("jpt_2_Down", &jpt_2_Down, "jpt_2_Down/F");
    
    T->Branch("jeta_2", &jeta_2, "jeta_2/F");
    T->Branch("jphi_2", &jphi_2, "jphi_2/F");
    T->Branch("jptraw_2", &jptraw_2, "jptraw_2/F");
    T->Branch("jptunc_2", &jptunc_2, "jptunc_2/F");
    T->Branch("jmva_2", &jmva_2, "jlrm_2/F");
    T->Branch("jlrm_2", &jlrm_2, "jlrm_2/F");
    T->Branch("jctm_2", &jctm_2, "jctm_2/I");
    
    T->Branch("mjj", &mjj, "mjj/F");
    T->Branch("mjj_Up", &mjj_Up, "mjj_Up/F");
    T->Branch("mjj_Down", &mjj_Down, "mjj_Down/F");
    
    T->Branch("jdeta", &jdeta, "jdeta/F");
    T->Branch("njetingap", &njetingap, "njetingap/I");
    
    T->Branch("nbtag", &nbtag, "nbtag/I");
    T->Branch("nbtag_Up", &nbtag_Up, "nbtag_Up/I");
    T->Branch("nbtag_Down", &nbtag_Down, "nbtag_Down/I");
    T->Branch("nbtag_noSF", &nbtag_noSF, "nbtag_noSF/I");
    T->Branch("bpt",   &bpt,   "bpt/F");
    T->Branch("beta",  &beta,  "beta/F");
    T->Branch("bphi",  &bphi,  "bphi/F");
    
    T->Branch("npartons",&npartons,"npartons/i");

    T->Branch("bdt", &bdt, "bdt/F");
    T->Branch("bdt_Up", &bdt_Up, "bdt_Up/F");
    T->Branch("bdt_Down", &bdt_Down, "bdt_Down/F");
    T->Branch("bdt_JES_Up", &bdt_JES_Up, "bdt_JES_Up/F");
    T->Branch("bdt_JES_Down", &bdt_JES_Down, "bdt_JES_Down/F");
    T->Branch("bdt_UnclusteredJES_Up", &bdt_UnclusteredJES_Up, "bdt_UnclusteredJES_Up/F");
    T->Branch("bdt_UnclusteredJES_Down", &bdt_UnclusteredJES_Down, "bdt_UnclusteredJES_Down/F");
    
    T->Branch("bdt_0jets", &bdt_0jets, "bdt_0jets/F");
    T->Branch("bdt_0jets_Up", &bdt_0jets_Up, "bdt_0jets_Up/F");
    T->Branch("bdt_0jets_Down", &bdt_0jets_Down, "bdt_0jets_Down/F");
    T->Branch("bdt_0jets_JES_Up", &bdt_0jets_JES_Up, "bdt_0jets_JES_Up/F");
    T->Branch("bdt_0jets_JES_Down", &bdt_0jets_JES_Down, "bdt_0jets_JES_Down/F");
    T->Branch("bdt_0jets_UnclusteredJES_Up", &bdt_0jets_UnclusteredJES_Up, "bdt_0jets_UnclusteredJES_Up/F");
    T->Branch("bdt_0jets_UnclusteredJES_Down", &bdt_0jets_UnclusteredJES_Down, "bdt_0jets_UnclusteredJES_Down/F");

    T->Branch("bdt_boosted", &bdt_boosted, "bdt_boosted/F");
    T->Branch("bdt_boosted_Up", &bdt_boosted_Up, "bdt_boosted_Up/F");
    T->Branch("bdt_boosted_Down", &bdt_boosted_Down, "bdt_boosted_Down/F");
    T->Branch("bdt_boosted_JES_Up", &bdt_boosted_JES_Up, "bdt_boosted_JES_Up/F");
    T->Branch("bdt_boosted_JES_Down", &bdt_boosted_JES_Down, "bdt_boosted_JES_Down/F");
    T->Branch("bdt_boosted_UnclusteredJES_Up", &bdt_boosted_UnclusteredJES_Up, "bdt_boosted_UnclusteredJES_Up/F");
    T->Branch("bdt_boosted_UnclusteredJES_Down", &bdt_boosted_UnclusteredJES_Down, "bdt_boosted_UnclusteredJES_Down/F");

    T->Branch("bdt_vbf", &bdt_vbf, "bdt_vbf/F");
    T->Branch("bdt_vbf_Up", &bdt_vbf_Up, "bdt_vbf_Up/F");
    T->Branch("bdt_vbf_Down", &bdt_vbf_Down, "bdt_vbf_Down/F");
    T->Branch("bdt_vbf_JES_Up", &bdt_vbf_JES_Up, "bdt_vbf_JES_Up/F");
    T->Branch("bdt_vbf_JES_Down", &bdt_vbf_JES_Down, "bdt_vbf_JES_Down/F");
    T->Branch("bdt_vbf_UnclusteredJES_Up", &bdt_vbf_UnclusteredJES_Up, "bdt_vbf_UnclusteredJES_Up/F");
    T->Branch("bdt_vbf_UnclusteredJES_Down", &bdt_vbf_UnclusteredJES_Down, "bdt_vbf_UnclusteredJES_Down/F");

    //---------------------------defining JES uncertainties------------------------->

    JESUncertainties * jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Summer16_UncertaintySources_AK4PFchs.txt");
    std::vector<std::string> uncertNames = jecUncertainties->getUncertNames();

    std::cout << "Number of uncertainties = " << uncertNames.size() << std::endl;

    int njetsUncUp[30];
    int njetsUncDown[30];
    float mjjUncUp[30];
    float mjjUncDown[30];

    int iUncert = 0;
    if (!isData) {
      for (auto const& Name : uncertNames) {
	TString name(Name);
	cout << name << endl;
	T->Branch("njets_"+name+"Up",&njetsUncUp[iUncert],"njets_"+name+"Up/I");
	T->Branch("njets_"+name+"Down",&njetsUncDown[iUncert],"njets_"+name+"Down/I");
	T->Branch("mjj_"+name+"Up",&mjjUncUp[iUncert],"njets_"+name+"Up/F");
	T->Branch("mjj_"+name+"Down",&mjjUncDown[iUncert],"njets_"+name+"Down/F");
	iUncert++;
      }
      cout << endl;
    }
    
    //-----------  Pileup reweighting official receipe------------------------>
    
    //Initialize Pile up object
    PileUp * PUofficial = new PileUp();
    
    if (applyPUreweighting_official) {
        TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpDataFile,"read");
        //std::cout<< "file PileUp data = " << filePUdistribution_data << std::endl;
        TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+PileUpMCFile, "read");
        //std::cout<<"file PileUp mc = " << filePUdistribution_MC <<std::endl;
        TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
        TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
        PUofficial->set_h_data(PU_data);
        PUofficial->set_h_MC(PU_mc);
    }

    //---------------- Lepton Scale Factors--------------------------------->
    
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso = new ScaleFactor();
    //std::cout<<"test1"<<std::endl;
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_muonTrig;
    SF_muonTrig = new ScaleFactor();
    //std::cout<<"test2"<<std::endl;
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTrigFile));
    //std::cout<<"test3"<<std::endl;
    
    //---------------- Loding recoil correction-------------------->
    RecoilCorrector recoilPFMetCorrector("HTT-utilities/RecoilCorrections/data/"+RecoilFileName);
    //RecoilCorrector recoilMvaMetCorrector("HTT-utilities/RecoilCorrections/data/"+RecoilMvaFileName);
    
    //--------------- ZPtMass Coorection--------------------------->
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
    
    //----------Event category weights for BDTs----------->
    
    TFile * fileCategoryWeights = new TFile(TString(cmsswBase)+"/src/"+CategoryWeightsFileName);
    
    std::cout<< "Filename is " << TString(cmsswBase) <<"/src/" << CategoryWeightsFileName << std::endl;
    if (fileCategoryWeights==NULL){
        std::cout << "File " << TString(cmsswBase) << "/src/" << CategoryWeightsFileName << "  does not exist!" << std::endl;
    }
    TH3D * jet0Hist = (TH3D*)fileCategoryWeights->Get(Jet0WeightsHist);
    TH3D * boostedHist = (TH3D*)fileCategoryWeights->Get(BoostedWeightsHist);
    TH3D * vbfHist = (TH3D*)fileCategoryWeights->Get(VBFWeightsHist);
    if (Jet0WeightsHist == NULL|| BoostedWeightsHist ==NULL || VBFWeightsHist == NULL){
        std::cout << "histogram "<< Jet0WeightsHist<<" and " << BoostedWeightsHist <<" and  " <<VBFWeightsHist  << " are not found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << CategoryWeightsFileName << std::endl;
        exit(-1);
    }
    
    //-----------SVFit----------------->
    edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
    TH1::AddDirectory(false);
    TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
    
    //-------------BTag scale factors--------------------->
    BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/CSVv2_ichep.csv");
    BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    reader_B.load(calib,BTagEntry::FLAV_B,"comb");
    reader_C.load(calib,BTagEntry::FLAV_C,"comb");
    reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");
    
    float etaBTAG[2] = {0.5,2.1};
    float ptBTAG[5] = {25.,35.,50.,100.,200.};
    
    std::cout << std::endl;
    for (int iEta=0; iEta<2; ++iEta) {
        for (int iPt=0; iPt<5; ++iPt) {
            float sfB = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
            float sfC = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
            float sfLight = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
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
    float MinBJetPt = 20.; // !!!!!
    
//    //-------------------BDTs ----------------------------------------->

//    TString jet0BDTweight = "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_0jets_May18_BDT.weights.xml";
//    TString boostBDTweight= "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_boosted_May25_BDT.weights.xml";
//    TString vbfBDTweight  = "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_vbf_May18_BDT.weights.xml";
    
//BDT weight for nomvis trainned samples
    TString jet0BDTweight = "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_0jets_nomvis_july7_BDT.weights.xml";
    TString boostBDTweight= "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_nomvis_boosted_july7_BDT.weights.xml";
    TString vbfBDTweight  = "/src/DesyTauAnalyses/NTupleMaker/data/mumu_BDTWeights/TMVA_RunBtoH_nomvis_vbf_July7_BDT.weights.xml";
      
    //This loads the library
    TMVA::Tools::Instance();
    
    //0jets
    TMVA::Reader *reader0jets = new TMVA::Reader("!V:!Color");
    reader0jets->AddVariable( "eta_tt",&eta_tt);
    reader0jets->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader0jets->AddVariable("dzeta",&dzeta);
    reader0jets->AddVariable( "met",&met);
    //    reader0jets->AddVariable("m_vis", &m_vis);
    reader0jets->AddVariable("costheta", &costheta);
    reader0jets->AddVariable("ptRatio", &ptRatio);
    reader0jets->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //boosted
    TMVA::Reader *readerboost = new TMVA::Reader("!V:!Color");
    readerboost->AddVariable( "eta_tt",&eta_tt);
    readerboost->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    readerboost->AddVariable("dzeta",&dzeta);
    readerboost->AddVariable( "met",&met);
    //    readerboost->AddVariable("m_vis", &m_vis);
    readerboost->AddVariable("costheta", &costheta);
    readerboost->AddVariable("ptRatio", &ptRatio);
    readerboost->AddVariable("pt_tot", &pt_tot);
    readerboost->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //vbf
    TMVA::Reader *readervbf = new TMVA::Reader("!V:!Color");
    readervbf->AddVariable("eta_tt",&eta_tt);
    readervbf->AddVariable("dphi_posmu_met",&dphi_posmu_met);
    readervbf->AddVariable("dzeta",&dzeta);
    readervbf->AddVariable("met",&met);
    //    readervbf->AddVariable("m_vis", &m_vis);
    readervbf->AddVariable("costheta", &costheta);
    readervbf->AddVariable("ptRatio", &ptRatio);
    readervbf->AddVariable("mjj", &mjj);
    //readervbf->AddVariable("jdeta", &jdeta);
    
    readervbf->BookMVA("BDT", cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for JESUp******************************************************************
    TMVA::Reader *reader_0jetsJESUp = new TMVA::Reader("!V:!Color");
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_0jetsJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Up);
    reader_0jetsJESUp->AddVariable("dzeta",&dzeta_JES_Up);
    reader_0jetsJESUp->AddVariable( "met",&met_JES_Up);
    //    reader_0jetsJESUp->AddVariable("m_vis", &m_vis);
    reader_0jetsJESUp->AddVariable("costheta", &costheta);
    reader_0jetsJESUp->AddVariable("ptRatio", &ptRatio);
    //BookMethod
    reader_0jetsJESUp->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    TMVA::Reader *reader_boostJESUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_boostJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Up);
    reader_boostJESUp->AddVariable("dzeta",&dzeta_JES_Up);
    reader_boostJESUp->AddVariable( "met",&met_JES_Up);
    //    reader_boostJESUp->AddVariable("m_vis", &m_vis);
    reader_boostJESUp->AddVariable("costheta", &costheta);
    reader_boostJESUp->AddVariable("ptRatio", &ptRatio);
    reader_boostJESUp->AddVariable("pt_tot", &pt_tot_JES_Up);
    //BookMethod
    reader_boostJESUp->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    TMVA::Reader *reader_vbfJESUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_vbfJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Up);
    reader_vbfJESUp->AddVariable("dzeta",&dzeta_JES_Up);
    reader_vbfJESUp->AddVariable( "met",&met_JES_Up);
    //    reader_vbfJESUp->AddVariable("m_vis", &m_vis);
    reader_vbfJESUp->AddVariable("costheta", &costheta);
    reader_vbfJESUp->AddVariable("ptRatio", &ptRatio);
    reader_vbfJESUp->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfJESUp->BookMVA("BDT", cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for JESDown****************************************************************
    TMVA::Reader *reader_0jetsJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_0jetsJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Down);
    reader_0jetsJESDown->AddVariable("dzeta",&dzeta_JES_Down);
    reader_0jetsJESDown->AddVariable( "met",&met_JES_Down);
    //    reader_0jetsJESDown->AddVariable("m_vis", &m_vis);
    reader_0jetsJESDown->AddVariable("costheta", &costheta);
    reader_0jetsJESDown->AddVariable("ptRatio", &ptRatio);
    //BookMethod
    reader_0jetsJESDown->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_boostJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_boostJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Down);
    reader_boostJESDown->AddVariable("dzeta",&dzeta_JES_Down);
    reader_boostJESDown->AddVariable( "met",&met_JES_Down);
    //    reader_boostJESDown->AddVariable("m_vis", &m_vis);
    reader_boostJESDown->AddVariable("costheta", &costheta);
    reader_boostJESDown->AddVariable("ptRatio", &ptRatio);
    reader_boostJESDown->AddVariable("pt_tot", &pt_tot_JES_Down);
    //BookMethod
    reader_boostJESDown->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_vbfJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_vbfJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_JES_Down);
    reader_vbfJESDown->AddVariable("dzeta",&dzeta_JES_Down);
    reader_vbfJESDown->AddVariable( "met",&met_JES_Down);
    //    reader_vbfJESDown->AddVariable("m_vis", &m_vis);
    reader_vbfJESDown->AddVariable("costheta", &costheta);
    reader_vbfJESDown->AddVariable("ptRatio", &ptRatio);
    reader_vbfJESDown->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfJESDown->BookMVA("BDT", cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for UnclustJESUp******************************************************************
    TMVA::Reader *reader_0jetsUnclustJESUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsUnclustJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_0jetsUnclustJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Up);
    reader_0jetsUnclustJESUp->AddVariable("dzeta",&dzeta_UnclusteredJES_Up);
    reader_0jetsUnclustJESUp->AddVariable( "met",&met_UnclusteredJES_Up);
    //   reader_0jetsUnclustJESUp->AddVariable("m_vis", &m_vis);
    reader_0jetsUnclustJESUp->AddVariable("costheta", &costheta);
    reader_0jetsUnclustJESUp->AddVariable("ptRatio", &ptRatio);
    //BookMethod
    reader_0jetsUnclustJESUp->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_boostUnclustJESUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostUnclustJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_boostUnclustJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Up);
    reader_boostUnclustJESUp->AddVariable("dzeta",&dzeta_UnclusteredJES_Up);
    reader_boostUnclustJESUp->AddVariable( "met",&met_UnclusteredJES_Up);
    //    reader_boostUnclustJESUp->AddVariable("m_vis", &m_vis);
    reader_boostUnclustJESUp->AddVariable("costheta", &costheta);
    reader_boostUnclustJESUp->AddVariable("ptRatio", &ptRatio);
    reader_boostUnclustJESUp->AddVariable("pt_tot", &pt_tot_UnclusteredJES_Up);
    //BookMethod
    reader_boostUnclustJESUp->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_vbfUnclustJESUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfUnclustJESUp->AddVariable( "eta_tt",&eta_tt);
    reader_vbfUnclustJESUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Up);
    reader_vbfUnclustJESUp->AddVariable("dzeta",&dzeta_UnclusteredJES_Up);
    reader_vbfUnclustJESUp->AddVariable( "met",&met_UnclusteredJES_Up);
    //    reader_vbfUnclustJESUp->AddVariable("m_vis", &m_vis);
    reader_vbfUnclustJESUp->AddVariable("costheta", &costheta);
    reader_vbfUnclustJESUp->AddVariable("ptRatio", &ptRatio);
    reader_vbfUnclustJESUp->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfUnclustJESUp->BookMVA("BDT",cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for UnclustJESDown******************************************************************
    TMVA::Reader *reader_0jetsUnclustJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsUnclustJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_0jetsUnclustJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Down);
    reader_0jetsUnclustJESDown->AddVariable("dzeta",&dzeta_UnclusteredJES_Down);
    reader_0jetsUnclustJESDown->AddVariable( "met",&met_UnclusteredJES_Down);
    //    reader_0jetsUnclustJESDown->AddVariable("m_vis", &m_vis);
    reader_0jetsUnclustJESDown->AddVariable("costheta", &costheta);
    reader_0jetsUnclustJESDown->AddVariable("ptRatio", &ptRatio);
    //BookMethod
    reader_0jetsUnclustJESDown->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_boostUnclustJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostUnclustJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_boostUnclustJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Down);
    reader_boostUnclustJESDown->AddVariable("dzeta",&dzeta_UnclusteredJES_Down);
    reader_boostUnclustJESDown->AddVariable( "met",&met_UnclusteredJES_Down);
    //    reader_boostUnclustJESDown->AddVariable("m_vis", &m_vis);
    reader_boostUnclustJESDown->AddVariable("costheta", &costheta);
    reader_boostUnclustJESDown->AddVariable("ptRatio", &ptRatio);
    reader_boostUnclustJESDown->AddVariable("pt_tot", &pt_tot_UnclusteredJES_Down);
    //BookMethod
    reader_boostUnclustJESDown->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_vbfUnclustJESDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfUnclustJESDown->AddVariable( "eta_tt",&eta_tt);
    reader_vbfUnclustJESDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met_UnclusteredJES_Down);
    reader_vbfUnclustJESDown->AddVariable("dzeta",&dzeta_UnclusteredJES_Down);
    reader_vbfUnclustJESDown->AddVariable( "met",&met_UnclusteredJES_Down);
    //    reader_vbfUnclustJESDown->AddVariable("m_vis", &m_vis);
    reader_vbfUnclustJESDown->AddVariable("costheta", &costheta);
    reader_vbfUnclustJESDown->AddVariable("ptRatio", &ptRatio);
    reader_vbfUnclustJESDown->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfUnclustJESDown->BookMVA("BDT", cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for MuScaleUp******************************************************************
    TMVA::Reader *reader_0jetsMuScaleUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsMuScaleUp->AddVariable( "eta_tt",&eta_tt_Up);
    reader_0jetsMuScaleUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_0jetsMuScaleUp->AddVariable("dzeta",&dzeta);
    reader_0jetsMuScaleUp->AddVariable( "met",&met);
    //   reader_0jetsMuScaleUp->AddVariable("m_vis", &m_vis_Up);
    reader_0jetsMuScaleUp->AddVariable("costheta", &costheta);
    reader_0jetsMuScaleUp->AddVariable("ptRatio", &ptRatio_Up);
    //BookMethod
    reader_0jetsMuScaleUp->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_boostMuScaleUp = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostMuScaleUp->AddVariable( "eta_tt",&eta_tt_Up);
    reader_boostMuScaleUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_boostMuScaleUp->AddVariable("dzeta",&dzeta);
    reader_boostMuScaleUp->AddVariable( "met",&met);
    //    reader_boostMuScaleUp->AddVariable("m_vis", &m_vis_Up);
    reader_boostMuScaleUp->AddVariable("costheta", &costheta);
    reader_boostMuScaleUp->AddVariable("ptRatio", &ptRatio_Up);
    reader_boostMuScaleUp->AddVariable("pt_tot", &pt_tot_Up);
    //BookMethod
    reader_boostMuScaleUp->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_vbfMuScaleUp = new TMVA::Reader("!V:!Color");
	
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfMuScaleUp->AddVariable( "eta_tt",&eta_tt_Up);
    reader_vbfMuScaleUp->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_vbfMuScaleUp->AddVariable("dzeta",&dzeta);
    reader_vbfMuScaleUp->AddVariable( "met",&met);
    //    reader_vbfMuScaleUp->AddVariable("m_vis", &m_vis_Up);
    reader_vbfMuScaleUp->AddVariable("costheta", &costheta);
    reader_vbfMuScaleUp->AddVariable("ptRatio", &ptRatio_Up);
    reader_vbfMuScaleUp->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfMuScaleUp->BookMVA("BDT", cmsswBase+vbfBDTweight);
    
    //Create TMVA reader for MuScaleDown******************************************************************
    TMVA::Reader *reader_0jetsMuScaleDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_0jetsMuScaleDown->AddVariable( "eta_tt",&eta_tt_Down);
    reader_0jetsMuScaleDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_0jetsMuScaleDown->AddVariable("dzeta",&dzeta);
    reader_0jetsMuScaleDown->AddVariable( "met",&met);
    //    reader_0jetsMuScaleDown->AddVariable("m_vis", &m_vis_Down);
    reader_0jetsMuScaleDown->AddVariable("costheta", &costheta);
    reader_0jetsMuScaleDown->AddVariable("ptRatio", &ptRatio_Down);
    //BookMethod
    reader_0jetsMuScaleDown->BookMVA("BDT", cmsswBase+jet0BDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_boostMuScaleDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_boostMuScaleDown->AddVariable( "eta_tt",&eta_tt_Down);
    reader_boostMuScaleDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_boostMuScaleDown->AddVariable("dzeta",&dzeta);
    reader_boostMuScaleDown->AddVariable( "met",&met);
    //    reader_boostMuScaleDown->AddVariable("m_vis", &m_vis_Down);
    reader_boostMuScaleDown->AddVariable("costheta", &costheta);
    reader_boostMuScaleDown->AddVariable("ptRatio", &ptRatio_Down);
    reader_boostMuScaleDown->AddVariable("pt_tot", &pt_tot_Down);
    //BookMethod
    reader_boostMuScaleDown->BookMVA("BDT", cmsswBase+boostBDTweight);
    
    //Create TMVA Reader Object
    TMVA::Reader *reader_vbfMuScaleDown = new TMVA::Reader("!V:!Color");
    
    //create set of variables as declared in weight file and declared them to reader
    reader_vbfMuScaleDown->AddVariable( "eta_tt",&eta_tt_Down);
    reader_vbfMuScaleDown->AddVariable( "dphi_posmu_met",&dphi_posmu_met);
    reader_vbfMuScaleDown->AddVariable("dzeta",&dzeta);
    reader_vbfMuScaleDown->AddVariable( "met",&met);
    //    reader_vbfMuScaleDown->AddVariable("m_vis", &m_vis_Down);
    reader_vbfMuScaleDown->AddVariable("costheta", &costheta);
    reader_vbfMuScaleDown->AddVariable("ptRatio", &ptRatio_Down);
    reader_vbfMuScaleDown->AddVariable("mjj", &mjj);
    //BookMethod
    reader_vbfMuScaleDown->BookMVA("BDT", cmsswBase+vbfBDTweight);

    //Met filters
    std::vector<TString> metFlags; metFlags.clear();
    metFlags.push_back("Flag_HBHENoiseFilter");
    metFlags.push_back("Flag_HBHENoiseIsoFilter");
    metFlags.push_back("Flag_globalTightHalo2016Filter");
    metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
    metFlags.push_back("Flag_goodVertices");
    if (isData)
      metFlags.push_back("Flag_eeBadScFilter");
    metFlags.push_back("Flag_BadPFMuonFilter");
    metFlags.push_back("Flag_BadChargedCandidateFilter");

    std::vector<TString> badChargedCandidateFlag; badChargedCandidateFlag.clear();
    badChargedCandidateFlag.push_back("Flag_BadChargedCandidateFilter");
    std::vector<TString> badPFMuonFlag; badPFMuonFlag.clear();
    badPFMuonFlag.push_back("Flag_BadPFMuonFilter");
    std::vector<TString> badGlobalMuonFlag; badGlobalMuonFlag.clear();
    badGlobalMuonFlag.push_back("Flag_BadGlobalMuonFilter");
    std::vector<TString> muonBadTrackFlag; muonBadTrackFlag.clear();
    muonBadTrackFlag.push_back("Flag_muonBadTrackFilter");
    std::vector<TString> chargedHadronTrackResolutionFlag; chargedHadronTrackResolutionFlag.clear();
    chargedHadronTrackResolutionFlag.push_back("Flag_chargedHadronTrackResolutionFilter");
    

    //------------------------------->
    
    int nFiles = 0;
    int nEvents = 0;
    int selEventsAllMuons = 0;
    int selEventsIdMuons = 0;
    int selEventsIsoMuons = 0;
    int nTotalFiles = 0;
    
    std::string dummy;
    // count number of files --->
    while (fileList0 >> dummy) nTotalFiles++;
    
    unsigned int RunMin = 9999999;
    unsigned int RunMax = 0;
    
    std::vector<unsigned int> allRuns; allRuns.clear();
    
    std::vector<Period> periods;
    //string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;
    
    if (isData) {
        std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
        if (inputFileStream.fail()) {
            std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
            std::cout << "please check" << std::endl;
            std::cout << "quitting program" << std::endl;
            exit(-1);
        }
        for(std::string s; std::getline(inputFileStream, s); )
        {
            periods.push_back(Period());
            std::stringstream ss(s);
            ss >> periods.back();
        }
    }
    
    //----Attention----//
    //--------file loop begins-------------->
    
    for (int iF=0; iF<nTotalFiles; ++iF) {
        std::string filen;
        fileList >> filen;
        
        std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
        TFile * file_ = TFile::Open(TString(filen));
        
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
        
        AC1B analysisTree(_tree);
        
        Long64_t numberOfEntries = analysisTree.GetEntries();
        
        std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
        
        //------------Looping over number of entries----------->
        
        for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) {
            
            //std::cout<<"enter into the loop"<<std::endl;
            analysisTree.GetEntry(iEntry);
            //std::cout << "number of events   = "<< nEvents<<std::endl;
            nEvents++;
            
            
            if (nEvents%10000==0)
                cout << "      processed " << nEvents << " events" << endl;
            
            weight = 1;
            genweight = 1;
            
            if (!isData) {
                genweight = 1;
                
                if (analysisTree.genweight<0) genweight = -1;
                
                //std::cout << "genweight = " << analysisTree.genweight << std::endl;
                weight *= genweight;
                histWeightsH->Fill(0.,weight);

            }
            else
                histWeightsH->Fill(0.,1.);
	    mcweight = 1;
	    if (!isData)
	      mcweight = analysisTree.genweight;
            puweight = 1;
            trigweight = 1;
            weightTrig_gen = 1;
            effweight_gen = 1;
            idweight_1 = 1;
            idweight_2 = 1;
            isoweight_1 = 1;
            isoweight_2 = 1;
            effweight = 1;
            topptweight = 1;
            zptmassweight = 1;
            trkeffweight = 1;
            jet0weight = 1;
            boostedweight = 1;
            vbfweight = 1;
            btag0weight =1;
            btag0weight_Up =1;
            btag0weight_Down =1;

	    metFilters_ = true;

	    badChargedCandidateFilter_ = true;
            badPFMuonFilter_ = true;
	    badGlobalMuonFilter_ = true;
	    muonBadTrackFilter_ = true;
	    chargedHadronTrackResolutionFilter_ = true;

	    badMuonFilter_ = true;
	    duplicateMuonFilter_ = true;
            
            TLorentzVector genZ; genZ.SetXYZT(0,0,0,0);
            TLorentzVector genZTT; genZTT.SetXYZT(0,0,0,0);
            TLorentzVector genV; genV.SetXYZT(0,0,0,0);
            TLorentzVector genL; genL.SetXYZT(0,0,0,0);
            std::vector<unsigned int> muTaus; muTaus.clear();
            std::vector<TLorentzVector> muTausLV; muTausLV.clear();
            //std::vector<TLorentzVector> promptTausLV; promptTausLV.clear();
            std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
            TLorentzVector gen_mu1; gen_mu1.SetXYZM(0.01,0,0,muonMass);
            TLorentzVector gen_mu2; gen_mu2.SetXYZM(0.01,0,0,muonMass);
            TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
            
            if (!isData){
                for (unsigned int igentau=0; igentau<analysisTree.gentau_count; ++igentau) {
                    TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
                                                        analysisTree.gentau_py[igentau],
                                                        analysisTree.gentau_pz[igentau],
                                                        analysisTree.gentau_e[igentau]);
                    if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
                        //promptTausLV.push_back(tauLV);
                        promptTausFirstCopy.push_back(tauLV);
                        promptTausLV += tauLV;
                        genV += tauLV; //added on Jan 25 this cause a problem in recoil correction as well as in zpt reweighting
                    }
                }
                
                for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
                    //	  cout << igen << "   pdgId = " << analysisTree.genparticles_pdgid[igen] << endl;
                    TLorentzVector genPart; genPart.SetXYZT(analysisTree.genparticles_px[igen],
                                                            analysisTree.genparticles_py[igen],
                                                            analysisTree.genparticles_pz[igen],
                                                            analysisTree.genparticles_e[igen]);
                    if (analysisTree.genparticles_pdgid[igen]==23) {
                        genZ.SetXYZT(analysisTree.genparticles_px[igen],
                                     analysisTree.genparticles_py[igen],
                                     analysisTree.genparticles_pz[igen],
                                     analysisTree.genparticles_e[igen]);
                    }
                    bool isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
                    bool isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
                    bool isLepton = isMuon || isElectron;
                    bool isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
                    fabs(analysisTree.genparticles_pdgid[igen])==14||
                    fabs(analysisTree.genparticles_pdgid[igen])==16;
                    bool isPrompt = analysisTree.genparticles_isPrompt[igen]||
                    analysisTree.genparticles_isPromptTauDecayProduct[igen];
                    
                    if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
                        if (isLepton) {
                            genV += genPart;
                            genL += genPart;
                        }
                        if (isMuon&&analysisTree.genparticles_isPromptTauDecayProduct[igen]) {
                            muTaus.push_back(igen);
                            muTausLV.push_back(genPart);
                        }
                        
                        if (isNeutrino)
                            genV += genPart;
                    }
                }
                
                //if (promptTausLV.size()==2)
                //  genZTT = promptTausLV[0] + promptTausLV[1];
                if (promptTausFirstCopy.size()==2)
                    genZTT = promptTausFirstCopy[0] + promptTausFirstCopy[1];
                
                // protection against vanishing
                if (genV.Pt()<0.005)   genV.SetXYZM(0.03,0.04,0.,91.2);
                if (genZTT.Pt()<0.005) genZTT.SetXYZM(0.03,0.04,0.,91.2);
                if (genZ.Pt()<0.005)   genZ.SetXYZM(0.03,0.04,0.,91.2);
                
                n_genZ_mass = genZ.M();
                n_genZ_Pt   = genZ.Pt();
                n_genZ_Eta  = genZ.Eta();
                n_genZ_Phi  = genZ.Phi();
                
                n_genV_mass = genV.M();
                n_genV_Pt   = genV.Pt();
                n_genV_Eta  = genV.Eta();
                n_genV_Phi  = genV.Phi();
                
                n_genZTT_mass = genZTT.M();
                n_genZTT_Pt   = genZTT.Pt();
                n_genZTT_Eta  = genZTT.Eta();
                n_genZTT_Phi  = genZTT.Phi();
                
                if (muTausLV.size()==2) {
                    if (muTausLV[0].Pt()>muTausLV[1].Pt()) {
                        gen_mu1 = muTausLV[0];
                        gen_mu2 = muTausLV[1];
                    }
                    else {
                        gen_mu1 = muTausLV[1];
                        gen_mu2 = muTausLV[0];
                    }
                }
                
                n_gen_mu1_Pt  = gen_mu1.Pt();
                n_gen_mu1_Eta = gen_mu1.Eta();
                n_gen_mu1_Phi = gen_mu1.Phi();
                n_gen_mu2_Pt  = gen_mu2.Pt();
                n_gen_mu2_Eta = gen_mu2.Eta();
                n_gen_mu2_Phi = gen_mu2.Phi();
                
                n_gen_taus = promptTausFirstCopy.size();
                n_gen_mutaus = muTausLV.size();
                
                float genVMass = -1.0;
                float genVEta;
                float genVPt = -1.0;
                
                genVMass = genV.M();
                genVEta  = genV.Eta();
                genVPt   = genV.Pt();
                
                //----------applying Zptmass weight------------
                if (applyZptmassCorr){
                    float zptmassweight = 1;
                    if (genVMass>50.0 && genVPt > 0.0) {
                        if (genVMass>1000.) genVMass = 1000.;
                        if (genVPt>1000.)   genVPt = 1000.;
                        zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(genVMass),histZMassPtWeights->GetYaxis()->FindBin(genVPt));
                        //std::cout << "genV Mass     =  " << genVMass <<     "          genVPt Pt   =      "<< genVPt <<std::endl;
                        //std::cout << "Mass bin  =  " <<histZMassPtWeights->GetXaxis()->FindBin(genVMass) << "       Pt bin  =  " << histZMassPtWeights->GetYaxis()->FindBin(genVPt) << std::endl;
                        //std::cout << "ztptmassweight is " << zptmassweight<< std::endl;
                        weight = weight*zptmassweight;
                    }
                }
                
                if (applyCategoryWeights){
                    
                    float genVAbsEta = abs(genVEta);
                    if (genVAbsEta >7.) genVAbsEta = 7.;
                    if (genVPt >1000.) genVPt = 1000.;
                    if (genVMass > 1000.) genVMass = 1000.;
                    float jet0Cat =1.;
                    float boosted = 1.;
                    float vbf = 1.;
                    if (genVPt>0.0){
                        //std::cout << "genV Mass     =  " << genVMass <<     "          genVPt Pt   =      "<< genVPt <<     "        genVAbsEta    =" << genVAbsEta << std::endl;
                        jet0Cat = jet0Hist->GetBinContent(jet0Hist->GetXaxis()->FindBin(genVAbsEta), jet0Hist->GetYaxis()->FindBin(genVPt), jet0Hist->GetZaxis()->FindBin(genVMass));
                        //cout << "  Eta Bin   = "<<jet0Hist->GetXaxis()->FindBin(genVAbsEta) <<
                        //"    Pt Bin  = "<<jet0Hist->GetYaxis()->FindBin(genVPt)<<
                        //"    Mass Bin = "<<jet0Hist->GetZaxis()->FindBin(genVMass)<<endl;
                        
                        //cout << " jet0Cat  =  " <<jet0Cat<<endl;
                        
                        if (jet0Cat<= 0.0) jet0weight = 1;
                        else jet0weight = jet0Cat;
                        
                        //cout<< " jet0weight = "<< jet0weight << endl;
                        
                        boosted = boostedHist->GetBinContent(boostedHist->GetXaxis()->FindBin(genVAbsEta), boostedHist->GetYaxis()->FindBin(genVPt),boostedHist->GetZaxis()->FindBin(genVMass));
//                        std::cout<< "   Eta bin     = " <<boostedHist->GetXaxis()->FindBin(genVAbsEta)<<
//                        "        pt bin  =   " <<  boostedHist->GetYaxis()->FindBin(genVPt)<<
//                        "        Mass bin =  " <<  boostedHist->GetZaxis()->FindBin(genVMass)<< std::endl;
                        
                        //std::cout << " boosted    =     " << boosted << std::endl;
                        
                        if (boosted <= 0.0) boostedweight = 1;
                        else boostedweight = boosted;
                        //std::cout << " boostedweight  =    " << boostedweight << std::endl;
                        
                        vbf = vbfHist->GetBinContent(vbfHist->GetXaxis()->FindBin(genVAbsEta), vbfHist->GetYaxis()->FindBin(genVPt), vbfHist->GetZaxis()->FindBin(genVMass));
                        //std::cout << "vbf    =   " << vbf << std::endl;
                        if (vbf <= 0.0) vbfweight = 1;
                        else vbfweight = vbf;
                        //std::cout << " vbfweight  =    " << vbfweight << std::endl;
                    }
                }
                
                if (n_gen_taus==2&&n_gen_mutaus==2&&n_genZTT_mass>60&&n_genZTT_mass<120) {
                    //			    std::cout << "Here we are" << std::endl;
                    
                    histZTTGenWeightsH->Fill(0.,genweight);
                    bool accept = n_gen_mu1_Pt>20 && n_gen_mu2_Pt>10 && fabs(n_gen_mu1_Eta)<2.4 && fabs(n_gen_mu2_Eta)<2.4;
                    
                    //fot single muon trigger initially we are considering either of the muon require to set the trigger and now we have decided to take only leading to fire the trigger and hence just consider the trigger weight for leading pt and eta.
                    if (accept) {
                        
                        float weight1 = (float)SF_muonIdIso->get_ScaleFactor(n_gen_mu1_Pt,fabs(n_gen_mu1_Eta));
                        float weight2 = (float)SF_muonIdIso->get_ScaleFactor(n_gen_mu2_Pt,fabs(n_gen_mu2_Eta));
                        
                        float weightTrig_gen = 1;
                        weightTrig_gen = SF_muonTrig->get_ScaleFactor(n_gen_mu1_Pt, fabs(n_gen_mu1_Eta));
                        
                        float effweight_gen = weight1*weight2*weightTrig_gen;
                        
                        histGenCutsGenWeightsH->Fill(0.,genweight);
                        histGenCutsWeightsH->Fill(0.,genweight*effweight_gen);
                    }
                }
                
                TW->Fill();
                
                //-------Appying Pileup correction-------->
                
                
                if (applyPUreweighting_official) {
                    nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
                    puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
                    
                    weight *= puweight;
                    PUweightsOfficialH->Fill(puweight);
                }

                //----------Applying Tau selection------------------>
                
                if (applyTauTauSelection) {
                    unsigned int nTaus = 0;
                    if (analysisTree.gentau_count>0) {
                        //	  cout << "Generated taus present" << endl;
                        for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++itau) {
                            // cout << itau << endl; "  : pt = "
                            //		 << analysisTree.gentau_visible_pt[itau]
                            //		 << "   eta = " <<  analysisTree.gentau_visible_eta[itau]
                            //		 << "   mother = " << int(analysisTree.gentau_mother[itau]) << endl;
                            if (int(analysisTree.gentau_mother[itau])==3) nTaus++;
                        }
                        // std::cout << "nTaus = " << nTaus << std::endl;//check
                    }
                    bool notTauTau = nTaus < 2;
                    
                    if (selectZToTauTauMuMu&&notTauTau) {
                        //	    std::cout << "Skipping event..." << std::endl;
                        //	    cout << endl;
                        continue;
                    }
                    
                    if (!selectZToTauTauMuMu&&!notTauTau) {
                        //	    std::cout << "Skipping event..." << std::endl;
                        //	    cout << endl;
                        continue;
                    }
                    //	  cout << endl;
                }
                
                if (applyTopPtReweighting) {
                    float topPt = -1;
                    float antitopPt = -1;
                    for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
                        
                        if (analysisTree.genparticles_pdgid[igen]==6)
                            topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                                analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
                        
                        if (analysisTree.genparticles_pdgid[igen]==-6)
                            antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
                                                    analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
                        
                        
                    }
                    if (topPt>0&&antitopPt>0) {
                        topptweight = topPtWeight(topPt,antitopPt);
                        //cout << "toppt = " << topPt << "   antitoppt = " << antitopPt << "   weight = " << topptweight << endl;
                        
                        weight *= topptweight;
                    }
                }

            }

            run = int(analysisTree.event_run);
            lumi = int(analysisTree.event_luminosityblock);
            evt = int(analysisTree.event_nr);
            
            if (isData) {
                bool lumi = false;
                int n=analysisTree.event_run;
                int lum = analysisTree.event_luminosityblock;
                
                std::string num = std::to_string(n);
                std::string lnum = std::to_string(lum);
                for(const auto& a : periods)
                {
                    
                    if ( num.c_str() ==  a.name ) {
                        for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
                            if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
                        }
                        auto last = std::prev(a.ranges.end());
                        if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
                    }
                    
                }
                if (!lumi) continue;
                
            }
            
            //std::cout << "passed lumi" << endl;
            
            npv = analysisTree.primvertex_count;
            npu = analysisTree.numtruepileupinteractions;
            rho = analysisTree.rho;
            
            npartons = analysisTree.genparticles_noutgoing;
            
            
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
            
            //--------------btag discriminant---------
            unsigned int nBTagDiscriminant = 0;
            for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
                TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
                if (discr=BTagDiscriminator)
                    nBTagDiscriminant = iBTag;
            }

            //----------trigger---------------
            bool isTriggerMuon = false;
            for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
                TString trigName(it->first);
                if (trigName.Contains(MuonTriggerName)) {
                    if (it->second==1)
                        isTriggerMuon = true;
                }
            }
            
            if (!isTriggerMuon) continue;
            
            unsigned int nMuonFilter = 0;
            bool isMuonFilter = false;
            
            unsigned int nSingleMuonFilter = 0;
            bool isSingleMuonFilter = false;
            
            unsigned int nfilters = analysisTree.run_hltfilters->size();
            
            for (unsigned int i=0; i<nfilters; ++i) {
                TString HLTFilter(analysisTree.run_hltfilters->at(i));
                if (HLTFilter==MuonFilterName) {
                    nMuonFilter = i;
                    isMuonFilter = true;
                }
                if (HLTFilter==SingleMuonFilterName) {
                    nSingleMuonFilter = i;
                    isSingleMuonFilter = true;
                }
            }
            if (!isMuonFilter) {
                cout << "Filter " << MuonFilterName << " not found " << endl;
                exit(-1);
            }
            if (!isSingleMuonFilter) {
                cout << "Filter " << SingleMuonFilterName << " not found " << endl;
                continue;
            }

	    // MET Filters
	    metFilters_ = metFiltersPasses(analysisTree,metFlags);
	    badChargedCandidateFilter_ = metFiltersPasses(analysisTree,badChargedCandidateFlag);
	    badPFMuonFilter_ = metFiltersPasses(analysisTree,badPFMuonFlag);
	    badGlobalMuonFilter_ = metFiltersPasses(analysisTree,badGlobalMuonFlag);
	    muonBadTrackFilter_ = metFiltersPasses(analysisTree,muonBadTrackFlag);
	    chargedHadronTrackResolutionFilter_ = metFiltersPasses(analysisTree,chargedHadronTrackResolutionFlag);
        

            //---------pfmet-----
            float met_ex = analysisTree.pfmetcorr_ex;
            float met_ey = analysisTree.pfmetcorr_ey;
            
            met = TMath::Sqrt(met_ex*met_ex+met_ey*met_ey);
            metphi = TMath::ATan2(met_ex,met_ey);
            
            //cout << "met before correction       = " << met << endl;
            
            float met_ex_uncorr = met_ex;
            float met_ey_uncorr = met_ey;
            
            metcov00 = analysisTree.pfmetcorr_sigxx;
            metcov01 = analysisTree.pfmetcorr_sigxy;
            metcov10 = analysisTree.pfmetcorr_sigyx;
            metcov11 = analysisTree.pfmetcorr_sigyy;
            
            //-------pfmet_JES
            float met_ex_JES_Up = analysisTree.pfmetcorr_ex_JetEnUp;
            float met_ey_JES_Up = analysisTree.pfmetcorr_ey_JetEnUp;
            
            met_JES_Up = TMath::Sqrt(met_ex_JES_Up*met_ex_JES_Up + met_ey_JES_Up*met_ey_JES_Up);
            //cout << " met_JES_Up      = "<< met_JES_Up<<endl;
            
            float met_ex_JES_Down = analysisTree.pfmetcorr_ex_JetEnDown;
            float met_ey_JES_Down = analysisTree.pfmetcorr_ey_JetEnDown;
            
            met_JES_Down = TMath::Sqrt(met_ex_JES_Down*met_ex_JES_Down + met_ey_JES_Down*met_ey_JES_Down);
            //cout << " met_JES_Down      = "<< met_JES_Down <<endl;
            //---------pfmet UnclusteredJES
            float met_ex_UnclusteredJES_Up = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            float met_ey_UnclusteredJES_Up =analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            
            met_UnclusteredJES_Up = TMath::Sqrt(met_ex_UnclusteredJES_Up*met_ex_UnclusteredJES_Up +met_ey_UnclusteredJES_Up*met_ey_UnclusteredJES_Up);
            //cout << "met_UnclusteredJES_Up    = " << met_UnclusteredJES_Up << endl;
            
            float met_ex_UnclusteredJES_Down = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            float met_ey_UnclusteredJES_Down =analysisTree.pfmetcorr_ey_UnclusteredEnDown;
            
            met_UnclusteredJES_Down = TMath::Sqrt(met_ex_UnclusteredJES_Down*met_ex_UnclusteredJES_Down +met_ey_UnclusteredJES_Down*met_ey_UnclusteredJES_Down);
            //cout << "met_UnclusteredJES_Down    = " << met_UnclusteredJES_Down << endl;
            
            //uncorrected pfmet------
            met_uncorr = met;
            metphi_uncorr = metphi;
            
            met_JES_Up_uncorr = met_JES_Up;
            met_JES_Down_uncorr = met_JES_Down;
            
            met_UnclusteredJES_Up_uncorr = met_UnclusteredJES_Up;
            met_UnclusteredJES_Down_uncorr = met_UnclusteredJES_Down;
            
            //input for recoil-----met----
            float metcorr_ex = met_ex;
            float metcorr_ey = met_ey;
            
            float metcorr_ex_JES_Up = met_ex_JES_Up;
            float metcorr_ey_JES_Up = met_ey_JES_Up;
            
            float metcorr_ex_JES_Down = met_ex_JES_Down;
            float metcorr_ey_JES_Down = met_ey_JES_Down;
            
            float metcorr_ex_UnclusteredJES_Up = met_ex_UnclusteredJES_Up;
            float metcorr_ey_UnclusteredJES_Up = met_ey_UnclusteredJES_Up;
            
            float metcorr_ex_UnclusteredJES_Down = met_ex_UnclusteredJES_Down;
            float metcorr_ey_UnclusteredJES_Down = met_ey_UnclusteredJES_Down;
            
            //-------------------------
            
            // muon selection
            vector<unsigned int> allMuons; allMuons.clear();
            vector<unsigned int> idMuons; idMuons.clear();
            vector<unsigned int> isoMuons; isoMuons.clear();
            vector<float> isoMuonsValue; isoMuonsValue.clear();
            vector<float> allMuonsIso; allMuonsIso.clear();
            vector<bool> isMuonPassedIdIso; isMuonPassedIdIso.clear();
            vector<bool> isMuonMatchedSingleMuFilter; isMuonMatchedSingleMuFilter.clear();
            
            for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
                allMuons.push_back(im);
                bool muPassed = true;
                bool muSingleMatched = false;
                
                if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
                if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) continue;
                if (fabs(analysisTree.muon_dxy[im])>dxyMuonLooseCut) continue;
                if (fabs(analysisTree.muon_dz[im])>dzMuonLooseCut) continue;
                if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) muPassed = false;
                if (fabs(analysisTree.muon_dz[im])>dzMuonCut) muPassed = false;

		
		if (analysisTree.muon_isBad[im]) badMuonFilter_ = false;//Updated on Jun5,2017
		if (analysisTree.muon_isDuplicate[im]) duplicateMuonFilter_ = false;//updated on Jun5, 2017
                
                bool goodGlob =
                analysisTree.muon_isGlobal[im] &&
                analysisTree.muon_normChi2[im] < 3 &&
                analysisTree.muon_combQ_chi2LocalPosition[im] < 12 &&
                analysisTree.muon_combQ_trkKink[im] < 20;
                
                bool ichepMed = analysisTree.muon_isLoose[im] &&
                analysisTree.muon_validFraction[im] > 0.49 &&
                analysisTree.muon_segmentComp[im] > (goodGlob ? 0.303 : 0.451);
                
                //std::cout << "data id check"<< std::endl;
                
                if (isData) {
                    //std::cout << "data id check"<< std::endl;
                    
                    if (analysisTree.event_run<=278808){
                        if (!ichepMed) muPassed = false;
                        
                        //if (muPassed) std::cout << "Applied ICHEP ID for run < 278808 " <<std::endl;
                    }
                    else {
                        if (!analysisTree.muon_isMedium[im]) muPassed = false;
                        //if (muPassed) std::cout << "Applied medium ID for run > 278808" <<std::endl;
                    }
                }
                else {
                if (!analysisTree.muon_isMedium[im]) muPassed = false;
		//if (muPassed) std::cout << " Applied medium muon ID for MC"<< std::endl;
                    
                }
                
                if (muPassed) idMuons.push_back(im);
                
                float absIso = 0;
                
                
                absIso = analysisTree.muon_r04_sumChargedHadronPt[im];
                float neutralIso = analysisTree.muon_r03_sumNeutralHadronEt[im] +
                analysisTree.muon_r04_sumPhotonEt[im] -
                0.5*analysisTree.muon_r04_sumPUPt[im];
                neutralIso = TMath::Max(float(0),neutralIso);
                absIso += neutralIso;
                
                
                float relIso = absIso/analysisTree.muon_pt[im];
                allMuonsIso.push_back(relIso);
                if (relIso>isoMuonCut) muPassed = false;
                if (muPassed) {
                    isoMuons.push_back(im);
                    isoMuonsValue.push_back(relIso);
                }
                
                isMuonPassedIdIso.push_back(muPassed);
                
                for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
                    float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
                                          analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                    if (dRtrig>DRTrigMatch) continue;
                    if (analysisTree.trigobject_filters[iT][nSingleMuonFilter]
                        && analysisTree.trigobject_pt[iT]>singleMuonTriggerPtCut
                        && fabs(analysisTree.trigobject_eta[iT])<singleMuonTriggerEtaCut)
                        muSingleMatched = true;
//                    because MC has filters
//                    if (!isData) {
//                        muSingleMatched = true;
//                    }
                }
                
                isMuonMatchedSingleMuFilter.push_back(muSingleMatched);
                
            }

	    metFilters_ = metFilters_ && badMuonFilter_ && duplicateMuonFilter_; //added on Jun 5, 2017
	    if (metFilters_<0.5) continue; //added on Jun 5, 2017
	    
            unsigned int indx1 = 0;
            unsigned int indx2 = 0;
            bool isIsoMuonsPair = false;
            float isoMin = 9999;
            
            if (isoMuons.size()>0) {
                for (unsigned int im1=0; im1<isoMuons.size(); ++im1) {
                    
                    unsigned int index1 = isoMuons[im1];
                    bool isMu1matched = false;
                    
                    if (analysisTree.muon_pt[index1]<ptMuonHighCut) continue;
                    for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
                        float dRtrig = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1], analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                        if (dRtrig>DRTrigMatch) continue;
                        if (analysisTree.trigobject_filters[iT][nMuonFilter] && analysisTree.muon_pt[index1] > ptMuonHighCut && fabs(analysisTree.muon_eta[index1]) < etaMuonHighCut)
                            isMu1matched = true;
                    }
                    //if(!isData) isMu1matched = true;
//                    if (isMu1matched) {
//                        for (unsigned int iMu=0; iMu<allMuons.size(); ++iMu) {
//                            unsigned int indexProbe = allMuons[iMu];
//                            if (index1==indexProbe) continue;
//                            float q_1 = analysisTree.muon_charge[index1];
//                            float q2 = analysisTree.muon_charge[indexProbe];
//                            if (q_1*q2>0) continue;
//                            float dR = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
//                                              analysisTree.muon_eta[indexProbe],analysisTree.muon_phi[indexProbe]);
//                            if (dR<dRleptonsCut) continue;
//                            
//                            TLorentzVector muon1; muon1.SetXYZM(analysisTree.muon_px[index1],
//                                                                analysisTree.muon_py[index1],
//                                                                analysisTree.muon_pz[index1],
//                                                                muonMass);
//                            TLorentzVector muon2; muon2.SetXYZM(analysisTree.muon_px[indexProbe],
//                                                                analysisTree.muon_py[indexProbe],
//                                                                analysisTree.muon_pz[indexProbe],
//                                                                muonMass);
//                           
//                        }
//                    }
                    for (unsigned int im2=im1+1; im2<isoMuons.size(); ++im2) {
                        unsigned int index2 = isoMuons[im2];
                        float q_1 = analysisTree.muon_charge[index1];
                        float q_2 = analysisTree.muon_charge[index2];
                        
                        //bool isMu2matched = false;
                        
                        if (analysisTree.muon_pt[index2]<ptMuonLowCut) continue;//added 26Aug
//                        for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
//                            float dRtrig = deltaR(analysisTree.muon_eta[index2],analysisTree.muon_phi[index2],
//                                                  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
//                            if (dRtrig>DRTrigMatch) continue;
//                            if (analysisTree.trigobject_filters[iT][nMuonFilter] &&
//                                analysisTree.muon_pt[index2] > ptMuonHighCut &&
//                                fabs(analysisTree.muon_eta[index2]) < etaMuonHighCut)
//                                isMu2matched = true;
//                            //std::cout << "trailing muon pt   =   " << analysisTree.muon_pt[index2]<< std::endl;
//                        }
                        //if (!isData) isMu2matched = true;
                        os = q_1*q_2 < 0;
                        //bool isPairSelected = q_1*q_2 > 0;
                        //if (oppositeSign) isPairSelected = q_1*q_2 < 0;
                        bool isTriggerMatch = isMu1matched;
                        
                        float dRmumu = deltaR(analysisTree.muon_eta[index1],analysisTree.muon_phi[index1],
                                              analysisTree.muon_eta[index2],analysisTree.muon_phi[index2]);
                        //if (isTriggerMatch && isPairSelected && dRmumu>dRleptonsCut) {
			if (isTriggerMatch && dRmumu>dRleptonsCut) {
                            bool sumIso = isoMuonsValue[im1]+isoMuonsValue[im2];
                            if (sumIso<isoMin) {
                                isIsoMuonsPair = true;
                                isoMin = sumIso;
                                if (analysisTree.muon_pt[index1]>analysisTree.muon_pt[index2]) {
                                    indx1 = index1;
                                    indx2 = index2;
                                }
                                else {
                                    indx2 = index1;
                                    indx1 = index2;
                                }
                            }
                        }
                    }
                }
            }

            if (isIsoMuonsPair) {
                
                
                
                TLorentzVector mu1; mu1.SetXYZM(analysisTree.muon_px[indx1],
                                                analysisTree.muon_py[indx1],
                                                analysisTree.muon_pz[indx1],
                                                muonMass);
                TLorentzVector mu1_Up; mu1_Up.SetXYZM((1.00 + muonScale)*analysisTree.muon_px[indx1],
                                                      (1.00 + muonScale)*analysisTree.muon_py[indx1],
                                                      (1.00 + muonScale)*analysisTree.muon_pz[indx1],
                                                      muonMass);
                TLorentzVector mu1_Down; mu1_Down.SetXYZM((1.00 - muonScale)*analysisTree.muon_px[indx1],
                                                          (1.00 - muonScale)*analysisTree.muon_py[indx1],
                                                          (1.00 - muonScale)*analysisTree.muon_pz[indx1],
                                                          muonMass);
                
                TLorentzVector mu2; mu2.SetXYZM(analysisTree.muon_px[indx2],
                                                analysisTree.muon_py[indx2],
                                                analysisTree.muon_pz[indx2],
                                                muonMass);
                TLorentzVector mu2_Up; mu2_Up.SetXYZM((1.00 + muonScale)*analysisTree.muon_px[indx2],
                                                      (1.00 + muonScale)*analysisTree.muon_py[indx2],
                                                      (1.00 + muonScale)*analysisTree.muon_pz[indx2],
                                                      muonMass);
                TLorentzVector mu2_Down; mu2_Down.SetXYZM((1.00 - muonScale)*analysisTree.muon_px[indx2],
                                                          (1.00 - muonScale)*analysisTree.muon_py[indx2],
                                                          (1.00 - muonScale)*analysisTree.muon_pz[indx2],
                                                          muonMass);
                genAccept = n_gen_taus==2 && n_gen_mutaus==2;
                genAccept = genAccept && genZTT.M()>60 && genZTT.M()<120;
                float ptLeadGen = TMath::Max(n_gen_mu1_Pt,n_gen_mu2_Pt);
                float ptTrailGen = TMath::Min(n_gen_mu1_Pt,n_gen_mu2_Pt);
                genAccept = genAccept && ptLeadGen>20;
                genAccept = genAccept && ptTrailGen>10;
                genAccept = genAccept && fabs(n_gen_mu1_Eta)<2.4;
                genAccept = genAccept && fabs(n_gen_mu2_Eta)<2.4;
                
                if(genAccept) histGenWeightH->Fill(1.,genweight);
                
                
                //-----------------filling muon variable----------->
                
                //Leading muon
                
                pt_1      = analysisTree.muon_pt[indx1];
                pt_1_Up   = (1.0+muonScale)* pt_1;
                pt_1_Down = (1.0-muonScale)* pt_1;
                q_1       = analysisTree.muon_charge[indx1];
                eta_1     = analysisTree.muon_eta[indx1];
                phi_1     = analysisTree.muon_phi[indx1];
                d0_1      = analysisTree.muon_dxy[indx1];
                dZ_1      = analysisTree.muon_dz[indx1];
                
                if (analysisTree.muon_dxyerr[indx1] != 0){
                    dcaSigd0_1 = log10(fabs(analysisTree.muon_dxy[indx1]/analysisTree.muon_dxyerr[indx1]));
                }
                if (analysisTree.muon_dzerr[indx1] != 0){
                    dcaSigdZ_1 = log10(fabs(analysisTree.muon_dz[indx1]/analysisTree.muon_dzerr[indx1]));
                }
                
                //Trailing muon
                
                pt_2      = analysisTree.muon_pt[indx2];
                pt_2_Up   = (1.0+muonScale)* pt_2;
                pt_2_Down = (1.0-muonScale)* pt_2;
                q_2       = analysisTree.muon_charge[indx2];
                eta_2     = analysisTree.muon_eta[indx2];
                phi_2     = analysisTree.muon_phi[indx2];
                d0_2      = analysisTree.muon_dxy[indx2];
                dZ_2      = analysisTree.muon_dz[indx2];
                
                if (analysisTree.muon_dxyerr[indx2] != 0){
                    dcaSigd0_2 = log10(fabs(analysisTree.muon_dxy[indx2]/analysisTree.muon_dxyerr[indx2]));
                }
                if (analysisTree.muon_dzerr[indx2] != 0){
                    dcaSigdZ_2 = log10(fabs(analysisTree.muon_dz[indx2]/analysisTree.muon_dzerr[indx2]));
                }
                
                //--------dimuon system-------->
                TLorentzVector dimuon = mu1 + mu2;
                TLorentzVector dimuon_Up = mu1_Up + mu2_Up;
                TLorentzVector dimuon_Down = mu1_Down + mu2_Down;
                float visZPx=dimuon.Px();
                float visZPy=dimuon.Py();
                
                m_vis = dimuon.M();
                m_vis_Up = dimuon_Up.M();
                m_vis_Down = dimuon_Down.M();
                m_vis_scaleUp = m_vis;
                m_vis_scaleDown = m_vis;
                m_vis_resoUp = m_vis;
                m_vis_resoDown = m_vis;
                
                pt_tt = dimuon.Pt();
                eta_tt = dimuon.Eta();
                dphi_tt = dPhiFrom2P(mu1.Px(), mu1.Py(), mu2.Px(), mu2.Py());
                dr_tt = deltaR(mu1.Eta(), mu1.Phi(), mu2.Eta(), mu2.Phi());

		pt_tt_Up = dimuon_Up.Pt();
                eta_tt_Up = dimuon_Up.Eta();
                dphi_tt_Up = dPhiFrom2P(mu1_Up.Px(), mu1_Up.Py(), mu2_Up.Px(), mu2_Up.Py());
                dr_tt_Up = deltaR(mu1_Up.Eta(), mu1_Up.Phi(), mu2_Up.Eta(), mu2_Up.Phi());

		pt_tt_Down = dimuon_Down.Pt();
                eta_tt_Down = dimuon_Down.Eta();
                dphi_tt_Down = dPhiFrom2P(mu1_Down.Px(), mu1_Down.Py(), mu2_Down.Px(), mu2_Down.Py());
                dr_tt_Down = deltaR(mu1_Down.Eta(), mu1_Down.Phi(), mu2_Down.Eta(), mu2_Down.Phi());
		
		float sumMuonPt = (mu1.Pt()+mu2.Pt());
		float sumMuonPt_Up = (mu1_Up.Pt()+mu2_Up.Pt());
		float sumMuonPt_Down = (mu1_Down.Pt()+mu2_Down.Pt());
		
                //float ptRatio = 0.0;
                if (sumMuonPt != 0)
                    ptRatio = (pt_tt/sumMuonPt);
		if (sumMuonPt_Up != 0)
                    ptRatio_Up = (pt_tt_Up/sumMuonPt_Up);
		if (sumMuonPt_Down != 0)
                    ptRatio_Down = (pt_tt_Down/sumMuonPt_Down);
            
                //dca of dimuon system
                for(unsigned int dimu=0; dimu<analysisTree.dimuon_count; ++dimu){
                    if (analysisTree.dimuon_dist2DE != 0){
                        sig2Mu2D = (analysisTree.dimuon_dist2D[dimu]/analysisTree.dimuon_dist2DE[dimu]);
                        dcasig2Mu2D = log10(fabs(analysisTree.dimuon_dist2D[dimu]/analysisTree.dimuon_dist2DE[dimu]));
                    }
                    
                    if (analysisTree.dimuon_dist3DE != 0){
                        sig2Mu3D = (analysisTree.dimuon_dist3D[dimu]/analysisTree.dimuon_dist3DE[dimu]);
                        dcasig2Mu3D = log10(fabs(analysisTree.dimuon_dist3D[dimu]/analysisTree.dimuon_dist3DE[dimu]));
                    }
                }

                //---------counting jets ---------->
                vector<unsigned int> jets; jets.clear();
                vector<unsigned int> jetsUp; jetsUp.clear();
                vector<unsigned int> jetsDown; jetsDown.clear();
                vector<unsigned int> jetspt20; jetspt20.clear();
                vector<unsigned int> bjets; bjets.clear();
                vector<unsigned int> bjets_Down; bjets_Down.clear();
                vector<unsigned int> bjets_Up; bjets_Up.clear();
                vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
                vector<unsigned int> bjetsRaw; bjetsRaw.clear();
                
                int indexLeadingJet = -1;
                float ptLeadingJet = -1;
                
                int indexSubLeadingJet = -1;
                float ptSubLeadingJet = -1;
                
                int indexLeadingBJet = -1;
                float ptLeadingBJet = -1;
                
                for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {
                    
                    float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                    float jetEta = analysisTree.pfjet_eta[jet];
                    if (absJetEta>jetEtaCut) continue;
                    
                    float jetPt = analysisTree.pfjet_pt[jet];
                    float jetPtDown = analysisTree.pfjet_pt[jet]*(1.0-analysisTree.pfjet_jecUncertainty[jet]);
                    float jetPtUp   = analysisTree.pfjet_pt[jet]*(1.0+analysisTree.pfjet_jecUncertainty[jet]);
                    //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
                    if (jetPtDown<jetPtLowCut) continue;
                    //std::cout << "jetpt =  " << jetPt << std::endl;
                    
                    //-----pfjetID------------
                    bool isPFJetId = looseJetiD(analysisTree,int(jet));
                    if (!isPFJetId) continue;
                    
                    bool cleanedJet = true;
                    
                    float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet], eta_1, phi_1);
                    
                    if (dR1<dRJetLeptonCut) cleanedJet = false;
                    
                    
                    float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet], eta_2, eta_2);
                    
                    if (dR2<dRJetLeptonCut) cleanedJet = false;
                    
                    //jet id
                    if (!cleanedJet) continue;
                    
                    if (jetPt>jetPtLowCut)
                        jetspt20.push_back(jet);
                    
                    if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance
                        
                        bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                        bool taggedUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                        bool taggedDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                        bool taggedRaw = tagged;
                        
                        if (!isData) {
                            int flavor = abs(analysisTree.pfjet_flavour[jet]);
                            
                            double jet_scalefactor = 1;
			    double jet_scalefactor_down = 1;
			    double jet_scalefactor_up = 1;
                            double JetPtForBTag = jetPt;
                            double tageff = 1;
                            
                            if (flavor==5) {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                                jet_scalefactor_down = reader_B.eval_auto_bounds("down",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                                jet_scalefactor_up = reader_B.eval_auto_bounds("up",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                                tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else if (flavor==4) {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                                jet_scalefactor_down = reader_C.eval_auto_bounds("down",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                                jet_scalefactor_up = reader_C.eval_auto_bounds("up",BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                                tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else {
                                if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                                if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                                jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                                jet_scalefactor_down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                                jet_scalefactor_up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                                tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                            }
                            
			    //			    std::cout << "SF : central = " << jet_scalefactor
			    //	      << "   up = " << jet_scalefactor_up
			    //	      << "   down = " << jet_scalefactor_down << std::endl;

                            if (tageff<1e-5)      tageff = 1e-5;
                            if (tageff>0.99999)   tageff = 0.99999;
                            rand.SetSeed((int)((jetEta+5)*100000));
                            double rannum = rand.Rndm();
                            
                            if (jet_scalefactor<1 && tagged) { // downgrade
                                double fraction = 1-jet_scalefactor;
                                if (rannum<fraction) {
                                    tagged = false;
                                    //		std::cout << "downgrading " << std::endl;
                                }
                            }
                            if (jet_scalefactor>1 && !tagged) { // upgrade
                                double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                                if (rannum<fraction) {
                                    tagged = true;
                                    //		std::cout << "upgrading " << std::endl;
                                }
                            }
                            if (jet_scalefactor_up<1 && taggedUp) { // downgrade
                                double fraction = 1-jet_scalefactor_up;
                                if (rannum<fraction) {
                                    taggedUp = false;
                                    //		std::cout << "downgrading " << std::endl;
                                }
                            }
                            if (jet_scalefactor_up>1 && !taggedUp) { // upgrade
                                double fraction = (jet_scalefactor_up-1.0)/(1.0/tageff-1.0);
                                if (rannum<fraction) {
                                    taggedUp = true;
                                    //		std::cout << "upgrading " << std::endl;
                                }
                            }
                            if (jet_scalefactor_down<1 && taggedDown) { // downgrade
                                double fraction = 1-jet_scalefactor_down;
                                if (rannum<fraction) {
                                    taggedDown = false;
                                    //		std::cout << "downgrading " << std::endl;
                                }
                            }
                            if (jet_scalefactor_down>1 && !taggedDown) { // upgrade
                                double fraction = (jet_scalefactor_down-1.0)/(1.0/tageff-1.0);
                                if (rannum<fraction) {
                                    taggedDown = true;
                                    //		std::cout << "upgrading " << std::endl;
                                }
                            }
                        }
                        if (taggedRaw)
                            bjetsRaw.push_back(jet);
                        
                        if (tagged) {
                            bjets.push_back(jet);
                            if (jetPt>ptLeadingBJet) {
                                    ptLeadingBJet = jetPt;
                                    indexLeadingBJet = jet;
                            }
                        }
                        
                        if (taggedUp) 
			  bjets_Up.push_back(jet);

			if (taggedDown)
			  bjets_Down.push_back(jet);

                    }
                    
                    if (jetPtUp>jetPtHighCut)
                        jetsUp.push_back(jet);
                    
                    if (jetPtDown>jetPtHighCut)
                        jetsDown.push_back(jet);
                    
                    if (jetPt>jetPtHighCut){
                        jets.push_back(jet);
                    
                        if (indexLeadingJet>=0) {
                            if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet) {
                                indexSubLeadingJet = jet;
                                ptSubLeadingJet = jetPt;
                            }
                        }
                        
                        if (jetPt>ptLeadingJet) {
                            //indexSubLeadingJet = indexLeadingJet;
                            //ptSubLeadingJet = ptLeadingJet;
                            indexLeadingJet = jet;
                            ptLeadingJet = jetPt;
                        }
                    }// This introduced to fix mjj
                }
                njets = jets.size();
                njets_Up = jetsUp.size();
                njets_Down = jetsDown.size();

		//-----jet uncertainty calculation-----

		int njetsMax = njets;
		
		if (!isData) {
		  jecUncertainties->runOnEvent(analysisTree,eta_1,phi_1,eta_2,phi_2);
		  
		  iUncert = 0;
		  for (auto const& Name : uncertNames) {
		    njetsUncUp[iUncert] = jecUncertainties->getNJets(Name,true);
		    njetsUncDown[iUncert] = jecUncertainties->getNJets(Name,false);
		    mjjUncUp[iUncert] = jecUncertainties->getMjj(Name,true);
		    mjjUncDown[iUncert] = jecUncertainties->getMjj(Name,false);
		    if (njetsUncUp[iUncert]>njetsMax) njetsMax = njetsUncUp[iUncert];
		    if (njetsUncDown[iUncert]>njetsMax) njetsMax = njetsUncDown[iUncert];
		    iUncert++;
		  }
		}

                int njetsCheckup = jecUncertainties->getNJets();
		
                njetspt20 = jetspt20.size();
                nbtag = bjets.size();
		nbtag_Up = bjets_Up.size();
		nbtag_Down = bjets_Down.size();
                nbtag_noSF = bjetsRaw.size();
                
                if (!isData) {
                    int nnbtag = nbtag_noSF;
                    btag0weight = 0;
                    btag0weight_Up = 0;
                    btag0weight_Down = 0;
                    
                    if (nnbtag<=2) {
                        
                        double b1Pt = 1;
                        double b2Pt = 1;
                        int b1Flav = 0;
                        int b2Flav = 0;
                        if (nnbtag>=1) {
                            int b1index = bjetsRaw.at(0);
                            b1Pt = analysisTree.pfjet_pt[b1index];
                            b1Flav = analysisTree.pfjet_flavour[b1index];
                        }
                        if (nnbtag==2) {
                            int b2index = bjetsRaw.at(1);
                            b2Pt = analysisTree.pfjet_pt[b2index];
                            b2Flav = analysisTree.pfjet_flavour[b2index];
                        }
                        
                        btag0weight = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,0,0));
                        btag0weight_Up = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,1,0));
                        btag0weight_Down = float(bTagEventWeight(nnbtag,b1Pt,b1Flav,b2Pt,b2Flav,1,-1,0));
                        
                    }
                    
                    //	      cout << "nbtag(raw) = " << nnbtag << "  weight(central) = " << btag0weight 
                    //		   << "   weight(up) = " << btag0weight_Up  
                    //		   << "   weight(down) = " << btag0weight_Down << endl;
                }
                
                bpt = -9999;
                beta = -9999;
                bphi = -9999;
                
                if (indexLeadingBJet>=0) {
                    bpt = analysisTree.pfjet_pt[indexLeadingBJet];
                    beta = analysisTree.pfjet_eta[indexLeadingBJet];
                    bphi = analysisTree.pfjet_phi[indexLeadingBJet];
                }
                
                jpt_1 = -9999;
                jpt_1_Up = -9999;
                jpt_1_Down = -9999;
                
                jeta_1 = -9999;
                jphi_1 = -9999;
                jptraw_1 = -9999;
                jptunc_1 = -9999;
                jmva_1 = -9999;
                jlrm_1 = -9999;
                jctm_1 = -9999;

                if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
                    cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
                
                if (indexLeadingJet>=0) {
                    jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
                    //std::cout << "jpt1 = "<< jpt_1 << std::endl;
                    jpt_1_Up = analysisTree.pfjet_pt[indexLeadingJet]*(1+analysisTree.pfjet_jecUncertainty[indexLeadingJet]);
                    jpt_1_Down = analysisTree.pfjet_pt[indexLeadingJet]*(1-analysisTree.pfjet_jecUncertainty[indexLeadingJet]);
                    jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
                    jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
                    jptraw_1 = analysisTree.pfjet_pt[indexLeadingJet]*analysisTree.pfjet_energycorr[indexLeadingJet];
                    jmva_1 = analysisTree.pfjet_pu_jet_full_mva[indexLeadingJet];
                }
                
                jpt_2 = -9999;
                jpt_2_Up = -9999;
                jpt_2_Down = -9999;
                
                jeta_2 = -9999;
                jphi_2 = -9999;
                jptraw_2 = -9999;
                jptunc_2 = -9999;
                jmva_2 = -9999;
                jlrm_2 = -9999;
                jctm_2 = -9999;
                
                if (indexSubLeadingJet>=0) {
                    jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
                    //std::cout << "jpt2 = "<< jpt_2 << std::endl;
                    jpt_2_Up = analysisTree.pfjet_pt[indexSubLeadingJet]*(1+analysisTree.pfjet_jecUncertainty[indexSubLeadingJet]);
                    jpt_2_Down = analysisTree.pfjet_pt[indexSubLeadingJet]*(1-analysisTree.pfjet_jecUncertainty[indexSubLeadingJet]);
                    jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
                    jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
                    jptraw_2 = analysisTree.pfjet_pt[indexSubLeadingJet]*analysisTree.pfjet_energycorr[indexSubLeadingJet];
                    jmva_2 = analysisTree.pfjet_pu_jet_full_mva[indexSubLeadingJet];
                }
                
                mjj =  -9999;
                jdeta =  -9999;
                njetingap = 0;
                
                if (indexLeadingJet>=0 && indexSubLeadingJet>=0){
                    
                    float unc1Up   = 1 + analysisTree.pfjet_jecUncertainty[indexLeadingJet];
                    float unc1Down = 1 - analysisTree.pfjet_jecUncertainty[indexLeadingJet];
                    
                    float unc2Up   = 1 + analysisTree.pfjet_jecUncertainty[indexSubLeadingJet];
                    float unc2Down = 1 - analysisTree.pfjet_jecUncertainty[indexSubLeadingJet];
                    
                    TLorentzVector jet1; jet1.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet],
                                                         analysisTree.pfjet_py[indexLeadingJet],
                                                         analysisTree.pfjet_pz[indexLeadingJet],
                                                         analysisTree.pfjet_e[indexLeadingJet]);
                    
                    TLorentzVector jet1Up; jet1Up.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet]*unc1Up,
                                                             analysisTree.pfjet_py[indexLeadingJet]*unc1Up,
                                                             analysisTree.pfjet_pz[indexLeadingJet]*unc1Up,
                                                             analysisTree.pfjet_e[indexLeadingJet]*unc1Up);
                    
                    TLorentzVector jet1Down; jet1Down.SetPxPyPzE(analysisTree.pfjet_px[indexLeadingJet]*unc1Down,
                                                                 analysisTree.pfjet_py[indexLeadingJet]*unc1Down,
                                                                 analysisTree.pfjet_pz[indexLeadingJet]*unc1Down,
                                                                 analysisTree.pfjet_e[indexLeadingJet]*unc1Down);
                    
                    TLorentzVector jet2; jet2.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet],
                                                         analysisTree.pfjet_py[indexSubLeadingJet],
                                                         analysisTree.pfjet_pz[indexSubLeadingJet],
                                                         analysisTree.pfjet_e[indexSubLeadingJet]);
                    
                    
                    TLorentzVector jet2Up; jet2Up.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet]*unc2Up,
                                                             analysisTree.pfjet_py[indexSubLeadingJet]*unc2Up,
                                                             analysisTree.pfjet_pz[indexSubLeadingJet]*unc2Up,
                                                             analysisTree.pfjet_e[indexSubLeadingJet]*unc2Up);
                    
                    TLorentzVector jet2Down; jet2Down.SetPxPyPzE(analysisTree.pfjet_px[indexSubLeadingJet]*unc2Down,
                                                                 analysisTree.pfjet_py[indexSubLeadingJet]*unc2Down,
                                                                 analysisTree.pfjet_pz[indexSubLeadingJet]*unc2Down,
                                                                 analysisTree.pfjet_e[indexSubLeadingJet]*unc2Down);
                    
                    mjj = (jet1+jet2).M();
                    mjj_Up = (jet1Up+jet2Up).M();
                    mjj_Down = (jet1Down+jet2Down).M();
                    
                    //	      std::cout << "mjj = " << mjj << " + " << mjj_Up << " - " << mjj_Down << std::endl;
                    
                    jdeta = abs(analysisTree.pfjet_eta[indexLeadingJet]- analysisTree.pfjet_eta[indexSubLeadingJet]);
                    
                    float etamax = analysisTree.pfjet_eta[indexLeadingJet];
                    float etamin = analysisTree.pfjet_eta[indexSubLeadingJet];
                    if (etamax<etamin) {
                        float tmp = etamax;
                        etamax = etamin;
                        etamin = tmp;
                    }
                    for (unsigned int jet=0; jet<jetspt20.size(); ++jet) {
                        int index = jetspt20.at(jet);
                        float etaX = analysisTree.pfjet_eta[index];
                        if (index!=indexLeadingJet&&index!=indexSubLeadingJet&&etaX>etamin&&etaX<etamax)
                            njetingap++;
                    }
                    
                }
                
                // Lepton SF, recoil, tracking eff----->
                
                if (!isData) {
                    
                    if (applyLeptonSF) {
                        //std::cout << "applying LeptonSF" << std::endl;
                        //Official Scale factor method
                        
                        isoweight_1=(float)SF_muonIdIso->get_ScaleFactor(double(pt_1), double(eta_1));
                        
                        isoweight_2=(float)SF_muonIdIso->get_ScaleFactor(double(pt_2),double(eta_2));
                        
                        MuSF_IdIso_Mu1H->Fill(isoweight_1);
                        MuSF_IdIso_Mu2H->Fill(isoweight_2);
                        
                        weight = weight*isoweight_1*isoweight_2;
                        
                        trigweight = (float)SF_muonTrig->get_ScaleFactor(double(pt_1), double(eta_1));
                        weight = weight* trigweight;
                        
                        effweight = isoweight_1*isoweight_2*trigweight;
                        
                        
                    }
                    
                    if (applyTrackEff){
                        double efficiency1 = eff(double(eta_1));
                        double efficiency2 = eff(double(eta_2));
                        double trackeff = 1 - (1-efficiency1)*(1-efficiency2);
                        if (trackeff > 0){
                            double trkeffweight = trackeff;
                            weight *=trkeffweight;
                        }
                    }
            
                    if (applyMEtRecoilCorrections) {
                        //-----pfmet------
		      recoilPFMetCorrector.CorrectByMeanResolution(met_ex, met_ey, genV.Px(), genV.Py(), genL.Px(), genL.Py(),
								   njets,metcorr_ex,metcorr_ey);
		      // std::cout << "PFMet : (" << met_ex << "," << met_ey << ")  " << "  (" << metcorr_ex << "," << metcorr_ey << ")" << std::endl;
		      recoilPFMetCorrector.CorrectByMeanResolution(met_ex_JES_Up, met_ey_JES_Up, genV.Px(), genV.Py(), genL.Px(), genL.Py(),
								   njets,metcorr_ex_JES_Up, metcorr_ey_JES_Up);
		      recoilPFMetCorrector.CorrectByMeanResolution(met_ex_JES_Down, met_ey_JES_Down, genV.Px(), genV.Py(), genL.Px(), genL.Py(),
								   njets,metcorr_ex_JES_Down, metcorr_ey_JES_Down);
		      recoilPFMetCorrector.CorrectByMeanResolution(met_ex_UnclusteredJES_Up, met_ey_UnclusteredJES_Up, genV.Px(), genV.Py(),
								   genL.Px(), genL.Py(), njets, metcorr_ex_UnclusteredJES_Up, metcorr_ey_UnclusteredJES_Up);
		      recoilPFMetCorrector.CorrectByMeanResolution(met_ex_UnclusteredJES_Down, met_ey_UnclusteredJES_Down, genV.Px(), genV.Py(),
								   genL.Px(), genL.Py(), njets, metcorr_ex_UnclusteredJES_Down, metcorr_ey_UnclusteredJES_Down);
		    }
		}
		
                met = TMath::Sqrt(metcorr_ex*metcorr_ex+metcorr_ey*metcorr_ey);
                //std::cout << "corrected met =  "<< met << std::endl;
                metphi = TMath::ATan2(metcorr_ex,metcorr_ey);
                
                met_JES_Up = TMath::Sqrt(metcorr_ex_JES_Up*metcorr_ex_JES_Up + metcorr_ey_JES_Up*metcorr_ey_JES_Up);
                //cout<< "corrected met_JES_Up  = "<< met_JES_Up <<endl;
                met_JES_Down = TMath::Sqrt(metcorr_ex_JES_Down*metcorr_ex_JES_Down + metcorr_ey_JES_Down*metcorr_ey_JES_Down);
                //cout<< "corrected met_JES_Down  = "<< met_JES_Down <<endl;
                met_UnclusteredJES_Up = TMath::Sqrt(metcorr_ex_UnclusteredJES_Up*metcorr_ex_UnclusteredJES_Up
						    + metcorr_ey_UnclusteredJES_Up*metcorr_ey_UnclusteredJES_Up);
                //cout << "corrected met_UnclusteredJES_Up = "<< met_UnclusteredJES_Up<<endl;
                met_UnclusteredJES_Down = TMath::Sqrt(metcorr_ex_UnclusteredJES_Down*metcorr_ex_UnclusteredJES_Down
						      + metcorr_ey_UnclusteredJES_Down*metcorr_ey_UnclusteredJES_Down);
                //cout << "corrected met_UnclusteredJES_Down = "<< met_UnclusteredJES_Down<<endl;
                
                float genmet_ex = analysisTree.genmet_ex;
                float genmet_ey = analysisTree.genmet_ey;
                
                genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
                genmetphi = TMath::ATan2(genmet_ey,genmet_ex);

		//defining metx_mu and mety_mu used in MSVvariation
		px_mu1     = pt_1 * TMath::Cos(phi_1);
		py_mu1     = pt_1 * TMath::Sin(phi_1);
		 
		px_mu2     = pt_2 * TMath::Cos(phi_2);
		py_mu2     = pt_2 * TMath::Sin(phi_2);
		
		float px_mu1Up   = pt_1_Up * TMath::Cos(phi_1);
		float px_mu1Down = pt_1_Down * TMath::Cos(phi_1);

		float py_mu1Up   = pt_1_Up * TMath::Sin(phi_1);
		float py_mu1Down = pt_1_Down * TMath::Sin(phi_1);

		float px_mu2Up   = pt_2_Up * TMath::Cos(phi_2);
		float px_mu2Down = pt_2_Down * TMath::Cos(phi_2);

		float py_mu2Up   = pt_2_Up * TMath::Sin(phi_2);
		float py_mu2Down = pt_2_Down * TMath::Sin(phi_2);

		met_x = metcorr_ex;
		met_y = metcorr_ey;
		
		float metx_muUp = met_x + (px_mu1 + px_mu2) - (px_mu1Up + px_mu2Up);
		float mety_muUp = met_y + (py_mu1 + py_mu2) - (py_mu1Up + py_mu2Up);

		float met_Up = TMath::Sqrt(metx_muUp*metx_muUp+mety_muUp*mety_muUp);

		float metx_muDown = met_x + (px_mu1 + px_mu2) - (px_mu1Down + px_mu2Down);
		float mety_muDown = met_y + (py_mu1 + py_mu2) - (py_mu1Down + py_mu2Down);

		float met_Down = TMath::Sqrt(metx_muDown*metx_muDown + mety_muDown*mety_muDown);

		//Total transverse momentum calculations for (boosted categoty BDT trainning)
		pt_tot = -9999;
		pt_tot_Up = -9999;
		pt_tot_Down = -9999;
		pt_tot_JES_Up =-9999;
		pt_tot_JES_Down =-9999;
		pt_tot_UnclusteredJES_Up = -9999;
		pt_tot_UnclusteredJES_Down =-9999;

		pt_tot_x  = met_x + px_mu1 + px_mu2;
		pt_tot_y  = met_y + py_mu1 + py_mu2;

		float pt_totx_JESUp =  metcorr_ex_JES_Up + px_mu1 + px_mu2;
		float pt_toty_JESUp =  metcorr_ey_JES_Up + py_mu1 + py_mu2;
		
		float pt_totx_JESDown =  metcorr_ex_JES_Down + px_mu1 + px_mu2;
		float pt_toty_JESDown =  metcorr_ey_JES_Down + py_mu1 + py_mu2;

		float pt_totx_UnclusteredJESUp =  metcorr_ex_UnclusteredJES_Up + px_mu1 + px_mu2;
		float pt_toty_UnclusteredJESUp =  metcorr_ey_UnclusteredJES_Up + py_mu1 + py_mu2;
		
		float pt_totx_UnclusteredJESDown =  metcorr_ex_UnclusteredJES_Down + px_mu1 + px_mu2;
		float pt_toty_UnclusteredJESDown =  metcorr_ey_UnclusteredJES_Down + py_mu1 + py_mu2;

		float pt_tot_xmuUp   = metx_muUp + px_mu1Up + px_mu2Up;
		float pt_tot_xmuDown = metx_muDown + px_mu1Down +  px_mu2Down;

		float pt_tot_ymuUp   = mety_muUp + py_mu1Up + py_mu2Up;
		float pt_tot_ymuDown = mety_muDown + py_mu1Down + py_mu2Down;

		pt_tot = TMath::Sqrt(pt_tot_x*pt_tot_x + pt_tot_y*pt_tot_y);
		pt_tot_Up   = TMath::Sqrt(pt_tot_xmuUp*pt_tot_xmuUp + pt_tot_ymuUp*pt_tot_ymuUp);
		pt_tot_Down = TMath::Sqrt(pt_tot_xmuDown*pt_tot_xmuDown + pt_tot_ymuDown*pt_tot_ymuDown);
		pt_tot_JES_Up   = TMath::Sqrt(pt_totx_JESUp*pt_totx_JESUp + pt_toty_JESUp*pt_toty_JESUp);
		pt_tot_JES_Down = TMath::Sqrt(pt_totx_JESDown*pt_totx_JESDown + pt_toty_JESDown*pt_toty_JESDown);
		pt_tot_UnclusteredJES_Up   = TMath::Sqrt(pt_totx_UnclusteredJESUp*pt_totx_UnclusteredJESUp + pt_toty_UnclusteredJESUp*pt_toty_UnclusteredJESUp);
		pt_tot_UnclusteredJES_Down = TMath::Sqrt(pt_totx_UnclusteredJESDown*pt_totx_UnclusteredJESDown + pt_toty_UnclusteredJESDown*pt_toty_UnclusteredJESDown);
		
                //----bisector of dimuon transerve momenta for defining the dzeta variables-------
                
                // bisector of electron and muon transverse momenta
                float mu1UnitX = mu1.Px()/mu1.Pt();
                float mu1UnitY = mu1.Py()/mu1.Pt();
                
                float mu2UnitX = mu2.Px()/mu2.Pt();
                float mu2UnitY = mu2.Py()/mu2.Pt();
                
                float zetaX = mu1UnitX + mu2UnitX;
                float zetaY = mu1UnitY + mu2UnitY;
                
                float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
                
                zetaX = zetaX/normZeta;
                zetaY = zetaY/normZeta;
                
                float vectorX = met_ex + mu2.Px() + mu1.Px();
                float vectorY = met_ey + mu2.Py() + mu1.Py();
                
                float vectorVisX = mu2.Px() + mu1.Px();
                float vectorVisY = mu2.Py() + mu1.Py();

		float vectorVisX_Up = mu2_Up.Px() + mu1_Up.Px();
		float vectorVisY_Up = mu2_Up.Py() +mu1_Up.Py();
		
                float vectorVisX_Down = mu2_Down.Px() + mu1_Down.Px();
		float vectorVisY_Down = mu2_Down.Py() +mu1_Down.Py();
		
                pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
		float pzetavis_Up = vectorVisX_Up*zetaX + vectorVisY_Up*zetaY;
		float pzetavis_Down = vectorVisX_Down*zetaX + vectorVisY_Down*zetaY;
                
                // computation of DZeta variable
                // pfmet
                computeDzeta(metcorr_ex,metcorr_ey,
                             zetaX,zetaY,pzetavis,pzetamiss,dzeta); // corrected
                //cout<<"dzeta  = " << dzeta<<endl;
                computeDzeta(met_ex_uncorr, met_ey_uncorr,
                             zetaX,zetaY,pzetavis,pzetamiss_uncorr, dzeta_uncorr); //uncorrected
                //cout << "dzeta  uncorrected = " << dzeta_uncorr<< endl;
		computeDzeta(metx_muUp, mety_muUp,
			     zetaX,zetaY,pzetavis_Up,pzetamiss_Up,dzeta_Up);
		computeDzeta(metx_muDown, mety_muDown,
			     zetaX,zetaY,pzetavis_Down,pzetamiss_Down,dzeta_Up);
                computeDzeta(metcorr_ex_JES_Up,metcorr_ey_JES_Up,
                             zetaX,zetaY,pzetavis,pzetamiss_JES_Up,dzeta_JES_Up); // JESUp
                
                computeDzeta(metcorr_ex_JES_Down,metcorr_ey_JES_Down,
                             zetaX,zetaY,pzetavis,pzetamiss_JES_Down,dzeta_JES_Down); // JESDown
                
                computeDzeta(metcorr_ex_UnclusteredJES_Up,metcorr_ey_UnclusteredJES_Up,
                             zetaX,zetaY,pzetavis,pzetamiss_UnclusteredJES_Up,dzeta_UnclusteredJES_Up); // UnclusteredUp
                
                computeDzeta(metcorr_ex_UnclusteredJES_Down,metcorr_ey_UnclusteredJES_Down,
                             zetaX,zetaY,pzetavis,pzetamiss_UnclusteredJES_Down,dzeta_UnclusteredJES_Down); // UnclusteredDown
                
                //genmet
                computeDzeta(genmet_ex,genmet_ey,
                             zetaX,zetaY,pzetavis,pzetamiss_genmet, dzeta_genmet);
                
                //-------------------cos(theta*) discriminator for higgs----->
                
                double E1 = EFromPandM0(muonMass, pt_1, eta_1);
                //std::cout << "E1 =" << E1 << std::endl;
                double E2 = EFromPandM0(muonMass, pt_2, eta_2);
                //std::cout << "E2 =" << E2 << std::endl;
                TLorentzVector mu1_CM; mu1_CM.SetPxPyPzE(analysisTree.muon_px[indx1],
                                                         analysisTree.muon_py[indx1],
                                                         analysisTree.muon_pz[indx1],
                                                         E1);
                TLorentzVector mu2_CM; mu2_CM.SetPxPyPzE(analysisTree.muon_px[indx2],
                                                         analysisTree.muon_py[indx2],
                                                         analysisTree.muon_pz[indx2],
                                                         E2);
                TLorentzVector dimuon_CM = mu1_CM + mu2_CM;
                
                costheta_1 = cosRestFrame(dimuon_CM, mu1_CM);
                costheta_2 = cosRestFrame(dimuon_CM, mu2_CM);
                
                if (q_1 > 0){
                    costheta = costheta_1;
                    //std::cout << " Its a leading Muon.\n"
                    //<< "cos(omega*) =  " <<cosd<< std::endl;
                }
                else {
                    if (q_2>0) {
                        costheta = costheta_2;
                        //std::cout << "Its a trailing muon. \n"
                        //<< "cos(omega*) =  " <<cosd<< std::endl;
                        
                    }
                }
                
                //new check implementation
                
                // check pos and neg muon pt and eta
                //if (!oppositeSign){
		if(!os){
                    pos_mupt1 = -9999;
                    pos_mupt2 = -9999;
                    neg_mupt1 = -9999;
                    neg_mupt2 = -9999;
                    
                    if (q_1 > 0)
                        pos_mupt1 = pt_1;
                        
                        
                    else neg_mupt1 = pt_1;
                    
                    if (q_2 > 0)
                        pos_mupt2 = pt_2;
                        
                    else
                        neg_mupt2 = pt_2;
                }

                //dphi_mumet and phi_twomu
                dphi_mumet_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metcorr_ex,metcorr_ey);
                //cout << "dphi_mumet1 = " << dphi_mumet_1 << endl;
                dphi_mumet_uncorr_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],met_ex_uncorr, met_ey_uncorr);
                //cout << "dphi_mumet1 uncorrected = " << dphi_mumet_uncorr_1 << endl;
		dphi_mumet_Up_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metx_muUp,mety_muUp);
		dphi_mumet_Down_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metx_muDown,mety_muDown);
                dphi_mumet_JES_Up_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metcorr_ex_JES_Up,metcorr_ey_JES_Up);
                dphi_mumet_JES_Down_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metcorr_ex_JES_Down, metcorr_ey_JES_Down);
                dphi_mumet_UnclusteredJES_Up_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metcorr_ex_UnclusteredJES_Up, metcorr_ey_UnclusteredJES_Up);
                dphi_mumet_UnclusteredJES_Down_1 = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1],metcorr_ex_UnclusteredJES_Down, metcorr_ey_UnclusteredJES_Down);
                
                dphi_mumet_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metcorr_ex,metcorr_ey);
                //cout << "dphi_mumet2 = " << dphi_mumet_2 << endl;
                dphi_mumet_uncorr_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],met_ex_uncorr, met_ey_uncorr);
                //cout << "dphi_mumet2 uncorrected = " << dphi_mumet_uncorr_2 << endl;
		dphi_mumet_Up_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metx_muUp,mety_muUp);
		dphi_mumet_Down_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metx_muDown,mety_muDown);
                dphi_mumet_JES_Up_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metcorr_ex_JES_Up,metcorr_ey_JES_Up);
                dphi_mumet_JES_Down_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metcorr_ex_JES_Down, metcorr_ey_JES_Down);
                dphi_mumet_UnclusteredJES_Up_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metcorr_ex_UnclusteredJES_Up, metcorr_ey_UnclusteredJES_Up);
                dphi_mumet_UnclusteredJES_Down_2 = dPhiFrom2P(analysisTree.muon_px[indx2],analysisTree.muon_py[indx2],metcorr_ex_UnclusteredJES_Down, metcorr_ey_UnclusteredJES_Down);
                
                if (q_1>0){
                    dphi_posmu_met = dphi_mumet_1;
                    dphi_posmu_met_uncorr = dphi_mumet_uncorr_1;
		    dphi_posmu_met_Up = dphi_mumet_Up_1;
		    dphi_posmu_met_Down = dphi_mumet_Down_1;
                    dphi_posmu_met_JES_Up = dphi_mumet_JES_Up_1;
                    dphi_posmu_met_JES_Down = dphi_mumet_JES_Down_1;
                    dphi_posmu_met_UnclusteredJES_Up = dphi_mumet_UnclusteredJES_Up_1;
                    dphi_posmu_met_UnclusteredJES_Down = dphi_mumet_UnclusteredJES_Down_1;
                }
                else{
                    dphi_posmu_met = dphi_mumet_2;
                    dphi_posmu_met_uncorr = dphi_mumet_uncorr_2;
		    dphi_posmu_met_Up = dphi_mumet_Up_2;
		    dphi_posmu_met_Down = dphi_mumet_Down_2;
                    dphi_posmu_met_JES_Up = dphi_mumet_JES_Up_2;
                    dphi_posmu_met_JES_Down = dphi_mumet_JES_Down_2;
                    dphi_posmu_met_UnclusteredJES_Up = dphi_mumet_UnclusteredJES_Up_2;
                    dphi_posmu_met_UnclusteredJES_Down = dphi_mumet_UnclusteredJES_Down_2;
                }
                
                dphi_twomu = dPhiFrom2P(analysisTree.muon_px[indx1],analysisTree.muon_py[indx1], analysisTree.muon_px[indx2],analysisTree.muon_py[indx2]);

		//-----------Evaluating BDT central values-----------------
                if (applyBDT){

		  bdt_0jets = reader0jets->EvaluateMVA("BDT");
		  bdt_boosted = readerboost->EvaluateMVA("BDT");
		  bdt_vbf = readervbf->EvaluateMVA("BDT");

		  //Evaluating BDT Variations
		  if (applyBDTvariations){
		    bdt_0jets_Up                    = reader_0jetsMuScaleUp->EvaluateMVA("BDT");
		    bdt_0jets_Down                  = reader_0jetsMuScaleDown->EvaluateMVA("BDT");
		    bdt_0jets_JES_Up                = reader_0jetsJESUp->EvaluateMVA("BDT");
		    bdt_0jets_JES_Down              = reader_0jetsJESDown->EvaluateMVA("BDT");
		    bdt_0jets_UnclusteredJES_Up     = reader_0jetsUnclustJESUp->EvaluateMVA("BDT");
		    bdt_0jets_UnclusteredJES_Down   = reader_0jetsUnclustJESDown->EvaluateMVA("BDT");
		    
		    bdt_boosted_Up                  = reader_boostMuScaleUp->EvaluateMVA("BDT");
		    bdt_boosted_Down                = reader_boostMuScaleDown->EvaluateMVA("BDT");
		    bdt_boosted_JES_Up              = reader_boostJESUp->EvaluateMVA("BDT");
		    bdt_boosted_JES_Down            = reader_boostJESDown->EvaluateMVA("BDT");
		    bdt_boosted_UnclusteredJES_Up   = reader_boostUnclustJESUp->EvaluateMVA("BDT");
		    bdt_boosted_UnclusteredJES_Down = reader_boostUnclustJESDown->EvaluateMVA("BDT");
		    
		    bdt_vbf_Up                      = reader_vbfMuScaleUp->EvaluateMVA("BDT");
		    bdt_vbf_Down                    = reader_vbfMuScaleDown->EvaluateMVA("BDT");
		    bdt_vbf_JES_Up                  = reader_vbfJESUp->EvaluateMVA("BDT");
		    bdt_vbf_JES_Down                = reader_vbfJESDown->EvaluateMVA("BDT");
		    bdt_vbf_UnclusteredJES_Up       = reader_vbfUnclustJESUp->EvaluateMVA("BDT");
		    bdt_vbf_UnclusteredJES_Down     = reader_vbfUnclustJESDown->EvaluateMVA("BDT");
		  }
		}
                                
		//-----------computingSVfitmass--------------
		
		m_sv           = -9999;
		m_sv_muUp      = -9999;
		m_sv_muDown    = -9999;
		m_sv_scaleUp   = -9999;
		m_sv_scaleDown = -9999;
		m_sv_resoUp    = -9999;
		m_sv_resoDown  = -9999;
		float px_sv    = -9999;
		float py_sv    = -9999;
		float msvmet_ex= -9999;
		float msvmet_ey= -9999;
                
		if (applySVFit){
			double measuredMETx = metcorr_ex;
			double measuredMETy = metcorr_ey;
			
			// define MET covariance
			TMatrixD covMET(2, 2);
			
			// std::cout << "covmetxx " << n_covmet_xx << "\t covmetxy " << n_covmet_xy  << "\t covmetyy " << n_covmet_yy <<std::endl;
			
			covMET[0][0] = metcov00;
			covMET[1][0] = metcov10;
			covMET[0][1] = metcov01;
		  covMET[1][1] = metcov11;
		  
		  /* ///Earlier definition
		  // define lepton four vectors
		  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
                  
		  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3));
		  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3));
                  
		  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
		  algo.addLogM(false);
		  algo.shiftVisPt(true, inputFile_visPtResolution);
		  algo.integrateMarkovChain();
		  m_sv = algo.getMass(); // return value is in units of GeV
		  
		  bool algoVerify = false;
		  algoVerify = algo.isValidSolution();
		  if ( algoVerify ) {
		  //std::cout << "SVfit mass (pfmet)    = " << m_sv << std::endl;
		  } else {
		  std::cout << "sorry -- status of NLL is not valid [" << algoVerify << "]" << std::endl;
		  }
		  
		  */
		  //Updated definition
		  //defining muon 4-vectors
		  svFitStandalone::MeasuredTauLepton svFitMu1(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu2(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu1Up(svFitStandalone::kTauToMuDecay, pt_1_Up, eta_1, phi_1, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu1Down(svFitStandalone::kTauToMuDecay, pt_1_Down, eta_1, phi_1, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu2Up(svFitStandalone::kTauToMuDecay, pt_2_Up, eta_2, phi_2, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu2Down(svFitStandalone::kTauToMuDecay, pt_2_Down, eta_2, phi_2, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu1scaleUp(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu1scaleDown(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu2scaleUp(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu2scaleDown(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu1resoUp(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu1resoDown(svFitStandalone::kTauToMuDecay, pt_1, eta_1, phi_1, 105.658e-3);
		  
		  svFitStandalone::MeasuredTauLepton svFitMu2resoUp(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3);
		  svFitStandalone::MeasuredTauLepton svFitMu2resoDown(svFitStandalone::kTauToMuDecay, pt_2, eta_2, phi_2, 105.658e-3);
		  
		  //central value
		  if (njets==0 && bdt_0jets>0.2){
		    SVfitStandaloneAlgorithm algo = SVFitMassComputation(svFitMu1, svFitMu2, measuredMETx, measuredMETy, covMET, inputFile_visPtResolution);
		    
		    m_sv = algo.getMass(); // return value is in units of GeV
		    if ( !algo.isValidSolution() ) 
		      std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl; 
		    pt_sv = algo.pt();
		    eta_sv = algo.eta();
		    phi_sv = algo.phi();
		    
		    px_sv = pt_sv*TMath::Cos(phi_sv);
		    py_sv = pt_sv*TMath::Sin(phi_sv);
		    
		    msvmet_ex = px_sv - dimuon.Px();
		    msvmet_ey = py_sv - dimuon.Py();
		    
		    msvmet = TMath::Sqrt(msvmet_ex*msvmet_ex+msvmet_ey*msvmet_ey);
		    msvmetphi = TMath::ATan2(msvmet_ey,msvmet_ex);
		    
		    m_sv_muUp      = m_sv;
		    m_sv_muDown    = m_sv;
		    m_sv_scaleUp   = m_sv;
		    m_sv_scaleDown = m_sv;
		    m_sv_resoUp    = m_sv;
		    m_sv_resoDown  = m_sv;
		    
		  }
		  
		  else if ((njets==1 || (njets==2 && mjj<300) || njets>2) && bdt_boosted>0.0){
		    SVfitStandaloneAlgorithm algo = SVFitMassComputation(svFitMu1, svFitMu2, measuredMETx, measuredMETy, covMET, inputFile_visPtResolution);
		    m_sv = algo.getMass(); // return value is in units of GeV
		    if ( !algo.isValidSolution() ) 
		      std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl; 
		    pt_sv = algo.pt();
		    eta_sv = algo.eta();
		    phi_sv = algo.phi();
		    
		    px_sv = pt_sv*TMath::Cos(phi_sv);
		    py_sv = pt_sv*TMath::Sin(phi_sv);
		    
		    msvmet_ex = px_sv - dimuon.Px();
		    msvmet_ey = py_sv - dimuon.Py();
		    
		    msvmet = TMath::Sqrt(msvmet_ex*msvmet_ex+msvmet_ey*msvmet_ey);
		    msvmetphi = TMath::ATan2(msvmet_ey,msvmet_ex);
		    m_sv_muUp      = m_sv;
		    m_sv_muDown    = m_sv;
		    m_sv_scaleUp   = m_sv;
		    m_sv_scaleDown = m_sv;
		    m_sv_resoUp    = m_sv;
		    m_sv_resoDown  = m_sv;
		  }
		  else if ((njets==2 && mjj>300) && bdt_vbf>0.8){
		    SVfitStandaloneAlgorithm algo = SVFitMassComputation(svFitMu1, svFitMu2, measuredMETx, measuredMETy, covMET, inputFile_visPtResolution);
		    m_sv = algo.getMass(); // return value is in units of GeV
		    if ( !algo.isValidSolution() ) 
		      std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl; 
		    pt_sv = algo.pt();
		    eta_sv = algo.eta();
		    phi_sv = algo.phi();
		    
		    px_sv = pt_sv*TMath::Cos(phi_sv);
		    py_sv = pt_sv*TMath::Sin(phi_sv);
		    
		    msvmet_ex = px_sv - dimuon.Px();
		    msvmet_ey = py_sv - dimuon.Py();
		    
		    msvmet = TMath::Sqrt(msvmet_ex*msvmet_ex+msvmet_ey*msvmet_ey);
		    msvmetphi = TMath::ATan2(msvmet_ey,msvmet_ex);
//		    m_sv_muUp      = m_sv;
//		    m_sv_muDown    = m_sv;
//		    m_sv_scaleUp   = m_sv;
//		    m_sv_scaleDown = m_sv;
//		    m_sv_resoUp    = m_sv;
//		    m_sv_resoDown  = m_sv;
		  }
		  
		  if (applyMSVvariation){
		    //SVfit scale variation with BDT cuts
		    if (njets==0 && bdt_0jets>0.2){
		      SVfitStandaloneAlgorithm algo_muUp = SVFitMassComputation(svFitMu1Up, svFitMu2Up, metx_muUp, mety_muUp,
										covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_muDown = SVFitMassComputation(svFitMu1Down, svFitMu2Down, metx_muDown, mety_muDown,
										  covMET, inputFile_visPtResolution);
		      
		      SVfitStandaloneAlgorithm algo_scaleUp = SVFitMassComputation(svFitMu1scaleUp, svFitMu2scaleUp, metcorr_ex_JES_Up, metcorr_ey_JES_Up,
										   covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_scaleDown = SVFitMassComputation(svFitMu1scaleDown, svFitMu2scaleDown, metcorr_ex_JES_Down,
										     metcorr_ey_JES_Down, covMET, inputFile_visPtResolution);
		      
		      SVfitStandaloneAlgorithm algo_resoUp = SVFitMassComputation(svFitMu1resoUp, svFitMu2resoUp, metcorr_ex_UnclusteredJES_Up,
										  metcorr_ey_UnclusteredJES_Up, covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_resoDown = SVFitMassComputation(svFitMu1resoDown, svFitMu2resoDown, metcorr_ex_UnclusteredJES_Down,
										    metcorr_ey_UnclusteredJES_Down, covMET, inputFile_visPtResolution);
					
		      m_sv_muUp      = algo_muUp.getMass();
		      m_sv_muDown    = algo_muDown.getMass();
		      
		      m_sv_scaleUp   = algo_scaleUp.getMass();
		      m_sv_scaleDown = algo_scaleDown.getMass();
		      
		      m_sv_resoUp    = algo_resoUp.getMass();
		      m_sv_resoDown  = algo_resoDown.getMass();
		      
		      pt_sv_muUp      = algo_muUp.pt();
		      eta_sv_muUp     = algo_muUp.eta();
		      phi_sv_muUp     = algo_muUp.phi();
		      
		      pt_sv_scaleUp   = algo_scaleUp.pt();
		      eta_sv_scaleUp  = algo_scaleUp.eta();
		      phi_sv_scaleUp  = algo_scaleUp.phi();
		      
		      pt_sv_resoUp   = algo_resoUp.pt();
		      eta_sv_resoUp  = algo_resoUp.eta();
		      phi_sv_resoUp  = algo_resoUp.phi();
		      
		      pt_sv_muDown  = algo_muDown.pt();
		      eta_sv_muDown = algo_muDown.eta();
		      phi_sv_muDown = algo_muDown.phi();
		      
		      pt_sv_scaleDown = algo_scaleDown.pt();
		      eta_sv_scaleDown= algo_scaleDown.eta();
		      phi_sv_scaleDown= algo_scaleDown.phi();
		      
		      pt_sv_resoDown  = algo_resoDown.pt();
		      eta_sv_resoDown = algo_resoDown.eta();
		      phi_sv_resoDown = algo_resoDown.phi();
		      
		    }
		    else if ((njets==1 || (njets==2 && mjj<300) || njets>2) && bdt_boosted>0.0){
		      SVfitStandaloneAlgorithm algo_muUp = SVFitMassComputation(svFitMu1Up, svFitMu2Up, metx_muUp, mety_muUp,
										covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_muDown = SVFitMassComputation(svFitMu1Down, svFitMu2Down, metx_muDown, mety_muDown,
										  covMET, inputFile_visPtResolution);
		      
		      SVfitStandaloneAlgorithm algo_scaleUp = SVFitMassComputation(svFitMu1scaleUp, svFitMu2scaleUp, metcorr_ex_JES_Up, metcorr_ey_JES_Up,
										   covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_scaleDown = SVFitMassComputation(svFitMu1scaleDown, svFitMu2scaleDown, metcorr_ex_JES_Down,
										     metcorr_ey_JES_Down, covMET, inputFile_visPtResolution);
		      
		      SVfitStandaloneAlgorithm algo_resoUp = SVFitMassComputation(svFitMu1resoUp, svFitMu2resoUp, metcorr_ex_UnclusteredJES_Up,
										  metcorr_ey_UnclusteredJES_Up, covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_resoDown = SVFitMassComputation(svFitMu1resoDown, svFitMu2resoDown, metcorr_ex_UnclusteredJES_Down,
										    metcorr_ey_UnclusteredJES_Down, covMET, inputFile_visPtResolution);
		      
		      m_sv_muUp   = algo_muUp.getMass();
		      m_sv_muDown = algo_muDown.getMass();
		      
		      m_sv_scaleUp   = algo_scaleUp.getMass();
		      m_sv_scaleDown = algo_scaleDown.getMass();
		      
		      m_sv_resoUp    = algo_resoUp.getMass();
		      m_sv_resoDown  = algo_resoDown.getMass();
		      
		      pt_sv_muUp  = algo_muUp.pt();
		      eta_sv_muUp = algo_muUp.eta();
		      phi_sv_muUp = algo_muUp.phi();
		      
		      pt_sv_scaleUp   = algo_scaleUp.pt();
		      eta_sv_scaleUp  = algo_scaleUp.eta();
		      phi_sv_scaleUp  = algo_scaleUp.phi();
		      
		      pt_sv_resoUp   = algo_resoUp.pt();
		      eta_sv_resoUp  = algo_resoUp.eta();
		      phi_sv_resoUp  = algo_resoUp.phi();
		      
		      pt_sv_muDown  = algo_muDown.pt();
		      eta_sv_muDown = algo_muDown.eta();
		      phi_sv_muDown = algo_muDown.phi();
		      
		      pt_sv_muDown  = algo_muDown.pt();
		      eta_sv_muDown = algo_muDown.eta();
		      phi_sv_muDown = algo_muDown.phi();
		      
		      pt_sv_scaleDown = algo_scaleDown.pt();
		      eta_sv_scaleDown= algo_scaleDown.eta();
		      phi_sv_scaleDown= algo_scaleDown.phi();
		      
		      pt_sv_resoDown  = algo_resoDown.pt();
		      eta_sv_resoDown = algo_resoDown.eta();
		      phi_sv_resoDown = algo_resoDown.phi();
		    }
		    else if ((njets==2 && mjj>300) && bdt_vbf>0.8){
		      SVfitStandaloneAlgorithm algo_muUp = SVFitMassComputation(svFitMu1Up, svFitMu2Up, metx_muUp, mety_muUp,
										covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_muDown = SVFitMassComputation(svFitMu1Down, svFitMu2Down, metx_muDown, mety_muDown,
										  covMET, inputFile_visPtResolution);
		      
		      SVfitStandaloneAlgorithm algo_scaleUp = SVFitMassComputation(svFitMu1scaleUp, svFitMu2scaleUp, metcorr_ex_JES_Up, metcorr_ey_JES_Up,
										   covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_scaleDown = SVFitMassComputation(svFitMu1scaleDown, svFitMu2scaleDown, metcorr_ex_JES_Down,
										     metcorr_ey_JES_Down, covMET, inputFile_visPtResolution);
					
		      SVfitStandaloneAlgorithm algo_resoUp = SVFitMassComputation(svFitMu1resoUp, svFitMu2resoUp, metcorr_ex_UnclusteredJES_Up,
										  metcorr_ey_UnclusteredJES_Up, covMET, inputFile_visPtResolution);
		      SVfitStandaloneAlgorithm algo_resoDown = SVFitMassComputation(svFitMu1resoDown, svFitMu2resoDown, metcorr_ex_UnclusteredJES_Down,
										    metcorr_ey_UnclusteredJES_Down, covMET, inputFile_visPtResolution);
		      
		      m_sv_muUp   = algo_muUp.getMass();
		      m_sv_muDown = algo_muDown.getMass();
		      
		      m_sv_scaleUp   = algo_scaleUp.getMass();
		      m_sv_scaleDown = algo_scaleDown.getMass();
		      
		      m_sv_resoUp    = algo_resoUp.getMass();
		      m_sv_resoDown  = algo_resoDown.getMass();
		      
		      pt_sv_muUp  = algo_muUp.pt();
		      eta_sv_muUp = algo_muUp.eta();
		      phi_sv_muUp = algo_muUp.phi();
		      
		      pt_sv_scaleUp   = algo_scaleUp.pt();
		      eta_sv_scaleUp  = algo_scaleUp.eta();
		      phi_sv_scaleUp  = algo_scaleUp.phi();
		      
		      pt_sv_resoUp   = algo_resoUp.pt();
		      eta_sv_resoUp  = algo_resoUp.eta();
		      phi_sv_resoUp  = algo_resoUp.phi();
		      
		      pt_sv_muDown  = algo_muDown.pt();
		      eta_sv_muDown = algo_muDown.eta();
		      phi_sv_muDown = algo_muDown.phi();
		      
		      pt_sv_scaleDown = algo_scaleDown.pt();
		      eta_sv_scaleDown= algo_scaleDown.eta();
		      phi_sv_scaleDown= algo_scaleDown.phi();
		      
		      pt_sv_resoDown  = algo_resoDown.pt();
		      eta_sv_resoDown = algo_resoDown.eta();
		      phi_sv_resoDown = algo_resoDown.phi();
		    }
		  }
		  
		}             
                bool recAccept = pt_1>20 && pt_2>10;
                recAccept = recAccept && fabs(eta_1)<2.4 && fabs(eta_2)<2.4;
                // filling ntuple
                if (genAccept) {
		  if (recAccept) {
                        histRecCutsWeightsH->Fill(0.,weight);
                        histRecCutsGenWeightsH->Fill(0.,genweight);
                        if (njets==0 && bdt_0jets>0.2) {
			  histBDTCutWeightsH->Fill(0.,weight);
			  histBDTCutGenWeightsH->Fill(0.,genweight);
                        }
			else if ((njets==1 || (njets==2 && mjj<300) || njets>2) && bdt_boosted>0.0){
			  histBDTCutWeightsH->Fill(0.,weight);
			  histBDTCutGenWeightsH->Fill(0.,genweight);
                        }
                        else if ((njets==2 && mjj>300) && bdt_vbf>0.8){
			  histBDTCutWeightsH->Fill(0.,weight);
			  histBDTCutGenWeightsH->Fill(0.,genweight);
			}
                    }
                }
                
                bool isZTTMM = n_gen_taus==2&&n_gen_mutaus==2;
                if (recAccept&&isZTTMM) {
                    noutDenRecCutsWeightsH->Fill(0.,weight);
                    noutDenRecCutsGenWeightsH->Fill(0.,genweight);
                    if (genAccept) {
                        noutNumRecCutsWeightsH->Fill(0.,weight);
                        noutNumRecCutsGenWeightsH->Fill(0.,genweight);
                    }
                    if (njets==0 && bdt_0jets>0.2) {
		      noutDenBDTCutWeightsH->Fill(0.,weight);
		      noutDenBDTCutGenWeightsH->Fill(0.,genweight);
		      if (genAccept) {
			noutNumBDTCutWeightsH->Fill(0.,weight);
			noutNumBDTCutGenWeightsH->Fill(0.,genweight);
		      }
                    }
		    else if ((njets==1 || (njets==2 && mjj<300) || njets>2) && bdt_boosted>0.0){
		      noutDenBDTCutWeightsH->Fill(0.,weight);
		      noutDenBDTCutGenWeightsH->Fill(0.,genweight);
		      if (genAccept) {
			noutNumBDTCutWeightsH->Fill(0.,weight);
			noutNumBDTCutGenWeightsH->Fill(0.,genweight);
		      }
		    }
		    else if ((njets==2 && mjj>300) && bdt_vbf>0.8){
		      noutDenBDTCutWeightsH->Fill(0.,weight);
		      noutDenBDTCutGenWeightsH->Fill(0.,genweight);
		      if (genAccept) {
			noutNumBDTCutWeightsH->Fill(0.,weight);
			noutNumBDTCutGenWeightsH->Fill(0.,genweight);
		      }
		    }
                }
                
                T->Fill();

            }
            
            if (isIsoMuonsPair) selEventsIsoMuons++;
            
        }// end of file processing (loop over events in one file)
        nFiles++;
        delete _tree;
        file_->Close();
        delete file_;
        
    }
    std::cout << std::endl;
    int allEvents = int(inputEventsH->GetEntries());
    std::cout << "Total number of input events                     = " << allEvents << std::endl;
    std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
    std::cout << "Total number of selected events (iso muon pairs) = " << selEventsIsoMuons << std::endl;
    std::cout << std::endl;
    std::cout << "RunMin = " << RunMin << std::endl;
    std::cout << "RunMax = " << RunMax << std::endl;
    
    //cout << "weight used:" << weight << std::endl;
    
    // using object as comp
    std::sort (allRuns.begin(), allRuns.end(), myobject);
    std::cout << "Runs : ";
    for (unsigned int iR=0; iR<allRuns.size(); ++iR)
        std::cout << " " << allRuns.at(iR);
    std::cout << std::endl;
    
    file->Write();
    file->Close();
    delete file;
    
}
