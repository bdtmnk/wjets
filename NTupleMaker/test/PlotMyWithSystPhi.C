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
//#include <cmath>
//WithoutNPV/
void PlotMyWithSyst(
            TString Variable = "mu_relIso[0]",
            TString Suffix = "_IsoWeight",
            TString xtitle = "#mu_{relIso} [GeV]",
            bool logY = false,
            bool syst = true,
            TString channel = "Wmnu"
          )
{

    //ModTDRStyle();
    TString Run = "ABCD";
    TString ytitle = "Events";
    TString suffix = "_Wmnu";
    TString directory = "/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/";

    bool blindData = false;
    bool MCDataDriven = true;
    int nbMin = 4;
    int nbMax = 11;
    bool plotLeg = true;
    int position = 0;
    // 0 - right, 1 - left, 2 - central
    bool showSignal = true;
    float QCDweight = 1;
    float EWweight = 1;
    double lumi = 0;
    
    if (Run=="ABC"){lumi=26.87;}
    if (Run=="D"){lumi=35.12;}  
    else lumi=26.87+35.12;
    TFile *file;

    //if (channel=="Wmnu")    
    file = new TFile(directory+"RootWithWeights"+"/"+Run+"/"+Variable+".root");

    cout<<directory+"RootWithWeights"+"/"+Run+"/"+Variable+".root"<<endl;	
    cout<<"Reading the histograms  :"<<endl;
    cout<<"data_obs_"+Variable+suffix<<endl;   
    TH1D *histData = (TH1D*)file->Get("data_obs_"+Variable+suffix);
      
    cout<<"DY :"<<endl;
    TH1D * DY = (TH1D*)file->Get("DY_"+Variable+suffix); 
    cout<<"QCD :"<<endl;
    TH1D * QCD = (TH1D*)file->Get("QCD_"+Variable+suffix); 
    cout<<"VV :"<<endl;
    TH1D * VV = (TH1D*)file->Get("VV_"+Variable+suffix); 
    cout<<"ST :"<<endl;
    TH1D * ST = (TH1D*)file->Get("ST_"+Variable+suffix); 
    cout<<"TT :"<<endl;
    TH1D * TT = (TH1D*)file->Get("TT_"+Variable+suffix); 
    cout<<"WJ :"<<endl;
    TH1D * WJ = (TH1D*)file->Get("WJ_"+Variable+suffix); 

    cout <<histData->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<DY->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<QCD->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<VV->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<ST->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<TT->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;
    cout <<WJ->GetBinContent(20) <<"  !!!!!!!!!!!!!!!1"<<endl;



	int nBins  =   histData->GetNbinsX();


  std::cout << "BKG : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/WJ->GetSumOfWeights() << "+/-" 
	    << TMath::Sqrt(histData->GetSumOfWeights())/WJ->GetSumOfWeights() << std::endl;

///////////////////////////////////////////////////////////////////////////////////////////////
    //ModTDRStyle();

    ModTDRStyle();

  
    TCanvas* canv1 = new TCanvas("c1", "c1");
    canv1->cd();
    std::vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
    pads[0]->SetLogy(logY);
    
    std::vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);\
	if (Variable.Contains("MT")) h = CreateAxisHists(2, histData, -3, 3-0.01);
	if (Variable.Contains("MT")) h = CreateAxisHists(2, histData, 0, 200-0.01);
	if (Variable.Contains("Iso")) h = CreateAxisHists(2, histData, 0, 0.3-0.01);
	if(!logY && Variable.Contains("pt")) h = CreateAxisHists(2, histData, 0, 200-0.01);
	if(!logY && Variable.Contains("pt[0]")) h = CreateAxisHists(2, histData, 25, 100-0.01);
	if(!logY && Variable.Contains("Ut")) h = CreateAxisHists(2, histData, 0, 200-0.01);
	if(!logY && Variable.Contains("Utr")) h = CreateAxisHists(2, histData, -100, 100-0.01);
	if(!logY && Variable.Contains("Ucol")) h = CreateAxisHists(2, histData, -150, 50-0.01);
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
    SetupTwoPadSplitAsRatio(pads, "Data/MC", true, 0.0, 2.1);
    StandardAxes(h[1]->GetXaxis(), h[0]->GetYaxis(),xtitle_ ,units);
    h[1]->GetYaxis()->SetNdivisions(4);
    h[1]->GetXaxis()->SetTitleOffset(1.2);
    h[1]->GetYaxis()->SetTitleOffset(2.0);
    pads[0]->cd();
    h[0]->GetYaxis()->SetTitleOffset(2.0);
    pads[1]->SetGrid(0,1);
    //it complains if the minimum is set to 0 and you try to set log y scale
    if(logY) h[0]->SetMinimum(1000);
    pads[0]->cd();
    
    // Setup legend
    TLegend *legend = PositionedLegend(0.55, 0.20, 3, 0.03);
    legend->SetTextFont(42);
    legend-> SetNColumns(2);
    histData->SetMarkerColor(1);
    histData->SetLineColor(1);
    histData->SetFillColor(1);
    histData->SetFillStyle(0);
    histData->SetLineWidth(2);
    histData->SetMarkerStyle(20);
    histData->SetMarkerSize(1.1);
    
    legend->AddEntry(histData, "Observed", "ple");
    
    InitHist(QCD,TColor::GetColor("#FFCCFF"));
//    InitHist(DY,TColor::GetColor("#DE5A6A"));
    InitHist(DY,kYellow);
//    InitHist(TT,TColor::GetColor("#9999CC"));
    InitHist(TT,kRed);
 //   InitHist(WJ,TColor::GetColor("#6F2D35"));
    InitHist(WJ,kBlue);
    InitHist(VV,TColor::GetColor("#4496C8"));
    InitHist(ST,TColor::GetColor("#FFCC66"));
/*
    WJ->SetFillColorAlpha(0.5);
    TT->SetFillColorAlpha(0.5);
    DY->SetFillColorAlpha(0.5);
    ST->SetFillColorAlpha(0.5);
    QCD->SetFillColorAlpha(0.5);  
    VV->SetFillColorAlpha(0.5);
*/    

    WJ->SetFillColorAlpha(kBlue,0.5);
    legend->AddEntry(WJ,"WJETS","f");
    legend->AddEntry(TT,"t#bar{t}","f");
    legend->AddEntry(DY,"DY","f");
    //  leg->AddEntry(VV,"VV+VVV","f");
    legend->AddEntry(ST,"ST","f");
    legend->AddEntry(QCD,"QCD","f");
    legend->AddEntry(VV,"VV","f");


    TH1D * DYn = (TH1D*)DY->Clone("DYn");
    TH1D * QCDn = (TH1D*)QCD->Clone("QCDn");
    TH1D * VVn = (TH1D*)VV->Clone("VVn");
    TH1D * STn = (TH1D*)ST->Clone("STn");
    TH1D * TTn = (TH1D*)TT->Clone("TTn");
    TH1D * WJn = (TH1D*)WJ->Clone("WJn");


	WJn->Add(WJ,TT,1,-1);
	TTn->Add(TT,DY,1,-1);
	DYn->Add(DY,ST,1,-1);
	STn->Add(ST,QCD,1,-1);
	QCDn->Add(QCD,VV,1,-1);

    WJn->Scale(EWweight);
    TTn->Scale(EWweight);
    DYn->Scale(EWweight);
    STn->Scale(EWweight); 
    QCDn->Scale(QCDweight);   
    VV->Scale(EWweight);

    QCD->Add(QCDn,VV);
    ST->Add(STn,QCD);
    DY->Add(DYn,ST);
    TT->Add(TTn,DY);
    WJ->Add(WJn,TT);


    WJ->SetMarkerStyle(0);
    TT->SetMarkerStyle(0);
    DY->SetMarkerStyle(0);
    ST->SetMarkerStyle(0); 
    QCD->SetMarkerStyle(0);   
    VV->SetMarkerStyle(0);

    WJ->SetMarkerSize(0);
    TT->SetMarkerSize(0);
    DY->SetMarkerSize(0);
    ST->SetMarkerSize(0); 
    QCD->SetMarkerSize(0);   
    VV->SetMarkerSize(0);


    WJ->SetMinimum(1000);
    WJ->Draw("hist sameh");
    TT->Draw("hist sameh");
    DY->Draw("hist sameh");
    ST->Draw("hist sameh"); 
    QCD->Draw("hist sameh");   
    VV->Draw("hist sameh");

    canv1->Update();
    
    if (blindData)
    {
        for (int iB=nbMin; iB<=nbMax; ++iB)
        {
            histData->SetBinContent(iB,-1);
            histData->SetBinError(iB,0);
        }
    }

  std::cout << "BKG : " << WJ->GetSumOfWeights() << " : " << WJ->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT : " << histData->GetSumOfWeights() << " : " << histData->Integral(0,nBins+1) << std::endl;
  std::cout << "DAT/BKG = " << histData->GetSumOfWeights()/WJ->GetSumOfWeights() << "+/-" 
	    << TMath::Sqrt(histData->GetSumOfWeights())/WJ->GetSumOfWeights() << std::endl;

///////////////// numerical values for MT comparison 
	float mean =0;
	float RMS =0;
	float chi2 =0;
	if (Variable.Contains("MT"))
	{

   	 TH1D * WWcut = (TH1D*)WJ->Clone("WWcut");
	 TH1D * DATAcut = (TH1D*)histData->Clone("DATAcut");
    	for (int iB=1; iB<=nBins; ++iB) 
		{
		if (iB < WWcut->GetXaxis()->FindBin(50) || iB > (WWcut->GetXaxis()->FindBin(120)-1)) {WWcut->SetBinContent(iB,0);DATAcut->SetBinContent(iB,0);/*cout<<WWcut->GetXaxis()->GetBinCenter(iB)<<endl;*/}
		}
	mean = WWcut->GetMean();
	RMS = WWcut->GetRMS();
	chi2 = DATAcut->Chi2Test(WWcut,"WWCHI2/NDF");
	cout <<"mean  "<<mean << endl;
	cout <<"RMS  "<<RMS << endl;
	cout <<"chi2  "<<chi2 << endl;
	//chi2 = DATAcut->Chi2Test(WWcut,"WWCHI2");
	//cout <<"chi2  "<<chi2 << endl;
	//cout <<"DATAcut->GetNDF()  "<<DATAcut->GetNDF() << endl;
	}

///////////////// numerical values for MT comparison end

    TH1D * bkgdErr = (TH1D*)WJ->Clone("bkgdErr");
    TH1D * bkgdErr1 = (TH1D*)WJ->Clone("bkgdErr1");
    TH1D * bkgdErr2 = (TH1D*)WJ->Clone("bkgdErr2");
    TH1D * bkgdErr3 = (TH1D*)WJ->Clone("bkgdErr3");
    TH1D * bkgdErr4 = (TH1D*)WJ->Clone("bkgdErr4");
    TH1D * bkgdErr5 = (TH1D*)WJ->Clone("bkgdErr5");

    TH1D * WJ_JetResUp;
    TH1D * WJ_JetResDown;
    TH1D * WJ_JetEnUp;
    TH1D * WJ_JetEnDown;
    TH1D * WJ_UnclusteredEnDown;
    TH1D * WJ_UnclusteredEnUp;

    TH1D * QCD_JetResUp;
    TH1D * QCD_JetResDown;
    TH1D * QCD_JetEnUp;
    TH1D * QCD_JetEnDown;
    TH1D * QCD_UnclusteredEnDown;
    TH1D * QCD_UnclusteredEnUp;

cout<<"Start Systematics"<<endl;

if (syst)


{
    
    TString TVariable = "pfmetcorr_pt";
    if(Variable.Contains("pfmet_pt") ){
         TVariable = "pfmetcorr_pt";
    }
    
    if(Variable.Contains("pfmet_pt") ){
         TVariable = "pfmetcorr_pt";
    }

    if(Variable.Contains("pfmet_pt") ){
         TVariable = "pfmetcorr_pt";
    }
    if (Variable.Contains("MTpuppi")){
         TVariable = "MTpuppi";
    }
    if (Variable.Contains("pfMT_corr")){
         TVariable = "MT";
    }
    if (Variable.Contains("puppi_pt")){
         TVariable = "puppi_pt";
    }
    if (Variable.Contains("pfMT_corr_XY")){
	TVariable = "MT_XY";
    }
    if (Variable.Contains("pfmet_pt_corr_XY")){
	TVariable = "pfmet_pt_XY";
    }

    TFile * file1 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_JetResUp.root");
    TFile * file2 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_JetResDown.root");
    TFile * file3 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_JetEnUp.root");
    TFile * file4 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_JetEnDown.root");
    TFile * file5 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_UnclusteredEnDown.root");
    TFile * file6 = new TFile(directory+"RootWithWeights/"+Run+"/"+TVariable+"_UnclusteredEnUp.root");


    WJ_JetResUp = (TH1D*)file1->Get("WJ_"+TVariable+"_JetResUp"+suffix); 
    WJ_JetResDown = (TH1D*)file2->Get("WJ_"+TVariable+"_JetResDown"+suffix); 

    WJ_JetEnUp = (TH1D*)file3->Get("WJ_"+TVariable+"_JetEnUp"+suffix); 
    WJ_JetEnDown = (TH1D*)file4->Get("WJ_"+TVariable+"_JetEnDown"+suffix); 
    
    WJ_UnclusteredEnDown = (TH1D*)file5->Get("WJ_"+TVariable+"_UnclusteredEnDown"+suffix); 
    WJ_UnclusteredEnUp = (TH1D*)file6->Get("WJ_"+TVariable+"_UnclusteredEnUp"+suffix); 

    QCD_JetResUp = (TH1D*)file1->Get("QCD_"+TVariable+"_JetResUp"+suffix); 
    QCD_JetResDown = (TH1D*)file2->Get("QCD_"+TVariable+"_JetResDown"+suffix); 
    QCD_JetEnUp = (TH1D*)file3->Get("QCD_"+TVariable+"_JetEnUp"+suffix); 
    QCD_JetEnDown = (TH1D*)file4->Get("QCD_"+TVariable+"_JetEnDown"+suffix); 
    QCD_UnclusteredEnDown = (TH1D*)file5->Get("QCD_"+TVariable+"_UnclusteredEnDown"+suffix); 
    QCD_UnclusteredEnUp = (TH1D*)file6->Get("QCD_"+TVariable+"_UnclusteredEnUp"+suffix); 




cout <<WJ_JetResUp->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout <<WJ_JetResDown->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout <<WJ_JetEnUp->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout <<WJ_JetEnDown->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout <<WJ_UnclusteredEnDown->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout <<WJ_UnclusteredEnUp->GetBinContent(3) <<"  !!!!!!!!!!!!!!!1"<<endl;
cout << "systematics!!!" <<endl;


    WJn->Scale(0);
    WJn->Add(WJ_JetResUp,QCD_JetResUp,1,-1);
    WJ_JetResUp->Add(WJn,QCD_JetResUp,EWweight,QCDweight);

    WJn->Scale(0);
    WJn->Add(WJ_JetResDown,QCD_JetResDown,1,-1);
    WJ_JetResDown->Add(WJn,QCD_JetResDown,EWweight,QCDweight);

    WJn->Scale(0);
    WJn->Add(WJ_JetEnUp,QCD_JetEnUp,1,-1);
    WJ_JetEnUp->Add(WJn,QCD_JetEnUp,EWweight,QCDweight);

    WJn->Scale(0);
    WJn->Add(WJ_JetEnDown,QCD_JetEnDown,1,-1);
    WJ_JetEnDown->Add(WJn,QCD_JetEnDown,EWweight,QCDweight);

    WJn->Scale(0);
    WJn->Add(WJ_UnclusteredEnDown,QCD_UnclusteredEnDown,1,-1);
    WJ_UnclusteredEnDown->Add(WJn,QCD_UnclusteredEnDown,EWweight,QCDweight);

    WJn->Scale(0);
    WJn->Add(WJ_UnclusteredEnUp,QCD_UnclusteredEnUp,1,-1);
    WJ_UnclusteredEnUp->Add(WJn,QCD_UnclusteredEnUp,EWweight,QCDweight);

}


    float errLumi = 0.06;
    float errMuon = 0.03;
    //float errElectron = 0.04;
    for (int iB=1; iB<=nBins; ++iB) {
        //QCD->SetBinError(iB,0);
        //VV->SetBinError(iB,0);
        //TT->SetBinError(iB,0);
        WJ->SetBinError(iB,0);
        //DY->SetBinError(iB,0);
        //ZTT->SetBinError(iB,0);
        //float eStat =  bkgdErr->GetBinError(iB);
        float eStat =  TMath::Sqrt(bkgdErr->GetBinContent(iB));
        float X = bkgdErr->GetBinContent(iB);
        float eLumi = errLumi * X;

	float eSyst_JetResUp = 0;
	float eSyst_JetResDown = 0;
	float eSyst_JetEnUp = 0;
	float eSyst_JetEnDown = 0;
	float eSyst_UnclusteredEnDown = 0;
	float eSyst_UnclusteredEnUp = 0;
	float eSyst_JetRes = 0;
	float eSyst_JetEn = 0;
	float eSyst_UnclusteredEn = 0;

if (syst)
{
	 eSyst_JetResUp = TMath::Abs(WJ_JetResUp->GetBinContent(iB)- WJ->GetBinContent(iB));
	 eSyst_JetResDown = TMath::Abs(WJ_JetResDown->GetBinContent(iB)- WJ->GetBinContent(iB));
	 eSyst_JetEnUp = TMath::Abs(WJ_JetEnUp->GetBinContent(iB)- WJ->GetBinContent(iB));
	 eSyst_JetEnDown = TMath::Abs(WJ_JetEnDown->GetBinContent(iB)- WJ->GetBinContent(iB));
	 eSyst_UnclusteredEnDown = TMath::Abs(WJ_UnclusteredEnDown->GetBinContent(iB)- WJ->GetBinContent(iB));
	 eSyst_UnclusteredEnUp = TMath::Abs(WJ_UnclusteredEnUp->GetBinContent(iB)- WJ->GetBinContent(iB));

	 eSyst_JetRes = TMath::Max(eSyst_JetResUp,eSyst_JetResDown);
	 eSyst_JetEn = TMath::Max(eSyst_JetEnUp,eSyst_JetEnDown);
	 eSyst_UnclusteredEn = TMath::Max(eSyst_UnclusteredEnUp,eSyst_UnclusteredEnDown);
}
        float eMuon = errMuon * X;
        //float eElectron = errElectron * X;
        //float eBkg = dummy->GetBinError(iB);
        //float Err = TMath::Sqrt(eStat*eStat);
	float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+
                            eSyst_JetRes*eSyst_JetRes+
                            eSyst_JetEn*eSyst_JetEn+
                            eSyst_UnclusteredEn*eSyst_UnclusteredEn+eMuon*eMuon);

	float Err1 = TMath::Sqrt(eStat*eStat+eMuon*eMuon+eLumi*eLumi);
	float Err2 = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eMuon*eMuon+eSyst_JetRes*eSyst_JetRes);
	float Err3 = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eSyst_JetRes*eSyst_JetRes+eMuon*eMuon+eSyst_JetEn*eSyst_JetEn);
//cout << eStat<< "  "<<eLumi <<"  "<< eMuon<<"  " << eSyst_JetRes<<"  " << eSyst_JetEn<<"  " << eSyst_UnclusteredEn<<"  " << iB<<endl;
//cout << Err<< "  "<<Err1 <<"  " << Err2<<"  " << Err3<<"  " <<endl;
        bkgdErr->SetBinError(iB,Err);
        bkgdErr1->SetBinError(iB,Err1);
        bkgdErr2->SetBinError(iB,Err2);
        bkgdErr3->SetBinError(iB,Err3);
    }

    
    bkgdErr->SetMarkerSize(0);    bkgdErr1->SetMarkerSize(0);    bkgdErr2->SetMarkerSize(0);    bkgdErr3->SetMarkerSize(0);    bkgdErr4->SetMarkerSize(0);  
/*    int new_idx = CreateTransparentColor(13,0.4);
    int new_idx1 = CreateTransparentColor(kGreen+2,1);
    int new_idx2 = CreateTransparentColor(kGreen,1);
    int new_idx3 = CreateTransparentColor(kYellow+1,1);*/

    int new_idx1 = CreateTransparentColor(kRed+1,1);
    int new_idx2 = CreateTransparentColor(kGreen,1);
    int new_idx3 = CreateTransparentColor(kBlue+2,1);
    int new_idx = CreateTransparentColor(kCyan+1,1);

    bkgdErr->SetFillColor(new_idx);
    bkgdErr1->SetFillColor(new_idx1);
    bkgdErr2->SetFillColor(new_idx2);
    bkgdErr3->SetFillColor(new_idx3);
    
    bkgdErr->SetFillStyle(3002);
    bkgdErr1->SetFillStyle(3002);
    bkgdErr2->SetFillStyle(3002);
    bkgdErr3->SetFillStyle(3002);


    bkgdErr->SetLineWidth(1);        
    bkgdErr1->SetLineWidth(1);        bkgdErr2->SetLineWidth(1);        bkgdErr3->SetLineWidth(1); 
    //bkgdErr->Draw("e2same");     //bkgdErr4->Draw("e2same");    bkgdErr3->Draw("e2same");    bkgdErr2->Draw("e2same");    bkgdErr1->Draw("e2same");


    legend->AddEntry(bkgdErr1, "Stat + Lumi + Muon uncertainty" , "F" );
if (syst)
{
    legend->AddEntry(bkgdErr2, "+ JER" , "F" );
    legend->AddEntry(bkgdErr3, "+ JES " , "F" );
    legend->AddEntry(bkgdErr, "+ Uncl" , "F" );
}

if (Variable.Contains("MT"))
{
    legend->AddEntry((TObject*)0, Form("Mean = %.1f", mean), "" );
    legend->AddEntry((TObject*)0, "", "" );
    legend->AddEntry((TObject*)0, Form("RMS = %.1f", RMS), "" );
    legend->AddEntry((TObject*)0, "", "" );
    legend->AddEntry((TObject*)0, Form("#chi_{2}/ndf = %.1f", chi2), "" );
}

    canv1->Update();

    TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
    TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
    
    TH1D * ratioErrH1 = (TH1D*)bkgdErr1->Clone("ratioErrH1");
    TH1D * ratioErrH2 = (TH1D*)bkgdErr2->Clone("ratioErrH2");
    TH1D * ratioErrH3 = (TH1D*)bkgdErr3->Clone("ratioErrH3");

    for (int iB=1; iB<=nBins; ++iB) {
        float x1 = histData->GetBinContent(iB);
        float x2 = WJ->GetBinContent(iB);
//cout <<iB<<"  "<< x1 <<"  "<< x2 << "   " << endl;
        ratioErrH->SetBinContent(iB,1.0);        ratioErrH1->SetBinContent(iB,1.0);        ratioErrH2->SetBinContent(iB,1.0);        ratioErrH3->SetBinContent(iB,1.0);    
        ratioErrH->SetBinError(iB,0.0);        ratioErrH1->SetBinError(iB,0.0);        ratioErrH2->SetBinError(iB,0.0);        ratioErrH3->SetBinError(iB,0.0);      
        float xBkg = bkgdErr->GetBinContent(iB);
        float errBkg = bkgdErr->GetBinError(iB);
        float errBkg1 = bkgdErr1->GetBinError(iB);
        float errBkg2 = bkgdErr2->GetBinError(iB);
        float errBkg3 = bkgdErr3->GetBinError(iB);
        if (xBkg>0) {
            float relErr = errBkg/xBkg;
            float relErr1 = errBkg1/xBkg;
            float relErr2 = errBkg2/xBkg;
            float relErr3 = errBkg3/xBkg;

            ratioErrH->SetBinError(iB,relErr);
            ratioErrH1->SetBinError(iB,relErr1);
            ratioErrH2->SetBinError(iB,relErr2);
            ratioErrH3->SetBinError(iB,relErr3);
            
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
    ratioErrH->GetYaxis()->SetRangeUser(0.0,2.0);
    //ratioErrH->GetXaxis()->SetRangeUser(50,120);
    //ratioErrH->Draw("e2same");

if (syst)
{
    ratioErrH->Draw("e2same");
    ratioErrH3->Draw("e2same");
    ratioErrH2->Draw("e2same");
}
    ratioErrH1->Draw("e2same");

    TLine *line = new TLine( h[1]->GetXaxis()->GetXmin(),1, h[1]->GetXaxis()->GetXmax()-0.01,1);
	if (Variable.Contains("MT")) line = new TLine( 0,1, 120-0.01,1);
	if (Variable.Contains("Iso")) line = new TLine( 0,1, 0.3-0.01,1);
	if(!logY && Variable.Contains("pt")) line = new TLine( 0,1, 200-0.01,1);
 	if(!logY && Variable.Contains("pt[0]")) line = new TLine( 25,1, 100-0.01,1);
	if(!logY && Variable.Contains("Ut")) line = new TLine( 0,1, 200-0.01,1); 
	if(!logY && Variable.Contains("Utr")) line = new TLine( -100,1, 100-0.01,1); 
	if(!logY && Variable.Contains("Ucol")) line = new TLine( -150,1, 50-0.01,1);
    line->SetLineColor(kBlack);
    line->Draw("same");

    //ratioH->SetMarkerStyle(2);
    ratioH->SetMarkerSize(0.5);
    ratioH->Draw("pe0same");


    pads[0]->cd();
    histData->Draw("pesame");
    
    FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);
    DrawCMSLogo(pads[0], "CMS", "", 11, 0.045, 0.035, 1.2);
    //DrawCMSLogo(pads[0], "CMS", "Preliminary", 11, 0.045, 0.001, 1.2);
    //DrawTitle(pads[0], "CMS Preliminary", 1);
    TString _lumi = "26.87";
    if (Run=="ABC") _lumi="26.87";
    if (Run=="D") _lumi = "35.12";
    _lumi = "58.82";

    DrawTitle(pads[0], _lumi+"fb^{-1} (13 TeV)", 3);
    if (channel=="Wmnu") DrawTitle(pads[0], "W to #mu#nu", 2);
       if (channel=="Wmnu") DrawTitle(pads[0], "W to #mu#nu", 2);
       if (channel=="Wenu") DrawTitle(pads[0], "W to e#nu", 2);
    FixBoxPadding(pads[0], legend, 0.05);
    legend->Draw();
    FixOverlay();
    canv1->Update();
    pads[0]->GetFrame()->Draw();
    if(logY) suffix=suffix+"_Log";
   
    //if (channel=="Wmnu") 
   canv1->Print("/nfs/dust/cms/user/dydukhle/STAU/METStudy/CMSSW_10_2_15_patch2/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/Final/Jan_Final_"+Variable+Suffix+suffix+".pdf");
    //canv1->Print("/nfs/dust/cms/user/bobovnii/new/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/test/WmnuPlots/Final/Final_"+Variable+Suffix+suffix+".png");


}
