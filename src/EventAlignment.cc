/*!
        \file                DeltaClusterAnalysis.cc
        \brief               Analysis to study the efficiency of the DUT 
        \author              S. Seif El Nasr-Storey 
        \date                02/11/16
        Support :            mail to : Sarah.Storey@cern.ch  
*/

#include "TROOT.h"
#include "TInterpreter.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include "EventAlignment.h"
using std::vector;
using std::map;

using namespace std;

const double pi = std::acos(-1);
const float maxChi2 = 1.0; 
const int minStrip = 0; //64; 
const int maxStrip = 254; //190;
const int sizeWindow = 50.0 ; 

EventAlignment::EventAlignment(int pNumberEvents_perGroup , int pMaxOffset ) :
  BeamAnaBase::BeamAnaBase()
{
  fNumberEvents_perGroup = pNumberEvents_perGroup;
  fMaxEventOffset = pMaxOffset;
  fL1_sigma=1;
}

void EventAlignment::eventLoop()
{

}
void EventAlignment::beginJob() 
{
 BeamAnaBase::beginJob();
  nEntries_ = analysisTree()->GetEntries();
  hist_ = outFile();
  setAddresses();
  histoBooker();
  analysisTree()->GetEntry(0);
  getCbcConfig(condEv()->cwd, condEv()->window);
  alParsLocal_ = aLparameteres();
  std::cout << "Correlation window = " << condEv()->cwd << std::endl ;
}
void EventAlignment::histoBooker()
{
  hist_->bookEventHistograms();
  hist_->bookTelescopeAnalysisHistograms();
  hist_->bookTrackMatchHistograms();
  hist_->bookFeI4Histograms();
  hist_->bookDUTHistograms("det0");
  hist_->bookDUTHistograms("det1");
  hist_->bookCBCHistograms("cbc0");
  hist_->bookCBCHistograms("cbc1");
  

} 
void EventAlignment::clearEvent() 
{
  BeamAnaBase::clearEvent();
}
void EventAlignment::endJob() 
{
  BeamAnaBase::endJob();
  hist_->closeFile();
}

EventAlignment::~EventAlignment()
{
  delete hist_;
}
void EventAlignment::clearTrackMatchhistograms()
{
  char buffer[250];
  
  std::vector<std::string> histoNames1{"deltaXPos_trkfei4", "deltaYPos_trkfei4"};//
  for(unsigned int i = 0 ; i < histoNames1.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch",histoNames1[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }

  std::vector<std::string> histoNames2{"deltaXPos_trkfei4_M", "deltaYPos_trkfei4_M"};//
  for(unsigned int i = 0 ; i < histoNames2.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch",histoNames2[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }
}
void EventAlignment::clearDUThistograms(const char* buffer)
{
  //std::string buffer =  pDetectorID;
  //sprintf(buffer,"det%d", pDetectorID);

  std::vector<std::string> histoNames{"Clusters_Width", "MatchedClusters_Width", "UnMatchedClusters_Width", "Stubs_Width", "MatchedStubs_Width" , "UnMatchedStubs_Width", "nClusters_noTrks", "nClusters_noFeI4", "nClusters_wTrks" , "nClusters_wFeI4" };
  for(unsigned int i = 0 ; i < histoNames.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName(buffer,histoNames[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }

  std::vector<std::string> histoNames2{"Clusters_noTrks", "Clusters_noFeI4", "Clusters_wTrks" , "Clusters_wTrks" , "hitRes" , "clusterRes" , "stubRes" };
  for(unsigned int i = 0 ; i < histoNames2.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName(buffer,histoNames2[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }


  std::vector<std::string> histoNames3{"deltaXhits", "deltaXhits_SingleHitEvents", "deltaXhits_SingleHitSingleTrackEvents" , "deltaXclusters" , "deltaXclusters_SingleClstrEvents" , "deltaXclusters_SingleClstrSingleTrkEvents" , "deltaXstubs"};
  for(unsigned int i = 0 ; i < histoNames3.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName(buffer,histoNames3[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }

  std::vector<std::string> histoNames4{"deltaXhits_Trks", "deltaXclusters_Trks" , "stubRes_Chi2" , "stubRes_TDC" ,"stubRes_CW"};
  for(unsigned int i = 0 ; i < histoNames4.size()  ; i++ )
  {
    TH2D* h0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName(buffer,histoNames4[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      for(unsigned int k = 0 ; k < h0->GetNbinsY() ; k++ ) 
      {
        h0->SetBinContent(j,k,0);
      }
    }
  }

  std::vector<std::string> histoNames5{"deltaX2D_SingleHits", "deltaX2D_Hits", "deltaX2D_SingleClusters" , "deltaX2D_Clusters" , "deltaX2D_Stubs" };
  for(unsigned int i = 0 ; i < histoNames5.size()  ; i++ )
  {
    TH2D* h0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName(buffer,histoNames5[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      for(unsigned int k = 0 ; k < h0->GetNbinsY() ; k++ ) 
      {
        h0->SetBinContent(j,k,0);
      }
    }
  }

  std::vector<std::string> histoNames6{"nFiducialTracks" , "nFiducialTracks_Both" };
  for(unsigned int i = 0 ; i < histoNames6.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName(buffer,histoNames6[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }

  std::vector<std::string> profileNames{"noiseOccupancy", "hitOccupancy" , "nTracks" , "noHitOccupancy" , "stubOccupancy" ,"noiseStubOccupancy" };
  for(unsigned int i = 0 ; i < profileNames.size()  ; i++ )
  {
    TProfile* h0 = dynamic_cast<TProfile*>(hist_->GetProfileByName(buffer,profileNames[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    { 
      h0->SetBinEntries(j,0); 
      h0->SetBinContent(j,0);
      h0->SetBinError(j,0);
    }
  }

  std::vector<std::string> profileNames2{"noiseOccupancy_xTrk", "hitOccupancy_xTrk" , "nTracks_xTrk" , "noHitOccupancy_xTrk" , "stubOccupancy_xTrk" ,"noiseStubOccupancy_xTrk" , "clusterEfficiency_perStrip" , "clusterEfficiency_c" , "clusterEfficiency",  "stubEfficiency_perStrip" , "stubEfficiency_c" , "clusterEfficiency"};
  for(unsigned int i = 0 ; i < profileNames2.size()  ; i++ )
  {
    TProfile* h0 = dynamic_cast<TProfile*>(hist_->GetProfileByName(buffer,profileNames2[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    { 
      h0->SetBinEntries(j,0); 
      h0->SetBinContent(j,0);
      h0->SetBinError(j,0);
    }
  }

  std::vector<std::string> histoNames7{"nMatchedTracks_TDC", "nMatchedClusters_TDC", "nUnMatchedClusters_TDC" , "nMatchedStubs_TDC" , "nUnMatchedStubs_TDC" , "deltaX2D_ClustersW" , "deltaX2D_StubsW"};
  for(unsigned int i = 0 ; i < histoNames7.size()  ; i++ )
  {
    TH2D* h0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName(buffer,histoNames7[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      for(unsigned int k = 0 ; k < h0->GetNbinsY() ; k++ ) 
      {
        h0->SetBinContent(j,k,0);
      }
    }
  }

  std::vector<std::string> histoNames8{"noisePerStrip"  };
  for(unsigned int i = 0 ; i < histoNames8.size()  ; i++ )
  {
    TH1D* h0 = dynamic_cast<TH1D*>(hist_->GetHistoByName(buffer,histoNames8[i]));
    for( unsigned int j = 0 ; j < h0->GetNbinsX() ; j++ )
    {  
      h0->SetBinContent(j,0);
    }
  }

}
void EventAlignment::clearHistograms()
{
  fFidTracks_Det0 = 0 ;
  fFidTracks_Det1 = 0 ; 
  fFidTracks_wHits_Det0 = 0 ;
  fFidTracks_wHits_Det1 = 0 ;
  fFidTracks_Single_Det0 = 0 ;
  fFidTracks_Single_Det1 = 0 ;
  fFidTracks_Single_wHits_Det0 = 0 ;
  fFidTracks_Single_wHits_Det1 = 0 ;
  fHits_Det0= 0 ; 
  fHits_Det1= 0 ; 
  fClusters_Det0 = 0; 
  fClusters_Det1 = 0 ;
  fFidTracks = 0 ;

  clearDUThistograms("det0");
  clearDUThistograms("det1");
  
  TH1D* hDet0_fTracks = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks_DUT"));  
  TH1D* hDet0_fTracks_SH_Events = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks_SHevents_DUT"));  
  TH1D* hDet0_fTracks_SC_Events = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks_SCevents_DUT"));  
  for(unsigned int i = 0 ; i < hDet0_fTracks->GetNbinsX() ; i++ )
  {
    hDet0_fTracks->SetBinContent(i,0);
    hDet0_fTracks_SH_Events->SetBinContent(i,0);
    hDet0_fTracks_SC_Events->SetBinContent(i,0);
  }

  

  TH2D* h2D_FeI4 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","HitMap")); 
  TH2D* h2D_Trks = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","TkHitMap")); 
  for( unsigned int i = 0 ; i < h2D_FeI4->GetNbinsX(); i++ )
  {
    for(unsigned int j = 0 ; j < h2D_FeI4->GetNbinsY() ; j++ )
    {
      h2D_FeI4->SetBinContent(i,j,0);
      h2D_Trks->SetBinContent(i,j,0);
    }
  }

  TH2D* hDet0_TDC2 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","ResTDC_Clusters_Det0"));  
  TH2D* hDet1_TDC2 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","ResTDC_Clusters_Det1"));  
  for(unsigned int i = 0 ; i < hDet0_TDC2->GetNbinsX() ; i++ )
  {
     for(unsigned int j = 0 ; j < hDet0_TDC2->GetNbinsY() ; j++ )
    {
      hDet0_TDC2->SetBinContent(i,j,0);
      hDet1_TDC2->SetBinContent(i,j,0);
    }
  }

  TH1D* hDet0_TDC = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","tdc_Cluster_dut0"));  
  TH1D* hDet1_TDC = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","tdc_Cluster_dut1"));  
  for(unsigned int i = 0 ; i < hDet0_TDC->GetNbinsX() ; i++ )
  {
    hDet0_TDC->SetBinContent(i,0);
    hDet1_TDC->SetBinContent(i,0);
  }

  TH1D* hDet0_Tks_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","TkImpact_Strip_Det0"));  
  TH1D* hDet1_Tks_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","TkImpact_Strip_Det1"));  

  TH1D* hDet0_STks_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleTkImpact_Strip_Det0"));  
  TH1D* hDet1_STks_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleTkImpact_Strip_Det1"));  

  TH1D* hDet0_Hits_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Strip_Det0"));  
  TH1D* hDet1_Hits_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Strip_Det1"));  
  
  TH1D* hDet0_SingleHits_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleHit_Strip_Det0"));  
  TH1D* hDet1_SingleHits_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleHit_Strip_Det1"));  
  
  TH1D* hDet0_Clusters_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Strip_Det0"));  
  TH1D* hDet1_Clusters_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Strip_Det1")); 

  
  TH1D* hDet0_SingleClusters_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleCluster_Strip_Det0"));  
  TH1D* hDet1_SingleClusters_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleCluster_Strip_Det1")); 

  TH1D* hStubs_Strips = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Stub_Strip"));  
  for( int i = 0 ; i< hDet0_Tks_Strips->GetNbinsX() ; i++ )
  {
    hDet0_Tks_Strips->SetBinContent(i,0);
    hDet1_Tks_Strips->SetBinContent(i,0);
    hDet0_Hits_Strips->SetBinContent(i,0);
    hDet1_Hits_Strips->SetBinContent(i,0);
    hDet0_Clusters_Strips->SetBinContent(i,0);
    hDet1_Clusters_Strips->SetBinContent(i,0);
    hStubs_Strips->SetBinContent(i,0);

    hDet0_SingleHits_Strips->SetBinContent(i,0);
    hDet1_SingleHits_Strips->SetBinContent(i,0);
    hDet0_SingleClusters_Strips->SetBinContent(i,0);
    hDet1_SingleClusters_Strips->SetBinContent(i,0);

    hDet0_STks_Strips->SetBinContent(i,0);
    hDet1_STks_Strips->SetBinContent(i,1);

  }

  TH1D* hDet0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXPos_dut0")); 
  TH1D* hDet1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXPos_dut1")); 
  TH1D* hSH_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXPos_dut0_SingleTks")); 
  TH1D* hSH_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXPos_dut1_SingleTks")); 

  TH1D* hCluster0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXCluster_dut0")); 
  TH1D* hCluster1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXCluster_dut1")); 
  TH1D* hSingleCluster0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaX_SingleCluster_dut0")); 
  TH1D* hSingleCluster1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaX_SingleCluster_dut1")); 
  

  TH1D* hStub = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXStub")); 
  for( int i = 0 ; i < hDet0->GetNbinsX() ; i++ )
  {
    hDet0->SetBinContent(i,0);
    hDet1->SetBinContent(i,0);
    hSH_Det0->SetBinContent(i,0);
    hSH_Det1->SetBinContent(i,0);
    hCluster0->SetBinContent(i,0);
    hCluster1->SetBinContent(i,0);
    hStub->SetBinContent(i,0);
    hSingleCluster0->SetBinContent(i,0);
    hSingleCluster1->SetBinContent(i,0);

  }
  hDet0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hDet1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hCluster0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hCluster0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hStub->GetXaxis()->SetRangeUser(-18.0, 18.0 );


  TH1D* hTkImpact_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","TkImpact_Det0")); 
  TH1D* hTkImpact_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","TkImpact_Det1"));
  TH1D* hHit_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Det0")); 
  TH1D* hHit_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Det1")); 
  TH1D* hCluster_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Det0")); 
  TH1D* hCluster_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Det1")); 
  TH1D* hStubs = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Stubs")); 

  TH1D* hSTkImpact_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleTkImpact_Det0")); 
  TH1D* hSTkImpact_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleTkImpact_Det1"));
  TH1D* hSHit_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleHit_Det0")); 
  TH1D* hSHit_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleHit_Det1")); 
  TH1D* hSCluster_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleCluster_Det0")); 
  TH1D* hSCluster_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleCluster_Det1")); 
  TH1D* hSStubs = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","SingleStubs")); 


  TH1D* hAllHit_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Strip_Det0_All")); 
  TH1D* hAllHit_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Hit_Strip_Det1_All")); 
  TH1D* hAllCluster_Det0 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Strip_Det0_All")); 
  TH1D* hAllCluster_Det1 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","Cluster_Strip_Det1_All"));

  for( int i = 0 ; i < hTkImpact_Det0->GetNbinsX() ; i++ )
  {
    hTkImpact_Det0->SetBinContent(i,0);
    hTkImpact_Det1->SetBinContent(i,0);
    hHit_Det0->SetBinContent(i,0);
    hHit_Det1->SetBinContent(i,0);
    hCluster_Det0->SetBinContent(i,0);
    hCluster_Det1->SetBinContent(i,0);
    hStubs->SetBinContent(i,0);

    hSTkImpact_Det0->SetBinContent(i,0);
    hSTkImpact_Det1->SetBinContent(i,0);
    hSHit_Det0->SetBinContent(i,0);
    hSHit_Det1->SetBinContent(i,0);
    hSCluster_Det0->SetBinContent(i,0);
    hSCluster_Det1->SetBinContent(i,0);
    hSStubs->SetBinContent(i,0);

    hAllHit_Det0->SetBinContent(i,0);
    hAllHit_Det1->SetBinContent(i,0);
    hAllCluster_Det0->SetBinContent(i,0);
    hAllCluster_Det1->SetBinContent(i,0);
  }
  hTkImpact_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hTkImpact_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hHit_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hHit_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hCluster_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hCluster_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hStubs->GetXaxis()->SetRangeUser(-18.0, 18.0 );


  hSTkImpact_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSTkImpact_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSHit_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSHit_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSCluster_Det0->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSCluster_Det1->GetXaxis()->SetRangeUser(-18.0, 18.0 );
  hSStubs->GetXaxis()->SetRangeUser(-18.0, 18.0 );

  TH2D* h2D_dut0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","deltaXPos_dut0_2D")); 
  TH2D* h2D_dut1 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","deltaXPos_dut1_2D")); 
  for( int i = 0 ; i < h2D_dut0->GetNbinsX(); i++ )
  {
    for( int j = 0 ; j < h2D_dut0->GetNbinsY() ; j++ )
    {
      h2D_dut0->SetBinContent(i,j,0);
      h2D_dut1->SetBinContent(i,j,0);
    }
  }
  h2D_dut0->GetXaxis()->SetRangeUser(-18.0,18.0);
  h2D_dut0->GetYaxis()->SetRangeUser(-18.0,18.0);
  h2D_dut1->GetXaxis()->SetRangeUser(-18.0,18.0);
  h2D_dut1->GetYaxis()->SetRangeUser(-18.0,18.0);

  TH2D* h2D_cluster0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","wClusterXPos_dut0")); 
  TH2D* h2D_cluster1 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","wClusterXPos_dut1")); 
  for( int i = 0 ; i < h2D_cluster0->GetNbinsX(); i++ )
  {
    for( int j = 0 ; j < h2D_cluster0->GetNbinsY() ; j++ )
    {
      h2D_cluster0->SetBinContent(i,j,0);
      h2D_cluster1->SetBinContent(i,j,0);
    }
  }

  TH2D* h2D_clusters0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","Clusters_Det0")); 
  TH2D* h2D_clusters1 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","Clusters_Det1")); 
  TH2D* h2D_Tks0 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","Tks_Det0")); 
  TH2D* h2D_Tks1 = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("TrackMatch","Tks_Det1")); 
  for( int i = 0 ; i < h2D_clusters0->GetNbinsX(); i++ )
  {
    for( int j = 0 ; j < h2D_clusters0->GetNbinsY() ; j++ )
    {
      h2D_clusters0->SetBinContent(i,j,0);
      h2D_clusters1->SetBinContent(i,j,0);
    
      h2D_Tks0->SetBinContent(i,j,0);
      h2D_Tks1->SetBinContent(i,j,0);
    }
  }

  


  TH1D* h_nTrack = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks"));
  TH1D* h_nTracks2 = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks2"));
  TH1D* h_nTracks_M = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","nTracks_M")); 
  TH1F* h_Chi2 = dynamic_cast<TH1F*>(hist_->GetHistoByName("TrackMatch","hTkChi2")); 
  TH1I* h_nClustersFei4 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","nClustersFei4"));
  TH1I* h_nClustersFei4_M = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","nClustersFei4_M")); 
  TH1I* h_lv1 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","lv1_Cluster"));  
  TH1I* h_nClusters_Fei4 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","nClusters_Fei4")); 
  TH1I* h_nClusters_Det0 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","nClusters_Det0")); 
  TH1I* h_nClusters_Det1 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","nClusters_Det1")); 
  TH1I* h_wClusters_Det0 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","wClusters_Det0")); 
  TH1I* h_wClusters_Det1 = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","wClusters_Det1")); 

  for( int i = 0 ; i < h_nTrack->GetNbinsX(); i++)
  {
    h_nTrack->SetBinContent(i,0);
    h_nTracks2->SetBinContent(i,0);
    h_nTracks_M->SetBinContent(i,0);
    h_nClustersFei4->SetBinContent(i,0);
    h_nClustersFei4_M->SetBinContent(i,0);
    h_lv1->SetBinContent(i,0);
    h_Chi2->SetBinContent(i,0);
    h_nClusters_Fei4->SetBinContent(i,0);
    h_nClusters_Det0->SetBinContent(i,0);
    h_nClusters_Det1->SetBinContent(i,0);
    h_wClusters_Det0->SetBinContent(i,0);
    h_wClusters_Det1->SetBinContent(i,0);
  }

  
  
  TH1I* h_ToT = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","hToT_UM")); 
  TH1I* h_ToT_UM = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatch","hToT_M"));  
  for( int i = 0 ; i < h_ToT->GetNbinsX() ;i++)
  {
    h_ToT->SetBinContent(i,0);
    h_ToT_UM->SetBinContent(i,0);
  }

  TH1D* h_DeltaCluster = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaXCluster"));  
  for( int i = 0 ;i  < h_DeltaCluster->GetNbinsX() ; i++)
  {
    h_DeltaCluster->SetBinContent(i,0);
  } 

}
void EventAlignment::lv1FeI4()
{
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 
  
  double cColPitch_fei4 = 250e-3;
  double cRowPitch_fei4 =  50e-3;

  int nEvents = 0; 
  for (Long64_t jentry=0; jentry<nEntries_;jentry++) 
  {
    clearEvent();
    Long64_t ientry = analysisTree()->GetEntry(jentry);

    if (ientry < 0) break;
    // if there are no hits in the Fei4 plane ... ignore  the event!
    if(telEv()->xPos->empty())    continue;
    // if there are no hits in the Fei4 plane ... ignore  the event!
    if(fei4Ev()->nPixHits==0  )    continue;
    
    if( fei4Ev()->nPixHits == 1 )
    {
       int cColumn = (int)(fei4Ev()->col->at(0)-1);
       int cRow  = (int)(fei4Ev()->row->at(0)-1); 
       double xval = x0_FEI4 + (cColumn)*(cColPitch_fei4); 
       double yval = y0_FEI4 + (cRow)*(cRowPitch_fei4);
       
       //double xval = x0_FEI4 + (cColumn)*(cColPitch_fei4); 
       //double yval = y0_FEI4 + (cRow)*(cRowPitch_fei4);
       hist_->fillHist2D("FeI4","HitMap_Raw", xval , yval );
    }

    std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
    int cMinWidth_FeI4Cluster = 3 ; 
    tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
    if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
    clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );

    for (unsigned int i= 0  ; i < xPos_FeI4Clusters.size(); i++) 
    {
      hist_->fillHist1D("FeI4","lv1", lv1_FeI4Clusters[i] );
      hist_->fillHist1D("FeI4","hTotFei4", tot_FeI4Clusters[i] );
      double Q = FeI4_ToT_to_Q(tot_FeI4Clusters[i])*1e-3; 
      hist_->fillHist1D("FeI4","hQFeI4", Q );
      hist_->fillHist2D("FeI4","HitMap", xPos_FeI4Clusters[i] , yPos_FeI4Clusters[i] );
       
    }
    hist_->fillHist1D("FeI4","nClusters", xPos_FeI4Clusters.size() );
      
    if( nEvents%100000 == 0 ) std::cout << "Finished FeI4 event " << nEvents << " ." << std::endl;
    nEvents++;
  }

  TH1D* htmp = dynamic_cast<TH1D*>(hist_->GetHistoByName("FeI4","lv1"));   
  double center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  fL1_mean = center;
  fL1_sigma = 1.0;

  double xMin = 0;
  double xMax = 100;
  htmp = dynamic_cast<TH1D*>(hist_->GetHistoByName("FeI4","hQFeI4"));  
  center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  TF1* f = new TF1("fToT", "landau", xMin , xMax );
  f->SetParameter(0, 1*htmp->GetMaximum() );
  f->SetParameter(1, 1*center );
  f->SetParameter(2, 2);
  f->SetLineColor(kGreen);
  htmp->Fit(f);
  fFeI4_Q_mpv = f->GetParameter(1);
  fFeI4_Q_sigma = f->GetParameter(2);
}
std::vector<double> EventAlignment::getRecoStubSeeds()
{
    std::vector<double> stubSeeds ; stubSeeds.clear();
    for(auto& stub : dutRecoStubmap()->at("C0"))  
    {
      auto& cluster = stub.seeding;
      double xStubPos = (cluster->x-nstrips()/2)*dutpitch(); 
      stubSeeds.push_back( xStubPos );
    }
    return stubSeeds;
}
std::vector<double> EventAlignment::getRecoStubMatches()
{
    std::vector<double> stubMatches ; stubMatches.clear();
    for(auto& stub : dutRecoStubmap()->at("C0"))  
    {
      auto& cluster = stub.matched;
      double xStubPos = (cluster->x-nstrips()/2)*dutpitch(); 
      stubMatches.push_back( xStubPos );
    }
    return stubMatches;
}
std::vector<double> EventAlignment::getRecoStubMatchWidths()
{
    std::vector<double> stubMatches ; stubMatches.clear();
    for(auto& stub : dutRecoStubmap()->at("C0"))  
    {
      auto& cluster = stub.matched;
      double xStubPos = (double)(cluster->size);
      stubMatches.push_back( xStubPos );
    }
    return stubMatches;
}
std::vector<double> EventAlignment::getRecoStubSeedWidths()
{
    std::vector<double> stubMatches ; stubMatches.clear();
    for(auto& stub : dutRecoStubmap()->at("C0"))  
    {
      auto& cluster = stub.seeding;
      double xStubPos = (double)(cluster->size);
      stubMatches.push_back( xStubPos );
    }
    return stubMatches;
}
std::vector<double> EventAlignment::getRecoStubPositions( )
{
    std::vector<double> stubPositions ; stubPositions.clear();
    for(auto& stub : dutRecoStubmap()->at("C0"))  
    {
      double xStubPos = (stub.x-nstrips()/2)*dutpitch(); 
#ifdef NOV_15  
        xStubPos = -1*xStubPos;
#endif
      stubPositions.push_back( xStubPos );
    }
    return stubPositions;
}
std::vector<double> EventAlignment::getCbcStubPositions( )
{
    std::vector<double> stubPositions ; stubPositions.clear();
    for(auto& stub : dutEv()->stubs)  
    {
      double xStubPos = (stub->x-nstrips()/2)*dutpitch(); 
#ifdef NOV_15 
        xStubPos = -1*xStubPos;
#endif      
      stubPositions.push_back( xStubPos );
    }
    return stubPositions;
}
std::vector<double> EventAlignment::getStubPositions()
{
  return getRecoStubPositions();
}

std::vector<double> EventAlignment::getClusterPositions( int pDetectorID , int maxClusterWidth )
{
   std::vector<double> clusterPositions ; clusterPositions.clear();
   std::string cDetectorName = "det" + std::to_string(pDetectorID) + "C0";
   for(auto& cl : dutRecoClmap()->at(cDetectorName) ) 
   {
      double xClusterPos = (cl.x-nstrips()/2)*dutpitch(); 
#ifdef NOV_15 
        xClusterPos = -1*xClusterPos;
#endif
      if( cl.size <= maxClusterWidth && cl.size > 0 ) clusterPositions.push_back( xClusterPos );
   }
   return clusterPositions;
}
std::vector<int> EventAlignment::getClusterWidths( int pDetectorID  , int maxClusterWidth )
{
   std::vector<int> clusterWidths ; clusterWidths.clear();
   std::string cDetectorName = "det" + std::to_string(pDetectorID) + "C0";
   std::string cDetectorName2 = "det" + std::to_string(pDetectorID) ;
   for(auto& cl : dutRecoClmap()->at(cDetectorName) ) 
   {  
      // I'm going to double count but whatever for now...
      hist_->fillHist1D(cDetectorName2,"clusterWidth_Raw" , cl.size );
      if( cl.size <= maxClusterWidth && cl.size > 0 ) clusterWidths.push_back( cl.size );
   }
   return clusterWidths;
}
std::vector<double> EventAlignment::getHitPositions( vector<int> detHits )
{
  std::vector<double> hitPositions ; hitPositions.clear();
  for( auto cHitDet  : detHits )
  {
    double xHitPos =  (cHitDet - nstrips()/2.0 )*dutpitch() ;
#ifdef NOV_15 
    xHitPos = (nstrips()/2.0 - cHitDet )*dutpitch() ;
#endif
    hitPositions.push_back( xHitPos );
  }
  return hitPositions;
}

tbeam::FeIFourEvent* EventAlignment::sortFei4Hits(tbeam::FeIFourEvent* pFei4Ev )
{

  int debug = 0; 
  tbeam::FeIFourEvent* pFei4Ev_sorted = new tbeam::FeIFourEvent();
  pFei4Ev_sorted->euEvt = pFei4Ev->euEvt;
  
  // column and row pitch of FeI4 plane      
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 
  
  // first sort Fei4 hits by their time of arrival 
  std::vector<int> tHit ; 
  tHit.clear(); 
  std::vector<double> dHit;
  dHit.clear(); 
  for(unsigned int i = 0 ; i < pFei4Ev->lv1->size() ; i++ )
  {
    double xval = x0_FEI4 + (pFei4Ev->col->at(i)-1)*(cColPitch_fei4*1e-3); 
    double yval = y0_FEI4 + (pFei4Ev->row->at(i)-1 )*(cRowPitch_fei4*1e-3); 
    double distance = std::sqrt(std::pow(xval,2.0) + std::pow(yval,2.0));  
    tHit.push_back( pFei4Ev->lv1->at(i) );
    dHit.push_back( distance );
  }
  // sort and remove duplicates 
  std::sort (tHit.begin(), tHit.end());  
  tHit.erase( std::unique( tHit.begin(), tHit.end() ), tHit.end() );


  std::sort (dHit.begin(), dHit.end());  
  dHit.erase( std::unique( dHit.begin(), dHit.end() ), dHit.end() );

  // now re-arrange hits in FeI4 event so that they are sorted by time of arrival
  // and for events with the same time of arrival by distance from the center of the FeI4 plane
  std::vector<int> indices; indices.clear(); 
  for(unsigned int i = 0 ; i < dHit.size() ; i++ )
  {
    for( unsigned int j =0 ; j < pFei4Ev->nPixHits ; j++)
    {
      double xval = x0_FEI4 + (pFei4Ev->col->at(j)-1)*(cColPitch_fei4*1e-3); 
      double yval = y0_FEI4 + (pFei4Ev->row->at(j)-1 )*(cRowPitch_fei4*1e-3); 
      double distance = std::sqrt(std::pow(xval,2.0) + std::pow(yval,2.0)); 
      if ( distance == dHit[i] )
      {
        indices.push_back(j);
      }
    }
  }

  for(unsigned int i = 0; i < tHit.size(); i++ )
  {
    for( unsigned j = 0 ; j < indices.size(); j++ )
    {
      int index = indices[j];
      if( pFei4Ev->lv1->at(index) == tHit[i] )
      {
        pFei4Ev_sorted->lv1->push_back(  pFei4Ev->lv1->at(index) );
        pFei4Ev_sorted->row->push_back(  pFei4Ev->row->at(index) );
        pFei4Ev_sorted->col->push_back(  pFei4Ev->col->at(index) );
        pFei4Ev_sorted->tot->push_back(  pFei4Ev->tot->at(index) );
      }
    }
  }

  // set number of pixels hit
  pFei4Ev_sorted->nPixHits = pFei4Ev_sorted->lv1->size() ;
  
  if( debug == 1 )
  {
    for(unsigned int j = 0 ; j < pFei4Ev_sorted->nPixHits ; j++ )
    {
      double xval = x0_FEI4 + (pFei4Ev_sorted->col->at(j)-1)*(cColPitch_fei4*1e-3); 
      double yval = y0_FEI4 + (pFei4Ev_sorted->row->at(j)-1 )*(cRowPitch_fei4*1e-3); 
      double distance = std::sqrt(std::pow(xval,2.0) + std::pow(yval,2.0)); 
      std::cout << "Hit time = " << pFei4Ev_sorted->lv1->at(j) << "\t" ;
      std::cout << "(" << pFei4Ev_sorted->row->at(j)-1<< "," << pFei4Ev_sorted->col->at(j)-1 << ") Original Coordinates = (" << xval << "," << yval << ")" ;
      std::cout << " ToT = " << (int)(pFei4Ev_sorted->tot->at(j)) ;
      std::cout << " Distance from center of FeI4 plane = " << distance << " ." << std::endl; 
      
     }
  }
  return pFei4Ev_sorted; 
} 
void EventAlignment::clusterFeI4Hits(tbeam::FeIFourEvent* pFei4Ev , std::vector<double>& xPosClusters,std::vector<double>& yPosClusters, std::vector<double>& sizeClusters, std::vector<double>&  totClusters, std::vector<double>&  lv1Clusters, int minWidthCluster )
{
  int debug = 0 ;
  //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 
  // prepare Fei4 "super" columns here 
  int cNumRows_perSuperPixel = minWidthCluster; 
  int cNumColumns_perSuperPixel = minWidthCluster; 
  
  std::vector<int> clusterLV1; 
  std::vector<int> clusterHits, clusterRows , clusterColumns; 
  std::vector<double> clusterCenterX , clusterCenterY;
  std::vector<int> clusterToT; 
  clusterHits.clear(); 
  clusterRows.clear(); clusterColumns.clear();
  clusterCenterX.clear(); clusterCenterY.clear();
  clusterToT.clear();
  clusterLV1.clear();

  int minRow = (int)(pFei4Ev->row->at(0)-1);
  int maxRow = (int)(pFei4Ev->row->at(0)-1);
  int minColumn = (int)(pFei4Ev->col->at(0)-1);
  int maxColumn = (int)(pFei4Ev->col->at(0)-1);
  bool clusterAdded = false ; 
      
  xPosClusters.clear(); yPosClusters.clear(); sizeClusters.clear(); totClusters.clear(); lv1Clusters.clear();
  int rowOffset = std::min( std::fabs(maxRow  - pFei4Ev->row->at(0)-1) , std::fabs(minRow  - pFei4Ev->row->at(0)-1 ));
  int colOffset = std::min( std::fabs(minColumn  - pFei4Ev->col->at(0)-1) , std::fabs(maxColumn  - pFei4Ev->col->at(0)-1 ));
  
  double xCenterCluster, yCenterCluster;
  double distance=0.;
  double maxDistance = std::sqrt( std::pow(minWidthCluster*cColPitch_fei4*1e-3,2.0) + std::pow(minWidthCluster*cRowPitch_fei4*1e-3,2.0) );
  std::vector<int> indices; indices.clear(); 
  do
  {
    int startIndex = 0; 
    clusterHits.clear(); 
    clusterRows.clear(); clusterColumns.clear();
    clusterCenterX.clear(); clusterCenterY.clear();
    clusterToT.clear();
    clusterLV1.clear(); 

    double xval = x0_FEI4 + (pFei4Ev->col->at(startIndex)-1)*(cColPitch_fei4*1e-3); 
    double yval = y0_FEI4 + (pFei4Ev->row->at(startIndex)-1 )*(cRowPitch_fei4*1e-3); 
  
    clusterCenterX.push_back( xval ); clusterCenterY.push_back( yval ); 
    clusterHits.push_back(1);
    clusterRows.push_back(pFei4Ev->row->at(startIndex)-1);
    clusterColumns.push_back(pFei4Ev->col->at(startIndex)-1);
    clusterToT.push_back( pFei4Ev->tot->at(startIndex) );
    clusterLV1.push_back( pFei4Ev->lv1->at(startIndex) );
    maxRow = (int)(pFei4Ev->row->at(startIndex)-1);
    minRow =(int)(pFei4Ev->row->at(startIndex)-1);
    minColumn = (int)(pFei4Ev->col->at(startIndex)-1);
    maxColumn = (int)(pFei4Ev->col->at(startIndex)-1);
    indices.clear(); //cFei4Ev->lv1->clear(); cFei4Ev->row->clear(); cFei4Ev->col->clear(); cFei4Ev->tot->clear();
    tbeam::FeIFourEvent* cFei4Ev = new tbeam::FeIFourEvent();
    cFei4Ev->lv1->clear();
    cFei4Ev->row->clear();
    cFei4Ev->col->clear();
    cFei4Ev->tot->clear();
    for (unsigned int i = startIndex+1 ; i < pFei4Ev->nPixHits; i++) 
    {
      int cToT = (int)(pFei4Ev->tot->at(i));
      int cColumn = (int)(pFei4Ev->col->at(i)-1);
      int cRow  = (int)(pFei4Ev->row->at(i)-1); 
      xval = x0_FEI4 + (cColumn)*(cColPitch_fei4*1e-3); 
      yval = y0_FEI4 + (cRow)*(cRowPitch_fei4*1e-3); 
          
      rowOffset = std::min( std::fabs(maxRow  - cRow) , std::fabs(minRow  - cRow ));
      colOffset = std::min( std::fabs(minColumn  - cColumn) , std::fabs(maxColumn  - cColumn ));
      xCenterCluster = std::accumulate(clusterCenterX.begin(), clusterCenterX.end(), 0.0)/clusterCenterX.size();
      yCenterCluster = std::accumulate(clusterCenterY.begin(), clusterCenterY.end(), 0.0)/clusterCenterY.size();
      distance = std::sqrt( std::pow(xCenterCluster-xval,2.0) + std::pow(yCenterCluster-yval,2.0) );
      
      if( (rowOffset <= cNumRows_perSuperPixel && colOffset <= cNumColumns_perSuperPixel) || distance <= maxDistance )
      {
        clusterHits.push_back(1);
        clusterRows.push_back(pFei4Ev->row->at(i)-1);
        clusterColumns.push_back(pFei4Ev->col->at(i)-1);
        clusterCenterX.push_back( xval ); clusterCenterY.push_back( yval ); 
        clusterToT.push_back( cToT );
        clusterLV1.push_back( pFei4Ev->lv1->at(i) );

        // find min/max in cluster rows
        auto minElementRow = *std::min_element(std::begin(clusterRows), std::end(clusterRows));
        auto maxElementRow = *std::max_element(std::begin(clusterRows), std::end(clusterRows));
        minRow = (int)minElementRow;
        maxRow = (int)maxElementRow;
        // find min/max in cluster columns
        auto minElementColumn = *std::min_element(std::begin(clusterColumns), std::end(clusterColumns));
        auto maxElementColumn = *std::max_element(std::begin(clusterColumns), std::end(clusterColumns));
        minColumn = (int)minElementColumn;
        maxColumn = (int)maxElementColumn;
      }
      else
      {
         indices.push_back( i );
         cFei4Ev->lv1->push_back(pFei4Ev->lv1->at(i) );
         cFei4Ev->row->push_back(pFei4Ev->row->at(i) );
         cFei4Ev->col->push_back(pFei4Ev->col->at(i) );
         cFei4Ev->tot->push_back(pFei4Ev->tot->at(i) );
      }
    }
    pFei4Ev = new tbeam::FeIFourEvent();
    pFei4Ev->lv1->clear();
    pFei4Ev->row->clear();
    pFei4Ev->col->clear();
    pFei4Ev->tot->clear();
    for( unsigned int i = 0 ; i < indices.size(); i++ )
    {
      pFei4Ev->lv1->push_back(cFei4Ev->lv1->at(i) );
      pFei4Ev->row->push_back(cFei4Ev->row->at(i) );
      pFei4Ev->col->push_back(cFei4Ev->col->at(i) );
      pFei4Ev->tot->push_back(cFei4Ev->tot->at(i) );
    }  
    pFei4Ev->nPixHits = indices.size();
 

    double meanXtot = 0. ; 
    double meanYtot = 0. ;
    double meanX = 0.0 ;
    double meanY = 0.0 ; 
    double meanLV1 = 0.0; 
    double sumToT = 0;
    double totalToT = 0; 
    int nHits = 0 ; 
    if( debug == 1 )  std::cout << "Cluster found! Contains " <<  clusterHits.size() << " Hits." << std::endl;


    for(unsigned int index = 0 ; index < clusterToT.size(); index++ )
    {
        
      if( debug == 1 ) std::cout << "\tOriginal hits : (" << clusterRows[index]-1<< "," << clusterColumns[index]-1 << ") Original Coordinates : (" << clusterCenterX[index] << "," << clusterCenterY[index] << ")" << std::endl ;
      
      if( clusterToT[index] <= fTot_Max )
      {
        double Q = FeI4_ToT_to_Q( clusterToT[index] ); 
        //Q = std::pow(clusterToT[index],2);

        meanXtot += clusterCenterX[index]*( Q );
        meanYtot += clusterCenterY[index]*( Q );
        meanX += clusterCenterX[index];
        meanY += clusterCenterY[index];
        sumToT += Q; //std::pow(clusterToT[index],2.0);
        totalToT += clusterToT[index]; 
        meanLV1 += (double)(clusterLV1[index]);
        nHits+=1;  
      }
    }

        

    if( sumToT != 0 )
    {
      meanXtot = meanXtot/sumToT;
      meanYtot = meanYtot/sumToT;
      meanX = meanX/nHits;
      meanY = meanY/nHits;
      meanLV1 = meanLV1/nHits;
    }
    
    if( debug == 1 )  std::cout << "\t\tGeometric Mean : Cluster center (X) = " <<  meanX << "\tCluster center (Y) = " <<  meanY <<  std::endl; 
    if( debug == 1 )  std::cout << "\t\tToT Mean : Cluster center (X) = " <<  meanXtot << "\tCluster center (Y) = " <<  meanYtot <<  std::endl; 

    if( sumToT != 0 )
    {
        int diffLV1 ;
        std::vector<int>::iterator minLocator = std::min_element(std::begin(clusterLV1), std::end(clusterLV1));
        int id_min = (int)(std::distance(std::begin(clusterLV1), minLocator));
        std::vector<int>::iterator maxLocator = std::max_element(std::begin(clusterLV1), std::end(clusterLV1));
        int id_max = (int)(std::distance(std::begin(clusterLV1), maxLocator));
        diffLV1 = std::fabs( clusterLV1[id_min] - clusterLV1[id_max]);

        if( diffLV1 == 0 )
        {
          xPosClusters.push_back( meanXtot );
          yPosClusters.push_back( meanYtot );
          sizeClusters.push_back( nHits );
          totClusters.push_back( totalToT );
          lv1Clusters.push_back( meanLV1 );
        }
    }
  }while( indices.size() != 0 ); //&& xPosClusters.size() < 1 );
}

void EventAlignment::removeDuplicateTracks(std::vector<int>& idTks)
{
  idTks.clear();
  // sort and remove duplicate tracks
  std::vector<int64_t>* xTracks = new std::vector<int64_t>();  // in um 
  std::vector<int64_t>* yTracks = new std::vector<int64_t>();  // in um
  
  //std::cout << "Number of tracks  = " << telEv()->xPos->size() << "\t" << std::endl;
  for( unsigned int i = 0 ; i < telEv()->xPos->size() ; i++)
  {
    int64_t xTk = (int64_t)(telEv()->xPos->at(i)*1e-3*1e6);
    int64_t yTk = (int64_t)(telEv()->yPos->at(i)*1e-3*1e6);
    if( (find (xTracks->begin(), xTracks->end(), xTk) == xTracks->end()) && (find (yTracks->begin(), yTracks->end(), yTk) == yTracks->end()) ) 
    {
       xTracks->push_back( xTk );
       yTracks->push_back( yTk );
       idTks.push_back( i ); 
    }
  }
}
void EventAlignment::FindDetFeI4CorrelationWindow(int nEventsToUse )
{
  int groupID = 0;
  int evOffset = 0; 

  double rotation = 0;
 //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;

  for (Long64_t jentry=(groupID*nEventsToUse) ; jentry<(groupID+1)*nEventsToUse ;jentry++) 
  {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      //if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      //if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;

      
      setDetChannelVectors();
      const auto& d0c0 = *det0C0();
      const auto& d0c1 = *det0C1();
      const auto& d1c0 = *det1C0();
      const auto& d1c1 = *det1C1();   

      std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
      std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
      std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
      std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
      std::vector<double> x_Stubs = getStubPositions();

      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;
      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;

      // remove track duplicates 
      std::vector<int> idTks;
      removeDuplicateTracks(idTks);
      
      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
      
      for( unsigned int i = 0 ; i < x_Clusters_Det0.size() ; i++ )
      {
        
        double xDet0 = x_Clusters_Det0[i] - alParsLocal_.d0Offset();

        std::vector<double> HitResidual; HitResidual.clear(); 
        std::vector<int> HitIndex; HitIndex.clear();
        for( int j = 0 ; j < xPos_FeI4Clusters.size()  ; j++ )
        {
          if( std::fabs(lv1_FeI4Clusters[j] - fL1_mean) > fL1_sigma ) continue;
          double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[j])*1e-3;
          //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 

          double xval = xPos_FeI4Clusters[j];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
          double yval = yPos_FeI4Clusters[j]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          yval = -1*yval; 
          HitResidual.push_back( std::fabs( yval - xDet0)  );
          HitIndex.push_back( j );
        }
        if( HitResidual.size() == 0 ) continue; 

        std::vector<double>::iterator minLocator = std::min_element(std::begin(HitResidual), std::end(HitResidual));
        int id = (int)(std::distance(std::begin(HitResidual), minLocator));
        id = HitIndex[id];

        double xval = xPos_FeI4Clusters[id];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[id]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
        yval = -1*yval; 
        double residual =  yval - xDet0 ;
        double residual_abs = std::fabs(residual);
        
        hist_->fillHist2D("FeI4", "x0_Correlation_Search", yval , xDet0 )  ; 
        hist_->fillHist1D("FeI4", "x0_Residual_Search", residual ) ; 
      }
  }

  TH2D* h = dynamic_cast<TH2D*>(hist_->Get2DHistoByName("FeI4","x0_Correlation_Search"));
  h->GetXaxis()->SetRangeUser(-10,10);
  h->GetYaxis()->SetRangeUser(-10,10);

  TF1* f_Correlation = new TF1("f_Correlation", Utility::correlationFunction , -10.0 , 10.0 , 2 );
  std::cout << "Finding correlation between FeI4 clusters and DUT clusters using the first " << nEventsToUse << " events\n.";
  f_Correlation->SetLineStyle(2);
  f_Correlation->SetParameter(1,1);
  h->Fit(f_Correlation);
  fFeI4_Correlation_y0 = f_Correlation->GetParameter(0);
  fFeI4_Correlation_m = f_Correlation->GetParameter(1);
  fFeI4_Correlation_y0_unc = f_Correlation->GetParError(0);
  fFeI4_Correlation_m_unc = f_Correlation->GetParError(1);
}
void EventAlignment::CorrelateDetFeI4(int groupID, int evOffset )
{
  double rotation = 0;
 //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;

  for (Long64_t jentry=(groupID*fNumberEvents_perGroup) ; jentry<(groupID+1)*fNumberEvents_perGroup ;jentry++) 
  {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      //if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      //if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;

      
      setDetChannelVectors();
      const auto& d0c0 = *det0C0();
      const auto& d0c1 = *det0C1();
      const auto& d1c0 = *det1C0();
      const auto& d1c1 = *det1C1();   

      std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
      std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
      std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
      std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
      std::vector<double> x_Stubs = getStubPositions();


    
      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;
      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;

      // remove track duplicates 
      std::vector<int> idTks;
      removeDuplicateTracks(idTks);
      
      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
      
      // "noise" DUT hits w.r.t. telescope tracks 
      if( idTks.size() == 0 )
      {
        hist_->fillHist1D("det0", "nClusters_noTrks", x_Clusters_Det0.size() );
        hist_->fillHist1D("det1", "nClusters_noTrks", x_Clusters_Det1.size() );
        for( unsigned int i = 0 ; i < std::max(x_Clusters_Det0.size(), x_Clusters_Det1.size()) ; i++ )
        {
          if( i < x_Clusters_Det0.size()) hist_->fillHist1D("det0", "Clusters_noTrks", x_Clusters_Det0[i] );
          if( i < x_Clusters_Det1.size()) hist_->fillHist1D("det1", "Clusters_noTrks", x_Clusters_Det1[i] );
        }
      }
      else
      {
        hist_->fillHist1D("det0", "nClusters_wTrks", x_Clusters_Det0.size() );
        hist_->fillHist1D("det1", "nClusters_wTrks", x_Clusters_Det1.size() );
        for( unsigned int i = 0 ; i < std::max(x_Clusters_Det0.size(), x_Clusters_Det1.size()) ; i++ )
        {
          if( i < x_Clusters_Det0.size()) hist_->fillHist1D("det0", "Clusters_wTrks", x_Clusters_Det0[i] );
          if( i < x_Clusters_Det1.size()) hist_->fillHist1D("det1", "Clusters_wTrks", x_Clusters_Det1[i] );
        } 
      }
      if( idTks.size() == 0 ) continue; 


      // "noise" DUT hits w.r.t. FeI4 
      if( xPos_FeI4Clusters.size() == 0 ) 
      {
        hist_->fillHist1D("det0", "nClusters_noFeI4", x_Clusters_Det0.size() );
        hist_->fillHist1D("det1", "nClusters_noFeI4", x_Clusters_Det1.size() );
        for( unsigned int i = 0 ; i < std::max(x_Clusters_Det0.size(), x_Clusters_Det1.size()) ; i++ )
        {
          if( i < x_Clusters_Det0.size()) hist_->fillHist1D("det0", "Clusters_noFeI4", x_Clusters_Det0[i] );
          if( i < x_Clusters_Det1.size()) hist_->fillHist1D("det1", "Clusters_noFeI4", x_Clusters_Det1[i] );
        }
      }
      else
      {
        hist_->fillHist1D("det0", "nClusters_wFeI4", x_Clusters_Det0.size() );
        hist_->fillHist1D("det1", "nClusters_wFeI4", x_Clusters_Det1.size() );
        for( unsigned int i = 0 ; i < std::max(x_Clusters_Det0.size(), x_Clusters_Det1.size()) ; i++ )
        {
          if( i < x_Clusters_Det0.size()) hist_->fillHist1D("det0", "Clusters_wFeI4", x_Clusters_Det0[i] );
          if( i < x_Clusters_Det1.size()) hist_->fillHist1D("det1", "Clusters_wFeI4", x_Clusters_Det1[i] );
        }
      }

      hist_->fillHist1D("TrackMatch", "nClusters_Det0", x_Clusters_Det0.size() );
      hist_->fillHist1D("TrackMatch", "nClusters_Det1", x_Clusters_Det1.size() );
      for( unsigned int i = 0 ; i < x_Clusters_Det0.size() ; i++ )
      {
        
        double xDet0 = x_Clusters_Det0[i] - alParsLocal_.d0Offset();

        std::vector<double> HitResidual; HitResidual.clear(); 
        std::vector<int> HitIndex; HitIndex.clear();
        for( int j = 0 ; j < xPos_FeI4Clusters.size()  ; j++ )
        {
          if( std::fabs(lv1_FeI4Clusters[j] - fL1_mean) > fL1_sigma ) continue;
          double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[j])*1e-3;
          //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 

          double xval = xPos_FeI4Clusters[j];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
          double yval = yPos_FeI4Clusters[j]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          yval = -1*yval; 
          double residual; 
          residual = yval - xDet0;
            
          HitResidual.push_back( std::fabs( residual)  );
          HitIndex.push_back( j );
        }
        if( HitResidual.size() == 0 ) continue; 
        hist_->fillHist1D("TrackMatch", "nClusters_Fei4", HitResidual.size() );
      
        std::vector<double>::iterator minLocator = std::min_element(std::begin(HitResidual), std::end(HitResidual));
        int id = (int)(std::distance(std::begin(HitResidual), minLocator));
        id = HitIndex[id];

        double xval = xPos_FeI4Clusters[id];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[id]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
        yval = -1*yval; 
        double residual; 
        residual =  yval - xDet0 ;
        double residual_abs = std::fabs(residual);
        
        // if the closest cluster in the event is 
        // if( residual >=  5.0*(cRowPitch_fei4/std::sqrt(12)) )  continue; 

        hist_->fillHist1D("FeI4", "lv1_M", (double)(lv1_FeI4Clusters[id]) );
        hist_->fillHist2D("FeI4", "x0_Correlation", yval , xDet0 )  ; 
        //std::cout << "FeI4 cluster at " << yval << " track impact at " << xDet0 << std::endl;
        hist_->fillHist1D("FeI4", "x0_Residual", residual ) ; 


      }
  }
} 
void EventAlignment::CorrelateFeI4Tracks(int groupID, int evOffset )
{
  double rotation = 0;
 //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;

  for (Long64_t jentry=(groupID*fNumberEvents_perGroup) ; jentry<(groupID+1)*fNumberEvents_perGroup ;jentry++) 
  {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;
      


      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;
      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;

      // remove track duplicates 
      std::vector<int> idTks;
      removeDuplicateTracks(idTks);
      if( idTks.size() == 0 ) continue; 
      
      for( int j = 0 ; j < idTks.size() ; j++ )
      {
          unsigned int itk = idTks[j];
          double dydz = telEv()->dydz->at(itk); 
          double dxdz = telEv()->dxdz->at(itk); 
          double chi2 = telEv()->chi2->at(itk); 
          double nDOF = telEv()->ndof->at(itk); 

          hist_->fillHist1D("TelescopeAnalysis","tk_dYdZ", dydz  );
          hist_->fillHist1D("TelescopeAnalysis","tk_dXdZ", dxdz );
          hist_->fillHist1D("TelescopeAnalysis","hTkChi2", chi2 );
          hist_->fillHist1D("TelescopeAnalysis","hTkNdof", nDOF );
          hist_->fillHist1D("TelescopeAnalysis","hTkChi2_Ndof", chi2/nDOF );
      
          if( idTks.size() == 1 ) 
          {
            hist_->fillHist1D("TelescopeAnalysis","hTkChi2_SingleTks", chi2/nDOF );
          }
          else
          {
            hist_->fillHist1D("TelescopeAnalysis","hTkChi2_MultiTks", chi2/nDOF );
          }
      }  

      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );

      int nClusters = 0 ;
      hist_->fillHist1D("TelescopeAnalysis","nTracks", idTks.size()  );
      for( int i = 0 ; i < xPos_FeI4Clusters.size() ; i++ )
      {
        if( std::fabs(lv1_FeI4Clusters[i] - fL1_mean) > fL1_sigma ) continue;
        double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[i])*1e-3;
        //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 

        nClusters++;
        hist_->fillHist1D("TelescopeAnalysis","hTotFei4", tot_FeI4Clusters[i]  );
        hist_->fillHist1D("TelescopeAnalysis","hLv1Fei4", lv1_FeI4Clusters[i]  );
       
        double xval = xPos_FeI4Clusters[i]*std::cos(rotation*pi/180) - yPos_FeI4Clusters[i]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[i]*std::cos(rotation*pi/180) + xPos_FeI4Clusters[i]*std::sin(rotation*pi/180); 
      
        std::vector<int> matchedTrackID;
        matchedTrackID.clear();
        std::vector<float> matchedTrackResidual, matchedTrackResidualX , matchedTrackResidualY , matchedTrackChi2 ;
        matchedTrackResidual.clear(); 
        matchedTrackResidualX.clear();
        matchedTrackResidualY.clear();
        for( int j = 0 ; j < idTks.size() ; j++ )
        {
          unsigned int itk = idTks[j];
          double tkX =  -1*telEv()->xPos->at(itk); //-1*telEv()->xPos->at(itk); 
          double tkY =  telEv()->yPos->at(itk); 
          double chi2 = telEv()->chi2->at(itk); 
          double distance = std::sqrt( std::pow(xval-tkX,2.0) + std::pow(yval - tkY,2.0) ) ;

          matchedTrackID.push_back( itk );
#ifdef NOV_15
            distance = std::sqrt( std::pow(-1*xval-tkY,2.0) + std::pow(yval - tkX,2.0) ) ;
            matchedTrackResidualX.push_back(-1*xval-tkY);
            //distance = std::sqrt( std::pow(-1*xval-tkY,2.0) + std::pow(yval - tkX,2.0) ) ;
            //matchedTrackResidualX.push_back(-1*xval-tkY);
            matchedTrackResidualY.push_back(yval-tkX);
#elif defined(MAY_16) || defined(OCT_16)
          //original
            distance = std::sqrt( std::pow(xval-tkX,2.0) + std::pow(yval - tkY,2.0) ) ;
            matchedTrackResidualX.push_back(xval-tkX);
            matchedTrackResidualY.push_back(yval-tkY);
#endif
          matchedTrackResidual.push_back(distance);
          matchedTrackChi2.push_back(chi2);

        }

        std::vector<float>::iterator minLocator = std::min_element(std::begin(matchedTrackResidual), std::end(matchedTrackResidual));
        int id = (int)(std::distance(std::begin(matchedTrackResidual), minLocator));
        unsigned int itk = (unsigned int)(matchedTrackID[id]);
        hist_->fillHist1D("TelescopeAnalysis","deltaXPos", matchedTrackResidualX[id]  );
        hist_->fillHist1D("TelescopeAnalysis","deltaYPos", matchedTrackResidualY[id]  );
        hist_->fillHist1D("TelescopeAnalysis","deltaDistance", matchedTrackResidual[id]  );
        hist_->fillHist2D("TelescopeAnalysis","tkXPosVsHtXPos",  xval ,  -1*telEv()->xPos->at(itk)   );
        hist_->fillHist2D("TelescopeAnalysis","tkYPosVsHtYPos",  yval ,     telEv()->yPos->at(itk)   );
        
        double x_Residual = xval - ( -1*telEv()->xPos->at(itk) );
        double y_Residual = yval - ( telEv()->yPos->at(itk) ); 
        hist_->fillHist2D("FeI4","Residuals",  x_Residual ,  y_Residual );
        

        matchedTrackID.clear();
        matchedTrackResidual.clear(); 
        matchedTrackResidualX.clear();
        matchedTrackResidualY.clear();
      }
      if( nClusters > 0 ) hist_->fillHist1D("TelescopeAnalysis","nClustersFei4", nClusters );
      
  }
}
int EventAlignment::MatchDetFeI4( int groupID , int prevOffset , double acceptanceWindow )
{
  double rotation = 0;
 //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;
  

  int totalNumClusters = 0; 
  std::vector<int> offsets{prevOffset-1,prevOffset,prevOffset+1};
  std::vector<int> correlated_hits; correlated_hits.clear();
  std::vector<double> fraction_correlated_hits ; fraction_correlated_hits.clear(); 
  for( int index = 0 ; index < offsets.size() ; index++)
  {
    int evOffset = offsets[index];
    int nCorrelatedEntries =0 ;
    int nEntries = 0 ; 
    for (Long64_t jentry=(groupID*fNumberEvents_perGroup + 1) ; jentry<(groupID+1)*fNumberEvents_perGroup + 1 ;jentry++) 
    {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;

      
      setDetChannelVectors();
      const auto& d0c0 = *det0C0();
      const auto& d0c1 = *det0C1();
      const auto& d1c0 = *det1C0();
      const auto& d1c1 = *det1C1();   

      std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
      std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
      std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
      std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
      std::vector<double> x_Stubs = getStubPositions();


      // remove track duplicates 
      std::vector<int> idTks;
      // can also use the built-in duplicate track removal function here 
      removeDuplicateTracks(idTks);
      if( idTks.size() == 0 ) continue; 
    
      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;

      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;
      
      // I group the FeI4 hits into clusters 
      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
      
      for( unsigned int i = 0 ; i < x_Clusters_Det0.size() ; i++ )
      {
        std::vector<double> HitResidual; HitResidual.clear(); 
        for( int j = 0 ; j < xPos_FeI4Clusters.size()  ; j++ )
        {
          if( std::fabs(lv1_FeI4Clusters[j] - fL1_mean) > fL1_sigma ) continue;
          //double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[j])*1e-3;
          //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 


          double yval = yPos_FeI4Clusters[j]*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          yval = -1*yval;
          HitResidual.push_back( std::fabs(yval - x_Clusters_Det0[i]) );
          
        }
        if( HitResidual.size() == 0 ) continue ; 

        std::vector<double>::iterator minLocator = std::min_element(std::begin(HitResidual), std::end(HitResidual));
        int id = (int)(std::distance(std::begin(HitResidual), minLocator));
        double yval = yPos_FeI4Clusters[id];
        yval = -1*yval;
        double residual_abs = std::fabs(yval - (x_Clusters_Det0[i]-alParsLocal_.d0Offset() ) ) ;
        double residual =  x_Clusters_Det0[i] - yval;
        
        if( std::fabs(residual-alParsLocal_.d0Offset()) > acceptanceWindow )
        { 
          hist_->fillHist1D("FeI4", "x0_Residual_Evnt", (int)jentry , 1 );
        }
        else
        {  
          hist_->fillHist1D("FeI4", "x0_Residual_Evnt", (int)jentry , 0 );
        }

        if( std::fabs(residual-alParsLocal_.d0Offset()) <= acceptanceWindow )
        { 
          nCorrelatedEntries++;
        }
        nEntries++; 
      }
      //if( nGoodEvents%10000 == 0 ) std::cout << "Events processed = " << nGoodEvents << std::endl; 
      nGoodEvents++;
    }
    //std::cout << "\tFor offset = " << offsets[index] << " , number of correlated clusters = " << nCorrelatedEntries << std::endl;
    //<< " which is " << fraction_correlated_hits[fraction_correlated_hits.size()-1] << " percent of the FeI4 clusters." << std::endl; 
    correlated_hits.push_back( nCorrelatedEntries );
  }
  std::vector<int>::iterator maxLocator = std::max_element(std::begin(correlated_hits), std::end(correlated_hits));
  int id = (int)(std::distance(std::begin(correlated_hits), maxLocator));
  int optimalOffset = offsets[id];
  return (int)(optimalOffset);
}
int EventAlignment::MatchDetFeI4_v2( int groupID , int prevOffset , double pCorrelation , double minCorrelation )
{
  double rotation = 0;
 //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;
  

  int totalNumClusters = 0; 
  std::vector<int> offsets{prevOffset-3,prevOffset-2,prevOffset-1,prevOffset,prevOffset+1,prevOffset+2,prevOffset+3};
  std::vector<int> correlated_hits; correlated_hits.clear();
  std::vector<double> fraction_correlated_hits ; fraction_correlated_hits.clear(); 
  for( int index = 0 ; index < offsets.size() ; index++)
  {
    int evOffset = offsets[index];
    int nCorrelatedEntries =0 ;
    int nEntries = 0 ; 
    for (Long64_t jentry=(groupID*fNumberEvents_perGroup + 1) ; jentry<(groupID+1)*fNumberEvents_perGroup + 1 ;jentry++) 
    {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;

      
      setDetChannelVectors();
      const auto& d0c0 = *det0C0();
      const auto& d0c1 = *det0C1();
      const auto& d1c0 = *det1C0();
      const auto& d1c1 = *det1C1();   

      std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
      std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
      std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
      std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
      std::vector<double> x_Stubs = getStubPositions();


      // remove track duplicates 
      std::vector<int> idTks;
      // can also use the built-in duplicate track removal function here 
      removeDuplicateTracks(idTks);
      if( idTks.size() == 0 ) continue; 
    
      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;

      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;
      
      // I group the FeI4 hits into clusters 
      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
      
      for( unsigned int i = 0 ; i < x_Clusters_Det0.size() ; i++ )
      {
        std::vector<double> HitResidual; HitResidual.clear(); 
        double xDet0 = x_Clusters_Det0[i] - alParsLocal_.d0Offset();
        for( int j = 0 ; j < xPos_FeI4Clusters.size()  ; j++ )
        {
          if( std::fabs(lv1_FeI4Clusters[j] - fL1_mean) > fL1_sigma ) continue;
          //double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[j])*1e-3;
          //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 
          
          double xval = xPos_FeI4Clusters[j];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
          double yval = yPos_FeI4Clusters[j]*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          yval = -1*yval;
          double expected_HitPosition = fFeI4_Correlation_y0 + fFeI4_Correlation_m*yval; 
          double residual = expected_HitPosition - xDet0 ; 
          double residual_abs = std::fabs(residual);
          HitResidual.push_back( residual_abs );
          
        }
        if( HitResidual.size() == 0 ) continue ; 

        std::vector<double>::iterator minLocator = std::min_element(std::begin(HitResidual), std::end(HitResidual));
        int id = (int)(std::distance(std::begin(HitResidual), minLocator));
        double xval = xPos_FeI4Clusters[id];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[id];
        yval = -1*yval; 
        double expected_HitPosition = fFeI4_Correlation_y0 + fFeI4_Correlation_m*yval; 
        double minExpected_HitPosition = (fFeI4_Correlation_y0-pCorrelation*fFeI4_Correlation_y0_unc) + (fFeI4_Correlation_m+pCorrelation*fFeI4_Correlation_m_unc)*yval; 
        double maxExpected_HitPosition = (fFeI4_Correlation_y0+pCorrelation*fFeI4_Correlation_y0_unc) + (fFeI4_Correlation_m-pCorrelation*fFeI4_Correlation_m_unc)*yval; 
        //double minExpected_HitPosition = (fFeI4_Correlation_y0-k*fFeI4_Correlation_y0_unc) + (fFeI4_Correlation_m-k*fFeI4_Correlation_m_unc)*yval; 
        //double maxExpected_HitPosition = (fFeI4_Correlation_y0+k*fFeI4_Correlation_y0_unc) + (fFeI4_Correlation_m+k*fFeI4_Correlation_m_unc)*yval; 
        double acceptanceWindow = std::fabs( maxExpected_HitPosition - minExpected_HitPosition);
        //std::cout << "DUT hit = " << xDet0 << " extrapolated from FeI4 hit " << expected_HitPosition << " -acceptance window = " << acceptanceWindow << std::endl;
          
        double residual = expected_HitPosition - xDet0 ; 
        double residual_abs = std::fabs(residual);
        
        if( residual_abs > acceptanceWindow )
        { 
          hist_->fillHist1D("FeI4", "x0_Residual_Evnt", (int)jentry , 1 );
        }
        else
        {  
          hist_->fillHist1D("FeI4", "x0_Residual_Evnt", (int)jentry , 0 );
        }

        if( residual_abs <= acceptanceWindow )
        { 
          nCorrelatedEntries++;
        }
        nEntries++; 
      }
      //if( nGoodEvents%10000 == 0 ) std::cout << "Events processed = " << nGoodEvents << std::endl; 
      nGoodEvents++;
    }
    if( nEntries == 0 ) continue ; 
    //std::cout << "\tFor offset = " << offsets[index] << " , number of correlated clusters = " << nCorrelatedEntries << std::endl;
    //std::cout  << " which is " << fraction_correlated_hits[fraction_correlated_hits.size()-1] << " percent of the FeI4 clusters." << std::endl; 
    correlated_hits.push_back( nCorrelatedEntries );
    fraction_correlated_hits.push_back( (double)(nCorrelatedEntries)/nEntries );
  }
  if( fraction_correlated_hits.size() > 0 )
  {
    std::vector<int>::iterator maxLocator = std::max_element(std::begin(correlated_hits), std::end(correlated_hits));
    int id = (int)(std::distance(std::begin(correlated_hits), maxLocator));
    int optimalOffset = offsets[id];
    double optimalCorrelation = fraction_correlated_hits[id];
    //std::cout << "For group = " << groupID << " offset =  " << offsets[index] << " and " << 100.0*optimalCorrelation << " percent of the FeI4 clusters correlated with DUT hits." << std::endl; 
    hist_->fillHistProfile("FeI4", "EventAlignment_Correlation", groupID, optimalCorrelation*100 ) ;
    hist_->fillHistProfile("FeI4", "EventAlignment_Offset", groupID, optimalOffset ) ; 
    sprintf(buffer,"For group %d : offset @ maximum correlation found to be %d , maximum correlation %.3f percent\n" , groupID , optimalOffset , 100*optimalCorrelation) ;
    if( optimalCorrelation == 0 || optimalCorrelation <  minCorrelation )
    { 
      sprintf(buffer,"\t!!! For group %d : maximum correlation %.3f percent !!!\n" , groupID  , 100*optimalCorrelation) ;
      std::cout << buffer ; 
      return -100;
    }
    else
    { 
      return (int)(optimalOffset);
    }
  }
  else
  {
    sprintf(buffer,"\t!!! For group %d : ZERO good events found !!!\n" , groupID  ) ;
    std::cout << buffer ; 
    return -100;
  } 
}
double EventAlignment::returnStrip( double xHit, double pitch )
{

  if( pitch == 0 ) pitch = dutpitch();
  double iStrip = xHit/pitch - nstrips()/2.0; 
#ifdef NOV_15
     iStrip = -1*xHit/pitch + nstrips()/2.0; 
#endif
  return iStrip;
}
double EventAlignment::returnHit( double iStrip, double pitch )
{

  if( pitch == 0 ) pitch = dutpitch();
  double xHit;
#if defined(MAY_16) || defined(OCT_16)
  xHit = (iStrip - nstrips()/2.0)*pitch; 
#elif NOV_15 
  xHit = (nstrips()/2.0 - iStrip)*pitch;
#endif
  return xHit;
}
void EventAlignment::clusterCbcHits(std::vector<double> xHits, std::vector<double>& xPosClusters, std::vector<double>& sizeClusters , double maxClusterWidth )
{
    int debug = 0;
    if( debug ) for( unsigned int i = 0 ; i < xHits.size() ; i++ ) std::cout << returnStrip(xHits[i]) << "\t";
    if( debug ) std::cout << std::endl;

    xPosClusters.clear();  sizeClusters.clear();
    int startIndex = 0 ; 
    std::vector<double> cluster; cluster.clear(); 
    cluster.push_back( returnStrip(xHits[startIndex]) );
    //if(debug) std::cout << "\t\t hit " << xHits[startIndex] << " strip = " << returnStrip(xHits[startIndex]) << "\n";
    startIndex++;
    while( startIndex < xHits.size() )
    {
        //if( debug) std::cout << "\t\t hit " << xHits[startIndex] << " strip = " << returnStrip(xHits[startIndex]) << "\n";
        double cDifference = std::fabs( returnStrip(xHits[startIndex]) - cluster[cluster.size()-1] )  ;
        //std::cout << "#" << returnStrip(xHits[startIndex]) << " - " << cluster[cluster.size()-1] << " : " << cDifference <<  std::endl ;
        //char buffer[256];
        //sprintf(buffer,"Difference = %.3f\n" , cDifference);
        //if( debug ) std::cout << buffer;
        if( cDifference <= 1 )
        {
            cluster.push_back(returnStrip(xHits[startIndex]));
        }
        else
        {
            //if( debug ) std::cout << "Difference == " << cDifference << " adding...." << returnStrip(xHits[startIndex]) << " to cluster " << std::endl;
            double cCenter = 0;
            // if odd
            if( cluster.size()%2 != 0 ) cCenter = cluster[ (cluster.size()/2) ];
            else cCenter = cluster[0] + (cluster[cluster.size()-1] - cluster[0])/2.0;
            int cCenter_hStrips = 2*cCenter ;
            
            if( cluster.size() <= maxClusterWidth )
            {
              xPosClusters.push_back( cCenter_hStrips);
              sizeClusters.push_back( cluster.size() );
            
              if( debug )
              {
                std::cout << "\t\t\t Found a cluster of size " << cluster.size() << " with center " << cCenter << "[" << cCenter_hStrips << " half strips] : "  ;
                for( unsigned int i = 0 ; i < cluster.size(); i++ )
                {
                    std::cout <<  cluster[i] << " " ;
                }
                std::cout << std::endl;
              }
            }

            cluster.clear();
            //if( debug ) std::cout << "\t\t\t adding " <<  returnStrip(xHits[startIndex]) << " to the cluster." << std::endl;
            cluster.push_back( returnStrip(xHits[startIndex]) );
            //if( debug ) std::cout << "\t\t\t checking .. just added " << cluster[ cluster.size() - 1]  << " to the cluster." << std::endl;
            
        }
        startIndex++;
    }
    if( cluster.size() > 0 )
    {
        double cCenter = 0;
        // if odd
        if( cluster.size()%2 != 0 ) cCenter = cluster[ (cluster.size()/2) ];
        else cCenter = cluster[0] + (cluster[cluster.size()-1] - cluster[0])/2.0;
        int cCenter_hStrips = 2*cCenter ;
        
        if( cluster.size() <= maxClusterWidth )
        {
          xPosClusters.push_back( cCenter_hStrips);
          sizeClusters.push_back( cluster.size() );
        
          if( debug )
          {
            std::cout << "\t\t\t Found a cluster of size " << cluster.size() << " with center " << cCenter << "[" << cCenter_hStrips << " half strips] : "  ;
            for( unsigned int i = 0 ; i < cluster.size(); i++ )
            {
                std::cout <<  cluster[i] << " " ;
            }
            std::cout << std::endl;
          }
        }

        // if( debug )
        // {
        //   std::cout << "\t\t\t Found a cluster of size " << cluster.size() << " with center " << cCenter << "[" << cCenter_hStrips << " half strips] : "  ;
        //   for( unsigned int i = 0 ; i < cluster.size(); i++ )
        //   {
        //       std::cout <<  cluster[i] << " " ;
        //   }
        //   std::cout << std::endl;
        // }

    }
}
uint8_t EventAlignment::findStubs(std::vector<FoundStubInfo>& xStubs_Seeds , std::vector<FoundStubInfo>& xStubs_Matched , std::vector<double> xHits_Det0 , std::vector<double> xHits_Det1 , int correlationWindow_Strips , bool halfStripRes )
{
    int debug = 0 ; 
    std::vector<double> xPosClusters_Det0; std::vector<double> sizeClusters_Det0;
    std::vector<double> xPosClusters_Det1; std::vector<double> sizeClusters_Det1;
    
    if( debug ) std::cout << "Det0 : \n";
    clusterCbcHits(xHits_Det0, xPosClusters_Det0 , sizeClusters_Det0);
    if( debug )  std::cout << "Det1 : \n";
    clusterCbcHits(xHits_Det1, xPosClusters_Det1 , sizeClusters_Det1);
    
    xStubs_Seeds.clear();
    xStubs_Matched.clear();
    
    FoundStub cSeed , cMatched;
    uint8_t stubWord = 0 ;
    uint8_t cStubCbc0 , cStubCbc1; 
    cStubCbc0 = 0 ; cStubCbc1 = 0;
    for( unsigned int j = 0 ; j < xPosClusters_Det1.size() ; j++)
    {
        for( unsigned int i = 0 ; i < xPosClusters_Det0.size(); i++ )
        {
            double cBend_hStrips = ( xPosClusters_Det1[j] - xPosClusters_Det0[i]) ;
            double cBend_fStrips = ( (int)(xPosClusters_Det1[j]/2.0) - (int)(xPosClusters_Det0[i]/2.0))  ;
            bool stubFindingCondition = halfStripRes ? ( std::fabs(cBend_hStrips) <= 2.0*correlationWindow_Strips) : ( std::fabs(cBend_fStrips) <= 1.0*correlationWindow_Strips);
            //if( cBend_hStrips <= 2.0*correlationWindow_Strips)
            //if( cBend_fStrips <= 1.0*correlationWindow_Strips)
            if( stubFindingCondition )
            {
                int chipID = (  xPosClusters_Det0[i] > 254 ) ? 1 : 0 ;
                //double cBend = cBend_hStrips/2.0;
                double cBend = halfStripRes ? cBend_hStrips/2.0 : cBend_fStrips; 
                // if( debug )
                // {
                //   std::cout << "Det0: Found a cluster of size " << sizeClusters_Det0[i] << " with center " << (int)(xPosClusters_Det0[i]/2) << "[" << xPosClusters_Det0[i] << " half strips] " << std::endl ;
                //   std::cout << "Det1: Found a cluster of size " << sizeClusters_Det1[j] << " with center " << (int)(xPosClusters_Det1[j]/2) << "[" << xPosClusters_Det1[j] << " half strips] " ;
                //   std::cout << "\n\t with a bend of " << cBend << "." ;
                //   std::cout << "\t....chip 0x" << std::hex << std::pow(2,chipID) << std::dec <<  "\n";
                // }
                
                cSeed = std::pair<int, int> (xPosClusters_Det1[j] , sizeClusters_Det1[j] );
                cMatched = std::pair<int, int> (xPosClusters_Det0[i] , sizeClusters_Det0[i] );

                if( cSeed.first < 254 ) cStubCbc0 = 1 ; 
                if( cSeed.first > 254 ) cStubCbc1 = 1 ; 

                FoundStubInfo cMatched_full = std::pair<FoundStub, double> (cMatched  , cBend);
                FoundStubInfo cSeed_full = std::pair<FoundStub, double> (cSeed  , cBend);
                
                xStubs_Seeds.push_back(cSeed_full);
                xStubs_Matched.push_back(cMatched_full);
            }
        }
    
    }
    //sortStubs( xStubs_Seeds, xStubs_Matched, correlationWindow_Strips);
    cStubCbc0 = 0 ; cStubCbc1 = 0;
    for( unsigned int i = 0 ; i < xStubs_Seeds.size(); i++ ) 
    {
      cSeed = xStubs_Seeds[i].first ; 
      cMatched = xStubs_Matched[i].first ; 

      if( cSeed.first < 254 ) cStubCbc0 = 1 ; 
      if( cSeed.first > 254 ) cStubCbc1 = 1 ; 
      if( debug ) std::cout << "\tSeed cluster " << (int)(cSeed.first/2.) << " [w=" << cSeed.second << "] matched cluster " << (int)(cMatched.first/2.0) << " [w=" << cMatched.second << "] bend = " << xStubs_Seeds[i].second << " strip(s).\n";
    }
    stubWord = (cStubCbc1 << 1) |  (cStubCbc0 << 0) ;
    return stubWord;
}
/*void EventAlignment::sortStubs(std::vector<FoundStubInfo>& xStubs_Seeds , std::vector<FoundStubInfo>& xStubs_Matched , int correlationWindow_Strips )
{
    std::vector<FoundStubInfo> xStubs_Seeds_Sorted , xStubs_Matched_Sorted; 
    xStubs_Seeds_Sorted.clear(); xStubs_Matched_Sorted.clear(); 

    int debug = 0 ; 
    std::cout << "Now sorting through all the founds stubs.\n";
    std::vector<int> seedCenters, matchedCenters , bends; 
    seedCenters.clear(); matchedCenters.clear(); bends.clear();

    for( unsigned int i = 0 ; i < xStubs_Seeds.size() ; i++ )
    {
      FoundStub seedCluster = xStubs_Seeds[i].first; 
      FoundStub matchedCluster = xStubs_Matched[i].first; 
      int cBend = xStubs_Seeds[i].second;

      int index = i ; 
      for(unsigned int j = i+1 ; j < xStubs_Seeds.size(); j++ )
      {
        FoundStub csCluster = xStubs_Seeds[j].first; 
        FoundStub cmCluster = xStubs_Matched[j].first; 
        
        if( std::fabs(csCluster.first - seedCluster.first) <= correlationWindow_Strips*2 )
        {
          if( csCluster.second <= cBend )
          {
            index = j ;
            cBend = csCluster.second;
          }
        }
      }

      seedCluster = xStubs_Seeds[index].first;
      matchedCluster = xStubs_Matched[index].first; 
      bool seedFound = find (seedCenters.begin(), seedCenters.end(), seedCluster.first) != seedCenters.end();
      //bool matchFound = find (matchedCenters.begin(), matchedCenters.end(), matchedCluster.first) != matchedCenters.end();
      //bool bendFound = find (bends.begin(), bends.end(), cBend) != bends.end();
      if( !seedFound  )
      {
        seedCenters.push_back( seedCluster.first );
        matchedCenters.push_back( seedCluster.first );
        bends.push_back( cBend );

        FoundStub cSeed = std::pair<int, int> (seedCluster.first, seedCluster.second  );
        FoundStub cMatched = std::pair<int, int> (matchedCluster.first, matchedCluster.second );
                
        FoundStubInfo cSeed_full = std::pair<FoundStub, int> (cSeed  , cBend);
        FoundStubInfo cMatched_full = std::pair<FoundStub, int> (cMatched  , cBend);
                
        xStubs_Seeds_Sorted.push_back(cSeed_full);
        xStubs_Matched_Sorted.push_back(cMatched_full);
        if( debug )
        {
          std::cout << "\tSeed cluster " << (int)(cSeed.first/2.) << " [w=" << cSeed.second << "] matched cluster " << (int)(cMatched.first/2.0) << " [w=" << cMatched.second << "] bend = " << cSeed_full.second/2 << " strips.\n";
        }
      }
  }
  xStubs_Seeds =  xStubs_Seeds_Sorted;
  xStubs_Matched = xStubs_Matched_Sorted;
}*/
int EventAlignment::findLongestCluster(std::vector<FoundStub> xStubs )
{
    if( xStubs.size() == 0 ) return 0;
    
    int cLength = xStubs[0].second;
    for( unsigned int i = 1 ; i < xStubs.size() ; i++ )
    {
        if( xStubs[i].second >= cLength ) cLength = xStubs[i].second;
    }
    return cLength;
}
void EventAlignment::dut_RecoStubEfficiency( int pGroupID , int cWindow )
{
  int debug = 0 ;
    
  int cEvCounter, cCounter; 
  cEvCounter = 0 ;
  cCounter = 0 ;
  for (Long64_t jentry=(pGroupID*fNumberEvents_perGroup ) ; jentry< (pGroupID+1)*fNumberEvents_perGroup ;jentry++)
  //for (Long64_t jentry=(pGroupID*fNumberEvents_perGroup ) ; jentry<(pGroupID+1)*fNumberEvents_perGroup + 1 ;jentry++)
  {
    clearEvent();
    Long64_t ientry = analysisTree()->GetEntry(jentry);

    if (ientry < 0) break;
    
    if(!isGoodEvent()) continue ;

    cEvCounter ++; 
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();   

    
    int tdcPhase  = condEv()->tdcPhase; 
    std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
    std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
    std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
    std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
    std::vector<int> w_Clusters_Det0 = getClusterWidths(0);
    std::vector<int> w_Clusters_Det1 = getClusterWidths(1);
    std::vector<double> x_CbcStubs = getCbcStubPositions();
    std::vector<double> x_RecoStubs = getRecoStubPositions();
    std::vector<double> x_RecoSeeds = getRecoStubSeeds();
    std::vector<double> x_RecoMatches = getRecoStubMatches();
    std::vector<double> w_RecoSeeds = getRecoStubSeedWidths();
    std::vector<double> w_RecoMatches = getRecoStubMatchWidths();


    // i can't reconstruct a stub if I only have hits in one layer
    // so ... have to ignore events with hits in only one of the two layers of the module
    if( xHits_Det0.size() == 0 || xHits_Det1.size() == 0) continue;
    
    uint8_t chipStub = (uint8_t)(dutEv()->stubWord);
    uint8_t recoStub = (uint8_t)(dutEv()->stubWordReco);
    
    std::vector<FoundStubInfo> xStubs_Seeds , xStubs_Matched;
    xStubs_Seeds.clear(); xStubs_Matched.clear();
    uint8_t analysisStub = findStubs( xStubs_Seeds ,  xStubs_Matched ,xHits_Det0 , xHits_Det1 , cWindow);
    //int cLength = ( xStubs_Seeds.size() == 0 ) ? 0 : findLongestCluster( xStubs_Seeds );
    int nRecoStubs =  x_RecoStubs.size();
    int nFoundStubs = xStubs_Seeds.size(); 
    if( debug  && nRecoStubs != nFoundStubs )
    //if( (debug || nRecoStubs != nFoundStubs) && recoStub!=0 )
    {
      if( chipStub == recoStub )
        std::cout << "Event " << jentry << ", TDC = " << tdcPhase << " window - " << cWindow << "\n-----\n";
      else
        std::cout << "!!Event " << jentry << ", TDC = " << tdcPhase << " window - " << cWindow << "\n-----\n";
    }
   
    for( unsigned int i = 0 ; i < x_RecoStubs.size(); i++ )
    {
      int iStrip = returnStrip(x_RecoStubs[i]);
      int iStrip_seed = returnStrip( x_RecoSeeds[i] );
      int iStrip_match = returnStrip( x_RecoMatches[i] );
      int cBend = ( iStrip_seed - iStrip_match ); 
      int cWidth_Seed = w_RecoSeeds[i]; 
      int cWidth_Match = w_RecoMatches[i]; 

       
      if( debug && nRecoStubs != nFoundStubs )
      //if( (debug || nRecoStubs != nFoundStubs) && recoStub!=0 )
        std::cout << "\tRECO stub @ seed == " <<  iStrip_seed << " match == " << iStrip_match << " bend == " << std::fabs(iStrip_seed-iStrip_match) << ". [ws= " << cWidth_Seed << " ,wm=" << cWidth_Match << "] \n";
      
      hist_->fillHist1D("det1","recoStubSeedWidth", cWidth_Seed );
      hist_->fillHist1D("det0","recoStubMatchWidth", cWidth_Match );
     
      hist_->fillHist1D("det1","recoStubSeeds" , iStrip );
      hist_->fillHist1D("det1","recoStubBend", cBend);


      hist_->fillHist2D("det1","recoStub_TDC", iStrip_seed , tdcPhase );
    }
    
    for( unsigned int i = 0 ; i < xStubs_Seeds.size() ; i++ )
    {
      FoundStub fCluster =  xStubs_Seeds[i].first; 
      int iStrip = (int)(fCluster.first/2.0); 
      int iStrip_match = (int)((xStubs_Matched[i].first).first/2.0);
      double cBend = (xStubs_Seeds[i].second);
      int cWidth_Seed = fCluster.second; 
      int cWidth_Match = (xStubs_Matched[i].first).second; 

      if( debug  && nRecoStubs != nFoundStubs)
      //if( (debug || nRecoStubs != nFoundStubs) && recoStub!=0 )
         std::cout << "\tFOUND stub @ seed == " <<  iStrip << " match == " << iStrip_match << " bend == " << cBend  << ". [ws= " << cWidth_Seed << " ,wm=" << cWidth_Match << "] \n";
      
      hist_->fillHist1D("det1","anaStubSeeds" , iStrip );
      hist_->fillHist1D("det1","anaStubSeedWidth", cWidth_Seed );
      hist_->fillHist1D("det0","anaStubMatchWidth", cWidth_Match );
      hist_->fillHist1D("det1","anaStubBend", cBend);


      hist_->fillHist2D("det1","anaStub_TDC", iStrip , tdcPhase );
      
    }
    if( debug  && nRecoStubs != nFoundStubs)
    //if( (debug || nRecoStubs != nFoundStubs) && recoStub!=0 )
      std::cout << "----- stub flag cbc = " << +chipStub << " --- stub flag RECO = " << +recoStub << ". Found : " << nFoundStubs << " stubs , reconstructed " << nRecoStubs << " stubs.\n";
    
    for( unsigned int i = 0 ; i < std::max(x_Clusters_Det1.size(),x_Clusters_Det0.size()) ; i++ )
    {
        if( i < x_Clusters_Det0.size() ) 
        {
          hist_->fillHist2D("det0","dutClusterWidth_TDC" , w_Clusters_Det0[i] , tdcPhase );
          hist_->fillHist1D("det0","clusterWidth" , w_Clusters_Det0[i]  );
        }
        if( i < x_Clusters_Det1.size() ) 
        {
          hist_->fillHist2D("det1","dutClusterWidth_TDC" , w_Clusters_Det1[i] , tdcPhase );
          hist_->fillHist1D("det1","clusterWidth" , w_Clusters_Det1[i]  );
        }
    }
  }
}
TF1*  EventAlignment::alignFeI4(TH1D* htmp , double pitch , double xExtra)
{
  xExtra = 100;
  htmp->GetXaxis()->SetRangeUser(-20,20);
  double center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  std::cout << "Center of FeI4 residual [X/Y] = " << center << std::endl;
  double xMin = std::max( htmp->GetXaxis()->GetXmin(), center-xExtra);
  double xMax = center-2*pitch;
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  double max0 = htmp->GetMaximum();

  xMin = center+2*pitch;
  xMax = std::min( htmp->GetXaxis()->GetXmax() , center+xExtra ); 
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  double max1 = htmp->GetMaximum();

  xMin = (max0 > max1) ? center-xExtra : center+ 2*pitch;
  xMax = (max0 > max1) ? center-2*pitch : center+xExtra;
  
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  std::cout << "Fitting background [FeI4 residual] between " << xMin << " and " << xMax << std::endl;
  TF1* f_b = new TF1("f_Background", Utility::normBackgroundFunction , xMin , xMax , 4 );
  f_b->FixParameter(0, 0.0*htmp->GetMaximum()); //f_b->SetParLimits(0, 0 , 100);
  f_b->SetParameter(1, htmp->GetMaximum());
  f_b->FixParameter(2, 0.0 );
  f_b->SetParameter(3, htmp->GetRMS() ); f_b->SetParLimits(3 , 0 , 10.0 );
  f_b->SetLineColor(kMagenta); f_b->SetLineStyle(2);
  htmp->Fit(f_b,"B");

  xMin = center - 2*pitch; 
  xMax = center + 2*pitch; 
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  //std::cout << "Fitting central region [FeI4 residual] between " << xMin << " and " << xMax << std::endl;
  //TF1* f_c = new TF1("f_Central", Utility::signalFunction  , xMin , xMax , 4 );
  //f_c->SetParameter(0, htmp->GetMaximum());
  //f_c->FixParameter(1, pitch); f_c->SetParLimits(1 , 0 , 2.0*pitch );
  //f_c->SetParameter(2, -1*center );
  //f_c->SetParameter(3, 20e-3 ); f_c->SetParLimits(3 , 0 , 1.0 );
  //f_c->SetLineColor(kRed); f_c->SetLineStyle(2);
  //htmp->Fit(f_c,"B");

  xMin = std::max(htmp->GetXaxis()->GetXmin(), center - 2.0) ;//center-xExtra;
  xMax = std::min(htmp->GetXaxis()->GetXmax(), center + 2.0) ;//center+xExtra;
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  
  std::cout << "Fitting FeI4 residual between " << xMin << " and " << xMax << std::endl;
  double a0 = 0.0;//f_c->GetParameter(0);
  double b0 = htmp->GetMaximum()*0.5;
  double b1 = f_b->GetParameter(1)*0.5;
  
  double xOffset_b = f_b->GetParameter(2);
  double xOffset_s = -1*center;

  double sigma_b = f_b->GetParameter(3);
  double sigma_s = 20e-3;

  //pitch = f_c->GetParameter(1);
  TF1* f = new TF1("f_Total", Utility::residualFeI4Function , xMin , xMax , 8);
  f->SetParameter(0 , b0) ; f->SetParName(0,"b0");
  f->SetParameter(1 , b1) ; f->SetParName(1,"b1");
  f->FixParameter(2 , a0) ; f->SetParName(2,"a0");

  f->SetParameter(3 , xOffset_b) ; f->SetParLimits(3 , 0.5*xOffset_b ,2.0*xOffset_b); f->SetParName(3,"offset_b");
  f->SetParameter(4 , xOffset_s) ; f->SetParLimits(4 , 0.5*xOffset_s ,2.0*xOffset_s); f->SetParName(4,"offset_s");

  f->SetParameter(5 , sigma_b) ;  f->SetParLimits(5 , 0 ,10.0); f->SetParName(5,"sigma_b");
  f->SetParameter(6 , sigma_s) ;  f->SetParLimits(6 , 0 ,1.0); f->SetParName(6,"sigma_s");

  f->SetParameter(7, pitch); f->SetParLimits(7 , 0 , 2.0*pitch); f->SetParName(7,"pitch");

  f->SetLineColor(kBlack); f->SetLineStyle(2);
  htmp->Fit(f,"SB");
  return f; 
}
TFitResultPtr EventAlignment::fitFeI4Residual( TH1D* htmp ,  TH1D* htmp_bkg ,  double pitch )
{
  double xMin = -12.0 ; double xMax = 12.0; 
  
  std::cout << "First trying to fit background.....\n";
  //TFitResultPtr result_Background =  fitFeI4_Background( htmp , pitch ) ;
  //TFitResult* f_result_B = result_Background.Get();
  // TF1* f_res_B = new TF1("f_b", Utility::normBackgroundFunction, xMin , xMax , 4);
  // f_res_B->SetParameters( f_result_B->GetParams() );
  // f_res_B->SetParName(0,"b0");
  // f_res_B->SetParName(1,"N_norm_b");
  // f_res_B->SetParName(2,"Offset_b");
  // f_res_B->SetParName(3,"sigma_b");
  // f_res_B->SetLineColor(kMagenta); //f_res_B->SetFillColor(kMagenta); f_res_B->SetFillStyle(3004);
  // f_res_B->Print();
  TF1* f_res_Bkg = new TF1("f_bkg", Utility::normBackgroundFunction, -10.0 , 10.0 , 4);
  f_res_Bkg->FixParameter(0,0.0); f_res_Bkg->SetParName(0,"b0"); 
  f_res_Bkg->SetParameter(1,htmp_bkg->GetMaximum()); f_res_Bkg->SetParName(1,"N_norm_b");
  f_res_Bkg->SetParameter(2,htmp_bkg->GetMean()); f_res_Bkg->SetParName(2,"Offset_b");
  f_res_Bkg->SetParameter(3,htmp_bkg->GetRMS()); f_res_Bkg->SetParName(3,"sigma_b");
  f_res_Bkg->SetLineColor(kMagenta); 
  TFitResultPtr result_Bfit =  htmp_bkg->Fit(f_res_Bkg,"qnSB");
  f_res_Bkg->Print();

  
  TFitResultPtr result_Central =  fitFeI4_Signal( htmp , pitch ) ;
  TFitResult* f_result_C = result_Central.Get();
  TF1* f_res_Central = new TF1("f_s", Utility::signalFunction , xMin , xMax , 4);
  f_res_Central->SetParameters( f_result_C->GetParams() );
  f_res_Central->SetParName(0,"Norm");
  f_res_Central->SetParName(1,"pitch");
  f_res_Central->SetParName(2,"offset");
  f_res_Central->SetParName(3,"sigma");
  f_res_Central->SetLineColor(kRed); f_res_Central->SetFillColor(kRed); f_res_Central->SetFillStyle(3004);
  f_res_Central->Print();


  double b0 = 0.15;
  double b1 = 0.01; 
  double xOffset_b = 0.0; //f_res_B->GetParameter( f_res_B->GetParNumber("Offset_b"));
  //double sigma_b = f_res_B->GetParameter( f_res_B->GetParNumber("sigma_b"));
  double sigma_b = f_res_Bkg->GetParameter( f_res_Bkg->GetParNumber("sigma_b")); //2.6
  double sigma_b_Error = f_res_Bkg->GetParError( f_res_Bkg->GetParNumber("sigma_b")); //2.6

  double N_norm = f_res_Central->GetParameter( f_res_Central->GetParNumber("Norm"));
  pitch = f_res_Central->GetParameter( f_res_Central->GetParNumber("pitch"));
  double xOffset_s = f_res_Central->GetParameter( f_res_Central->GetParNumber("offset"));
  double sigma_s0 = f_res_Central->GetParameter( f_res_Central->GetParNumber("sigma"));
  double sigma_s1 = 2.0*sigma_s0;
  
  xMin = -100*pitch; //-10*pitch;
  xMax =  100*pitch; //10*pitch;
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  TF1* f = new TF1("f_Total", Utility::residualFeI4_total , xMin , xMax , 9);
  f->SetLineColor(kCyan); f->SetLineStyle(2);
  // f->SetParameter(0,N_norm); f->SetParName(0,"N_norm");
  // f->SetParameter(1,b0); f->SetParName(1,"b0"); f->SetParLimits(1,0,1);
  // f->SetParameter(2,b1); f->SetParName(2,"b1"); f->SetParLimits(2,0,1);
  // f->SetParameter(3,xOffset_b); f->SetParName(3,"xOffset_b");
  // f->FixParameter(4,xOffset_s); f->SetParName(4,"xOffset_s");
  // f->FixParameter(5,sigma_b); f->SetParName(5,"sigma_b"); f->SetParError(5, sigma_b_Error);  //f->SetParLimits(1,0,1);
  // f->SetParameter(6,sigma_s0); f->SetParName(6,"sigma_s0"); //f->SetParLimits(6, 0 , 10.0*sigma_s0 );
  // f->SetParameter(7,sigma_s1); f->SetParName(7,"sigma_s1"); f->SetParLimits(7, 0 , 1.0 );
  // f->SetParameter(8,pitch); f->SetParName(8,"pitch"); //f->SetParLimits(2,0,1);
  // TFitResultPtr result =  htmp->Fit(f,"nSB");
  
  f->SetParameter(0,N_norm); f->SetParName(0,"N_norm");
  f->SetParameter(1,b0); f->SetParName(1,"b0"); f->SetParLimits(1,0,1);
  f->SetParameter(2,b1); f->SetParName(2,"b1"); f->SetParLimits(2,0,1);
  f->FixParameter(3,xOffset_b); f->SetParName(3,"xOffset_b");
  f->FixParameter(4,xOffset_s); f->SetParName(4,"xOffset_s");
  f->SetParameter(5,sigma_b); f->SetParName(5,"sigma_b"); f->SetParLimits(5, 0.75*sigma_b , 1.25*sigma_b);  //f->SetParLimits(1,0,1);
  f->SetParameter(6,sigma_s0); f->SetParName(6,"sigma_s0"); f->SetParLimits(6, 0.75*sigma_s0 , 1.25*sigma_s0 );
  f->SetParameter(7,sigma_s1); f->SetParName(7,"sigma_s1"); f->SetParLimits(7, 0 , 1.0 );
  f->SetParameter(8,pitch); f->SetParName(8,"pitch"); //f->SetParLimits(2,0,1);
  TFitResultPtr result = htmp->Fit(f,"SB");

  xMin = -12; 
  xMax =  12;
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  f->SetRange(xMin,xMax);
  return result;
}
TFitResultPtr EventAlignment::fitFeI4_Signal( TH1D* htmp , double pitch  )
{
  htmp->GetXaxis()->SetRangeUser(-20,20);
  double center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  
  double xMin = center - 5*pitch; 
  double xMax = center + 5*pitch; 
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  TF1* f_c = new TF1("f_Central", Utility::signalFunction , xMin , xMax , 4 );
  f_c->SetParameter(0, htmp->GetMaximum()); f_c->SetParName(0, "Norm"); 
  f_c->SetParameter(1, pitch); f_c->SetParLimits(1 , 0 , 2.0*pitch); f_c->SetParName(1, "pitch"); 
  f_c->SetParameter(2, -1*center ); f_c->SetParName(2, "offset"); 
  f_c->SetParameter(3, 20e-3 ); f_c->SetParLimits(3 , 0 , 100.0 ); f_c->SetParName(3, "sigma"); 
  
  f_c->SetLineColor(kRed); f_c->SetLineStyle(2);
  TFitResultPtr result = htmp->Fit(f_c,"qnSB");
  return result;
}
TFitResultPtr EventAlignment::fitFeI4_Background( TH1D* htmp , double pitch )
{
  
  htmp->GetXaxis()->SetRangeUser(-20,20);
  double center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  double max = htmp->GetMaximum();
  double xMin = -10.0 ;//center - 150*pitch;//(max0 > max1) ? center-xExtra : center+ 1*pitch;
  double xMax =  10.0 ;//center + 150*pitch;//(max0 > max1) ? center-1*pitch : center+xExtra;
  htmp->GetXaxis()->SetRangeUser(xMin,xMax);
  htmp->GetYaxis()->SetRangeUser(0.1,200 );
  
  center = ((double)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((double)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();
  TF1* f_b = new TF1("f_Background", Utility::normBackgroundFunction , xMin , xMax , 4 );
  f_b->FixParameter(0, 0.0);  f_b->SetParName(0,"b0"); 
  f_b->SetParameter(1, htmp->GetMaximum()); f_b->SetParName(1,"N_norm_b");
  f_b->SetParameter(2,  center); f_b->SetParName(2,"Offset_b");
  f_b->SetParameter(3, 2*htmp->GetRMS() ); f_b->SetParLimits(3 , 0 , 10.0 ); f_b->SetParName(3,"sigma_b");
  f_b->SetLineColor(kMagenta); f_b->SetLineStyle(2);
  
  htmp->GetYaxis()->SetRangeUser(0.1, 1.2*max );
  TFitResultPtr result = htmp->Fit(f_b,"nSB");
  return result;
}
void EventAlignment::doFeI4TrackFit_HighIntensity()
{
  double yMin , yMax; 
  double cColPitch_fei4 = 250e-3; 
  double cRowPitch_fei4 = 50e-3; 

  double fitRange = 1.5*std::max(cColPitch_fei4,cRowPitch_fei4);

  std::cout << "=========================================================" << std::endl;
  std::cout << "Fitting FeI4 residual Y" << std::endl;
  std::cout << "=========================================================" << std::endl;
  TH1D* htmpY = dynamic_cast<TH1D*>(hist_->GetHistoByName("TelescopeAnalysis","deltaYPos"));   
  TF1* f_resY = alignFeI4( htmpY , cRowPitch_fei4 , 20*cRowPitch_fei4 ) ;
    
  std::cout << "=========================================================" << std::endl;
  std::cout << "Fitting FeI4 residual X" << std::endl;
  std::cout << "=========================================================" << std::endl;
  TH1D* htmpX = dynamic_cast<TH1D*>(hist_->GetHistoByName("TelescopeAnalysis","deltaXPos"));   
  TF1* f_resX = alignFeI4( htmpX , cColPitch_fei4 , 20*cColPitch_fei4) ;

  double xMin = std::min( -1*f_resY->GetParameter(4) - 5*f_resY->GetParameter(5) , -1*f_resX->GetParameter(4) - 5*f_resX->GetParameter(5) ); 
  double xMax = std::max( -1*f_resY->GetParameter(4) + 5*f_resY->GetParameter(5) , -1*f_resX->GetParameter(4) + 5*f_resX->GetParameter(5) );
  
  //double xExtra = xMax - (-1*f_resY->GetParameter(4) ); 
  //f_resY = alignFeI4( htmpY , cRowPitch_fei4 , xExtra ) ;
  //xExtra = xMax - (-1*f_resX->GetParameter(4) ); 
  //f_resX = alignFeI4( htmpX , cColPitch_fei4 , xExtra ) ;
  htmpX->GetXaxis()->SetRangeUser(xMin,xMax);
  htmpY->GetXaxis()->SetRangeUser(xMin,xMax);
  f_resY->SetRange(xMin,xMax);
  f_resX->SetRange(xMin,xMax);
  
  fResidualY = std::pair<double, double>(-1*f_resY->GetParameter(4)  , f_resY->GetParameter(7)  ) ;
  fTkPointing_y = f_resY->GetParameter(6);
  fFeI4_Window_y = cRowPitch_fei4;//fResidualY.second +2*fTkPointing_y;

  fResidualX = std::pair<double, double>(-1*f_resX->GetParameter(4)  , f_resX->GetParameter(7)  ) ;
  fTkPointing_x = f_resX->GetParameter(6);
  fFeI4_Window_x = cColPitch_fei4 ;//fResidualX.second + 2*fTkPointing_x;

  TF1* f_SignalX = new TF1("f_SgnlX" , Utility::signalFunction , xMin , xMax , 4 );
  f_SignalX->SetParameter(0, f_resX->GetParameter(2) );
  f_SignalX->SetParameter(1, f_resX->GetParameter(7) );
  f_SignalX->SetParameter(2, f_resX->GetParameter(4) );
  f_SignalX->SetParameter(3, f_resX->GetParameter(6) );
  f_SignalX->SetLineColor(kRed); f_SignalX->SetFillColor(kRed); f_SignalX->SetLineStyle(2); f_SignalX->SetFillStyle(3002);
  f_SignalX->Draw("same");



  TF1* f_BkgroundX = new TF1("f_BkgX" , Utility::normBackgroundFunction , xMin , xMax , 4 );
  f_BkgroundX->SetParameter(0, f_resX->GetParameter(0) ) ; 
  f_BkgroundX->SetParameter(1, f_resX->GetParameter(1) ) ; 
  f_BkgroundX->SetParameter(2, f_resX->GetParameter(3) ) ; 
  f_BkgroundX->SetParameter(3, f_resX->GetParameter(5) ) ; 
  f_BkgroundX->SetLineColor(kGreen); f_BkgroundX->SetFillColor(kGreen); f_BkgroundX->SetLineStyle(2); f_BkgroundX->SetFillStyle(3004);
  f_BkgroundX->Draw("SAMe");



  TF1* f_SignalY = new TF1("f_SgnlY" , Utility::signalFunction , xMin , xMax , 4 );
  f_SignalY->SetParameter(0, f_resY->GetParameter(2) );
  f_SignalY->SetParameter(1, f_resY->GetParameter(7) );
  f_SignalY->SetParameter(2, f_resY->GetParameter(4) );
  f_SignalY->SetParameter(3, f_resY->GetParameter(6) );
  f_SignalY->SetLineColor(kRed); f_SignalY->SetFillColor(kRed); f_SignalY->SetLineStyle(2); f_SignalY->SetFillStyle(3002);
  f_SignalY->Draw("SAMe");



  TF1* f_BkgroundY = new TF1("f_BkgY" , Utility::normBackgroundFunction , xMin , xMax , 4 );
  f_BkgroundY->SetParameter(0, f_resX->GetParameter(0) ) ; 
  f_BkgroundY->SetParameter(1, f_resX->GetParameter(1) ) ; 
  f_BkgroundY->SetParameter(2, f_resX->GetParameter(3) ) ; 
  f_BkgroundY->SetParameter(3, f_resX->GetParameter(5) ) ; 
  f_BkgroundY->SetLineColor(kBlue); f_BkgroundY->SetFillColor(kBlue); f_BkgroundY->SetLineStyle(2); f_BkgroundY->SetFillStyle(3005);
  f_BkgroundY->Draw("SAMe");


  xMin = (-1*f_resX->GetParameter(4)) - f_resX->GetParameter(7);
  xMax = (-1*f_resX->GetParameter(4)) + f_resX->GetParameter(7);
  double binWidth = htmpX->GetXaxis()->GetBinWidth(1);
  double nInTimeTracks_X = f_SignalX->Integral(xMin,xMax)/binWidth;
  double nOutOfTimeTracks_X = f_BkgroundX->Integral(xMin,xMax)/binWidth;
  double nTracks_X = f_resX->Integral(xMin,xMax)/binWidth;
  double fractionBackground_X = nOutOfTimeTracks_X/nTracks_X;
  
  xMin = (-1*f_resY->GetParameter(4)) - f_resY->GetParameter(7);
  xMax = (-1*f_resY->GetParameter(4)) + f_resY->GetParameter(7);
  binWidth = htmpY->GetXaxis()->GetBinWidth(1);
  double nInTimeTracks_Y = f_SignalY->Integral(xMin,xMax)/binWidth;
  double nOutOfTimeTracks_Y = f_BkgroundY->Integral(xMin,xMax)/binWidth;
  double nTracks_Y = f_resY->Integral(xMin,xMax)/binWidth;
  double fractionBackground_Y = nOutOfTimeTracks_Y/nTracks_Y;
  

  xMin = fResidualX.first - fResidualX.second;
  xMax = fResidualX.first + fResidualX.second;
  yMin = fResidualY.first - fResidualY.second;
  yMax = fResidualY.first + fResidualY.second;
  
  double int_b0 = f_BkgroundX->Integral(xMin,xMax); 
  double int_b1 = f_BkgroundY->Integral(yMin,yMax); 

  double int_s0 = f_SignalX->Integral(xMin,xMax); 
  double int_s1 = f_SignalY->Integral(yMin,yMax); 

  double b = int_b0*int_b1; 
  double s = int_s0*int_s1;

 
  char buffer[250];
  std::cout << "Summary of the Fits to the residuals:\n";
  sprintf(buffer,"\t\t  Residual X : Offset = %.3f mm , Sigma = %.3f mm , Smear = %.3f mm\n" , fResidualX.first , fResidualX.second , fTkPointing_x);
  std::cout << buffer;
  sprintf(buffer,"\t\t  Residual Y : Offset = %.3f mm , Sigma = %.3f mm , Smear = %.3f mm\n" , fResidualY.first , fResidualY.second , fTkPointing_y);
  std::cout << buffer;
  // std::cout << "Summary of expected number of background (out of sync) tracks in matching window :\n";
  // sprintf(buffer,"\t\t  In X : P(Tk is noise)  = %.3f percent. Size of matching window = %.3f um \n" , 100*fractionBackground_X, 1e3*fFeI4_Window_x );
  // std::cout << buffer;
  // sprintf(buffer,"\t\t  In Y : P(Tk is noise)  = %.3f percent. Size of matching window = %.3f um \n" , 100*fractionBackground_Y, 1e3*fFeI4_Window_y );
  // std::cout << buffer;
  // std::cout << "Summary of expected number of background (out of sync) tracks in matching window :\n";
  // TString out;
  // out.Form("\t\t Background = %.3e\n" , b); std::cout << out.Data();
  // out.Form("\t\t Signal = %.3e\n" , s); std::cout << out.Data();
  // out.Form("\t\t Purity = %.3e\n" , s/(b+s)); std::cout << out.Data();

}
double EventAlignment::ExtrapolateFeI4Hit_DUT( int itrk ,int pDetectorID )
{
  double xOffset = (pDetectorID == 0 ) ? alParsLocal_.d0Offset() : alParsLocal_.d1Offset();
  double zOffset = (pDetectorID == 0 ) ? alParsLocal_.d0Z() : alParsLocal_.d1Z();
  char detName[250]; sprintf(detName,"det%d", pDetectorID);

  double xExtrapolatedTkdut0;
#if defined(MAY_16) || defined(OCT_16)
  xExtrapolatedTkdut0 = telEv()->yPos->at(itrk) + ( zOffset- alParsLocal_.FEI4z() )*telEv()->dydz->at(itrk);
#elif NOV_15
  xExtrapolatedTkdut0 = telEv()->xPos->at(itrk) + ( zOffset- alParsLocal_.FEI4z() )*telEv()->dxdz->at(itrk);
#endif
  double xTk0 =  xExtrapolatedTkdut0 + xOffset; 
  return xTk0; 
}
bool EventAlignment::checkHits(int pDetectorID , vector<int> detHits, double xTrk , int& pNumHits , std::vector<double>& residualsAbs , std::vector<double>& residuals )
{
  double Offset = (pDetectorID == 0 ) ? alParsLocal_.d0Offset() : alParsLocal_.d1Offset() ;
  std::string cHistoName = "deltaXPos_dut" + std::to_string(pDetectorID) ; 
  std::string cEffHistoName = "hHitEfficiency_dut" + std::to_string(pDetectorID) ; 
  std::string cPosHistoName = "hposx_dut" +  std::to_string(pDetectorID) ; 
  for( auto cHitDet  : detHits )
  {
      
    double xHitPos =  (cHitDet - nstrips()/2.0 )*dutpitch() ;
    double xHitResidual = xHitPos - xTrk ; 
    hist_->fillHist1D("TrackMatch",cPosHistoName, xHitPos);
    hist_->fillHist1D("Efficiency",cEffHistoName, cHitDet );
    //if( std::fabs(xHitPos-Offset) <= 10*dutpitch() )
    //{
    residuals.push_back( xHitResidual - Offset );
    residualsAbs.push_back( std::abs(xHitResidual-Offset) );
    //}
    pNumHits++;
    //if( std::fabs(xHitPos-xTrk) < pAcceptanceRegion ) pNumHits++;
  }
  return (pNumHits > 0) ? true : false; 
}
bool EventAlignment::checkHitIsFid( double xHit, int pDetectorID )
{
  char detName[250]; sprintf(detName,"det%d", pDetectorID);

  int cStrip = 999.;
  bool hit_fiducial = CheckFiducial( xHit , 0 , nstrips() , detName );
  //bool hit_fiducial = isTrkfiducial(xHit ,cStrip,detName) ; 
  return hit_fiducial;
}
void EventAlignment::CalculateResolution(int pDetectorID)
{
  char buffer[250];
  sprintf(buffer,"deltaXPos_dut%d", pDetectorID);
  TH1D* res = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch",buffer));   
  res->GetXaxis()->SetRangeUser(-18,18);
  double center = ((double)res->GetMaximumBin())*(res->GetXaxis()->GetXmax()-res->GetXaxis()->GetXmin())/((double)res->GetNbinsX()) + res->GetXaxis()->GetXmin();
  //res->SetAxisRange(center-5*dutpitch(), center+5*dutpitch(), "X");
  //sprintf(buffer,"fResidual_dut%d", pDetectorID);
  
  // TF1* fResidual_Dut = new TF1(buffer, "gaus",center-5*dutpitch(), center+5*dutpitch() ); 
  // fResidual_Dut->SetParameter(0, htmp->GetMaximum() );
  // fResidual_Dut->SetParameter(1, center );
  // fResidual_Dut->SetParameter(2, dutpitch() ); 
  // htmp->Fit(fResidual_Dut,"q");

  double windowSize = 5; 
  double xMin = center -0.1 ; double xMax = center + 0.1;
  double tkPointing = 0.01; 
  res->GetXaxis()->SetRangeUser(xMin,xMax);
  TF1* f_Signal = new TF1("f_Signal" , Utility::stepFunction , xMin , xMax , 4);
  f_Signal->SetParameter(0,res->GetMaximum());
  f_Signal->SetParameter(1,dutpitch()*0.5);
  f_Signal->SetParameter(2,tkPointing);
  f_Signal->SetParameter(3,-1*center);
  res->Fit(f_Signal,"q");
  
  // now given the first guess at the "signal" parameters fit the whole curve 
  xMin = center-windowSize*dutpitch(); xMax = center + windowSize*dutpitch(); 
  res->GetXaxis()->SetRangeUser(xMin,xMax);
  TF1* f_Total = new TF1("f_Total",Utility::residualFunction, xMin , xMax  , 8 );
  f_Total->SetParameter(0,f_Signal->GetParameter(0));
  f_Total->SetParameter(1,f_Signal->GetParameter(1));
  f_Total->SetParameter(2,f_Signal->GetParameter(2));
  f_Total->SetParameter(3,f_Signal->GetParameter(3));
  
  f_Total->SetParameter(4,1);
  f_Total->SetParameter(5,0.0*f_Signal->GetParameter(0));
  f_Total->SetParameter(6,f_Signal->GetParameter(3));
  f_Total->SetParameter(7,2.0); f_Total->SetParLimits(7, 0 , 10) ;
  f_Total->SetLineStyle(2);
  res->Fit(f_Total,"Bq");


  if( pDetectorID == 0 )
  { 
    fDet0_Offset = -1*f_Total->GetParameter(3);
    fDet0_Resolution = res->GetRMS();
    fDet0_Resolution_Unc = res->GetRMSError();

    //fDet0_Offset = fResidual_Dut->GetParameter(1);
    //fDet0_Resolution = fResidual_Dut->GetParameter(2) ;
    //fDet0_Resolution_Unc =  fResidual_Dut->GetParError(2);
  }
  else
  {

    fDet1_Offset = -1*f_Total->GetParameter(3);
    fDet1_Resolution = res->GetRMS();
    fDet1_Resolution_Unc = res->GetRMSError();
    
    // fDet1_Offset = fResidual_Dut->GetParameter(1);
    // fDet1_Resolution = fResidual_Dut->GetParameter(2) ;
    // fDet1_Resolution_Unc =  fResidual_Dut->GetParError(2);

  }
  res->GetXaxis()->SetRangeUser(-18,18.);
  f_Total->SetRange(-18,18.);
  
}
void EventAlignment::CalculateEfficiency_Hits( AllExtrapolatedHits pExtrapolatedHits , double acceptanceWindow )
{
  int cFidTrack_det0 = 0 ;
  int cFidTrack_det1 = 0 ;
  int cFidTracks_both = 0 ; 
  int cFidTracks_either = 0 ; 
  
  int cNumStubs = 0 ;
  int cNumHitsDet0 = 0; int cHitDet0 = 0 ;  int cNumClustersDet0 = 0 ;
  int cNumHitsDet1 = 0; int cHitDet1 = 0 ;  int cNumClustersDet1 = 0 ;
  int nGoodEvents = 0 ;
  char buffer[250];
  //double acceptanceWindow =  nSigma*dutpitch(); 
  for( unsigned int i = 0 ; i < (int)(pExtrapolatedHits.size()) ; i++ )
  {
    ExtrapolatedHits pExtrapolatedHits_forEntry = pExtrapolatedHits[i];
    Long64_t jentry = (Long64_t)(pExtrapolatedHits_forEntry.first);
    std::vector<DutHitsTracking> hits = pExtrapolatedHits_forEntry.second;

    int nTracks = hits.size(); 
    // only accept events with exactly one matched track
    //if( nTracks > 1 ) continue;

    clearEvent();
    
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();   

    int tdcPhase  = condEv()->tdcPhase; 
    std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
    std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
    std::vector<double> x_Clusters_Det0 = getClusterPositions(0);
    std::vector<double> x_Clusters_Det1 = getClusterPositions(1);
    std::vector<double> x_Stubs = getStubPositions();
    std::vector<int> w_Clusters_Det0 = getClusterWidths(0);
    std::vector<int> w_Clusters_Det1 = getClusterWidths(1);

  
    // identify the fiducial hits in detector 0 and detector 1 
    // only consider events where exactly one hit was found 
    //if( hits.size() != 1 ) continue; 
    int nFidHits_Det0 = 0;
    cFidTrack_det0 = 0;
    cFidTrack_det1 = 0;
    for( unsigned int j = 0 ; j < hits.size() ; j++ )
    {
      DutHitsTracking hitInfo_wTracking = hits[j]; 

      double chi2_Tk = hitInfo_wTracking.first; 
      //only consider tracks with a "good" chi2
      //if( chi2_Tk > maxChi2 ) continue ; 

      DutHits tkImpact = hitInfo_wTracking.second; 

      double xHit_Extrapolated_Det0 = tkImpact.first; 
      double xHit_Extrapolated_Det1 = tkImpact.second;

      bool hitFiducial_Det0 =  checkHitIsFid( xHit_Extrapolated_Det0 , 0 ); 
      bool hitFiducial_Det1 =  checkHitIsFid( xHit_Extrapolated_Det1 , 1 );

      if( hitFiducial_Det0  )
      { 
          nFidHits_Det0++;
          double xTrkImpact = xHit_Extrapolated_Det0; 
          int iStrip = xHit_Extrapolated_Det0/dutpitch() + nstrips()/2.0; 
          hist_->fillHist2D("TrackMatch","Tks_Det0", iStrip , tdcPhase ) ;

          fFidTracks_Single_Det0+=1; 
          cFidTrack_det0+=1 ; 
          
          hist_->fillHist1D("TrackMatch","TkImpact_Det0", xTrkImpact ) ;
          hist_->fillHist1D("TrackMatch","TkImpact_Strip_Det0", iStrip ) ;
          
          if( hits.size() == 1 ) hist_->fillHist1D("TrackMatch","SingleTkImpact_Det0", xTrkImpact ) ;
          
          for( unsigned int k = 0 ; k < xHits_Det0.size() ; k++)
          {
            iStrip = xHits_Det0[k]/dutpitch() + nstrips()/2.0; 
            if( checkHitIsFid(xHits_Det0[k] ,0) )
            {
              
              double residual = xHits_Det0[k] -  xTrkImpact ;
              hist_->fillHist1D("TrackMatch","deltaXPos_dut0", residual ) ;
              
              hist_->fillHist1D("det0","deltaXhits", residual ) ;
              hist_->fillHist2D("det0","deltaX2D_Hits" , residual , iStrip );
              hist_->fillHist2D("det0","deltaXhits_Trks", residual , chi2_Tk) ;
              if( xHits_Det0.size() == 1 )
              { 
                hist_->fillHist1D("TrackMatch","deltaXPos_dut0_SingleTks", residual ) ;
                hist_->fillHist1D("det0","deltaXhits_SingleHitEvents", residual ) ;
                hist_->fillHist2D("det0","deltaX2D_SingleHits" , residual , iStrip );
                if( nTracks == 1 )
                { 
                  hist_->fillHist1D("det0","deltaXhits_SingleHitSingleTrackEvents", residual ) ;
                }
              }
              
              if( std::fabs(residual) <= acceptanceWindow  )
              {
                  fHits_Det0 += 1 ; 
                  hist_->fillHist1D("TrackMatch","Hit_Det0", xHits_Det0[k] ) ;
                  hist_->fillHist1D("TrackMatch","Hit_Strip_Det0", iStrip ) ;
                  if( xHits_Det0.size() == 1 )
                  { 
                    hist_->fillHist1D("TrackMatch","SingleHit_Det0", xHits_Det0[k] ) ;
                    hist_->fillHist1D("TrackMatch","SingleHit_Strip_Det0", iStrip ) ;
                  }
              }
              hist_->fillHist1D("TrackMatch","Hit_Strip_Det0_All", iStrip ) ;
              hist_->fillHist2D("TrackMatch","deltaXPos_dut0_2D", residual , xHits_Det0[k]) ;
            }
          }

          hist_->fillHistProfile("det0","nTracks", iStrip , 1.0 ) ;
          for( unsigned int k = 0 ; k < x_Clusters_Det0.size() ; k++)
          {
            iStrip = x_Clusters_Det0[k]/dutpitch() + nstrips()/2.0; 
            if( checkHitIsFid(x_Clusters_Det0[k] ,0) )
            {
              double residual = x_Clusters_Det0[k] -  xTrkImpact ;
              hist_->fillHist1D("det0","deltaXclusters", residual ) ;
              hist_->fillHist1D("TrackMatch","deltaXCluster_dut0", residual ) ;

              hist_->fillHist2D("det0","deltaX2D_Clusters" , residual , iStrip );
              hist_->fillHist2D("det0","deltaXclusters_Trks", residual , chi2_Tk) ;
              if( nTracks == 1 )
              {
                hist_->fillHist1D("det0","deltaXclusters_SingleTrkEvents", residual ) ;
              }

              if( w_Clusters_Det0[k] == 1 )
              { 
                hist_->fillHist1D("TrackMatch","SingleCluster_Det0", x_Clusters_Det0[k] ) ;
                hist_->fillHist1D("TrackMatch","SingleTkImpact_Strip_Det0", iStrip ) ;
              
                hist_->fillHist1D("TrackMatch","deltaX_SingleCluster_dut0", residual ) ;
                hist_->fillHist1D("TrackMatch","nTracks_SCevents_DUT", 1 ) ;
                
                hist_->fillHist1D("det0","deltaXclusters_SingleClstrEvents", residual ) ;
                hist_->fillHist2D("det0","deltaX2D_SingleClusters" , residual , iStrip );
                if( nTracks == 1 ) hist_->fillHist1D("det0","deltaXclusters_SingleClstrSingleTrkEvents", residual ) ;
              }

              hist_->fillHist2D("TrackMatch","wClusterXPos_dut0", x_Clusters_Det0[k] ,  w_Clusters_Det0[k] ) ;
              hist_->fillHist2D("TrackMatch","ResTDC_Clusters_Det0", residual,  tdcPhase ) ;
              if( std::fabs(residual) <= acceptanceWindow )
              {
                fClusters_Det0 += 1 ; 
                hist_->fillHist1D("TrackMatch","Cluster_Det0", x_Clusters_Det0[k] ) ;
                hist_->fillHist1D("TrackMatch","Cluster_Strip_Det0", iStrip ) ;
                hist_->fillHist1D("TrackMatch","deltaXCluster_dut0_M", residual ) ;
                hist_->fillHist1D("TrackMatch","wClusters_Det0", w_Clusters_Det0[k] ) ;
                hist_->fillHist1D("TrackMatch","tdc_Cluster_dut0", tdcPhase ) ;
                hist_->fillHist2D("TrackMatch","Clusters_Det0", iStrip , tdcPhase ) ;

                hist_->fillHistProfile("det0","hitOccupancy", iStrip , 1.0 ) ;
                if( w_Clusters_Det0[k] == 1 ) hist_->fillHist1D("TrackMatch","SingleCluster_Strip_Det0", iStrip ) ;
              }
              else
              {
                hist_->fillHistProfile("det0","noHitOccupancy", iStrip , 1.0 ) ;
                hist_->fillHist1D("TrackMatch","deltaXCluster_dut0_NM", residual ) ;
              }
              
              residual = x_Clusters_Det0[k] -  (xTrkImpact-50*dutpitch() ); 
              if( std::fabs(residual) <= acceptanceWindow  ) 
              {
                hist_->fillHistProfile("det0","noiseOccupancy", iStrip , 1.0/(255 - acceptanceWindow/dutpitch() ) ) ;
              }
              residual = x_Clusters_Det0[k] -  (xTrkImpact+50*dutpitch() ); 
              if( std::fabs(residual) <= acceptanceWindow  ) 
              {
                hist_->fillHistProfile("det0","noiseOccupancy", iStrip , 1.0/(255 - acceptanceWindow/dutpitch() ) ) ;
              }

              hist_->fillHist1D("TrackMatch","Cluster_Strip_Det0_All", iStrip ) ;
            }
          }

      }
      if( hitFiducial_Det1 )
      { 
          double xHit = xHit_Extrapolated_Det1 ; 
          int iStrip = xHit/dutpitch() + nstrips()/2.0; 
          hist_->fillHist2D("TrackMatch","Tks_Det1", iStrip , tdcPhase ) ;

          fFidTracks_Single_Det1+=1; 
          cFidTrack_det1+=1 ; 
          
          hist_->fillHist1D("TrackMatch","TkImpact_Det1", xHit ) ;
          hist_->fillHist1D("TrackMatch","TkImpact_Strip_Det1", iStrip ) ;
          
          if( hits.size() == 1 ) hist_->fillHist1D("TrackMatch","SingleTkImpact_Det1", xHit ) ;

          for( unsigned int k = 0 ; k < xHits_Det1.size() ; k++)
          {
            iStrip = xHits_Det1[k]/dutpitch() + nstrips()/2.0; 
            if( checkHitIsFid(xHits_Det1[k] ,0) )
            {
              
              double residual = xHits_Det1[k] -  xHit ;
              hist_->fillHist1D("det1","deltaXhits", residual ) ;
              hist_->fillHist1D("TrackMatch","deltaXPos_dut1", residual ) ;
              
              hist_->fillHist2D("det1","deltaX2D_Hits" , residual , iStrip );
              if( hits.size() == 1 )
              { 
                hist_->fillHist2D("det1","deltaX2D_SingleHits" , residual , iStrip );
                hist_->fillHist1D("TrackMatch","deltaXPos_dut1_SingleTks", residual ) ;
                hist_->fillHist1D("det1","deltaXhits_SingleHitEvents", residual ) ;
                if( nTracks == 1 ) hist_->fillHist1D("det1","deltaXhits_SingleHitSingleTrackEvents", residual ) ;
              }
              
              if( std::fabs(residual) <= acceptanceWindow  )
              {
                fHits_Det1 += 1 ; 
                hist_->fillHist1D("TrackMatch","Hit_Det1", xHits_Det1[k] ) ;
                hist_->fillHist1D("TrackMatch","Hit_Strip_Det1", iStrip ) ;
                if( xHits_Det1.size() == 1 )
                { 
                  hist_->fillHist1D("TrackMatch","SingleHit_Det1", xHits_Det1[k] ) ;
                  hist_->fillHist1D("TrackMatch","SingleHit_Strip_Det1", iStrip ) ;
                }
              }
              hist_->fillHist1D("TrackMatch","Hit_Strip_Det1_All", iStrip ) ;
              hist_->fillHist2D("TrackMatch","deltaXPos_dut1_2D", residual , xHits_Det1[k]) ;
              hist_->fillHist2D("det1","deltaXhits_Trks", residual , chi2_Tk) ;
            }
          }
          for( unsigned int k = 0 ; k < x_Clusters_Det1.size() ; k++)
          {
            iStrip = x_Clusters_Det1[k]/dutpitch() + nstrips()/2.0; 
            if( checkHitIsFid(x_Clusters_Det1[k] ,0) )
            {
              double residual = x_Clusters_Det1[k] -  xHit ;
              hist_->fillHist1D("det1","deltaXclusters", residual ) ;
              hist_->fillHist2D("det1","deltaX2D_Clusters" , residual , iStrip );
              hist_->fillHist1D("TrackMatch","deltaXCluster_dut1", residual ) ;
              if( w_Clusters_Det1[k] == 1 )
              { 

                hist_->fillHist1D("TrackMatch","deltaX_SingleCluster_dut1", residual ) ;
                hist_->fillHist1D("TrackMatch","SingleTkImpact_Strip_Det1", iStrip ) ;

                hist_->fillHist2D("det1","deltaX2D_SingleClusters" , residual , iStrip );
                hist_->fillHist1D("det1","deltaXclusters_SingleClstrEvents", residual ) ;
                if( nTracks == 1 ) hist_->fillHist1D("det1","deltaXclusters_SingleClstrSingleTrkEvents", residual ) ;
              }
              hist_->fillHist2D("TrackMatch","wClusterXPos_dut1", x_Clusters_Det1[k] ,  w_Clusters_Det1[k] ) ;
              hist_->fillHist2D("TrackMatch","ResTDC_Clusters_Det1", residual,  tdcPhase ) ;
              hist_->fillHist2D("det1","deltaXclusters_Trks", residual , chi2_Tk) ;
              
              if( std::fabs(residual) <= acceptanceWindow )
              {
                fClusters_Det1 += 1 ; 
                hist_->fillHist1D("TrackMatch","Cluster_Det1", x_Clusters_Det1[k] ) ;
                hist_->fillHist1D("TrackMatch","Cluster_Strip_Det1", iStrip ) ;
                hist_->fillHist1D("TrackMatch","deltaXCluster_dut1_M", residual ) ;
                hist_->fillHist1D("TrackMatch","wClusters_Det1", w_Clusters_Det1[k] ) ; 
                hist_->fillHist1D("TrackMatch","tdc_Cluster_dut1", tdcPhase ) ;
                hist_->fillHist2D("TrackMatch","Clusters_Det1", iStrip , tdcPhase ) ;
                if( w_Clusters_Det1[k] == 1 ) hist_->fillHist1D("TrackMatch","SingleCluster_Strip_Det1", iStrip ) ;
              }
              else
              {
                hist_->fillHist1D("TrackMatch","deltaXCluster_dut1_NM", residual ) ;
              }
              hist_->fillHist1D("TrackMatch","Cluster_Strip_Det1_All", iStrip ) ;
            }
          }

      }
      if( hitFiducial_Det0 && hitFiducial_Det1 )
      {
        fFidTracks+=1; 
        for( unsigned int k = 0 ; k < x_Stubs.size() ; k++)
        {
          int iStrip = x_Stubs[k]/dutpitch() + nstrips()/2.0; 

          if( checkHitIsFid(x_Stubs[k] ,0) && checkHitIsFid(x_Stubs[k] ,1) )
          {
            hist_->fillHist1D("TrackMatch","Stubs", x_Stubs[k] ) ;
            if( x_Stubs.size() == 1 ) hist_->fillHist1D("TrackMatch","SingleStubs", x_Stubs[k] ) ;

            double residual = x_Stubs[k] -  xHit_Extrapolated_Det1  ;
            hist_->fillHist1D("TrackMatch","deltaXStub", residual ) ;
            hist_->fillHist2D("det1","deltaX2D_Stubs" , residual , iStrip );
            hist_->fillHist1D("det1", "deltaXstubs" , residual ) ;

            if ( std::fabs(residual) <= acceptanceWindow )
            {
              fStubs += 1;
              hist_->fillHist1D("TrackMatch","Stub_Strip", iStrip ) ;
              hist_->fillHistProfile("det1","stubOccupancy", iStrip , 1.0 ) ;
            }
          }
        }

        if( x_Clusters_Det0.size() == x_Clusters_Det1.size() && x_Clusters_Det0.size() == 1 )
        {
          hist_->fillHist1D("TrackMatch","deltaXCluster", x_Clusters_Det1[0]- x_Clusters_Det0[0] ) ;
        }
      }
    }
    if( nFidHits_Det0!= 0 ) hist_->fillHist1D("TrackMatch","nTracks_DUT", nFidHits_Det0 ) ;
    if( xHits_Det0.size() == 1 && nFidHits_Det0!=0 ) hist_->fillHist1D("TrackMatch","nTracks_SHevents_DUT", nFidHits_Det0 ) ;

    hist_->fillHist1D("det0","nFiducialTracks", cFidTrack_det0 ) ;
    hist_->fillHist1D("det1","nFiducialTracks", cFidTrack_det1 ) ;
  }
  
}
AllExtrapolatedHits EventAlignment::MatchFeI4Hits(int groupID, int evOffset, double windowSize  )
{
  AllExtrapolatedHits cExtrapolatedHits; cExtrapolatedHits.clear();
  double rotation = 0;
  //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  double cScale = windowSize; 
  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;
  

  for (Long64_t jentry=(groupID*fNumberEvents_perGroup ) ; jentry<(groupID+1)*fNumberEvents_perGroup ;jentry++) 
  {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;
      int eventID = jentry; 
      ExtrapolatedHits pMatch; 
      pMatch.first = eventID; 

    
      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;
      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;

      // remove track duplicates 
      std::vector<int> idTks;
      removeDuplicateTracks(idTks);
      if( idTks.size() == 0 ) continue; 

      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
    

      //use only events with one cluster in the FeI4 
      //if( xPos_FeI4Clusters.size() != 1 ) continue; 

      std::vector<int> matchedTkIDs; matchedTkIDs.clear();
      hist_->fillHist1D("TrackMatch","nTracks", idTks.size() );
      AllMatchedTracks clusters; clusters.clear();
      int nClusters= 0; 
      int mMatchedClusters=0;
      for( unsigned int iCluster = 0 ; iCluster < xPos_FeI4Clusters.size() ; iCluster++ )
      {
        if( std::fabs(lv1_FeI4Clusters[iCluster] - fL1_mean) > 1*fL1_sigma ) continue; 
        double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[iCluster])*1e-3;
        //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 

        hist_->fillHist2D("TrackMatch","HitMap", xPos_FeI4Clusters[iCluster] , yPos_FeI4Clusters[iCluster] );
        nClusters++;
        double xval = xPos_FeI4Clusters[iCluster];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[iCluster]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          
        std::vector<int> tkIDs; tkIDs.clear();
        std::vector<double> tkDistance; tkDistance.clear();
        for(unsigned int iTrk = 0 ; iTrk < idTks.size(); iTrk++ )
        {
           unsigned int itk = idTks[iTrk];
           double tkX =  -1*telEv()->xPos->at(itk); 
           double tkY =  telEv()->yPos->at(itk); 
           double chi2 = telEv()->chi2->at(itk); 

           double residualX = xval - tkX - fResidualX.first;
           double residualY = yval - tkY - fResidualY.first;

           double distance = std::sqrt( std::pow(residualX,2) + std::pow(residualY,2));
           if( std::fabs(residualX) <= cScale*fFeI4_Window_x &&  std::fabs(residualY) <= cScale*fFeI4_Window_y )
           {
              tkIDs.push_back( itk );
              tkDistance.push_back( distance );
           }
        }
        
        if( tkIDs.size() > 0 )
        { 
          mMatchedClusters++;

          hist_->fillHist1D("TrackMatch","hToT_M",  tot_FeI4Clusters[iCluster] );
          hist_->fillHist1D("TrackMatch","lv1_Cluster", lv1_FeI4Clusters[iCluster] );
        }
        else
        {
          hist_->fillHist1D("TrackMatch","hToT_UM",  tot_FeI4Clusters[iCluster] );
        }
        if( tkIDs.size() == 0 ) continue;
        
        hist_->fillHist1D("TrackMatch","nTracks_M", tkIDs.size() );
        std::vector<double>::iterator minLocator = std::min_element(std::begin(tkDistance), std::end(tkDistance));
        int id = (int)(std::distance(std::begin(tkDistance), minLocator));
        int itk = tkIDs[id];
        if( (find (matchedTkIDs.begin(), matchedTkIDs.end(), itk) == matchedTkIDs.end()) ) 
        {
          
          xval =  -1*telEv()->xPos->at(itk); 
          yval = telEv()->yPos->at(itk);
          double chi2 = telEv()->chi2->at(itk);
          double nDOF = telEv()->ndof->at(itk); 



          hist_->fillHist2D("TrackMatch","TkHitMap", xval , yval );

          DutHits hits;
          hits.first = ExtrapolateFeI4Hit_DUT( itk  , 0 );
          hits.second = ExtrapolateFeI4Hit_DUT( itk , 1 );

          DutHitsTracking hits_wTracking; 
          hits_wTracking.first = chi2/nDOF; 
          hits_wTracking.second = hits;
          
          std::vector<DutHitsTracking> hitContainter = pMatch.second; 
          hitContainter.push_back(hits_wTracking );
          pMatch.second = hitContainter;
          
          matchedTkIDs.push_back( itk);
          if( idTks.size() == 1 ) hist_->fillHist1D("TrackMatch","hTkChi2_SingleTks",  chi2/nDOF ) ;
          else  hist_->fillHist1D("TrackMatch","hTkChi2_MultiTks",  chi2/nDOF ) ;

          hist_->fillHist1D("TrackMatch","hTkChi2",  chi2/nDOF ) ;
          hist_->fillHist1D("TrackMatch","hTkNdof",  telEv()->ndof->at(itk));
        }
      }
      hist_->fillHist2D("TrackMatch","nTrkCluster_FeI4", nClusters , idTks.size() );
      if( nClusters > 0 ) hist_->fillHist1D("TrackMatch","nClustersFei4", nClusters );
      if( nClusters == 1 )
      {
        hist_->fillHist1D("TrackMatch","nMactchedClusters_1", mMatchedClusters );
      }     
      if( nClusters == 2 )
      {
        hist_->fillHist1D("TrackMatch","nMactchedClusters_2", mMatchedClusters );
      }
      if( nClusters == 3 )
      {
        hist_->fillHist1D("TrackMatch","nMactchedClusters_3", mMatchedClusters );
      }

      if( matchedTkIDs.size() == 0 ) continue; 
      hist_->fillHist1D("TrackMatch","nClustersFei4_M", mMatchedClusters);
      hist_->fillHist1D("TrackMatch","nTracks2", matchedTkIDs.size() );
      cExtrapolatedHits.push_back( pMatch );
  }
  return cExtrapolatedHits;
}
void EventAlignment::CalculateHitResidual( AllExtrapolatedHits pExtrapolatedHits  )
{
  char buffer[250];
  double chi2_cut  = 1.0; //0.20 ; 
  for( unsigned int i = 0 ; i < (int)(pExtrapolatedHits.size()) ; i++ )
  {
    ExtrapolatedHits pExtrapolatedHits_forEntry = pExtrapolatedHits[i];
   
    Long64_t jentry = (Long64_t)(pExtrapolatedHits_forEntry.first);
    std::vector<DutHitsTracking> hits = pExtrapolatedHits_forEntry.second;
    
    int nTracks = hits.size(); 
    // only accept events with exactly one matched track
    //if( nTracks > 3 ) continue;

    clearEvent();
    
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();   

    int tdcPhase  = condEv()->tdcPhase; 
    std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
    std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
    std::vector<double> x_Clusters_Det0 = getClusterPositions(0,3);
    std::vector<double> x_Clusters_Det1 = getClusterPositions(1,3);
    std::vector<int> w_Clusters_Det0 = getClusterWidths(0,3);
    std::vector<int> w_Clusters_Det1 = getClusterWidths(1,3);
    std::vector<double> x_CbcStubs = getCbcStubPositions();
    std::vector<double> x_RecoStubs = getRecoStubPositions();
    std::vector<double> x_RecoSeeds = getRecoStubSeeds();
    std::vector<double> x_RecoMatches = getRecoStubMatches();
    std::vector<double> w_RecoSeeds = getRecoStubSeedWidths();
    std::vector<double> w_RecoMatches = getRecoStubMatchWidths();

    // identify the fiducial hits in detector 0 and detector 1 
    // only consider events where exactly one hit was found 
    //if( hits.size() != 1 ) continue; 
    for( unsigned int j = 0 ; j < hits.size() ; j++ )
    {
      DutHitsTracking hitInfo_wTracking = hits[j]; 

      double chi2_Tk = hitInfo_wTracking.first; 
      //only consider tracks with a "good" chi2
      if( chi2_Tk > chi2_cut ) continue ; 

      DutHits tkImpact = hitInfo_wTracking.second; 

      double xHit_Extrapolated_Det0 = tkImpact.first; 
      double xHit_Extrapolated_Det1 = tkImpact.second;
      
      double trackHit_SeedLayer = returnStrip( xHit_Extrapolated_Det1 );
      double trackHit_CorrelationLayer = returnStrip( xHit_Extrapolated_Det0 );

      for( unsigned int k = 0 ; k < xHits_Det0.size(); k++ )
      {
        double xHit =  xHits_Det0[k]; 
        double hit_nStrips = returnStrip( xHit );
        //std::cout << "Xhit det = " << xHit <<  " impact " << xHit_Extrapolated_Det0 << std::endl ;
        //std::cout << "\t\t Strip det =  "  << hit_nStrips << " impact " << trackHit_CorrelationLayer << std::endl ;
        double residual = xHit - xHit_Extrapolated_Det0;
        hist_->fillHist1D("det0","hitRes", residual ) ;
      }

      for( unsigned int k = 0 ; k < xHits_Det1.size(); k++ )
      {
        double xHit =  xHits_Det1[k]; 
        double hit_nStrips = returnStrip( xHit );
        //std::cout << "Xhit det = " << xHit <<  " impact " << xHit_Extrapolated_Det0 << std::endl ;
        //std::cout << "\t\t Strip det =  "  << hit_nStrips << " impact " << trackHit_CorrelationLayer << std::endl ;
        double residual = xHit - xHit_Extrapolated_Det0;
        hist_->fillHist1D("det1","hitRes", residual ) ;
      }

    }
  }
}
void EventAlignment::alignX(std::vector<SingleGroup_Offset> gOffsets , double windowSize )
{
  double d0Offset_original = alParsLocal_.d0Offset();
  double d1Offset_original = alParsLocal_.d1Offset();

  int nGroups = 5; 
  char buffer[256];
  int nStrips = 4 ;
  fDet0_Offset=0;
  fDet1_Offset=1;
  alParsLocal_.d0Offset(fDet0_Offset);
  alParsLocal_.d1Offset(fDet1_Offset);  
  
  clearHistograms();
  clearTrackMatchhistograms();
  double matchingWindowSize_X = windowSize; 
  double matchingWindowSize_Y = windowSize;
  double area = matchingWindowSize_X*matchingWindowSize_X*(fResidualX.second)*(fResidualY.second);
  double nHits_Det0 , nHits_Det1 , nClusters_Det0 , nClusters_Det1 , nStubs_Det; 
  double nFidTracks_Det0 , nFidTracks_Det1 , nFidStubs ; 
  for(unsigned int i = 0 ; i < nGroups ; i++ )
  {
    int groupID = gOffsets[i].first ;
    int evOffset = gOffsets[i].second ;
    
    AllExtrapolatedHits cMatchedHits = CorrelateFeI4Tracks(groupID, evOffset, matchingWindowSize_X , matchingWindowSize_Y );
    CalculateHitResidual(cMatchedHits);
  }

  int pDetectorID = 0 ; 
  TH1D* res = dynamic_cast<TH1D*>(hist_->GetHistoByName("det0","hitRes"));
  double center = ((double)res->GetMaximumBin())*(res->GetXaxis()->GetXmax()-res->GetXaxis()->GetXmin())/((double)res->GetNbinsX()) + res->GetXaxis()->GetXmin();
  res->SetAxisRange(center-5*dutpitch(), center+5*dutpitch(), "X");
  sprintf(buffer,"fResidual_dut%d", pDetectorID);
  TF1* fResidual_Dut0 = new TF1(buffer, "gaus",center-5*dutpitch(), center+5*dutpitch() ); 
  fResidual_Dut0->SetParameter(0, res->GetMaximum() );
  fResidual_Dut0->SetParameter(1, center );
  fResidual_Dut0->SetParameter(2, dutpitch() ); 
  res->Fit(fResidual_Dut0,"q");
  fDet0_Offset = fResidual_Dut0->GetParameter(1);

  pDetectorID = 1 ; 
  res = dynamic_cast<TH1D*>(hist_->GetHistoByName("det1","hitRes"));
  center = ((double)res->GetMaximumBin())*(res->GetXaxis()->GetXmax()-res->GetXaxis()->GetXmin())/((double)res->GetNbinsX()) + res->GetXaxis()->GetXmin();
  res->SetAxisRange(center-5*dutpitch(), center+5*dutpitch(), "X");
  sprintf(buffer,"fResidual_dut%d", pDetectorID);
  TF1* fResidual_Dut1 = new TF1(buffer, "gaus",center-5*dutpitch(), center+5*dutpitch() ); 
  fResidual_Dut1->SetParameter(0, res->GetMaximum() );
  fResidual_Dut1->SetParameter(1, center );
  fResidual_Dut1->SetParameter(2, dutpitch() ); 
  res->Fit(fResidual_Dut0,"q");
  fDet1_Offset = fResidual_Dut1->GetParameter(1);

  sprintf(buffer, "Offset (X) :  Det0 = %.3f mm [originally %.3f mm], Det1 = %.3f mm  [originally %.3f mm]\n", fDet0_Offset , d0Offset_original ,  fDet1_Offset , d1Offset_original);
  std::cout << buffer ; 
  alParsLocal_.d0Offset(fDet0_Offset);
  alParsLocal_.d1Offset(fDet1_Offset);  
  clearHistograms();

  //setd0XOffset(fDet0_Offset);
  //setd1XOffset(fDet1_Offset);

  // int nGroups_Alignment = 100; 
  // alParsLocal_.d0Offset(0.);
  // alParsLocal_.d1Offset(0.);
  // //align detectors in X 
  // //setd0XOffset(0);
  // //setd1XOffset(0);
  // std::cout << "Aligning in X :" << std::endl; 
  // for( int i = 0 ; i < std::min( (int)(nGroups_Alignment), (int)(gOffsets.size()) ) ; i++ )
  // {
  //   int groupID = gOffsets[i].first;
  //   int evOffset = gOffsets[i].second;
  //   AllExtrapolatedHits cMatchedHits = MatchFeI4Hits(groupID, evOffset, windowSize );

  //   CalculateEfficiency_Hits(cMatchedHits);
  //   if( i%20 == 0 ) std::cout << "...." ; 
  //   if( i%100 == 0 ) std::cout << std::endl; 
  // }
  // std::cout << std::endl;
  // CalculateResolution(0);
  // CalculateResolution(1);
  // sprintf(buffer, "Offset (X) :  Det0 = %.3f mm , Det1 = %.3f mm  \n", fDet0_Offset , fDet1_Offset );
  // std::cout << buffer ; 
  // sprintf(buffer, "Resolution :  Det0 = %.3f mm , Det1 = %.3f mm  \n", fDet0_Resolution , fDet1_Resolution );
  // std::cout << buffer;
  // alParsLocal_.d0Offset(fDet0_Offset);
  // alParsLocal_.d1Offset(fDet1_Offset);  
  // //setd0XOffset(fDet0_Offset);
  // //setd1XOffset(fDet1_Offset);
  // clearHistograms();
}
AllExtrapolatedHits EventAlignment::CorrelateFeI4Tracks(int groupID, int evOffset, double windowSizeX , double windowSizeY )
{
  AllExtrapolatedHits cExtrapolatedHits; cExtrapolatedHits.clear();
  double rotation = 0;
  //default pitch and dimensions of fei4 plane
  double cColPitch_fei4 = 250.0 ; // um 
  double cRowPitch_fei4 =  50.0 ; // um 
  int nPixels = 2; 
  double minDistance = std::sqrt( std::pow(nPixels*cColPitch_fei4*1e-3 ,2.0) + std::pow(nPixels*cRowPitch_fei4*1e-3, 2.0 ) ) ;
  // where did these come from?!
  double x0_FEI4 = -9.875; // mm
  double y0_FEI4 = -8.375; // mm 

  int nEvents = 0 ;
  int nClusters = 0; 
  char buffer[256];
  int nGoodEvents = 0 ;
  
  std::vector<Long64_t> goodEntries; goodEntries.clear();
  std::vector<double> matchTrks_ResidualX; matchTrks_ResidualX.clear();
  for (Long64_t jentry=(groupID*fNumberEvents_perGroup ) ; jentry<(groupID+1)*fNumberEvents_perGroup ;jentry++) 
  {
      clearEvent();
      Long64_t ientry = analysisTree()->GetEntry(jentry);

      if (ientry < 0) break;
      
      // if there are no tracks in the telecscope ... ignore the event!
      if(telEv()->xPos->empty())    continue;
      // if there are no hits in the Fei4 plane ... ignore  the event!
      if(fei4Ev()->nPixHits==0  )    continue;

      if(!isGoodEvent()) continue ;
      int eventID = jentry; 

      clearEvent();
      Long64_t fei4_Entry = jentry;
      fei4_Entry =  jentry + (Long64_t)evOffset;
      ientry = analysisTree()->GetEntry(fei4_Entry);
      if (ientry < 0) break;

      // remove track duplicates 
      std::vector<int> idTks;
      removeDuplicateTracks(idTks);
      if( idTks.size() == 0 ) continue; 

      std::vector<double> xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ;
      int cMinWidth_FeI4Cluster = 3 ; 
      tbeam::FeIFourEvent* cFei4Ev_sorted = sortFei4Hits(fei4Ev());
      if( cFei4Ev_sorted->row->size() == 0 ) continue ; 
      clusterFeI4Hits(cFei4Ev_sorted , xPos_FeI4Clusters, yPos_FeI4Clusters , size_FeI4Clusters , tot_FeI4Clusters , lv1_FeI4Clusters ,   cMinWidth_FeI4Cluster );
      if( xPos_FeI4Clusters.size() == 0 ) continue; 

      ExtrapolatedHits pMatch; 
      pMatch.first = eventID; 

      
      std::vector<int> matchedTkIDs_Y; matchedTkIDs_Y.clear();
      std::vector<int> matchedTkIDs_X; matchedTkIDs_X.clear();
      std::vector<int> matchedTkIDs_X_else; matchedTkIDs_X_else.clear();
      AllMatchedTracks clusters; clusters.clear();
      int nClusters= 0; 
      int mMatchedClusters=0;
      for( unsigned int iCluster = 0 ; iCluster < xPos_FeI4Clusters.size() ; iCluster++ )
      {
        if( std::fabs(lv1_FeI4Clusters[iCluster] - fL1_mean) > 1*fL1_sigma ) continue; 
        double Q = FeI4_ToT_to_Q( tot_FeI4Clusters[iCluster])*1e-3;
        //if( std::fabs(Q - fFeI4_Q_mpv) > 3*fFeI4_Q_sigma) continue ; 

        nClusters++;
        double xval = xPos_FeI4Clusters[iCluster];//*std::cos(rotation*pi/180) - yPos_FeI4Clusters[j]*std::sin(rotation*pi/180) ;
        double yval = yPos_FeI4Clusters[iCluster]; //*std::cos(rotation*pi/180) + xPos_FeI4Clusters[j]*std::sin(rotation*pi/180); 
          
        std::vector<int> tkIDs; tkIDs.clear();
        std::vector<double> tkDistance; tkDistance.clear();
        std::vector<double> tkResidualY; tkResidualY.clear();
        std::vector<double> tkResidualX; tkResidualX.clear();


        std::vector<int> tkIDs_else; tkIDs_else.clear();
        std::vector<double> tkDistance_else; tkDistance_else.clear();
        std::vector<double> tkResidualY_else; tkResidualY_else.clear();
        std::vector<double> tkResidualX_else; tkResidualX_else.clear();
        for(unsigned int iTrk = 0 ; iTrk < idTks.size(); iTrk++ )
        {
           unsigned int itk = idTks[iTrk];
           double tkX =  -1*telEv()->xPos->at(itk); 
           double tkY =  telEv()->yPos->at(itk); 
           double chi2 = telEv()->chi2->at(itk); 

           double residualX, residualY; 
           
#if defined(MAY_16) || defined(OCT_16)
            residualX = (xval - fResidualX.first) - tkX ;
            residualY = (yval - fResidualY.first) - tkY ;
#elif NOV_15
            residualX = (-1*xval - fResidualX.first) - tkY ;
            residualY = (yval - fResidualY.first) - tkX ;
#endif
           double distance = std::sqrt( std::pow(residualX,2) + std::pow(residualY,2)); //std::fabs(residualX) ;
           
           if( std::fabs(residualX) < windowSizeX*fFeI4_Window_x )
           {
              tkIDs.push_back( itk );
              tkDistance.push_back( std::fabs(residualX) );
              tkResidualX.push_back( residualX );
              tkResidualY.push_back( residualY );

           }

           double residualX_else, residualY_else;
#if defined(MAY_16) || defined(OCT_16)
            residualX_else = (xval - (fResidualX.first-10*fResidualX.second)) - tkX ;
            residualY_else = (yval - fResidualY.first) - tkY ;
#elif NOV_15
           residualX_else = (-1*xval - (fResidualX.first-10*fResidualX.second)) - tkY ;
           residualY_else = (yval - fResidualY.first) - tkX ;
#endif
           if( std::fabs(residualX_else) < windowSizeX*fFeI4_Window_x )
           {
              tkIDs_else.push_back( itk );
              tkDistance_else.push_back( std::fabs(residualX_else) );
              tkResidualX_else.push_back( residualX_else );
              tkResidualY_else.push_back( residualY_else );
          }
           
        }
        if( tkIDs.size() == 0 ) continue;
        
        std::vector<double>::iterator minLocator = std::min_element(std::begin(tkDistance), std::end(tkDistance));
        int id = (int)(std::distance(std::begin(tkDistance), minLocator));
        int itk = tkIDs[id];
        if( (find (matchedTkIDs_X.begin(), matchedTkIDs_X.end(), itk) == matchedTkIDs_X.end()) ) 
        {
          
          xval =  -1*telEv()->xPos->at(itk); 
          yval = telEv()->yPos->at(itk);
          double chi2 = telEv()->chi2->at(itk);
          double nDOF = telEv()->ndof->at(itk); 

          //if( chi2/nDOF <= maxChi2 )
          //{
            hist_->fillHist1D("TrackMatch","deltaXPos_trkfei4", tkResidualX[id] ) ;
            hist_->fillHist1D("TrackMatch","deltaYPos_trkfei4", tkResidualY[id] ) ;

            matchedTkIDs_X.push_back( itk );
            if( std::fabs(tkResidualY[id]) <= windowSizeY*fFeI4_Window_y )
            { 
              if( (find (matchedTkIDs_Y.begin(), matchedTkIDs_Y.end(), itk) == matchedTkIDs_Y.end()) ) 
              {
                DutHits hits;
                hits.first = ExtrapolateFeI4Hit_DUT( itk  , 0 );
                hits.second = ExtrapolateFeI4Hit_DUT( itk , 1 );
                //hits.first = ExtrapolateFeI4Hit_DUT( yval , itk  , 0 );
                //hits.second = ExtrapolateFeI4Hit_DUT( yval , itk , 1 );
                bool hitFiducial_Det0 = checkHitIsFid( hits.first , 0 );
                bool hitFiducial_Det1 = checkHitIsFid( hits.second , 0 );
                // only consider tracks that are in the fiducial region of Det0 and Det1 
                if( hitFiducial_Det0 && hitFiducial_Det1 )
                {
                  double iStrip_Det0 = returnStrip( hits.first );
                  double iStrip_Det1 = returnStrip( hits.second );

                  hist_->fillHist1D("TrackMatch","trackImpact_Strips_Det0", (int)iStrip_Det0 ) ;
                  hist_->fillHist1D("TrackMatch","trackImpact_Strips_Det1", (int)iStrip_Det1 ) ;
    
                  DutHitsTracking hits_wTracking; 
                  hits_wTracking.first = chi2/nDOF; 
                  hits_wTracking.second = hits;
                    
                  std::vector<DutHitsTracking> hitContainter = pMatch.second; 
                  hitContainter.push_back(hits_wTracking );
                  pMatch.second = hitContainter;
              
                  matchedTkIDs_Y.push_back( itk );
                  
                  hist_->fillHist1D("TrackMatch","deltaXPos_trkfei4_M", tkResidualX[id] ) ;
                  hist_->fillHist1D("TrackMatch","deltaYPos_trkfei4_M", tkResidualY[id] ) ;
                }
              }
            }
          //}
          //std::cout << "Event " << jentry << " - Hit found : (" << hits.first << "," << hits.second << ")\n";
        }

        if( tkIDs_else.size() != 0 )
        {
          minLocator = std::min_element(std::begin(tkDistance_else), std::end(tkDistance_else));
          id = (int)(std::distance(std::begin(tkDistance_else), minLocator));
          itk = tkIDs_else[id];
          if( (find (matchedTkIDs_X_else.begin(), matchedTkIDs_X_else.end(), itk) == matchedTkIDs_X_else.end()) ) 
          {
            hist_->fillHist1D("TrackMatch","deltaXPos_trkfei4_UM", tkResidualX_else[id] ) ;
            hist_->fillHist1D("TrackMatch","deltaYPos_trkfei4_UM", tkResidualY_else[id] ) ;
            matchedTkIDs_X_else.push_back( itk );
          }
        }

      }
      if( matchedTkIDs_Y.size() == 0 ) continue; 
      goodEntries.push_back( jentry );   
      cExtrapolatedHits.push_back( pMatch );
      int nFidTracks = pMatch.second.size();
      if( nFidTracks > 0 ) 
      {
        hist_->fillHist1D("det0","nFiducialTracks", nFidTracks ) ;
        hist_->fillHist1D("det1","nFiducialTracks", nFidTracks ) ;
      }    
  }
  return cExtrapolatedHits;
}
void EventAlignment::removeDuplicateHits( AllExtrapolatedHits& pExtrapolatedHits )
{
  AllExtrapolatedHits cExtrapolatedHits; cExtrapolatedHits.clear();

  for( unsigned int i = 0 ; i < (int)(pExtrapolatedHits.size()) ; i++ )
  {
    ExtrapolatedHits pExtrapolatedHits_forEntry = pExtrapolatedHits[i];
    Long64_t jentry = (Long64_t)(pExtrapolatedHits_forEntry.first);
    std::vector<DutHitsTracking> hits = pExtrapolatedHits_forEntry.second;
    std::vector<DutHitsTracking> cHits ; cHits.clear();
    for( unsigned int j = 0 ; j < hits.size() ; j++ )
    {
      DutHitsTracking hitInfo_wTracking = hits[j]; 

      double chi2_Tk = hitInfo_wTracking.first; 
      
      DutHits tkImpact = hitInfo_wTracking.second; 

      double xHit_Extrapolated_Det0 = tkImpact.first; 
      double xHit_Extrapolated_Det1 = tkImpact.second;

      if( j == 0 )  cHits.push_back( hitInfo_wTracking );
      else
      {
        DutHitsTracking pHitInfo_wTracking = hits[j-1];
        if(  !((pHitInfo_wTracking.second).first == xHit_Extrapolated_Det0 &&  (pHitInfo_wTracking.second).second == xHit_Extrapolated_Det1 ) )
        {
          cHits.push_back( pHitInfo_wTracking );
        }
        else
        {
          std::cout << "Removed a duplicate!\n";
        }
      }
    
    }
  }
}
void EventAlignment::CalculateDUTefficiency( AllExtrapolatedHits pExtrapolatedHits , double acceptanceWindow , double correlationWindow )
{
  double acceptanceWindow_Strips = 2.0 ; 
  double cbcBoundary_Strips = nstrips()/2.0 ; 
  if( acceptanceWindow == 0 ) acceptanceWindow = acceptanceWindow_Strips*dutpitch();

  int cFidTrack_det0 = 0 ;
  int cFidTrack_det1 = 0 ;
  int cFidTracks_both = 0 ; 
  int cFidTracks_either = 0 ; 
  
  int cNumStubs = 0 ;
  int cNumHitsDet0 = 0; int cHitDet0 = 0 ;  int cNumClustersDet0 = 0 ;
  int cNumHitsDet1 = 0; int cHitDet1 = 0 ;  int cNumClustersDet1 = 0 ;
  int nGoodEvents = 0 ;
  char buffer[250];
  double chi2_cut  = 1.0; //0.20 ; 
  for( unsigned int i = 0 ; i < (int)(pExtrapolatedHits.size()) ; i++ )
  {
    ExtrapolatedHits pExtrapolatedHits_forEntry = pExtrapolatedHits[i];
    std::vector<double> evOffsets{ -2 , -1 , 0 , 1 , 2 };
    
    Long64_t jentry = (Long64_t)(pExtrapolatedHits_forEntry.first);
    std::vector<DutHitsTracking> hits = pExtrapolatedHits_forEntry.second;
    
    int nTracks = hits.size(); 

    // only accept events with exactly one matched track
    if( nTracks > 1 ) continue;

    clearEvent();
    
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();   

    int tdcPhase  = condEv()->tdcPhase; 
    std::vector<double> xHits_Det0 = getHitPositions( d0c0 );
    std::vector<double> xHits_Det1 = getHitPositions( d1c0 );
    std::vector<double> x_Clusters_Det0 = getClusterPositions(0,3);
    std::vector<double> x_Clusters_Det1 = getClusterPositions(1,3);
    std::vector<int> w_Clusters_Det0 = getClusterWidths(0,3);
    std::vector<int> w_Clusters_Det1 = getClusterWidths(1,3);
    std::vector<double> x_CbcStubs = getCbcStubPositions();
    std::vector<double> x_RecoStubs = getRecoStubPositions();
    std::vector<double> x_RecoSeeds = getRecoStubSeeds();
    std::vector<double> x_RecoMatches = getRecoStubMatches();
    std::vector<double> w_RecoSeeds = getRecoStubSeedWidths();
    std::vector<double> w_RecoMatches = getRecoStubMatchWidths();

    // identify the fiducial hits in detector 0 and detector 1 
    // only consider events where exactly one hit was found 
    //if( hits.size() != 1 ) continue; 
    
    int nFidHits_Det0 = 0;
    int nFidHits_Det1 = 0;
    cFidTrack_det0 = 0;
    cFidTrack_det1 = 0;

    int cFidTracks_cbc0 = 0 ;
    int cFidTracks_cbc1 = 0 ;
    int nEvents_funny = 0; 
    for( unsigned int j = 0 ; j < hits.size() ; j++ )
    {
      DutHitsTracking hitInfo_wTracking = hits[j]; 

      double chi2_Tk = hitInfo_wTracking.first; 
      //only consider tracks with a "good" chi2
      if( chi2_Tk > chi2_cut ) continue ; 

      DutHits tkImpact = hitInfo_wTracking.second; 

      double xHit_Extrapolated_Det0 = tkImpact.first; 
      double xHit_Extrapolated_Det1 = tkImpact.second;
      
      bool hitFiducial_Det0 =  checkHitIsFid( xHit_Extrapolated_Det0 , 0 ); 
      bool hitFiducial_Det1 =  checkHitIsFid( xHit_Extrapolated_Det1 , 1 );
      cFidTrack_det0 += hitFiducial_Det0 ? 1 : 0 ;
      cFidTrack_det1 += hitFiducial_Det1 ? 1 : 0 ;
      
      double trackHit_SeedLayer = returnStrip( xHit_Extrapolated_Det1 );
      double trackHit_CorrelationLayer = returnStrip( xHit_Extrapolated_Det0 );

      bool interestingEvent = false;
      if( hitFiducial_Det0 )
      {
        fFidTracks_Det0++; 
        hist_->fillHist1D("det0","fidTracks_chi2",chi2_Tk);
        hist_->fillHist1D("det0","fidTracks_impact_mm", xHit_Extrapolated_Det0);
        hist_->fillHist1D("det0","fidTracks_impact_strips", trackHit_CorrelationLayer);
        
        //cFidTracks_cbc0 += ( (int)(trackHit_CorrelationLayer) <= 254/2.0 ) ? 1 : 0 ;
        //cFidTracks_cbc1 += ( (int)(trackHit_CorrelationLayer) > 254/2.0 ) ? 1 : 0 ;
        if( trackHit_CorrelationLayer < cbcBoundary_Strips )
        { 
          cFidTracks_cbc0++; 
          hist_->fillHist1D("cbc0","fidTracks_impact_mm",xHit_Extrapolated_Det0 );
          hist_->fillHist1D("cbc0","nCandidateClusters", tdcPhase);
        }
        else
        { 
          cFidTracks_cbc1++; 
          hist_->fillHist1D("cbc1","fidTracks_impact_mm",xHit_Extrapolated_Det0 );
          hist_->fillHist1D("cbc1","nCandidateClusters", tdcPhase);
        }

        bool hitFound = false; 
        for( unsigned int k = 0 ; k < xHits_Det0.size(); k++ )
        {
          double xHit =  xHits_Det0[k]; 
          double hit_nStrips = returnStrip( xHit );
          //std::cout << "Xhit det = " << xHit <<  " impact " << xHit_Extrapolated_Det0 << std::endl ;
          //std::cout << "\t\t Strip det =  "  << hit_nStrips << " impact " << trackHit_CorrelationLayer << std::endl ;
          double residual = xHit - xHit_Extrapolated_Det0;
          hist_->fillHist1D("det0","hitRes", residual ) ;
          if( std::fabs( residual ) <= acceptanceWindow_Strips*dutpitch() )
          { 
            hitFound = true; 
          } 
        }
        if( hitFound ) fHits_Det0++;

        bool clusterFound = false; 
        double expected_Strip = (int)( returnStrip( xHit_Extrapolated_Det0 ) );
        hist_->fillHist2D("det0", "fidTracks_impact_mm_2D" , xHit_Extrapolated_Det0 , expected_Strip ) ;   
        for( unsigned int k = 0 ; k < x_Clusters_Det0.size(); k++ )
        {
          double xHit =  x_Clusters_Det0[k]; 
          double hit_nStrips = returnStrip( xHit );
          double residual = xHit - xHit_Extrapolated_Det0; 
          //std::cout << "Cluster residual = " << residual << std::endl;
          hist_->fillHist1D("det0","clusterRes", residual ) ;
           
          hist_->fillHist2D("det0","clusterRes_Chi2", residual , chi2_Tk ) ;
          hist_->fillHist2D("det0","clusterRes_TDC", residual , tdcPhase ) ;
          hist_->fillHist2D("det0","clusterRes_CW", residual , w_Clusters_Det0[k]) ;
           
          if( std::fabs( residual ) <= acceptanceWindow_Strips*dutpitch() )
          { 
              hist_->fillHist2D("det0","matchedCluster_CW", hit_nStrips , w_Clusters_Det0[k]);
              hist_->fillHist1D("det0","matchedCluster_mm", xHit_Extrapolated_Det0 );
              hist_->fillHistProfile("det0","matchedClusterWidth_mm" , xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
              hist_->fillHist1D("det0","matchedCluster_strips", hit_nStrips );
              
              if( hit_nStrips < cbcBoundary_Strips )
              { 
                hist_->fillHist1D("cbc0","matchedCluster_mm", xHit_Extrapolated_Det0 );
                hist_->fillHist2D("cbc0","matchedCluster_CW_mm", xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                hist_->fillHistProfile("cbc0","matchedClusterWidth_mm" , xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                hist_->fillHist1D("cbc0","nClusters_ANA", tdcPhase);
              }
              else
              {
                hist_->fillHist1D("cbc1","matchedCluster_mm", xHit_Extrapolated_Det0 );
                hist_->fillHist2D("cbc1","matchedCluster_CW_mm", xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                hist_->fillHistProfile("cbc1","matchedClusterWidth_mm" , xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                hist_->fillHist1D("cbc1","nClusters_ANA", tdcPhase);
              }
              
              hist_->fillHist2D("det0", "matchedCluster_mm_All_2D" , xHit_Extrapolated_Det0 , expected_Strip ) ;   
              if(w_Clusters_Det0[k] == 1 ) 
              { 
                hist_->fillHist2D("det0", "fidTracksMatched_impact_mm_2D" , xHit_Extrapolated_Det0 , expected_Strip ) ;   
                if( expected_Strip == hit_nStrips ) hist_->fillHist2D("det0", "matchedCluster_mm_2D" , xHit_Extrapolated_Det0, hit_nStrips ) ; 
                if( (hit_nStrips-1)<=  expected_Strip && expected_Strip < (hit_nStrips) ) hist_->fillHist2D("det0", "matchedCluster_nextStrip_mm_2D" , xHit_Extrapolated_Det0, hit_nStrips ) ; 
              }
              if(w_Clusters_Det0[k] == 2 ) 
              { 
                hist_->fillHist2D("det0", "fidTracksMatched_impact_mm_2D_CW2" , xHit_Extrapolated_Det0 , expected_Strip ) ;   
                if( expected_Strip == hit_nStrips ) hist_->fillHist2D("det0", "matchedCluster_mm_2D_CW2" , xHit_Extrapolated_Det0, hit_nStrips ) ; 
                if( (hit_nStrips-1)<=  expected_Strip && expected_Strip < (hit_nStrips) ) hist_->fillHist2D("det0", "matchedCluster_nextStrip_mm_2D_CW2" , xHit_Extrapolated_Det0, hit_nStrips ) ; 
              }
              

              if( xHit_Extrapolated_Det0 >= -1*0.09 && xHit_Extrapolated_Det0 <= 0.0*0.09 )
              {
                //std::cout << "Det0 : xTrack = " << xHit_Extrapolated_Det0 << " Center cluster @ strip " << hit_nStrips << " width cluster = " << w_Clusters_Det0[k] << std::endl;
                hist_->fillHist1D("cbc0","fidTracks_impact_mm_CentralRegion", xHit_Extrapolated_Det0 );
                hist_->fillHist1D("cbc1","fidTracks_impact_mm_CentralRegion", xHit_Extrapolated_Det0 );
                
                hist_->fillHist2D("cbc0","fidTracks_impact_mm_CentralRegion_2D", xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                hist_->fillHist2D("cbc1","fidTracks_impact_mm_CentralRegion_2D", xHit_Extrapolated_Det0 , w_Clusters_Det0[k]);
                if( hit_nStrips < cbcBoundary_Strips)
                { 
                 hist_->fillHist1D("cbc0","matchedCluster_mm_CentralRegion", xHit_Extrapolated_Det0 );
                 hist_->fillHist2D("cbc0","matchedCluster_mm_CentralRegion_2D", xHit_Extrapolated_Det0  , w_Clusters_Det0[k]);
                }
                else
                { 
                  hist_->fillHist1D("cbc1","matchedCluster_mm_CentralRegion", xHit_Extrapolated_Det0 );
                  hist_->fillHist2D("cbc1","matchedCluster_mm_CentralRegion_2D", xHit_Extrapolated_Det0  , w_Clusters_Det0[k]);
                }
              }
                

              clusterFound = true; 
           } 
        }
        if( clusterFound ) fClusters_Det0++;

      }
      
      if( hitFiducial_Det1 )
      {
        fFidTracks_Det1++;
        hist_->fillHist1D("det1","fidTracks_chi2",chi2_Tk);
        hist_->fillHist1D("det1","fidTracks_impact_mm", xHit_Extrapolated_Det1);
        hist_->fillHist1D("det1","fidTracks_impact_strips", trackHit_SeedLayer);

        //cFidTracks_cbc0 += ( (int)(trackHit_SeedLayer) <= 254/2.0 ) ? 1 : 0 ;
        //cFidTracks_cbc1 += ( (int)(trackHit_SeedLayer) > 254/2.0 ) ? 1 : 0 ;
        if( trackHit_SeedLayer < cbcBoundary_Strips)
        { 
          cFidTracks_cbc0++; 
          hist_->fillHist1D("cbc0","fidTracks_impact_mm",xHit_Extrapolated_Det1 );
          hist_->fillHist1D("cbc0","nCandidateClusters", tdcPhase);
        }
        else
        { 
          cFidTracks_cbc1++; 
          hist_->fillHist1D("cbc1","fidTracks_impact_mm",xHit_Extrapolated_Det1 );
          hist_->fillHist1D("cbc1","nCandidateClusters", tdcPhase);
        }

        bool hitFound = false; 
        for( unsigned int k = 0 ; k < xHits_Det1.size(); k++ )
        {
           double xHit =  xHits_Det1[k]; 
           double hit_nStrips = returnStrip( xHit );
           double residual = xHit - xHit_Extrapolated_Det1;
           hist_->fillHist1D("det1","hitRes", residual ) ;
           if( std::fabs( residual ) <= acceptanceWindow_Strips*dutpitch() )
           { 
              hitFound = true; 
           } 
        }
        if( hitFound ) fHits_Det1++;

        bool clusterFound = false; 
        double expected_Strip = (int)( returnStrip( xHit_Extrapolated_Det1 ) );
        hist_->fillHist2D("det1", "fidTracks_impact_mm_2D" , xHit_Extrapolated_Det1 , expected_Strip ) ;   
        for( unsigned int k = 0 ; k < x_Clusters_Det1.size(); k++ )
        {
           double xHit =  x_Clusters_Det1[k]; 
           double hit_nStrips = returnStrip( xHit );
           double residual = xHit - xHit_Extrapolated_Det1;
           //double hit_nStrips = returnStrip( x_Clusters_Det1[k] );
           //double residual = x_Clusters_Det1[k] - xHit_Extrapolated_Det1;
           hist_->fillHist1D("det1","clusterRes", residual ) ;
           hist_->fillHist2D("det1","clusterRes_Chi2", residual , chi2_Tk ) ;
           hist_->fillHist2D("det1","clusterRes_TDC", residual , tdcPhase ) ;
           hist_->fillHist2D("det1","clusterRes_CW", residual , w_Clusters_Det1[k]) ;

           if( std::fabs( residual ) <= acceptanceWindow_Strips*dutpitch() )
           {
              hist_->fillHist2D("det1","matchedCluster_CW", hit_nStrips , w_Clusters_Det1[k]);
              hist_->fillHist1D("det1","matchedCluster_mm", xHit_Extrapolated_Det1 );
              hist_->fillHist1D("det1","matchedCluster_strips", hit_nStrips );
              hist_->fillHistProfile("det1","matchedClusterWidth_mm" , xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);

              if( hit_nStrips < cbcBoundary_Strips )
              { 
                hist_->fillHist1D("cbc0","matchedCluster_mm", xHit_Extrapolated_Det1 );
                hist_->fillHist2D("cbc0","matchedCluster_CW_mm", xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                hist_->fillHistProfile("cbc0","matchedClusterWidth_mm" , xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                hist_->fillHist1D("cbc0","nClusters_ANA", tdcPhase);
              }
              else
              { 
                hist_->fillHist1D("cbc1","matchedCluster_mm", xHit_Extrapolated_Det1 );
                hist_->fillHist2D("cbc1","matchedCluster_CW_mm", xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                hist_->fillHistProfile("cbc1","matchedClusterWidth_mm" , xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                hist_->fillHist1D("cbc1","nClusters_ANA", tdcPhase);
              }
              clusterFound = true; 
              
              hist_->fillHist2D("det1", "matchedCluster_mm_All_2D" , xHit_Extrapolated_Det1 , expected_Strip ) ;   
              if(w_Clusters_Det1[k] == 1 ) 
              { 
                hist_->fillHist2D("det1", "fidTracksMatched_impact_mm_2D" , xHit_Extrapolated_Det1 , expected_Strip ) ;   
                if( expected_Strip == hit_nStrips ) hist_->fillHist2D("det1", "matchedCluster_mm_2D" , xHit_Extrapolated_Det1, hit_nStrips ) ; 
                if( (hit_nStrips-1)<=  expected_Strip && expected_Strip < (hit_nStrips) ) hist_->fillHist2D("det1", "matchedCluster_nextStrip_mm_2D" , xHit_Extrapolated_Det1, hit_nStrips ) ; 
              }
              if(w_Clusters_Det1[k] == 2 ) 
              { 
                hist_->fillHist2D("det1", "fidTracksMatched_impact_mm_2D_CW2" , xHit_Extrapolated_Det1 , expected_Strip ) ;   
                if( expected_Strip == hit_nStrips ) hist_->fillHist2D("det1", "matchedCluster_mm_2D_CW2" , xHit_Extrapolated_Det1, hit_nStrips ) ; 
                if( (hit_nStrips-1)<=  expected_Strip && expected_Strip < (hit_nStrips) ) hist_->fillHist2D("det1", "matchedCluster_nextStrip_mm_2D_CW2" , xHit_Extrapolated_Det1, hit_nStrips ) ; 
              }

              if( xHit_Extrapolated_Det1 >= -1*0.09 && xHit_Extrapolated_Det1 <= 0.0*0.09 )
              {
                //std::cout << "Det0 : xTrack = " << xHit_Extrapolated_Det1 << " Center cluster @ strip " << hit_nStrips << " width cluster = " << w_Clusters_Det1[k] << std::endl;
                hist_->fillHist1D("cbc0","fidTracks_impact_mm_CentralRegion", xHit_Extrapolated_Det1 );
                hist_->fillHist1D("cbc1","fidTracks_impact_mm_CentralRegion", xHit_Extrapolated_Det1 );
                
                hist_->fillHist2D("cbc0","fidTracks_impact_mm_CentralRegion_2D", xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                hist_->fillHist2D("cbc1","fidTracks_impact_mm_CentralRegion_2D", xHit_Extrapolated_Det1 , w_Clusters_Det1[k]);
                if( hit_nStrips < cbcBoundary_Strips)
                { 
                 hist_->fillHist1D("cbc0","matchedCluster_mm_CentralRegion", xHit_Extrapolated_Det1 );
                 hist_->fillHist2D("cbc0","matchedCluster_mm_CentralRegion_2D", xHit_Extrapolated_Det1  , w_Clusters_Det1[k]);
                }
                else
                { 
                  hist_->fillHist1D("cbc1","matchedCluster_mm_CentralRegion", xHit_Extrapolated_Det1 );
                  hist_->fillHist2D("cbc1","matchedCluster_mm_CentralRegion_2D", xHit_Extrapolated_Det1  , w_Clusters_Det1[k]);
                }
              }
           } 
        }
        if( clusterFound ) fClusters_Det1++;

      }
      if( hitFiducial_Det0 && hitFiducial_Det1 )
      {
        cFidTracks_both += 1; 
        bool trackSub_cbc0;
        bool trackSub_cbc1;
        double trackBend_Strips = trackHit_SeedLayer-trackHit_CorrelationLayer;
        hist_->fillHist1D("det1","fidTracks_Bends", trackBend_Strips);

        if( std::fabs(trackBend_Strips) <=  correlationWindow )
        {
          bool stubFound  = false;
          hist_->fillHist1D("det1","candidateStubs_impact_mm", xHit_Extrapolated_Det1);
          hist_->fillHist2D("det1","candidateStubs_impact_mm_2D", xHit_Extrapolated_Det1, tdcPhase);
          hist_->fillHist1D("det1","candidateStubs_impact_strips", trackHit_SeedLayer);
          
           // if track looks like it should generate a stub then fill candidate stub histogram
          if( trackHit_SeedLayer < cbcBoundary_Strips )// && (x_Clusters_Det0.size() == 1  && x_Clusters_Det1.size() == 1 ) )
          { 
            hist_->fillHist1D("cbc0","candidateStubs_impact_mm", xHit_Extrapolated_Det1 );
            hist_->fillHist2D("cbc0","candidateStubs_impact_mm_2D", xHit_Extrapolated_Det1 , tdcPhase );

            hist_->fillHist1D("cbc0","nCandidateStubs", tdcPhase);
            hist_->fillHist1D("cbc0","stubMatching", 0);
            hist_->fillHist2D("cbc0","stubMatching_2D", 0 , xHit_Extrapolated_Det1);
            hist_->fillHist2D("cbc0","stubMatching_TDC", 0 , tdcPhase);
          }
          else
          //else if  (x_Clusters_Det0.size() == 1  && x_Clusters_Det1.size() == 1 ) 
          {
            hist_->fillHist1D("cbc1","candidateStubs_impact_mm", xHit_Extrapolated_Det1 );
            hist_->fillHist2D("cbc1","candidateStubs_impact_mm_2D", xHit_Extrapolated_Det1 , tdcPhase );

            hist_->fillHist1D("cbc1","nCandidateStubs", tdcPhase);
            hist_->fillHist1D("cbc1","stubMatching", 0);
            hist_->fillHist2D("cbc1","stubMatching_2D", 0 , xHit_Extrapolated_Det1);
            hist_->fillHist2D("cbc1","stubMatching_TDC", 0 , tdcPhase);
          }
          
          if( x_Clusters_Det0.size() != 0 && x_Clusters_Det1.size() != 0)
          {
            uint8_t chipStub = (uint8_t)(dutEv()->stubWord);
            if( chipStub != 0 ) hist_->fillHist1D("det1","cbcStub", tdcPhase ) ;
            uint8_t recoStub = (uint8_t)(dutEv()->stubWordReco);
            if( recoStub != 0 ) hist_->fillHist1D("det1","recoStub", tdcPhase ) ;
            
            std::vector<FoundStubInfo> xStubs_Seeds , xStubs_Matched;
            xStubs_Seeds.clear(); xStubs_Matched.clear();
            uint8_t analysisStub = findStubs( xStubs_Seeds ,  xStubs_Matched ,xHits_Det0 , xHits_Det1 , correlationWindow, true );

            bool anaStub_1Strip_Found = false; 
            bool anaStub_2Strip_Found = false; 
            bool anaStub_3Strip_Found = false; 
            
            //if( x_Clusters_Det0.size() == 1  && x_Clusters_Det1.size() == 1 )
            //{
            if( chipStub == 3 ){ hist_->fillHist1D("cbc0","nStubs_CHIP", tdcPhase ); hist_->fillHist1D("cbc1","nStubs_CHIP", tdcPhase ); }
            if( recoStub == 3 ){ hist_->fillHist1D("cbc0","nStubs_RECO", tdcPhase ); hist_->fillHist1D("cbc1","nStubs_RECO", tdcPhase ); }
            if( chipStub == 1 ) hist_->fillHist1D("cbc0","nStubs_CHIP", tdcPhase ); 
            if( recoStub == 1 ) hist_->fillHist1D("cbc0","nStubs_RECO", tdcPhase ); 
            if( chipStub == 2 ) hist_->fillHist1D("cbc1","nStubs_CHIP", tdcPhase );
            if( recoStub == 2 ) hist_->fillHist1D("cbc1","nStubs_RECO", tdcPhase );
            //}

            bool stubFound_cbc0 = false;
            bool stubFound_cbc1 = false;
            std::vector<double> stubRes; stubRes.clear();  
            std::vector<int> stubID; stubID.clear();

            std::vector<double> stubRes_cbc0; stubRes_cbc0.clear();  
            std::vector<int> stubID_cbc0; stubID_cbc0.clear();

            std::vector<double> stubRes_cbc1; stubRes_cbc1.clear();  
            std::vector<int> stubID_cbc1; stubID_cbc1.clear();
            for( unsigned int k = 0; k < xStubs_Seeds.size(); k++ )
            {

              double stub_SeedCluster = (xStubs_Seeds[k].first).first/2.0; 
              double stub_SeedCluster_pos = returnHit( stub_SeedCluster );

              double stub_MatchedCluster = (xStubs_Matched[k].first).first/2.0; 
              double stub_MatchedCluster_pos = returnHit( stub_MatchedCluster );

              double residual_Seed = stub_SeedCluster_pos - xHit_Extrapolated_Det1;
              double residual_Matched = stub_MatchedCluster_pos - xHit_Extrapolated_Det0;

              hist_->fillHist2D("det1","exTrack_Corr" , xHit_Extrapolated_Det1  , stub_SeedCluster);
              
              hist_->fillHist2D("det0","stubRes_CW" , residual_Matched , (xStubs_Matched[k].first).second );
              hist_->fillHist2D("det1","stubRes_CW" , residual_Seed , (xStubs_Seeds[k].first).second );
              
              hist_->fillHist2D("det0","stubRes_TDC" , residual_Matched , tdcPhase );
              hist_->fillHist2D("det1","stubRes_TDC" , residual_Seed , tdcPhase);
              
              hist_->fillHist2D("det0","stubRes_Chi2" , residual_Matched , chi2_Tk );
              hist_->fillHist2D("det1","stubRes_Chi2" , residual_Seed , chi2_Tk);

              hist_->fillHist1D("det0","stubRes", residual_Matched ) ;
              hist_->fillHist1D("det1","stubRes", residual_Seed ) ;
              
              if( std::fabs( residual_Seed) <= 1.0)// && std::fabs( trackHit_CorrelationLayer - stub_MatchedCluster) <= 1.0 ) 
                anaStub_1Strip_Found = true; 
              if( std::fabs( residual_Seed) <= 2.0)// && std::fabs( trackHit_CorrelationLayer - stub_MatchedCluster) <= 2.0 )
              {
                anaStub_2Strip_Found = true; 
              }
              if( std::fabs( residual_Seed) <= 3.0)// && std::fabs( trackHit_CorrelationLayer - stub_MatchedCluster) <= 3.0 ) 
              { 
                anaStub_3Strip_Found = true; 
              }
              if( std::fabs(residual_Seed) <=  acceptanceWindow_Strips*dutpitch() ) 
              //&& std::fabs(residual_Matched) <= acceptanceWindow_Strips*dutpitch() )// && std::fabs( trackHit_CorrelationLayer - stub_MatchedCluster) <= 3.0 ) 
              { 
                stubFound = true; 
                stubRes.push_back( std::fabs(residual_Seed) + std::fabs(residual_Matched) );
                stubID.push_back(k);
                
                if( stub_SeedCluster < cbcBoundary_Strips )
                { 
                  stubFound_cbc0 = true;
                  //if( xHit_Extrapolated_Det1 > 0 ) printf("Found a stub in CBC0 although the track points to CBC1! Residual  = %.3f mm\n" , residual_Seed );
                  stubRes_cbc0.push_back( std::fabs(residual_Seed) + std::fabs(residual_Matched) );
                  stubID_cbc0.push_back(k);
                }
                else
                { 
                  stubFound_cbc1 = true;
                  //if( xHit_Extrapolated_Det1 < 0 ) printf("Found a stub in CBC1 although the track points to CBC0! Residual  = %.3f mm\n" , residual_Seed );
                  stubRes_cbc1.push_back( std::fabs(residual_Seed) + std::fabs(residual_Matched) );
                  stubID_cbc1.push_back(k);
                }
              }
            }
            double width = 0;

            if( stubFound )
            {
              std::vector<double>::iterator minLocator = std::min_element(std::begin(stubRes), std::end(stubRes));
              int id = (int)(std::distance(std::begin(stubRes), minLocator));
              double stub_SeedCluster = (xStubs_Seeds[id].first).first/2.0; 
              double stub_SeedCluster_pos = returnHit( stub_SeedCluster );
              width = (xStubs_Seeds[id].first).second;
              
              hist_->fillHist1D("det1","matchedStub_mm", xHit_Extrapolated_Det1 );
              hist_->fillHist2D("det1","matchedStub_mm_2D", xHit_Extrapolated_Det1 , tdcPhase );
              //use closest cluster
              hist_->fillHist1D("det1","matchedStub_strips", stub_SeedCluster );

              if( stubFound_cbc0 )
              { 
                
                if( stub_SeedCluster == 126 )
                {
                  //misMatchedStub_mm

                }
                hist_->fillHist1D("cbc0","matchedStub_mm", xHit_Extrapolated_Det1 );
                hist_->fillHist2D("cbc0","matchedStub_mm_2D", xHit_Extrapolated_Det1 , tdcPhase );
              }
              else if( stubFound_cbc1 )
              { 
                hist_->fillHist1D("cbc1","matchedStub_mm", xHit_Extrapolated_Det1 );
                hist_->fillHist2D("cbc1","matchedStub_mm_2D", xHit_Extrapolated_Det1 , tdcPhase );
              }
            }
            if( anaStub_1Strip_Found ) hist_->fillHist1D("det1","matchedStubs_1strip", tdcPhase ) ;
            if( anaStub_2Strip_Found ) hist_->fillHist1D("det1","matchedStubs_2strip", tdcPhase ) ;
            if( anaStub_3Strip_Found ) hist_->fillHist1D("det1","matchedStubs_3strip", tdcPhase ) ;

            if( trackHit_SeedLayer < cbcBoundary_Strips )// && (x_Clusters_Det0.size() == 1  && x_Clusters_Det1.size() == 1) )
            { 
              if( stubFound_cbc0 ) hist_->fillHist1D("cbc0","nStubs_ANA",tdcPhase);

              int iBin = 0;
              if( chipStub == 0 ) iBin = 1; 
              if( chipStub == 1 || chipStub == 3 ) iBin =2 ; 
              if( chipStub == 2 )iBin = 3; 
              hist_->fillHist1D("cbc0","stubMatching", iBin);
              hist_->fillHist2D("cbc0","stubMatching_2D", iBin , xHit_Extrapolated_Det1);
              hist_->fillHist2D("cbc0","stubMatching_CW", iBin , width);
              hist_->fillHist2D("cbc0","stubMatching_TDC", iBin , tdcPhase);
            }
            else
            //else if (x_Clusters_Det0.size() == 1  && x_Clusters_Det1.size() == 1)
            { 
              
              if( stubFound_cbc1 ) hist_->fillHist1D("cbc1","nStubs_ANA",tdcPhase);
              
              int iBin = 0;
              if( chipStub == 0 )iBin = 1; 
              if( chipStub == 2 || chipStub == 3 ) iBin = 2 ; 
              if( chipStub == 1 )iBin = 3; 
              hist_->fillHist1D("cbc1","stubMatching", iBin);
              hist_->fillHist2D("cbc1","stubMatching_2D", iBin , xHit_Extrapolated_Det1);
              hist_->fillHist2D("cbc1","stubMatching_CW", iBin , width);
              hist_->fillHist2D("cbc1","stubMatching_TDC", iBin , tdcPhase);
            }

          }
          else
          {
            if( trackHit_SeedLayer < cbcBoundary_Strips )
            { 
              hist_->fillHist1D("cbc0","stubMatching", 1);
              hist_->fillHist2D("cbc0","stubMatching_2D", 1 , xHit_Extrapolated_Det1);
              hist_->fillHist2D("cbc0","stubMatching_TDC", 1 , tdcPhase);
            }
            else
            {
              hist_->fillHist1D("cbc1","stubMatching", 1);
              hist_->fillHist2D("cbc1","stubMatching_2D", 1 , xHit_Extrapolated_Det1);
              hist_->fillHist2D("cbc1","stubMatching_TDC", 1 , tdcPhase);
            }
          }
          fStubs++;
          hist_->fillHist1D("det1","tracks_CW", tdcPhase ) ;
          if( stubFound ) fStubsMatchedToTracks++; 
        } 
      }

    }
    hist_->fillHist1D("det0","fidTracks_num",cFidTrack_det0);
    hist_->fillHist1D("det1","fidTracks_num",cFidTrack_det1);

    if( cFidTracks_cbc0 != 0 ) hist_->fillHist1D("cbc0","fidTracks_num",cFidTracks_cbc0);
    if( cFidTracks_cbc1 != 0 ) hist_->fillHist1D("cbc1","fidTracks_num",cFidTracks_cbc1);
        
    //if( cFidTrack_det0 > 0 ) hist_->fillHist1D("det0","nFiducialTracks", cFidTrack_det0 ) ;
    //if( cFidTrack_det1 > 0 ) hist_->fillHist1D("det1","nFiducialTracks", cFidTrack_det1 ) ;
    //hist_->fillHist1D("det1","nFiducialTracks_Both", cFidTracks_both ) ;
    
  }
}
void EventAlignment::estimateBackground(TFitResult* fitResult , double xMin , double xMax , double& meanIntegral , double& errorIntegral )
{
  const double* pars = fitResult->GetParams();
  const double* pars_Errors = fitResult->GetErrors();

  double chi2 = fitResult->Chi2()/fitResult->Ndf() ;

  
  double N_norm = pars[0];  
  double b0 = pars[1];  
  double xOffset_b = pars[3]; 
  double sigma_b = pars[5]; 
  
  double err_N_norm = pars_Errors[0]*std::sqrt(chi2);
  double err_b0 = pars_Errors[1]*std::sqrt(chi2);
  double err_xOffset = pars_Errors[3]*std::sqrt(chi2);
  double err_sigma_b = pars_Errors[5]*std::sqrt(chi2);
  
  const int nPars = std::pow(3.,4.);
  double* integrals = new double[nPars];
  unsigned int counter =0; 
  for( int i = 0 ; i < 3 ; i++)
  {
    double N = N_norm + (1-i)*err_N_norm;
    for(  int j = 0 ; j < 3 ; j++ )
    {
      double b = b0 + (1-j)*err_b0;
      for(  int k = 0 ; k < 3 ; k++ )
      {
        double xOffset = xOffset_b + (1-k)*err_xOffset;
        for(  int l = 0 ; l < 3 ; l++ )
        {
          double sigma = sigma_b + (1-l)*err_sigma_b;
          double pars_Background[4] = { 0.0 , b*N , xOffset , sigma};
          TF1 f_res_B("f_b", Utility::normBackgroundFunction, xMin , xMax , 4);
          f_res_B.SetParameters( pars_Background );
          TString out;
          double b_integral = f_res_B.Integral(xMin,xMax);
          integrals[counter] = b_integral;
          counter++;
        }
      }
    }
  }
  meanIntegral = TMath::Mean( nPars , integrals);
  errorIntegral = TMath::RMS( nPars , integrals);
}
void EventAlignment::estimateTrackContamination(double windowSize)
{
  TString out;
  double xMin , xMax , yMin , yMax; 
  double cColPitch_fei4 = 250e-3; 
  double cRowPitch_fei4 = 50e-3; 
 
  TH1D* res_UM = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaYPos_trkfei4_UM"));   
  TH1D* res_Aligned = dynamic_cast<TH1D*>(hist_->GetHistoByName("TrackMatch","deltaYPos_trkfei4"));   
  TFitResultPtr f_Result = fitFeI4Residual( res_Aligned , res_UM , cRowPitch_fei4 ) ;
  TFitResult* fitResult = f_Result.Get();
  
  double chi2 = fitResult->Chi2()/fitResult->Ndf() ;
  hist_->fillHist1D("TrackMatch","chi2_BackgroundFit", windowSize , chi2 ) ;

  double x0 = 0.0 - windowSize*fFeI4_Window_y;
  double x1 = 0.0 + windowSize*fFeI4_Window_y;
  double meanIntegral_total = res_Aligned->Integral(res_Aligned->GetXaxis()->FindBin( x0), res_Aligned->GetXaxis()->FindBin( x1 ));
  double errIntegral_total = std::sqrt(meanIntegral_total);

  double binWidth = res_Aligned->GetXaxis()->GetBinWidth(1);
  double meanIntegral_background, errIntegral_background ;
  estimateBackground(fitResult ,  x0 , x1 , meanIntegral_background, errIntegral_background ) ;
  fBackground = meanIntegral_background/binWidth;
  fBackground_Err = errIntegral_background/binWidth;
  double relErr_background = fBackground_Err/fBackground;
  
  double contamination = fBackground/meanIntegral_total;
  double error_contamination = relErr_background*contamination;



  // const double* pars = f_Result->GetParams();
  // double N_norm = pars[0];  
  // double b0 = pars[1];  
  // double b1 = pars[2];
  // double xOffset_b = pars[3]; 
  // double xOffset_s = pars[4];
  // double sigma_b = pars[5]; 
  // double sigma_s0 = pars[6];  
  // double sigma_s1 = pars[7];  
  // double pitch = pars[8];

  // xMin = -0.5;
  // xMax =  0.5;
  // TF1 f("f_Total", Utility::residualFeI4_total , xMin , xMax , 9);
  // f.SetParameters(pars);
  // f.SetLineColor(kBlack); f.SetLineStyle(2);

  // double pars_Background[4] = { 0.0 , b0*N_norm , xOffset_b , sigma_b};
  // TF1 f_res_B("f_b", Utility::normBackgroundFunction, xMin , xMax , 4);
  // f_res_B.SetParameters( pars_Background );
  // f_res_B.SetLineColor(kCyan); f_res_B.SetFillStyle(3001); f_res_B.SetFillColor(kCyan);

  // double pars_b[8] = { N_norm , b0 , b1 , xOffset_b , xOffset_s , sigma_b , sigma_s1 , pitch };
  // TF1 f_res_B1("f_b1", Utility::residualFeI4_background, xMin , xMax , 8);
  // f_res_B1.SetParameters( pars_b );
  // f_res_B1.SetLineColor(kMagenta); f_res_B1.SetFillStyle(3003); f_res_B1.SetFillColor(kMagenta);

  // double pars_Signal_1[4] = { N_norm*( 1- (b0+b1)) , pitch , xOffset_s , sigma_s0 };
  // TF1 f_res_S1("f_s1", Utility::signalFunction, xMin , xMax , 4);
  // f_res_S1.SetParameters( pars_Signal_1 );
  // f_res_S1.SetLineColor(kRed); f_res_S1.SetFillStyle(3001); f_res_S1.SetFillColor(kRed);


  // double x0 = 0.0 - windowSize*fFeI4_Window_y;
  // double x1 = 0.0 + windowSize*fFeI4_Window_y; 
 
  // double integral_b = f_res_B.Integral(x0,x1);
  // double integral_b1 =  f_res_B1.Integral(x0,x1);
  // double integral_t = f.Integral(x0,x1);
  // double error = f.IntegralError(x0,x1, fitResult->GetParams(), (fitResult->GetCovarianceMatrix()).GetMatrixArray() );
  // double relError = std::sqrt( std::pow(error/integral_t,2.0)*2.0 );
  
    
  // double Contamination_0 = integral_b/integral_t;
  // double Contamination_1 = integral_b1/integral_t;
  // double Error_Contamination_0 = std::sqrt(2.0)*relError;
  
  fTrackContamination = contamination;
  fTrackContamination_Err = error_contamination;
  // fTrackContamination_1 = Contamination_1;
  out.Form("Contamination = %.3f , error =  %.3f perc. \n" , 100*fTrackContamination , 100*relErr_background ); 
  std::cout << out.Data();
   
}
void EventAlignment::CalculateEfficiency(double matchingWindowSize)
{
  estimateTrackContamination(matchingWindowSize);

  double hitEfficiency_det0 = (double)(fHits_Det0)/(double)(fFidTracks_Det0);
  double err_hitEff_det0 = std::sqrt( std::fabs(1.0 - hitEfficiency_det0)/fFidTracks_Det0 );
  double hitEff_det0_corr = hitEfficiency_det0/(1.0-fTrackContamination*1e-2);
  double err_hitEfficiency_det0_corr = std::sqrt( std::pow(err_hitEff_det0/hitEfficiency_det0,2.0) + std::pow(fTrackContamination_Err/fTrackContamination,2.0) );

  double clusEfficiency_det0 = (double)(fClusters_Det0)/(double)(fFidTracks_Det0);
  double err_clusEff_det0 = std::sqrt( std::fabs(1.0 - clusEfficiency_det0)/fFidTracks_Det0 );
  double clusEff_det0_corr = clusEfficiency_det0/(1.0-fTrackContamination*1e-2);
  double err_clusEfficiency_det0_corr = std::sqrt( std::pow(err_clusEff_det0/clusEfficiency_det0,2.0) + std::pow(fTrackContamination_Err/fTrackContamination,2.0) );

  double hitEfficiency_det1 = (double)(fHits_Det1)/(double)(fFidTracks_Det1);
  double err_hitEff_det1 = std::sqrt( std::fabs(1.0 - hitEfficiency_det1)/fFidTracks_Det1 );
  double hitEff_det1_corr = hitEfficiency_det1/(1.0-fTrackContamination*1e-2);
  double err_hitEfficiency_det1_corr = std::sqrt( std::pow(err_hitEff_det1/hitEfficiency_det1,2.0) + std::pow(fTrackContamination_Err/fTrackContamination,2.0) );
  
  double clusEfficiency_det1 = (double)(fClusters_Det1)/(double)(fFidTracks_Det1);
  double err_clusEff_det1 = std::sqrt( std::fabs(1.0 - clusEfficiency_det1)/fFidTracks_Det1 );
  double clusEff_det1_corr = clusEfficiency_det1/(1.0-fTrackContamination*1e-2);
  double err_clusEfficiency_det1_corr = std::sqrt( std::pow(err_clusEff_det1/clusEfficiency_det1,2.0) + std::pow(fTrackContamination_Err/fTrackContamination,2.0) );

  double stubEfficiency = (double)(fStubsMatchedToTracks)/(double)(fStubs);
  double err_stubEff = std::sqrt( std::fabs(1.0 - stubEfficiency)/fStubs );
  double stubEff_corr = stubEfficiency/(1.0-fTrackContamination*1e-2);
  double err_stubEff_corr = std::sqrt( std::pow(err_stubEff/stubEfficiency,2.0) + std::pow(fTrackContamination_Err/fTrackContamination,2.0) );
  
  char buffer[250];
  sprintf(buffer,"Found %d hits  that match %d tracks in Det0 : hit efficiency = %.3f ± %.3f [raw] , %.3f ± %.3f [corrected] \n" , fHits_Det0 ,  fFidTracks_Det0 , 100.0*hitEfficiency_det0 , 100*err_hitEff_det0 ,  100*hitEff_det0_corr , 100 *err_hitEfficiency_det0_corr );
  std::cout << buffer;
  sprintf(buffer,"Found %d hits  that match %d tracks in Det1 : hit efficiency = %.3f ± %.3f [raw] , %.3f ± %.3f [corrected] \n" , fHits_Det1 ,  fFidTracks_Det1 , 100.0*hitEfficiency_det1 , 100*err_hitEff_det1 ,  100*hitEff_det1_corr , 100 *err_hitEfficiency_det1_corr );
  std::cout << buffer;
  sprintf(buffer,"Found %d clusters  that match %d tracks in Det0 : cluster efficiency = %.3f ± %.3f [raw] , %.3f ± %.3f [corrected] \n" , fClusters_Det0 ,  fFidTracks_Det0 , 100.0*clusEfficiency_det0 , 100*err_clusEff_det0 ,  100*clusEff_det0_corr , 100 *err_clusEfficiency_det0_corr );
  std::cout << buffer;
  sprintf(buffer,"Found %d clusters  that match %d tracks in Det1 : cluster efficiency = %.3f ± %.3f [raw] , %.3f ± %.3f [corrected] \n" , fClusters_Det1 ,  fFidTracks_Det1 , 100.0*clusEfficiency_det1 , 100*err_clusEff_det1 ,  100*clusEff_det1_corr , 100 *err_clusEfficiency_det1_corr );
  std::cout << buffer;
  sprintf(buffer,"Found %d stubs  that match %d tracks in the CBC : stub efficiency = %.3f ± %.3f [raw] , %.3f ± %.3f [corrected] \n" , fStubsMatchedToTracks ,  fStubs , 100.0*stubEfficiency , 100*err_stubEff ,  100*stubEff_corr , 100 *err_stubEff_corr );
  std::cout << buffer;

  double vcth = condEv()->vcth;
  sprintf(buffer,"./Efficiency.dat");   
  std::ofstream fileStream(buffer, std::ios_base::app | std::ios_base::out);
  sprintf(buffer,"%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", fL1_mean , vcth, 100*clusEfficiency_det0 , 100*err_clusEff_det0 , 100*clusEff_det0_corr , 100*err_clusEfficiency_det0_corr, 100*clusEfficiency_det1 , 100*err_clusEff_det1 , 100*clusEff_det1_corr , 100*err_clusEfficiency_det1_corr , 100*stubEfficiency , 100*err_stubEff, 100*stubEff_corr, 100*err_stubEff_corr, 100*fTrackContamination , 100*fTrackContamination_Err); 
  fileStream << buffer; 
  fileStream.close();

}
std::vector<SingleGroup_Offset> EventAlignment::alignEvents(std::string eAlignmentFile )
{
  std::vector<SingleGroup_Offset> gOffsets; gOffsets.clear();
  //double vcth = condEv()->vcth;
  //double cwd = condEv()->cwd;
  double cwd = 5;
  double vcth = 106.;

  // find MPV for trigger latency in FeI4
  //lv1FeI4();
  fL1_mean = 13.5; fL1_sigma=1;
  std::cout << "L1 cut around " << fL1_mean << std::endl;
  //fL1_mean = 13.5; fL1_sigma = 1;
  int numberEvents_perGroup = fNumberEvents_perGroup;
  
  char buffer[250];
  int NumEntries = nEntries_ ; //nEntries_;
  std::cout << "Vcth =  " << vcth << " Correlation Window = " << cwd << std::endl;
  std::cout << "Run has " << NumEntries << " triggers." << std::endl;
  std::vector<int> groupIDs; groupIDs.clear();
  std::vector<int> Offsets ; Offsets.clear(); 
  int nGroups = 1 + nEntries_/numberEvents_perGroup; 
  double cCorrelation=10; //10
  double cMinCorrelation = 0.0;
  if( eAlignmentFile.empty() )
  {
    // now find correlation window between FeI4 and DUT 
    std::cout << "Finding correct correlation window between FeI4 and DUT\n";
    FindDetFeI4CorrelationWindow();
  
    sprintf(buffer,"FeI4_EventAlignment.dat");   
    ofstream fileStream( buffer );
    std::cout << "Have to analyse : " << nGroups << " groups." << std::endl;
    int prevOptimalOffset = 0 ; 
    double acceptanceWindow = 1.00;  // for  0° 0.25 mm is ok
    //find correct event alignment between DUT and telescope 
    for( unsigned int i = 0 ; i < nGroups ; i++ )
    {
      if( (i*numberEvents_perGroup)%10000 == 0 ) std::cout << "Event re-alignment completed for " << (i*numberEvents_perGroup) << " events." << std::endl;
      int optimalOffset = 0;
      optimalOffset = MatchDetFeI4_v2(i,prevOptimalOffset, cCorrelation , cMinCorrelation );// MatchDetFeI4(i, prevOptimalOffset , acceptanceWindow);
      if( optimalOffset == -100 || std::fabs(optimalOffset) >= fMaxEventOffset ) continue ; 

      sprintf(buffer , "%d\t%d\n" , i , optimalOffset );
      fileStream << buffer;
      if( i == 0 )
      {  
        std::cout << buffer; 
      }
      if( optimalOffset != prevOptimalOffset )
      { 
        std::cout << buffer; 
      }
      prevOptimalOffset = optimalOffset; 
      Offsets.push_back( optimalOffset ); 
      groupIDs.push_back( i );
    }
    fileStream.close();
  }
  else
  { 
    std::cout << "Reading event alignment file\n";
    std::ifstream input(eAlignmentFile);
    while (input) 
    {
      int groupID ; 
      int optimalOffset; 
      input  >> groupID >> optimalOffset; 
      if( std::fabs(optimalOffset) < fMaxEventOffset ) 
      {
        Offsets.push_back( optimalOffset ); 
        groupIDs.push_back( groupID );
      }
    }
  }
  std::cout << "Have " << groupIDs.size() << " groups after event re-alignment. [ Looked at " << nGroups << " groups]." <<  std::endl;
  //return gOffsets;

  //now, just to be safe, ignore groups between offset skips 
  SingleGroup_Offset gOffset; 
  std::cout << "Discarding groups on boundaries of missed triggers.....\n";
  nGroups = groupIDs.size(); 
  for( unsigned int i = 0 ; i < nGroups ; i++)
  {
    int groupID = groupIDs[i];
    int optimalOffset = Offsets[i];
    //if( std::fabs(optimalOffset) < 12 || std::fabs(optimalOffset) > 15 ) continue;
    if( i == 0 )
    {
      if(Offsets[i] == Offsets[i+1] )
      {
      	gOffset.first = i; 
      	gOffset.second = Offsets[i];
      	gOffsets.push_back( gOffset ) ;
        // check correlation between module and FeI4
        CorrelateDetFeI4( i, Offsets[i] );
        CorrelateFeI4Tracks(i, Offsets[i] );
      	dut_RecoStubEfficiency(i);
      }
    }
    else if( i < groupIDs.size() - 1 )
    {
      if( Offsets[i-1] == Offsets[i] && Offsets[i] == Offsets[i+1] )
      {
        gOffset.first = i; 
      	gOffset.second = Offsets[i];
      	gOffsets.push_back( gOffset ) ;
        // check correlation between module and FeI4
        CorrelateDetFeI4( i, Offsets[i] );
        CorrelateFeI4Tracks(i, Offsets[i] );
      	dut_RecoStubEfficiency(i);
      }
    }
    else
    {
      if( Offsets[i-1] == Offsets[i])
      {
        gOffset.first = i; 
      	gOffset.second = Offsets[i];
      	gOffsets.push_back( gOffset ) ;
        // check correlation between module and FeI4
        CorrelateDetFeI4( i, Offsets[i] );
        CorrelateFeI4Tracks(i, Offsets[i] );
      	dut_RecoStubEfficiency(i);
      }
    }
    int nEvents = gOffsets.size()*numberEvents_perGroup;
    if( nEvents%10000 == 0 ) std::cout << "Finished re-aligning " << nEvents << " events with telescope data.\n";
  }
  std::cout << "Originally had : " << groupIDs.size() << " groups. Now have : " << gOffsets.size() << " groups." << std::endl;
  std::cout << "Now aligning FeI4 plane and extracting track matching parameters.\n";
  doFeI4TrackFit_HighIntensity();
  //alignX(gOffsets , 1.0 );
  //return gOffsets;
  
  // align 
  double matchingWindowSize =1.0; 
  int nStrips = 4 ;
  sprintf(buffer,"Mean L1 for FeI4 = %.3f\n" , fL1_mean);
  std::cout << buffer;    
  std::vector<double> l1Cut{ fL1_mean };
  std::vector<double> Steps_Window_Y{1.0}; 
  std::vector<double> Steps_Window_X{1.0};
  std::vector<int> nStrips_eff{1};
  int Nsteps_Y = Steps_Window_Y.size();
  int Nsteps_X = Steps_Window_X.size();
  int totalStubs_found = 0; 
  int totalStubs_matched = 0; 
  for( unsigned int index_Y = 0 ; index_Y < Nsteps_Y; index_Y++)
  {
    double matchingWindowSize_Y = Steps_Window_Y[index_Y]; 
    for( unsigned int index_X = 0 ; index_X < Nsteps_X ; index_X ++)
    {
      sprintf(buffer,"./stubEfficiency.dat");   
      std::ofstream fileStream1(buffer, std::ios_base::app | std::ios_base::out);
      for( unsigned int index_L1 = 0 ; index_L1 < l1Cut.size() ; index_L1++)
      {
        
        sprintf(buffer,"--------- Using mean L1 for FeI4 = %.3f\n" , l1Cut[index_L1] );
        std::cout << buffer;    
        fL1_mean = l1Cut[index_L1]; fL1_sigma = 1;
        clearHistograms();
        clearTrackMatchhistograms();
        double matchingWindowSize_X = Steps_Window_X[index_X]; 
        double area = matchingWindowSize_X*matchingWindowSize_X*(fResidualX.second)*(fResidualY.second);
        double nHits_Det0 , nHits_Det1 , nClusters_Det0 , nClusters_Det1 , nStubs_Det; 
        double nFidTracks_Det0 , nFidTracks_Det1 , nFidStubs ; 
        nGroups = gOffsets.size(); 
        for(unsigned int i = 0 ; i < nGroups ; i++ )
        {
          int groupID = gOffsets[i].first ;
          int evOffset = gOffsets[i].second ;
          
          AllExtrapolatedHits cMatchedHits = CorrelateFeI4Tracks(groupID, evOffset, matchingWindowSize_X , matchingWindowSize_Y );
          sprintf(buffer,"Group%d : %d matched hits.\n" , i , cMatchedHits.size() );
          //std::cout << buffer ; 
          // removeDuplicateHits( cMatchedHits );

          CalculateDUTefficiency(cMatchedHits,nStrips*dutpitch() );
          int nEvents =  i*numberEvents_perGroup;
          // TODO - count number of events properly 
          if( nEvents%20000 == 0 && nEvents > 100)
          { 
            nFidStubs = fStubs - totalStubs_found ; 
            nStubs_Det = fStubsMatchedToTracks - totalStubs_matched ; 
            
            double eff = nStubs_Det/nFidStubs; 
            double err_eff = std::sqrt( std::fabs(1.0 - eff)/nFidStubs );

            sprintf(buffer,"Finished looking for DUT hits in %d events : stub efficiency = %.3f ± %.3f [raw , %d stubs from telescope ]\n" , nEvents , 100.0*eff , 100*err_eff , (int)(nFidStubs) );
            std::cout << buffer;
            sprintf(buffer,"%d\t%d\t%.3f\t%d\t%.3e\t%.3e\t%.3f\t%.3f\n", (int)(vcth) , evOffset, l1Cut[index_L1] , (int)(nEvents) , (double)(nFidStubs), (double)nStubs_Det , 100*eff , 100*err_eff );
            fileStream1 << buffer ; 

            totalStubs_found = fStubs; 
            totalStubs_matched = fStubsMatchedToTracks;
            hist_->fillHistProfile("det1","stubEff_time", (int)(nEvents/20000) , eff ) ;
          }
          sprintf(buffer,"Calculated efficiency for %d events\n" , nEvents);
          if( nEvents%1000 == 0 ) std::cout << buffer;

        }
        CalculateEfficiency(matchingWindowSize_Y);
      }
      fileStream1.close();
        
    }
  }

  return gOffsets;
}
