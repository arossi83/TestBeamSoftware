/*!
        \file                Histogrammer.cc
        \brief               Class provides user to book, fill histograms
        \author              Suvankar Roy Chowdhury
        \date                05/07/16
        Support :            mail to : suvankar.roy.chowdhury@cern.ch
*/
#include "Histogrammer.h"
#include <climits>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "TLorentzVector.h"
#include "TFile.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

Histogrammer::Histogrammer(std::string& outFile) {
  fout_ = new TFile(TString(outFile),"RECREATE");
  isFileopen_ = true;
}

void Histogrammer::bookEventHistograms() {
  fout_->cd();
  fout_->mkdir("EventInfo");
  fout_->cd("EventInfo");
  new TH1I("nevents", "#Events", 10000001, -0.5, 10000000.5);
  new TH1I("dutAngle", "DUT Angle;DUTAngle;#Events", 1805, -90.25, 90.25);  
  new TH1I("alignAngle", "Alignment Angle;DUTAngle;#Events", 18050, -90.25, 90.25);  
  new TH1I("hvSettings", "High Voltage settings;HV;#Events", 1000,-0.5,999.5);
  new TH1I("vcth", "Vcth value;vcth;#Events", 200,-0.5,199.5);
  new TH1I("offset", ";offset;#Events", 200,-0.5,199.5);
  new TH1I("window", ";window;#Events", 200,-0.5,199.5);
  new TH1I("tilt", ";tilt;#Events", 200,-0.5,199.5);
  new TH1I("condData", ";condData;#Events", 20,-0.5,19.5);
  new TH1I("tdcPhase",";tdc;#Events",17,-0.5,16.5);
  new TH1I("isPeriodic",";isPeriodic;#Events",2,-0.5,1.5);
  new TH1I("isGoodFlag",";isGood;#Events",2,-0.5,1.5);
}

void Histogrammer::bookDUTHistograms(std::string det) 
{
  int nStrips = 254; 

#ifdef NOV_15
  nStrips = nStrips*4;
#endif

  TString d(det);
  fout_->cd();
  fout_->mkdir(d);
  fout_->cd(d);
  new TH2I("hitmapfull",d + " hitmap;strip no.;#Events",1016,-0.5,1015.5,2,-0.5,1.5);
  bookDUTHistoForColumn(d,"C0");
  //bookDUTHistoForColumn(d,"C1");

  new TH1D("nFiducialTracks","Number of Fiducial Tracks; Number of Fiducial Tracks; Number of Events", 100 ,-0.5,99.5); 
  new TH1D("nFiducialTracks_Both","Number of Fiducial Tracks; Number of Fiducial Tracks; Number of Events", 100 ,-0.5,99.5); 
 
  double xMin = -1.0*0.09*1.5*nStrips; 
  double xMax = 1.0*0.09*1.5*nStrips;
  double xBin = 0.09/1.0;  
  double xOffset = 0.0;//xBin/2.0;
 

  // DUT histograms 
  new TH1D("stubsFoundCbc" , "Number of Events" , 4 , -0.5 , 3.5 );
  new TH1D("stubsFoundReco" , "Number of Events" , 4 , -0.5 , 3.5 );
  TH1D* h = dynamic_cast<TH1D*>(Utility::getHist1D("stubsFoundCbc"));
  h->GetXaxis()->SetBinLabel(1,"none");
  h->GetXaxis()->SetBinLabel(2,"CBC_0");
  h->GetXaxis()->SetBinLabel(3,"CBC_1");
  h->GetXaxis()->SetBinLabel(4,"both");
  // h = dynamic_cast<TH1D*>(Utility::getHist1D("stubsFoundReco"));
  // h->GetXaxis()->SetBinLabel(1,"CBC_0");
  // h->GetXaxis()->SetBinLabel(2,"CBC_1");
  // h->GetXaxis()->SetBinLabel(3,"Both");
  new TProfile("stubEff_time" , "Stub Efficiency as a function of Time Elapsed; Number of Events [mod 20k]; Stub Efficiency" , 100 , -0.5 , 99.5 );
  
  new TH2D("checkEvAlignment","; Stub ID; Offset", 100000 ,-0.5, 100000-0.5 , 5.0 , -2.0-0.5 , 3.0 - 0.5 );
  new TH2D("exTrack_Corr","; Seed from extrapolated track; Seed from CBC stub [as per ANALYSIS]", nStrips+1,-0.5, nStrips - 0.5 , nStrips+1,-0.5, nStrips - 0.5 );
  
  new TH1D("cbcStub" , "Number of stubs [as per stub OR from CBC] that match tracks [ 1 strip window around impact point in each layer]; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );
  new TH1D("recoStub" , "Number of stubs [as per stub OR from RECO] that match tracks [ 1 strip window around impact point in each layer]; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );

  new TH1D("tracks_CW" , "Number of tracks with hits in the DUT that match the correlation window; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );
  new TH1D("matchedStubs_1strip" , "Number of stubs [as per ANALYSIS] that match tracks [ 1 strip window around impact point in each layer]; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );
  new TH1D("matchedStubs_2strip" , "Number of stubs [as per ANALYSIS] that match tracks [ 2 strip window around impact point in each layer]; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );
  new TH1D("matchedStubs_3strip" , "Number of stubs [as per ANALYSIS] that match tracks [ 3 strip window around impact point in each layer]; TDC Phase ; Number of Events" , 17, -0.5, 16.5 );
 
  new TH1D("recoStubSeedWidth" , "Seed Cluster Width as per RECO; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  new TH1D("anaStubSeedWidth" , "Seed Cluster Width as per ANALYSIS; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  new TH1D("recoStubMatchWidth" , "Matched Cluster Width as per RECO; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  new TH1D("anaStubMatchWidth" , "Matched Cluster Width as per ANALYSIS; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  
  new TH1D("recoStubBend" , "Bend information as per RECO; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  new TH1D("anaStubBend" , "Bend information as per ANALYSIS; Stub Bend [strips]; Number of Events" , 15 , -0.5 , 15 - 0.5);
  
  new TH1D("recoStubSeeds","Stub Seeds as per RECO; Seed Cluster Center [strips]; Number of Events", nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("anaStubSeeds","Stub Seeds as per ANALYSIS; Seed Cluster Center [strips]; Number of Events",  nStrips+1,-0.5, nStrips - 0.5);
  
  new TH2D("dutClusterWidth_TDC","Cluster width on DUT with TDC Phase [x]; Cluster Width [strips]; TDC Phase", 15 , -0.5 , 15 - 0.5 , 17, -0.5, 16.5 );
  new TH2D("anaStub_TDC","Stub Seeds with TDC Phase as per RECO [x]; Strip Number; TDC Phase", nStrips+1,-0.5, nStrips - 0.5, 17, -0.5, 16.5 );
  new TH2D("recoStub_TDC","Stub Seeds with TDC Phase as per ANALYSIS [x];  Strip Number; TDC Phase", nStrips+1,-0.5, nStrips - 0.5 , 17, -0.5, 16.5 );
  
  int nBins = (xMax-xMin)/((0.09/(20.0))) ; 
  new TH1D("fidTracks_Bends", ";Track Bend [Strips]; Number of Events" , 22/0.05 , -11.0-0.5, 11.0-0.5 );

  new TH1D("fidTracks_num", ";# Fiducial Tracks Per Event; Number of Events" , 100 , 0.0-0.5, 100.0-0.5 );
  new TH1D("fidTracks_chi2", ";#Chi^{2}_{Tracks}/nDOF_{Tracks}; Number of Tracks" , (10.0)/0.05 , 0.0-0.1, 10.0-0.1 );
  
  new TH2D("fidTracks_impact_mm_2D", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("fidTracksMatched_impact_mm_2D", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_mm_All_2D", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax, nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_mm_2D", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_nextStrip_mm_2D", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  
  new TH2D("fidTracks_impact_mm_2D_CW2", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("fidTracksMatched_impact_mm_2D_CW2", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_mm_2D_CW2", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_nextStrip_mm_2D_CW2", ";x_{Track} [mm]; n_{Strip}" , (xMax-xMin)/(0.1*xBin) , xMin , xMax,  nStrips+1,-0.5, nStrips - 0.5);
  

 
  new TH1D("fidTracks_impact_mm", ";x_{Track} [mm]; Number of Events" ,(xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset);
  new TH1D("fidTracks_impact_strips", ";x_{Track} [strips]; Number of Events" , nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("matchedCluster_mm", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset);
  new TH1D("matchedCluster_strips", ";x_{Cluster} [strips]; Number of Events" , nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("matchedCluster_CW", ";x_{Track} [mm]; Cluster Width [Strips]" ,  nStrips+1,-0.5, nStrips - 0.5 , 10 , -0.5 , 9.5);
  new TProfile("matchedClusterWidth_mm", ";x_{Track} [mm]; Average Cluster Width" ,(xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset);
  new TH1D("candidateStubs_impact_mm", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset);
  new TH2D("candidateStubs_impact_mm_2D", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset, 20 , -0.5 , 19.5);
  new TH1D("candidateStubs_impact_strips", ";x_{Track} [strips]; Number of Events" , nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("matchedStub_mm", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin-xOffset , xMax-xOffset);
  new TH2D("matchedStub_mm_2D", ";x_{Track} [mm]; TDC Phase" , (xMax-xMin)/xBin, xMin , xMax , 20 , -0.5 , 19.5);
  
  new TH1D("matchedStub_strips", ";x_{Cluster} [strips]; Number of Events" , nStrips+1,-0.5, nStrips - 0.5);


  new TH1D("hitRes","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]", (xMax-xMin)/(0.25*xBin), xMin , xMax );
  new TH1D("clusterRes","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]", (xMax-xMin)/(0.25*xBin), xMin , xMax );
  
  new TH2D("clusterRes_Chi2","Stub Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", (xMax-xMin)/(0.1*xBin), xMin , xMax  , (10.0)/0.05 , 0.0-0.1, 10.0-0.1 );
  new TH2D("clusterRes_TDC","Stub Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; TDC Phase", (xMax-xMin)/(0.1*xBin), xMin , xMax , 17, -0.5, 16.5 );
  new TH2D("clusterRes_CW","Stub Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Cluster Width [strips]", (xMax-xMin)/(0.1*xBin), xMin , xMax  , 17, -0.5, 16.5 );
  
  new TH1D("stubRes","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]", (xMax-xMin)/(0.1*xBin), xMin , xMax );
  new TH2D("stubRes_Chi2","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", (xMax-xMin)/(0.1*xBin), xMin , xMax , (10.0)/0.05 , 0.0-0.1, 10.0-0.1 );
  new TH2D("stubRes_TDC","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]; TDC Phase", (xMax-xMin)/(0.1*xBin), xMin , xMax , 17, -0.5, 16.5 );
  new TH2D("stubRes_CW","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]; Cluster Width [strips]", (xMax-xMin)/(0.1*xBin), xMin , xMax , 17, -0.5, 16.5 );
 

  new TH1D("stubRecoMismatch","Percentage of Events with a Stub Reco mismatch; Percentage of Events; Counts",  10/0.1, -0.5 , 10-0.5);
  new TH1D("stubMismatches_n","Number of Stubs for Mis-matched Events; Number of Stubs; Number of Events",  10 , -0.5 , 10-0.5);
  new TH1D("clusterWidth","Cluster Width for Det; Cluster Width [Strips]; Number of Events",  20, -0.5 ,20-0.5);
  new TH1D("clusterWidth_Raw","Cluster Width for Det; Cluster Width [Strips]; Number of Events",  20, -0.5 ,20-0.5);
  new TH1D("seedClusterWidth","Cluster Width for Stub Seed; Cluster Width [Strips]; Number of Events",  20, -0.5 ,20-0.5);
  
  new TH1D("stubSeeds","Stub Seeds as per analysis; Seed Cluster Center [half strips]; Number of Events",  510, -0.5 ,510-0.5);
  new TProfile("stubSeedsProf","Stub Seeds as per analysis; Seed Cluster Center [half strips]; Average Seed Cluster Width [strips]", 510, -0.5 ,510-0.5);
  
  
  new TH1D("stubMismatches","Number of Stub Mismatches [cbc/reco] per TDC phase; TDC Phase; Number of Mismatched",  17, -0.5, 16.5 );
  new TH1D("stubMatches","Number of Stub Matches [cbc/reco] per TDC phase; TDC Phase; Number of Matched Events",  17, -0.5, 16.5 );
  
  new TH1D("stubMismatches_w","Number of Stub Mismatches [cbc/reco] per maximum cluster width in event; Cluster Width; Number of Mismatched Events",  20, -0.5, 19.5 );
  new TH1D("stubMatches_w","Number of Stub Matches [cbc/reco] per maximum cluster width in event; Cluster Width; Number of Matched Events",  20, -0.5, 19.5 );
  
  new TH1D("dutClusters_cbc","Cluster Occupancy on DUT;Strip Number; Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("dutClusters_reco","Cluster Occupancy on DUT;Strip Number; Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("dutHits","Hit Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("cbcStubs","CBC Stub Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TH1D("recoStubs","RECO Stub Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);

  new TH2D("dutHits_TDC","Hit Occupancy on DUT with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5, 17, -0.5, 16.5 );
  new TH2D("cbcStubs_TDC","CBC Stub Occupancy on DUT with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5 , 17, -0.5, 16.5 );
    

  new TH1D("Clusters_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);
  new TH1D("MatchedClusters_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);
  new TH1D("UnMatchedClusters_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);

  new TH1D("Stubs_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);
  new TH1D("MatchedStubs_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);
  new TH1D("UnMatchedStubs_Width","Distribution of cluster widths; Cluster Width; Number of Events;" , 100 , -0.5 ,99.5);
  
  new TH1D("nClusters_noTrks","Number of DUT custers per event with no tracks; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
  new TH1D("nClusters_noFeI4","Number of DUT custers per event with no FeI4 hit; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
  new TH1D("nClusters_wTrks","Number of DUT custers per event with reconstructed tracks; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
  new TH1D("nClusters_wFeI4","Number of DUT custers per event with an FeI4 hit; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
 
  nBins = (xMax-xMin)/((0.09/1)) ; 
  new TH1D("Clusters_noTrks","Cluster Position (no reconstructed track); x_{DUT}; Number of Events [a.u]", nBins, xMin , xMax );
  new TH1D("Clusters_noFeI4","Cluster Position (no FeI4 hit); x_{DUT}; Number of Events [a.u]",nBins, xMin , xMax );
  new TH1D("Clusters_wTrks","Cluster Position (with reconstructed tracks); x_{DUT}; Number of Events [a.u]",nBins,  xMin , xMax );
  new TH1D("Clusters_wFeI4","Cluster Position (with FeI4 hit); x_{DUT}; Number of Events [a.u]",nBins,  xMin , xMax );
 
  nBins = (xMax-xMin)/((0.09/(5.0*1.5))) ; 
  new TH1D("deltaXhits","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  new TH1D("deltaXhits_SingleHitEvents","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  new TH1D("deltaXhits_SingleHitSingleTrackEvents","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; Number of Events}",nBins, xMin , xMax );
  
  new TH1D("deltaXclusters","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  new TH1D("deltaXclusters_SingleTrkEvents","Cluster Residuals on DUT [Single Matched Track]; x_{Cluster,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  new TH1D("deltaXclusters_SingleClstrEvents","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  new TH1D("deltaXclusters_SingleClstrSingleTrkEvents","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
 
  new TH1D("deltaXstubs","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]; Number of Events}",nBins,  xMin , xMax );
  
  new TH2D("deltaX2D_SingleHits","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; Strip number", nBins, xMin , xMax , nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("deltaX2D_Hits","Hit Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; Strip number", nBins, xMin , xMax , nStrips+1,-0.5, nStrips - 0.5);
  
  new TH2D("deltaX2D_SingleClusters","Cluster Residuals on DUT [w_{Cluster} == 1]; x_{Cluster,DUT} - x_{Track} [mm]; Strip number", nBins, xMin , xMax , nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("deltaX2D_Clusters","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Strip number", nBins,xMin , xMax , nStrips+1,-0.5, nStrips - 0.5);
  new TH2D("deltaX2D_Stubs","Stub Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Strip number", nBins,xMin , xMax , nStrips+1,-0.5, nStrips - 0.5);
  
  new TH2D("deltaX2D_ClustersW","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; Cluster Width", nBins,xMin , xMax , 100 , -0.5 ,99.5);
  new TH2D("deltaX2D_StubsW","Stub Residuals on DUT; x_{Stub,DUT} - x_{Track} [mm]; Largest Cluster in Stub Width", nBins,xMin , xMax , 100 , -0.5 ,99.5);
  

  new TH2D("deltaXhits_Trks","Cluster Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", nBins, xMin , xMax , (5.0)/0.05 , 0.0-0.1, 5.0-0.1); 
  new TH2D("deltaXclusters_Trks","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", nBins, xMin , xMax , (5.0)/0.05 , 0.0-0.1, 5.0-0.1); 
  
  new TH2D("nMatchedTracks_TDC","Number of Tracks for events with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5 , 17, -0.5, 16.5 );
  new TH2D("nMatchedClusters_TDC","Number of Clusters with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5, 17, -0.5, 16.5 );
  new TH2D("nUnMatchedClusters_TDC","Number of Clusters with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5 , 17, -0.5, 16.5 );
  new TH2D("nMatchedStubs_TDC","Number of Matched Stubs with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5, 17, -0.5, 16.5 );
  new TH2D("nUnMatchedStubs_TDC","Number of Un-matched Stubs with TDC Phase [x]; Strip number ; TDC Phase", nStrips+1,-0.5, nStrips - 0.5 , 17, -0.5, 16.5 );

  
  //new TH2D("deltaX_Single_hits_Trks","Cluster Residuals on DUT; x_{Hit,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", nBins, xMin , xMax , (5.0)/0.1 , 0.0-0.5, 5.0-0.5); 
  //new TH2D("deltaX_Single_clusters_Trks","Cluster Residuals on DUT; x_{Cluster,DUT} - x_{Track} [mm]; #Chi^{2}_{Tracks}/nDOF_{Tracks}", nBins, xMin , xMax , (5.0)/0.1 , 0.0-0.5, 5.0-0.5); 
  new TProfile("trkOccupancy","Noise Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("noiseOccupancy","Noise Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("hitOccupancy","Cluster Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("noHitOccupancy","Cluster Occupancy on DUT;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("nTracks","Number of Fiducial Tracks;Strip Number;Noise Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("stubOccupancy","Stub Occupancy;Strip Number;Stub Occupancy", nStrips+1,-0.5, nStrips - 0.5);
  new TProfile("noiseStubOccupancy","Stub Occupancy;Strip Number;Stub Occupancy", nStrips+1,-0.5, nStrips - 0.5);
 
  nBins = (xMax-xMin)/((0.09/(5.0*1.0))) ; 
  new TProfile("trkOccupancy_xTrk","Noise Occupancy on DUT; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("noiseOccupancy_xTrk","Noise Occupancy on DUT; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("hitOccupancy_xTrk","Cluster Occupancy on DUT; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("noHitOccupancy_xTrk","Cluster Occupancy on DUT; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("nTracks_xTrk","Number of Fiducial Tracks; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("stubOccupancy_xTrk","Stub Occupancy;S X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  new TProfile("noiseStubOccupancy_xTrk","Stub Occupancy; X_{Track Impact} [mm]; Occupancy", nBins,  xMin , xMax );
  
  new TProfile("clusterEfficiency","Measured cluster Efficiency; X_{Track Impact} [mm]; Measured Cluster Efficiency [%]", nBins,  xMin , xMax );
  new TProfile("clusterEfficiency_c","Corrected cluster Efficiency; X_{Track Impact} [mm]; Measured Cluster Efficiency [%]", nBins,  xMin , xMax );
  new TProfile("clusterEfficiency_perStrip","Cluster Efficiency per Strip; Strip Number; Cluster Efficiency [%]",  nStrips+1,-0.5, nStrips - 0.5 );
  
  new TProfile("stubEfficiency","Measured stub Efficiency; X_{Track Impact} [mm]; Measured Cluster Efficiency [%]", nBins,  xMin , xMax );
  new TProfile("stubEfficiency_c","Corrected stub Efficiency; X_{Track Impact} [mm]; Measured Cluster Efficiency [%]", nBins,  xMin , xMax );
  new TProfile("stubEfficiency_perStrip","Stub Efficiency per Strip; Strip Number; Cluster Efficiency [%]",  nStrips+1,-0.5, nStrips - 0.5 );
  

  new TH1D("noisePerStrip","Noise Occupancy on DUT; Strip Number; Number of Events" , nStrips+1,-0.5, nStrips - 0.5 );
}

void Histogrammer::bookDUTHistoForColumn(TString& d, TString c) 
{
  fout_->cd(d);    
  new TH1I("chsize" + c ,"dut0 channel occupancy per event" + c + ";#Channels;#Events",51,-0.5,50.5);
  new TH1I("hitmap" + c,"dut0 hitmap " + c + ";strip no.;#Events",1016,-0.5,1015.5);
  new TH1D("ncluster" + c,"#cluster dut0 " + c + ";#Clusters;#Events",51,-0.5,50.5);
  new TH1I("clusterWidth" + c,"dut0 clusterWidth " + c +";#ClusterWidth;#Events",128,-0.5,127.5);
  new TH1D("clusterPos" + c,"dut0 clusterPos " + c + ";Strip Number;#Events",1016,-0.5,1015.5);
  new TProfile("clusterWidthVsPosProf" + c,"dut0 clusterWidthVsPos " + c + ";Strip Number;Cluster Width",1016,-0.5,1015.5);
  new TH2D("clusterWidthVsPos2D" + c ,"dut0 clusterWidthVsPos " + c + ";Strip Number;Cluster Width",1016,-0.5,1015.5, 20,-0.5,19.5);
   new TH2D("nhitvsnclus" + c,"#Clusters vs #Hits;#Hits;#Clusters",50,-0.5,49.5, 50,-0.5,49.5);
   new TH2D("nhitvsHitClusPosDiff" + c,"Cluster-Hit MinPosDiff vs #Hits;#Hits;#Pos Diff",50,-0.5,49.5, 256,-0.5,255.5);
   new TH2D("propertyVsTDC2D" + c, "Hit Property vs TDC " + d + c + ";TDC;",17,-0.5,16.5, 10, 0.5, 10.5);
}

void Histogrammer::bookStubHistograms() 
{
  fout_->cd();
  fout_->mkdir("StubInfo");
  fout_->cd("StubInfo");

  new TH1I("cbcStubWord","Stub Bit from CBC",16,-0.5,15.5);
  new TH1I("recoStubWord","Stub Bit from offline CBC logic emulation",16,-0.5,15.5);
  new TH1I("nstubsFromCBCSword","Total number of stubs from CBC stub word",20,-.5,19.5);
  new TH1I("nstubsFromRecoSword","Total number of stubs from Reco Stub word",20,-.5,19.5); 
  new TH1I("nstubsFromReco","Total number of stubs from Reconstruction",20,-.5,19.5); 
  new TH1I("stubMatch","Matching between CBC Stub and Reco Stub",4,0.5,4.5);
  new TH1I("nstubsdiffSword","#StubsRecoStubword - #StubsfromStubWord",20,-0.5,19.5);
  new TH1I("nstubsdiff","#StubsReco - #StubsfromStubWord",20,-0.5,19.5);
  bookStubHistoForColumn("C0");
  //bookStubHistoForColumn("C1");
}

void Histogrammer::bookStubHistoForColumn(TString c) {
  fout_->cd("StubInfo");
  new TH1I("nstubReco" + c,"Number of stubs for " + c + " from offline reconstruction;#stubs;Events",20,-0.5,19.5);
}

void Histogrammer::bookCorrelationHistograms() {
  fout_->mkdir("Correlation");
  fout_->cd("Correlation");
  bookCorrelationHistoForColumn("C0");
  //bookCorrelationHistoForColumn("C1");
}
    
void Histogrammer::bookCorrelationHistoForColumn(TString c) {
  new TH1D("cor_hit" + c,"Sensor Hit Correlation " + c,4,0.5,4.5);
  new TH1I("nclusterdiff" + c,"Difference in #clusters between dut0 and dut1() for " + c + ";#cluster_{det0} - #cluster_  {det1_};Events",20,-0.5,19.5);
  //new TH1D("nhitdiff" + c,"#Hits det0 - #Hits det1 " + c,100,-.5,99.5);
}


void Histogrammer::bookTrackMatchHistograms() 
{
  int nStrips = 254; 
#ifdef NOV_15 
  nStrips = nStrips*4;
#endif
  fout_->mkdir("TrackMatch");
  fout_->cd("TrackMatch");
  new TH1D("nTrackParams","#Tracks Telescope;#tracks;#events",30,-0.5,29.5);
  new TH1D("nTracks","Number Tracks within Correlation Window Per Cluster (after duplicate removal);#tracks;#events",100,-0.5,100.5);
  new TH1D("nTracks_M","Number Tracks within Correlation Window Per Cluster (after duplicate removal);#tracks;#events",100,-0.5,100.5);
  new TH1D("nTracks2","Number Tracks Per;#tracks;#events",100,-0.5,100.5);
  new TH1D("nUniqueTracks","Number of Tracks Per Event;Number of Tracks; Counts",100,-0.5,100.5);
  new TH1D("nTrackParamsNodupl","#Tracks Telescope after duplicate removal;#tracks;#events",30,-0.5,29.5);

  new TH1F("hTkChi2", "Chi Squared Tracks", (100.0)/0.05 , -0.5, 99.5);
  new TH1F("hTkChi2_SingleTks", "Chi Squared Tracks", (100.0)/0.05 , -0.5, 99.5);
  new TH1F("hTkChi2_MultiTks", "Chi Squared Tracks", (100.0)/0.05 , -0.5, 99.5);
  new TH1F("hTkNdof", "NDOF Tracks", (100.0)/1.0 , -0.5, 99.5);
  

  new TH1I("nPxlClusterFei4","Number of Pixels per cluster in Fei4 plane;Pixels in Cluster;Clusters",30,-0.5,29.5);
  new TH1I("nClustersFei4","Number of FeI4 Clusters per event;Percentage;Count",(100) ,-0.5,99.5);
  new TH1I("nClustersFei4_M","Number of FeI4 Clusters per event;Percentage;Count",(100) ,-0.5,99.5);
  new TH1I("hToT_M","Total Time over Threshold of Matched clusters;Tot of Matched Cluster;Count", 200 ,-0.5,199.5);
  new TH1I("hToT_UM","Total Time over Threshold of Un-Matched clusters;Tot of Matched Cluster;Count", 200 ,-0.5,200.5);
  new TH1I("lv1_Cluster","Number of FeI4 Clusters per event;Percentage;Count",(100) ,-0.5,99.5);
  new TH1D("nTracks_DUT","Number of Fiducial Tracks;Number of Tracks; Number of Events",(100) ,-0.5,99.5);
  new TH1D("nTracks_SHevents_DUT","Number of Fiducial Tracks;Number of Tracks; Number of Events",(100) ,-0.5,99.5);
  new TH1D("nTracks_SCevents_DUT","Number of Fiducial Tracks;Number of Tracks; Number of Events",(100) ,-0.5,99.5);
  

  new TH1D("trackImpact_Strips_Det0","; Track Impact [strips]; Number of Events", nStrips +1 , -0.5 , nStrips -0.5);
  new TH1D("trackImpact_Strips_Det1","; Track Impact [strips]; Number of Events", nStrips +1 , -0.5 , nStrips -0.5);

  new TH1I("nMactchedClusters_1","Number of Matched Clusters [<= number of clusters in event];Number of Clusters Matched; Number of Events;",(5) ,-0.5,4.5);
  new TH1I("nMactchedClusters_2","Number of Matched Clusters [<= number of clusters in event];Number of Clusters Matched; Number of Events;",(5) ,-0.5,4.5);
  new TH1I("nMactchedClusters_3","Number of Matched Clusters [<= number of clusters in event];Number of Clusters Matched; Number of Events;",(5) ,-0.5,4.5);

  new TH2D("nTrkCluster_FeI4","Trk/Cluster Correlation; Number of FeI4 Clusters; Number of Tracks",(100) ,-0.5,99.5 , (100) ,-0.5,99.5);
  
   
  new TH1I("nClusters_Fei4","Number of FeI4 Clusters per event;Number of Clusters;Count",(100) ,-0.5,99.5);
  new TH1I("nClusters_Det0","Number of Det0 Clusters per event;Number of Clusters;Count",(100) ,-0.5,99.5);
  new TH1I("nClusters_Det1","Number of Det1 Clusters per event;Number of Clusters;Count",(100) ,-0.5,99.5);
  
  new TH1I("wClusters_Det0","Width of Det0 Clusters;Cluster Width;Count",(100) ,-0.5,99.5);
  new TH1I("wClusters_Det1","Width of Det1 Clusters;Cluster Width;Count",(100) ,-0.5,99.5);
  

  new TH1D("hposx","Xpos of Matched Cluster on FeI4 plane; x_{FeI4} [mm]; Events [a.u]",40/(250e-3),-20.,20.);
  new TH1D("hposy","Ypos of Matched Cluster on FeI4 plane; y_{FeI4} [mm]; Events [a.u]",40/(50e-3),-20.,20.);
  new TH1D("hposx_dut0","Xpos of hit on DUT0 plane; x_{DUT} [mm]; Events [a.u]",90/(90e-3),-45.,45.);
  new TH1D("hposx_dut1","Xpos of hit on DUT1 plane; x_{DUT} [mm]; Events [a.u]",90/(90e-3),-45.,45.);
  new TH1D("hposCluster_dut0","Xpos of Cluster in Det0; x_{Cluster} [mm]; Count [a.u]",(36/(0.09/10)), -18.0, 18.0);
  new TH1D("hposCluster_dut1","Xpos of Cluster in Det1; x_{Cluster} [mm]; Count [a.u]",(36/(0.09/10)), -18.0, 18.0);

  TH1D* h1d = dynamic_cast<TH1D*>(Utility::getHist1D("hposx_dut0"));
  h1d->SetLineColor(kRed);
  h1d->GetXaxis()->SetRangeUser(-12,12);
  h1d = dynamic_cast<TH1D*>(Utility::getHist1D("hposx_dut1"));
  h1d->SetLineColor(kRed);
  h1d->GetXaxis()->SetRangeUser(-12,12);

  new TH1D("hposxTk","Xpos of matched track on FeI4 plane; x_{Track} [mm]; Events [a.u]",40/(250e-3),-20.,20.);
  new TH1D("hposyTk","Ypos of matched track on FeI4 plane; y_{Track} [mm]; Events [a.u]",40/(50e-3),-20.,20.);
  new TH1D("hposxTk_A","Xpos of matched track on FeI4 plane (after alignment); x_{Track} [mm]; Events [a.u]",40/(250e-3),-20.,20.);
  new TH1D("hposyTk_A","Ypos of matched track on FeI4 plane (after alignment); y_{Track} [mm]; Events [a.u]",40/(50e-3),-20.,20.);
  new TH1D("hposxTk_dut0","Xpos of extrapolated track at DUT0 plane before alignment(#fid trk=1); x_{Track} [mm]; Events [a.u]",40/(50e-3),-20.,20.);
  new TH1D("hposxTk_dut1","Xpos of extrapolated track at DUT0 plane before alignment(#fid trk=1); x_{Track} [mm]; Events [a.u]",40/(50e-3),-20.,20.);
  new TH1D("hminposClsDUT0","Xpos of the cluster hit at DUT0 plane with min dist from track(#fid trk=1)",100,-20.,20.);
  new TH1D("hminposClsDUT1","Xpos of the cluster hit at DUT1 plane with min dist from track(#fid trk=1)",100,-20.,20.);
  new TH1D("hminposStub","Xpos of the stub hit with min dist from trk(#fid trk=1)",100,-20.,20.);
   h1d = dynamic_cast<TH1D*>(Utility::getHist1D("hposxTk_dut0"));
  h1d->GetXaxis()->SetRangeUser(-12,12);
  h1d->SetLineColor(kBlue);
  h1d->SetFillColor(kBlue);
  h1d->SetFillStyle(3002);
  h1d->GetXaxis()->SetRangeUser(-12,12);
  h1d = dynamic_cast<TH1D*>(Utility::getHist1D("hposxTk_dut1"));
  h1d->GetXaxis()->SetRangeUser(-12,12);
  h1d->SetLineColor(kBlue);
  h1d->SetFillColor(kBlue);
  h1d->SetFillStyle(3002);
  

  new TH2D("TrkCluster_Det0","Trk/Cluster Correlation;x_{Track} [mm]; x_{Cluster} [mm]",(36/(0.09/10)), -18.0, 18.0, (36/(0.09/10)), -18.0, 18.0 );
  new TH2D("TrkCluster_Det1","Trk/Cluster Correlation;x_{Track} [mm]; x_{Cluster} [mm]",(36/(0.09/10)), -18.0, 18.0, (36/(0.09/10)), -18.0, 18.0 );
  new TH2D("minclsTrkPoscorrD0","Cluster xTrk Pos Correlation;trk;cluster",nStrips+1,-0.5,nStrips - 0.5, nStrips+1,-0.5,nStrips - 0.5 );
  new TH2D("minclsTrkPoscorrD1","Cluster xTrk Pos Correlation;trk;cluster",nStrips+1,-0.5,nStrips - 0.5, nStrips+1,-0.5,nStrips - 0.5 );
  TH2D* h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minclsTrkPoscorrD0"));
  h2d->SetOption("colz");
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minclsTrkPoscorrD1"));
  h2d->SetOption("colz");


  //new TH1D("residualDUT0multitrkfidNodupl","ClusterResidual at DUT0 plane(fiducial)(#trk>1, no duplicate tracks)",400,-20.,20.);
  //new TH1D("residualDUT1multitrkfidNodupl","ClusterResidual at DUT1 plane(fiducial)(#trk>1, no duplicate tracks)",400,-20.,20.);
  new TH1D("minresidualDUT0_1trkfid","ClusterResidual at DUT0 plane(fiducial)(#trk=1)",1000,-10.,10.);
  new TH1D("minresidualDUT1_1trkfid","ClusterResidual at DUT1 plane(fiducial)(#trk=1)",1000,-10.,10.);
  new TH1D("clswidthDUT0_1trkfid","ClusterWidth(cluster matched to track) at DUT0 plane(fiducial)(#trk=1)",50,-0.5,49.5);
  new TH1D("clswidthDUT1_1trkfid","ClusterWidth(cluster matched to track) at DUT1 plane(fiducial)(#trk=1)",50,-0.5,49.5);


  new TH2D("minstubTrkPoscorrD1_all","Closest-Stub xTrk Pos Correlation;trk;stub",1016,-0.5,1015.5, 1016,-0.5,1015.5);
  new TH2D("minstubTrkPoscorrD1_matched","Closest-Stub xTrk Pos Correlation(matched);trk;stub",1016,-0.5,1015.5, 1016,-0.5,1015.5);
  new TH2D("minstubTrkPoscorrD1_CBC","Closest-Stub xTrk Pos Correlation(matched);trk;stub",1016,-0.5,1015.5, 1016,-0.5,1015.5);
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minstubTrkPoscorrD1_all"));
  //new TH1D("sresidualC0multitrkfidNodupl","Stub Residual at DUT0 plane(fiducial)(#trk>1, no duplicate tracks)",800,-20.,20.);
  new TH2D("minstubTrkPoscorrD1","Cluster xTrk Pos Correlation;trk;stub",nStrips+1,-0.5,nStrips - 0.5, nStrips+1,-0.5,nStrips - 0.5 );
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minstubTrkPoscorrD1"));
  h2d->SetOption("colz");

  new TH1D("sminresidualC0_1trkfid","Stub Residual at DUT1 plane(fiducial)(#trk=1)",1000,-10.,10.);

  new TH1I("trkcluseff","",12,-0.5,11.5);

  TH1I* h = dynamic_cast<TH1I*>(Utility::getHist1D("trkcluseff"));
  h->GetXaxis()->SetBinLabel(1,"xtkNodupl(=1)");
  h->GetXaxis()->SetBinLabel(2,"xtkFidD0");
  h->GetXaxis()->SetBinLabel(3,"xtkFidD1");
  h->GetXaxis()->SetBinLabel(4,"xtkFidD0 && xtkFidD1");
  h->GetXaxis()->SetBinLabel(5,"xtkClsMatchD0");
  h->GetXaxis()->SetBinLabel(6,"xtkClsMatchD1");
  h->GetXaxis()->SetBinLabel(7,"xtkClsMatchD0_&&_D1");
  h->GetXaxis()->SetBinLabel(8,"no-match_D0&&_D1");
  h->GetXaxis()->SetBinLabel(9,"xtkStubMatchC0");
  h->GetXaxis()->SetBinLabel(10,"xtkCBCStubMatchC0");
  h->GetXaxis()->SetBinLabel(11,"xtkCBCStubWrongC0");
  h->GetXaxis()->SetBinLabel(12,"xtkCBCStub>1");

  new TH1I("trkCBCeff","",7,-0.5,6.5);
  TH1I* hc = dynamic_cast<TH1I*>(Utility::getHist1D("trkCBCeff"));
  hc->GetXaxis()->SetBinLabel(1,"xtkFidChip0");
  hc->GetXaxis()->SetBinLabel(2,"StubMatchChip0");
  hc->GetXaxis()->SetBinLabel(3,"xtkFidChip1");
  hc->GetXaxis()->SetBinLabel(4,"StubMatchChip1");
  hc->GetXaxis()->SetBinLabel(5,"xtkFidChip0 && StubChip1");
  hc->GetXaxis()->SetBinLabel(6,"xtkFidChip1 && StubChip0");
  hc->GetXaxis()->SetBinLabel(7,"CBCStub>1");


  new TH1I("effVtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effCBCVtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effChip0Vtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effChip1Vtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effRECOVtdc_num",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effVtdc_den",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effChip0Vtdc_den",";TDC;#Events",17,-0.5,16.5);
  new TH1I("effChip1Vtdc_den",";TDC;#Events",17,-0.5,16.5);
  
  //new TH1I("effVtdc_den",";TDC;#Count",17,-0.5,16.5);
  new TH1I("effVtdc_den",";TDC;#Count",17,-0.5,16.5);

  new TH1D("deltaXPos_trkfei4", "Difference in Track impact and Hit X Position after alignment; x_{FeI4} - (x_{Track} + x_{Offset}); Counts;", (24/(250e-3/10.0)), -12.0, 12.0);
  new TH1D("deltaYPos_trkfei4", "Difference in Track Impact and Hit Y Position after alignment; y_{FeI4} - (y_{Track} + y_{Offset}); Counts;", (24/(50e-3/10.0)), -12.0, 12.0);
 
  new TH1D("deltaXPos_trkfei4_UM", "Difference in Track impact and Hit X Position after alignment; x_{FeI4} - (x_{Track} + x_{Offset}); Counts;", (24/(250e-3/10.0)), -12.0, 12.0);
  new TH1D("deltaYPos_trkfei4_UM", "Difference in Track Impact and Hit Y Position after alignment; y_{FeI4} - (y_{Track} + y_{Offset}); Counts;", (24/(50e-3/10.0)), -12.0, 12.0);
  
  new TH1D("deltaXPos_trkfei4_M", "Difference in Track impact and Hit X Position after alignment and matching; x_{FeI4} - (x_{Track} + x_{Offset}); Counts;", (40/0.5e-3), -20.0, 20.0);
  new TH1D("deltaYPos_trkfei4_M", "Difference in Track Impact and Hit Y Position after alignment and matching; y_{FeI4} - (y_{Track} + y_{Offset}); Counts;", (40/0.5e-3), -20.0, 20.0);
  
  int nBins = (36/(0.09/200)); //(36/(0.09/50));
  new TH1F("deltaXPos_trkDUT","Difference in Track Impact and Hit X Position in Det0; x_{Track} - x_{DUT} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1F("deltaXPos_trkfei4_SingleCluster", "Difference in Track impact and Hit X Position after alignment; x_{FeI4} - (x_{Track} + x_{Offset}); Counts;", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4_SingleCluster", "Difference in Track Impact and Hit Y Position after alignment; y_{FeI4} - (y_{Track} + y_{Offset}); Counts;", 40000, -20.0, 20.0);
  
  new TH1D("tdc_Cluster_dut0","TDC phase of Cluster in Det0; TDC Phase; Number of Events;", 20,-0.5,20.-0.5);
  new TH1D("tdc_Cluster_dut1","TDC phase of Cluster in Det1; TDC Phase; Number of Events;", 20,-0.5,20.-0.5);
  
  
  new TH2D("Tks_Det0" , "Tk Impact Det0; Strip Number; TDC Phase", nStrips+1 , -0.5 , nStrips-0.5, 20,-0.5,20.-0.5);
  new TH2D("Tks_Det1" , "Tk Impact Det1; Strip Number; TDC Phase", nStrips+1 , -0.5 , nStrips-0.5, 20,-0.5,20.-0.5);
  new TH2D("Clusters_Det0" , "Cluster Position Det0; Strip Number; TDC Phase", nStrips+1 , -0.5 , nStrips-0.5, 20,-0.5,20.-0.5);
  new TH2D("Clusters_Det1" , "Cluster Position Det1; Strip Number; TDC Phase", nStrips+1 , -0.5 , nStrips-0.5, 20,-0.5,20.-0.5);

  new TH2D("ResTDC_Clusters_Det0" , "Det0 ; x_{Cluster} - x_{Track} [mm]; TDC Phase", (36/(0.09/5)), -18.0, 18.0, 20,-0.5,20.-0.5);
  new TH2D("ResTDC_Clusters_Det1" , "Det1 ; x_{Cluster} - x_{Track} [mm]; TDC Phase", (36/(0.09/5)), -18.0, 18.0, 20,-0.5,20.-0.5);



  new TH1I("trkFidStability","",300,0,300);
  new TH1I("clsD0Stability","",300,0,300);
  new TH1I("clsD1Stability","",300,0,300);
  new TH1I("clsBothStability","",300,0,300);
  new TH1I("stubStability","",300,0,300);
  new TH1I("cbcStability","",300,0,300);

  fout_->mkdir("TrackMatch/CBCcheck");
  fout_->cd("TrackMatch/CBCcheck");
  new TH1I("EventN","EventN;Event Number",300000,0,300000);
  new TH1D("NClusterDUT0","NClusterDUT0;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("NClusterDUT1","NClusterDUT1;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("ClusterPosDUT0","ClusterPosDUT0;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("ClusterPosDUT1","ClusterPosDUT1;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubPosDUT1","RECOStubPos;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubResDUT1","RECOStubRes;Strip Number;#Events",5000,-8,8);
  new TH1D("ClusterDiff","ClusterDiff;Strip Number;#Events",20,0,20);
  new TH2D("ClusDifvsRes","ClusDifvsRes;Strip Number;Res. (mm)",40,-20,20,500,-4,4);

  fout_->mkdir("TrackMatch/RECOcheck");
  fout_->cd("TrackMatch/RECOcheck");
  new TH1I("EventN","EventN;Event Number",300000,0,300000);
  new TH1D("NClusterDUT0","NClusterDUT0;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("NClusterDUT1","NClusterDUT1;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("ClusterPosDUT0","ClusterPosDUT0;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("ClusterPosDUT1","ClusterPosDUT1;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubPosDUT1","RECOStubPos;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubResDUT1","RECOStubRes;Strip Number;#Events",5000,-8,8);
  new TH1D("ClusterDiff","ClusterDiff;Strip Number;#Events",20,0,20);
  new TH2D("ClusDifvsRes","ClusDifvsRes;Strip Number;Res. (mm)",40,-20,20,500,-4,4);

  fout_->mkdir("TrackMatch/Bothcheck");
  fout_->cd("TrackMatch/Bothcheck");
  new TH1I("EventN","EventN;Event Number",300000,0,300000);
  new TH1D("NClusterDUT0","NClusterDUT0;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("NClusterDUT1","NClusterDUT1;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("ClusterPosDUT0","ClusterPosDUT0;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("ClusterPosDUT1","ClusterPosDUT1;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubPosDUT1","RECOStubPos;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubResDUT1","RECOStubRes;Strip Number;#Events",5000,-8,8);
  new TH1D("ClusterDiff","ClusterDiff;Strip Number;#Events",20,0,20);
  new TH2D("ClusDifvsRes","ClusDifvsRes;Strip Number;Res. (mm)",40,-20,20,500,-4,4);

  fout_->mkdir("TrackMatch/CBCOnlycheck");
  fout_->cd("TrackMatch/CBCOnlycheck");
  new TH1I("EventN","EventN;Event Number",300000,0,300000);
  new TH1D("NClusterDUT0","NClusterDUT0;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("NClusterDUT1","NClusterDUT1;#Clusters;#Events",51,-0.5,50.5);
  new TH1D("ClusterPosDUT0","ClusterPosDUT0;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("ClusterPosDUT1","ClusterPosDUT1;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubPosDUT1","RECOStubPos;Strip Number;#Events",1016,-0.5,1015.5);
  new TH1D("RECOStubResDUT1","RECOStubRes;Strip Number;#Events",5000,-8,8);
  new TH1D("ClusterDiff","ClusterDiff;Strip Number;#Events",20,0,20);
  new TH2D("ClusDifvsRes","ClusDifvsRes;Strip Number;Res. (mm)",40,-20,20,500,-4,4);

  fout_->cd("TrackMatch");
  new TProfile("nfidtrk_1k","#Events with 1 Fiducial Tracks;Event Number (#times 1000);Entries for 1000 events",10000,0.5,10000.5);
  new TProfile("nmatchedStub_1k","#Events with Stub Matched;Event Number (#times 1000);Entries for 1000 events",10000,-0.5,10000.5);

  new TH1D("deltaXPos_dut0","Difference in Extrapolated Track and Hit X Position in Det0; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXPos_dut1","Difference in Extrapolated Track and Hit X Position in Det1; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXPos_dut0_SingleTks","Difference in Extrapolated Track and Hit X Position in Det0; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXPos_dut1_SingleTks","Difference in Extrapolated Track and Hit X Position in Det1; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
 
  new TH1D("deltaXCluster_dut0","Difference in Extrapolated Track and Cluster X Position in Det0; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXCluster_dut1","Difference in Extrapolated Track and Cluster X Position in Det1; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXCluster_dut0_M","Difference in Extrapolated Track and Hit X Position in Det0; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXCluster_dut1_M","Difference in Extrapolated Track and Hit X Position in Det1; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXCluster_dut0_NM","Difference in Extrapolated Track and Hit X Position in Det0; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXCluster_dut1_NM","Difference in Extrapolated Track and Hit X Position in Det1; x_{DUT} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  
  new TH1D("deltaX_SingleCluster_dut0","Difference in Extrapolated Track and Cluster X Position in Det0; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaX_SingleCluster_dut1","Difference in Extrapolated Track and Cluster X Position in Det1; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXCluster_dut0_SingleTks","Difference in Extrapolated Track and Cluster X Position in Det0; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  new TH1D("deltaXCluster_dut1_SingleTks","Difference in Extrapolated Track and Cluster X Position in Det1; x_{Cluster} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXStub","Difference in Extrapolated Track and Stub X Position in; x_{Stub} - x_{Track} [mm]; Count [a.u]",nBins, -18.0, 18.0);
  

  new TH1D("nFidTks_dut0","Number of Fiducial Tracks per Event in Det0; #Tracks ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("nFidTks_dut1","Number of Fiducial Tracks per Event in Det0; #Tracks ; #Count [a.u]",100,-0.5,100.5);

  new TH1D("nFidTks_dut0_SingleTks","Number of Fiducial Tracks per Event in Det0; #Tracks ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("nFidTks_dut1_SingleTks","Number of Fiducial Tracks per Event in Det0; #Tracks ; #Count [a.u]",100,-0.5,100.5);

  new TH1D("nHits_dut0","Number of Hits per Event in Det0; #Hits ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("nHits_dut1","Number of Hits per Event in Det1; #Hits ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("nClusters_dut0","Number of Clusters per Event in Det0; #Clusters ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("nClusters_dut1","Number of Clusters per Event in Det1; #Clusters ; #Count [a.u]",100,-0.5,100.5);
  new TH1D("hClusterWidth_dut0","Cluster Width in Det0; Cluster Width (#Strips) ; Count [a.u]",50,-0.5,49.5);
  new TH1D("hClusterWidth_dut1","Cluster Width in Det1; Cluster Width (#Strips) ; Count [a.u]",50,-0.5,49.5);
  new TH1D("nStubs","Number of Stubs per Event; #Stubs ; #Count [a.u]",100,-0.5,100.5);
  
  //new TH1D("sresidualC0multitrkfidNodupl","Stub Residual at DUT0 plane(fiducial)(#trk>1, no duplicate tracks)",800,-20.,20.);
  new TH2D("correlationTksFeI4Hits","Matched Trk and Fei4 Hit Correlation;# FeI4 Hits; # Matched Tracks; # Fei4 Hits",30,-0.5,29.5,30,-0.5,29.5);
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("correlationTksFeI4Hits"));
  h2d->SetOption("colz");



  new TH2F("hHitTkX", "Correlation between (X) FeI4 hit and matched track; x_{FeI4} [mm]; x_{Track} [mm]",(int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/250e-3) , -20.0, 20.0); 
  new TH2F("hHitTkY", "Correlation between (Y) FeI4 hit and matched track; y_{FeI4} [mm]; y_{Track} [mm]",(int)(40.0/50e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0); 
  new TH2F("TkMapM", "Track Impact on FeI4 Plane (Matched);x_{Track} [mm]; y_{Track} [mm]",(int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0); 
  new TH2F("HtMapUM","FEI4 Hit Map (Un-matched); x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  new TH2F("HtMapM","FEI4 Hit Map (Matched); x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  
  new TH1D("Det0_ZAlignment", "Resolution of Det0 for different Z_{offsets}; #Delta Z_{offset} [mm] ; #sigma [nm]", (1000)/(0.5) , -0.5 , 999.5 );
  new TH1D("Det1_ZAlignment", "Resolution of Det1 for different Z_{offsets}; #Delta Z_{offset} [mm] ; #sigma [nm]", (1000)/(0.5) ,  -0.5 , 999.5 );


  nBins = (36/(0.09/1)) ; 
  new TH1D("TkImpact_Det0","Track Impact Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("TkImpact_Det1","Track Impact Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("Hit_Det0","Hit Position Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("Hit_Det1","Hit Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("Cluster_Det0","Cluster Position Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("Cluster_Det1","Cluster Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("Stubs","Stub Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("deltaXCluster","#Delta x_{Clusters}; x_{Cluster}; Number of Events [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("SingleTkImpact_Det0","Track Impact Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleTkImpact_Det1","Track Impact Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  
  new TH1D("SingleClusterTrksImpact_Det0","Track Impact Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleClusterTrksImpact_Det1","Track Impact Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);

  new TH1D("SingleHit_Det0","Hit Position Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleHit_Det1","Hit Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleCluster_Det0","Cluster Position Det0; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleCluster_Det1","Cluster Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  new TH1D("SingleStubs","Stub Position Det1; x_{DUT}; Number of Events [a.u]",nBins, -18.0, 18.0);
  
  new TH2D("deltaXPos_dut0_2D","Test ; x_{DUT} - x_{Track} [mm]; x_{DUT} [mm]",(36/(0.09/5)), -18.0, 18.0, (36/(0.09/1)), -18.0, 18.0 );
  new TH2D("deltaXPos_dut1_2D","Test ; x_{DUT} - x_{Track} [mm]; x_{DUT} [mm]",(36/(0.09/5)), -18.0, 18.0, (36/(0.09/1)), -18.0, 18.0 );
  
  new TH2D("wClusterXPos_dut0","Test ; x_{Cluster}  [mm]; Cluster Width",(36/(0.09/1)), -18.0, 18.0, (100) ,-0.5,99.5 );
  new TH2D("wClusterXPos_dut1","Test ; x_{Cluster}  [mm]; Cluster Width",(36/(0.09/1)), -18.0, 18.0, (100) ,-0.5,99.5 );
 
  new TH1D("TkImpact_Strip_Det0" , "Track Impact Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("TkImpact_Strip_Det1" , "Track Impact Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);

  new TH1D("SingleTkImpact_Strip_Det0","Track Impact Det0; Strip Number; Number of Events [a.u]",254 , -0.5 ,254-0.5);
  new TH1D("SingleTkImpact_Strip_Det1","Track Impact Det1; Strip Number; Number of Events [a.u]",254 , -0.5 ,254-0.5);
  
  new TH1D("Hit_Strip_Det0" , "Hit Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("Hit_Strip_Det1" , "Hit Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);

  new TH1D("SingleHit_Strip_Det0" , "Hit Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("SingleHit_Strip_Det1" , "Hit Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);

  new TH1D("Hit_Strip_Det0_All" , "Hit Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("Hit_Strip_Det1_All" , "Hit Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);

  new TH1D("Cluster_Strip_Det0" , "Cluster Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("Cluster_Strip_Det1" , "Cluster Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);


  new TH1D("SingleCluster_Strip_Det0" , "Cluster Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("SingleCluster_Strip_Det1" , "Cluster Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);

  new TH1D("Cluster_Strip_Det0_All" , "Cluster Position Det0; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH1D("Cluster_Strip_Det1_All" , "Cluster Position Det1; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);


  new TH1D("Stub_Strip" , "Stub Position; Strip Number; Number of Events", 254 , -0.5 ,254-0.5);
  new TH2D("HitMap","FEI4 Hit Map; x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  new TH2D("TkHitMap","Tk Hit Map; x_{Track} [mm]; y_{Track} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);

  new TH1D("chi2_BackgroundFit","#Chi^{2}/NdOF Background Fit; Size of matching window [1/p]; #Chi^{2}/NdOF", 100/0.05 ,-0.5, 100.-0.5);
  new TH1D("det0_Corrected_ClusterEfficiency","Corrected Det0 Cluster Efficiency; Size of matching window [1/p]; Efficiency ", 100/0.05 ,-0.5, 100.-0.5);
  


  new TH1D("hposyTkDUT0","Ypos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",500,-100.,100.);
  new TH1D("hposyTkDUT1","Ypos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",500,-100.,100.);

  new TH1D("hposxTkDUT0","Xpos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",500,-100.,100.);
  new TH1D("hposxTkDUT1","Xpos of extrapolated track at DUT0 plane after alignment(#fid trk=1)",500,-100.,100.);


  new TH2D("minstubTrkPoscorrD1_all","Closest-Stub xTrk Pos Correlation;trk;stub",1016,-0.5,1015.5, 1016,-0.5,1015.5);
  new TH2D("minstubTrkPoscorrD1_matched","Closest-Stub xTrk Pos Correlation(matched);trk;stub",1016,-0.5,1015.5, 1016,-0.5,1015.5);
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minstubTrkPoscorrD1_all"));
  h2d->SetOption("colz");
  h2d = dynamic_cast<TH2D*>(Utility::getHist2D("minstubTrkPoscorrD1_matched"));
  h2d->SetOption("colz");

  new TH1D("stubPos_all","Stub Pos;stub pos (strip)",1016,-0.5,1015.5);
  new TH1D("stubPos_matched","Stub Pos;stub pos (strip)",1016,-0.5,1015.5);

  new TH1D("trackPos_all","Stub Pos;stub pos (strip)",1016,-0.5,1015.5);
  new TH1D("trackPos_matched","Stub Pos;stub pos (strip)",1016,-0.5,1015.5);


  //new TH1F("deltaXPos_trkfei4", "Difference in matched Track impact and Hit X Position", 40000, -100.0, 100.0);
  //new TH1F("deltaYPos_trkfei4", "Difference in matched Track Impact and Hit Y Position", 40000, -100.0, 100.0);

  new TProfile("nfidtrk_1k","#Events with 1 Fiducial Tracks;Event Number (#times 1000);Entries for 1000 events",10000,0.5,10000.5);
  new TProfile("nmatchedStub_1k","#Events with Stub Matched;Event Number (#times 1000);Entries for 1000 events",10000,-0.5,10000.5);
}

void Histogrammer::bookTelescopeAnalysisHistograms() {
  fout_->cd();
  fout_->mkdir("TelescopeAnalysis");
  fout_->cd("TelescopeAnalysis");
  new TH1D("nTrack_noduplicate","#Tracks Telescope(after duplicate removal);#tracks;#events",30,-0.5,29.5);
  new TH1D("nTrack_fiducial","#Tracks Telescope(fiducial);#tracks;#events",30,-0.5,29.5);
 
  new TH1F("deltaXPos_fit", "Difference in Track impact and Hit X Position", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_fit", "Difference in Track Impact and Hit Y Position", 40000, -20.0, 20.0);

  new TH1F("hMatchingEfficiency_Angle","Cluster Matching Efficiency;#theta;#Matching Efficiency", (30/1.0) , -15.0 , 15.0 );
  new TH1D("hMatchingEfficiency_YOffset","Cluster Matching Efficiency;Y_{Offset} [#mum]; Matching Efficiency", (int)(200/10.0) , -200.5 , 199.5 );
 
  new TH1D("tk_dYdZ","dydz of tracks;dydz; [a.u.]; Count [a.u]",(20e-4/1e-5),-10e-4,10e-4);
  new TH1D("tk_dXdZ","dxdz of tracks;dxdz; [a.u.]; Count [a.u]",(20e-4/1e-5),-10e-4,10e-4);


  new TH1I("nTrack","Number Telescope Tracks Per Event;#tracks;#events",30,-0.5,100.5);
  new TH1I("nTracks","Number Telescope Tracks Per Event (after duplicate removal);#tracks;#events",100,-0.5,100.5);
  new TH1I("nTrackM","#Tracks Telescope (Matched);#tracks;#events",30,-0.5,29.5);
  new TH1I("nhitsFei4","#Hits in FeI4;#hitsFEI4;#Events",100,-0.5,99.5);
  new TH1D("nPxlClusterFei4","Number of Pixels per cluster in Fei4 plane;Pixels in Cluster;Clusters",100,-0.5,99.5);
  new TH1D("nClustersFei4","Number of Clusters per Event;Number of Clusters Per Event;Count",100,-0.5,99.5);
  new TH1F("hTotFei4","Time Over Threshold FeI4 Hits;Time Over Threshold [a.u];Number of Events [a.u]",100,-0.5,99.5);
  new TH1F("hLv1Fei4","L1 FeI4 Hits;L1 FeI4 Cluster [a.u];Number of Events [a.u]",100,-0.5,99.5);
 
  new TH1F("TkXPos", "Track Impact on Fei4 Plane; FEI4 X (mm); Number of Events",(int)(40.0/(250e-3)) , -20.0, 20.0 );
  new TH1F("TkYPos", "Track Impact on Fei4 Plane; FEI4 Y (mm); Number of Events",(int)(40.0/(50e-3)) , -20.0, 20.0 );
  new TH1F("HtColumn", "Hit Column", 85, -0.5, 84.5);
  new TH1F("HtRow", "Hit Row", 340, -0.5, 339.5);
  new TH1F("HtXPos", "FEI4 X; x_{FeI4} (mm); Number of Events",(int)(40.0/(250e-3)) , -20.0, 20.0 );
  new TH1F("HtYPos", "FEI4 Y; y_{FeI4} (mm); Number of Events",(int)(40.0/(50e-3)) , -20.0, 20.0 );
  new TH1F("TkXPosDUT0", "Track Impact on DUT0 (extrapolated); FEI4 X (mm); Number of Events",(int)(40.0/(250e-3)) , -20.0, 20.0 );
  new TH1F("TkYPosDUT1", "Track Impact on DUT1 (extrapolated); FEI4 Y (mm); Number of Events",(int)(40.0/(50e-3)) , -20.0, 20.0 );
  new TH1F("hTkChi2", "Chi Squared Tracks", (100.0)/0.05 , -0.5, 99.5);
  new TH1F("hTkNdof", "NDOF Tracks", (100.0)/1.0 , -0.5, 99.5);
  new TH1F("hTkChi2_Ndof", "NDOF Tracks", (10.0)/0.01 , -0.5, 10.0 - 0.5 );
  new TH1F("hTkChi2_SingleTks", "Chi Squared Tracks", (10.0)/0.01 , -0.5, 10.0 - 0.5 );
  new TH1F("hTkChi2_MultiTks", "Chi Squared Tracks", (10.0)/0.01 , -0.5, 10.0 - 0.5 );
  

  int nBins = (20e-3/20);
  new TH1D("deltaXPos", "Difference in Track impact and FeI4 Cluster X Position; x_{FeI4} - x_{Track}; Counts", (40.0)/5e-3 , -20.0, 20.0);
  new TH1D("deltaYPos", "Difference in Track Impact and FeI4 Cluster Y Position; y_{FeI4} - y_{Track}; Counts", (40.0)/5e-3 , -20.0, 20.0);
#ifdef NOV_15 
    new TH1D("deltaXPos", "Difference in Track Impact and FeI4 Cluster X Position; -x_{FeI4} - y_{Track}; Counts", (40.0)/5e-3 , -20.0, 20.0);
    new TH1D("deltaYPos", "Difference in Track impact and FeI4 Cluster Y Position; y_{FeI4} - x_{Track}; Counts", (40.0)/5e-3 , -20.0, 20.0);
#endif
  new TH1F("deltaDistance", "Distance between track and hit",  40000, -20.0, 20.0);
  new TH1F("deltaXPos_SinglePixel", "Difference in Track impact and Hit X Position (Single Pixel Hit)", 40000 , -20.0, 20.0);
  new TH1F("deltaXPos_SinglePixel_BSub", "Difference in Track impact and Hit X Position (Single Pixel Hit)", 40000 , -20.0, 20.0);
  new TH1F("deltaYPos_SinglePixel_BSub", "Difference in Track impact and Hit X Position (Single Pixel Hit)", 40000 , -20.0, 20.0);
  new TH1F("deltaYPos_SinglePixel", "Difference in Track Impact and Hit Y Position (Single Pixel Hit)", 40000, -20.0, 20.0);
  new TH1F("deltaXPos_TwoPixel", "Difference in Track impact and Hit X Position (Two Pixel Hit)", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_TwoPixel", "Difference in Track Impact and Hit Y Position (Two Pixel Hit)", 40000, -20.0, 20.0);
  new TH1F("deltaXPos_MultiPixel", "Difference in Track impact and Hit X Position (Multi Pixel Hit)", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_MultiPixel", "Difference in Track Impact and Hit Y Position (Multi Pixel Hit)", 40000, -20.0, 20.0);

  new TH2F("TkMap", " Track Impact on FeI4 Plane;x_{Track} [mm]; y_{Track} [mm]",(int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0); 
  new TH2F("HtMap","FEI4 Hit Map; x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  new TH2F("HtMapSingle","FEI4 Hit Map (Single Pixel Hits); x_{FeI4} [mm]; y_{FeI4} [mm" , (int)(20.0/250e-3) , -10.0, 10.0,  (int)(20.0/50e-3) , -10.0, 10.0); 
  new TH2F("HtMapTwo","FEI4 Hit Map (Two Pixel Hits);x_{FeI4} [mm]; y_{FeI4} [mm" ,  (int)(20.0/250e-3) , -10.0, 10.0,  (int)(20.0/50e-3) , -10.0, 10.0); 
  new TH2F("HtMapMulti","FEI4 Hit Map (Multi Pixel Hits);x_{FeI4} [mm]; y_{FeI4} [mm" ,  (int)(20.0/250e-3) , -10.0, 10.0,  (int)(20.0/50e-3) , -10.0, 10.0); 
  new TH2F("tkXPosVsHtXPos", "All Telescope Tracks;x_{FeI4} [mm]; x_{Track} [mm]",(int)(20.0/250e-3) , -10.0, 10.0,  (int)(20.0/50e-3) , -10.0, 10.0); 
  new TH2F("tkYPosVsHtYPos", "All Telescope Tracks;y_{FeI4} [mm]; y_{Track} [mm]",(int)(20.0/250e-3) , -10.0, 10.0,  (int)(20.0/50e-3) , -10.0, 10.0); 
  new TH2F("hHitTkX", "Correlation between (X) FeI4 hit telescope tracks; x_{FeI4} [mm]; x_{Track} [mm]",(int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/250e-3) , -20.0, 20.0); 
  new TH2F("hHitTkY", "Correlation between (Y) FeI4 hit telescope tracks; y_{FeI4} [mm]; y_{Track} [mm]",(int)(40.0/50e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0); 
  
  new TH2F("tkXPosVsHtXPosM", "Matched Telescope Tracks ( Xpos Track Impact vs. Xpos FeI4-Hit );Xpos of FeI4-Hit(mm);Xpos of Track Impact(mm)", 96, -12.0, 12.0, 96, -12.0, 12.0);
  new TH2F("tkYPosVsHtYPosM", "Matched Telescope Tracks ( Ypos Track Impact vs. Ypos FeI4-Hit );Ypos of FeI4-Hit(mm);Ypos of Track Impact(mm)", 240, -12.0, 12.0, 240, -12.0, 12.0);

  new TH2F("checkImpact","checkImpact;Track # (Xmin);Track # (Ymin)",30,-0.5,29.5,30,-0.5,29.5 );

  TH2F* h = dynamic_cast<TH2F*>(Utility::getHist2D("tkXPosVsHtXPos"));
  h->SetOption("colz");
  h = dynamic_cast<TH2F*>(Utility::getHist2D("tkYPosVsHtYPos"));
  h->SetOption("colz");
  h = dynamic_cast<TH2F*>(Utility::getHist2D("HtMap"));
  h->SetOption("colz");
  h = dynamic_cast<TH2F*>(Utility::getHist2D("tkXPosVsHtXPosM"));
  h->SetOption("colz");
  h = dynamic_cast<TH2F*>(Utility::getHist2D("tkYPosVsHtYPosM"));
  h->SetOption("colz");
  
  new TH1F("deltaXPos_trkfei4", "Difference in Track impact and Hit X Position after alignment; x_{FeI4} - (x_{Track} + x_{Offset}); Counts;", (40/(10e-3)), -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4", "Difference in Track Impact and Hit Y Position after alignment; y_{FeI4} - (y_{Track} + y_{Offset}); Counts;", (40/(10e-3)), -20.0, 20.0);
  new TH1F("deltaXPos_trkfei4M", "Difference in matched Track impact and Hit X Position", 40000, -20.0, 20.0);
  new TH1F("deltaYPos_trkfei4M", "Difference in matched Track Impact and Hit Y Position", 40000, -20.0, 20.0);

  TH1F* h1 = dynamic_cast<TH1F*>(Utility::getHist1D("deltaXPos_trkfei4"));
  h1->SetLineColor(kRed);
  h1 = dynamic_cast<TH1F*>(Utility::getHist1D("deltaYPos_trkfei4"));
  h1->SetLineColor(kRed);
  h1 = dynamic_cast<TH1F*>(Utility::getHist1D("deltaXPos_trkfei4M"));
  h1->SetLineColor(kGreen);
  h1->SetFillColor(kGreen);
  h1->SetFillStyle(3001);
  h1 = dynamic_cast<TH1F*>(Utility::getHist1D("deltaYPos_trkfei4M"));
  h1->SetLineColor(kGreen);
  h1->SetFillColor(kGreen);
  h1->SetFillStyle(3001);
  
  TH1D* h1f = dynamic_cast<TH1D*>(Utility::getHist1D("deltaXPos"));
  h1f->SetLineColor(kGreen);
  h1f = dynamic_cast<TH1D*>(Utility::getHist1D("deltaYPos"));
  h1f->SetLineColor(kBlue);

}

void Histogrammer::bookTrackFitHistograms(float zMin, float zStep, int zNsteps){
  fout_->cd();
  fout_->mkdir("TrackFit");
  fout_->cd("TrackFit");

  new TH1I("d0_1tk1Hit_diffX","X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1Hit_diffX","X_{TkAtDUT}-X_{DUT}, d1",100000,-100,100);

  new TH1I("d0_1tk1Hit_diffX_bis","X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1Hit_diffX_bis","X_{TkAtDUT}-X_{DUT}, d1",100000,-100,100);

  new TH1I("d0_1tk1Hit_diffX_ter","X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1Hit_diffX_ter","X_{TkAtDUT}-X_{DUT}, d1",100000,-100,100);

  for (int iz=0; iz<zNsteps; iz++){
    new TH1I(Form("d0_1tk1Hit_diffX_iz%i", iz),"X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
    new TH1I(Form("d1_1tk1Hit_diffX_iz%i", iz),"X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
  }

  float zMax = zMin + ((float)zNsteps) * zStep;
  float shift = zStep/2.;

  new TH1F("d0_offsetVsZ", "x_{DUT} offset vs injected z_{DUT}, d0", zNsteps, zMin-shift, zMax-shift);
  new TH1F("d1_offsetVsZ", "x_{DUT} offset vs injected z_{DUT},d1", zNsteps, zMin-shift, zMax-shift);

  new TH1F("d0_chi2VsZ","chi2 vs injected z_{DUT}, d0", zNsteps, zMin-shift, zMax-shift);
  new TH1F("d1_chi2VsZ","chi2 vs injected z_{DUT}, d1", zNsteps, zMin-shift, zMax-shift);

  new TH1I("d0_1tk1Hit_diffX_aligned","X_{TkAtDUT}-X_{DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1Hit_diffX_aligned","X_{TkAtDUT}-X_{DUT}, d1",100000,-100,100);
  new TH1I("d0_1tk1ClusterBothPlanes_diffX_aligned","X_{TkAtDUT}-X_{Cls,DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1ClusterBothPlanes_diffX_aligned","X_{TkAtDUT}-X_{Cls,DUT}, d1",100000,-100,100);
  new TH1I("d0_1tk1ClusterBothPlanesConstraint_diffX_aligned","X_{TkAtDUT}-X_{Cls,DUT}, d0",100000,-100,100);
  new TH1I("d1_1tk1ClusterBothPlanesConstraint_diffX_aligned","X_{TkAtDUT}-X_{Cls,DUT}, d1",100000,-100,100);

  new TH1F("bothPlanes_chi2VsTheta","chi2 vs injected #theta", 41, -20.-0.5, 21.-0.5);
  new TH1F("bothPlanesConstraint_chi2VsTheta","chi2 vs injected #theta", 41, -20.-0.5, 21.-0.5);
  new TH1F("bothPlanesConstraint_chi2VsDeltaZ","chi2 vs injected #deltaZ", 41, 0.-0.125, 10.25-0.125);

}

void Histogrammer::bookFeI4Histograms()
{
  fout_->mkdir("FeI4");
  fout_->cd("FeI4");
  new TH1D("x0_Residual_Evnt",";Event number; Number of Mismatched Clusters",(300e3/5) , -0.5, 300e3 -0.5);
  new TH1D("lv1","Number of level1 accept;Level1 Accept; [a.u.]; Count [a.u]",(100/1.0),-0.5,99.5);
  new TH1D("lv1_M","Number of level1 accept;Level1 Accept; [a.u.]; Count [a.u]",(100/1.0),-0.5,99.5);
  
  new TH1I("nHits","Number of FeI4 hits per event;Pixels in Cluster; [Hits]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nClusters","Number of FeI4 custers per event; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
  new TH1I("hTotFei4","Total ToT FeI4 Cluster; [Clusters];Count [a.u]", 100 ,-0.5,99.5);
  new TH1D("hQFeI4","Charge Deposited in FeI4 per Event; Q [kElectrons]; Number of Events",(1000/1.5) , -0.5,  1000-0.5);
  new TH1I("nPixels","Number of Pixels per cluster in Fei4 plane;Pixels in Cluster; Count[a.u]",100,-0.5,99.5);
  
  new TH1I("nTracks","Number of Tracks per event; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_event","Number of Tracks within acceptance window per event; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_cluster","Number of Tracks within acceptance window per cluster in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  
  new TH1I("nTracks_1Cluster","Number of Tracks for events with 1 cluster in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_2Cluster","Number of Tracks for events with 2 clusters in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_3Cluster","Number of Tracks for events with 3 clusters in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_4Cluster","Number of Tracks for events with 4 clusters in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_5Cluster","Number of Tracks for events with 5 clusters in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  new TH1I("nTracks_6Cluster","Number of Tracks for events with 6 clusters in Fei4 plane; [Tracks]; Count [a.u]",100,-0.5,99.5);
  
  double Min_FeI4 = -15.0;
  double Max_FeI4 =  15.0; 
  double Range_FeI4 = (Max_FeI4-Min_FeI4);
  int nResX_binSize = (int)(Range_FeI4/(25e-3));
  int nResY_binSize = (int)(Range_FeI4/(25e-3));
  double Min_DUT = -100.0 ; 
  double Max_DUT =  100.0 ; 
  double Range_DUT = (Max_DUT-Min_DUT);
  new TH1D("tk_dYdZ","dydz of tracks matched to FeI4 clusters;dydz; [a.u.]; Count [a.u]",(20e-4/1e-5),-10e-4,10e-4);
  new TH1D("tk_dXdZ","dxdz of tracks matched to FeI4 clusters;dxdz; [a.u.]; Count [a.u]",(20e-4/1e-5),-10e-4,10e-4);
  
  new TH2D("x0_Correlation","; x_{FeI4,Cluster} [mm]; x_{DUT0,Cluster}", Range_FeI4/(50e-3/5), Min_FeI4, Max_FeI4 , Range_DUT/(0.09), Min_DUT, Max_DUT );
  new TH2D("x1_Correlation","; x_{FeI4,Cluster} [mm]; x_{DUT1,Cluster}", Range_FeI4/(50e-3/5), Min_FeI4, Max_FeI4 , Range_DUT/(0.09), Min_DUT, Max_DUT );
  
  new TH2D("nClusters0_Correlation","; n_{FeI4,Cluster} [mm]; n_{DUT0,Cluster}",100,-0.5,99.5 , 100,-0.5,99.5);
  new TH2D("nClusters1_Correlation","; n_{FeI4,Cluster} [mm]; n_{DUT1,Cluster}",100,-0.5,99.5 , 100,-0.5,99.5);
  new TH2F("HitMap","FEI4 Hit Map; x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  
  new TH2F("HitMap_Raw","FEI4 Hit Map; x_{FeI4} [mm]; y_{FeI4} [mm]" , (int)(40.0/250e-3) , -20.0, 20.0,  (int)(40.0/50e-3) , -20.0, 20.0);
  
  
  new TH2D("Residuals", "Difference in Track Impact and FeI4 Cluster; x_{FeI4} - x_{Track} [mm]; y_{FeI4} - y_{Track} [mm]", nResX_binSize, Min_FeI4 , Max_FeI4,  nResY_binSize,  Min_FeI4 , Max_FeI4);
  
  new TH1D("x0_Residual","; x_{FeI4,Cluster} - x_{DUT0,Cluster} [mm]; Counts [a.u]", std::max(Range_FeI4,Range_DUT)/(250e-3/1.25), std::min(Min_FeI4,Min_DUT), std::max(Max_FeI4,Max_DUT));  
  new TH1D("x1_Residual","; x_{FeI4,Cluster} - x_{DUT1,Cluster} [mm]; Counts [a.u]", std::max(Range_FeI4,Range_DUT)/(250e-3/1.25), std::min(Min_FeI4,Min_DUT), std::max(Max_FeI4,Max_DUT));  

  new TH2D("x0_Correlation_Search","; x_{FeI4,Cluster} [mm]; x_{DUT0,Cluster}", Range_FeI4/(50e-3/5), Min_FeI4, Max_FeI4 , Range_DUT/(0.09), Min_DUT, Max_DUT );
  new TH1D("x0_Residual_Search","; x_{FeI4,Cluster} - x_{DUT0,Cluster} [mm]; Counts [a.u]", std::max(Range_FeI4,Range_DUT)/(250e-3/1.25), std::min(Min_FeI4,Min_DUT), std::max(Max_FeI4,Max_DUT));  

  int groupID_min = 0 - 0.5; 
  int groupID_max = 1e4 - 0.5 ; 
  int nGroups = groupID_max - groupID_min ; 
  new TProfile("EventAlignment_Correlation", "; Group ID; Percentage of Correlated Events;" , nGroups , groupID_min , groupID_max);
  new TProfile("EventAlignment_Offset", "; Group ID; Percentage of Correlated Events;" , nGroups , groupID_min , groupID_max);

} 

void Histogrammer::bookCBCHistograms(std::string cbc)
{
  TString d(cbc);
  fout_->cd();
  fout_->mkdir(d);
  fout_->cd(d);
  int nStrips = 254; 
#ifdef NOV_15 
  nStrips = nStrips*4;
#endif
  double xMin = -1.0*0.09*1.5*nStrips; 
  double xMax = 1.0*0.09*1.5*nStrips;
  double xBin = 0.09/1.0;  
  double xOffset = 0.0;//xBin/2.0;
  new TH1D("fidTracks_num", ";# Fiducial Tracks Per Event; Number of Events" , 100 , 0.0-0.5, 100.0-0.5 );
  new TH1D("fidTracks_impact_mm", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin , xMax );
  new TH1D("fidTracks_impact_strips", ";x_{Track} [mm]; Number of Events" , nStrips+1, -0.5, nStrips-0.5);
  new TH1D("matchedCluster_mm", ";x_{Cluster} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin , xMax );
  
  new TH1D("fidTracks_impact_mm_CentralRegion", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/(0.1*xBin) , xMin , xMax );
  new TH1D("matchedCluster_mm_CentralRegion", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/(0.1*xBin) , xMin , xMax );
  new TH2D("fidTracks_impact_mm_CentralRegion_2D", ";x_{Track} [mm]; Cluster Width [Strips]" , (xMax-xMin)/(0.1*xBin) , xMin , xMax , 5 , -0.5 , 5-0.5) ;
  new TH2D("matchedCluster_mm_CentralRegion_2D", ";x_{Track} [mm]; Cluster Width [Strips]" , (xMax-xMin)/(0.1*xBin) , xMin , xMax ,5  , -0.5 , 5-0.5 );
  
  new TH1D("matchedCluster_strips", ";x_{Cluster} [mm]; Number of Events" , nStrips+1, -0.5, nStrips-0.5);
  new TH2D("matchedCluster_CW", ";x_{Cluster} [mm]; Cluster Width [Strips]" ,  nStrips+1, -0.5, nStrips-0.5, 10 , -0.5 , 9.5);
  new TH2D("matchedCluster_CW_mm", ";x_{Track} [mm]; Cluster Width [Strips]" , (xMax-xMin)/xBin, xMin , xMax, 10 , -0.5 , 9.5);
  new TProfile("matchedClusterWidth_mm", ";x_{Track} [mm]; Average Cluster Width [Strips]" , (xMax-xMin)/xBin, xMin , xMax );
  
  new TH1D("matchedStub_mm", ";x_{Stub} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin , xMax );
  new TH2D("matchedStub_mm_2D", ";x_{Track} [mm]; TDC Phase" , (xMax-xMin)/xBin, xMin , xMax , 20 , -0.5 , 19.5);
  new TH1D("candidateStubs_impact_mm", ";x_{Track} [mm]; Number of Events" , (xMax-xMin)/xBin, xMin , xMax );
  new TH2D("candidateStubs_impact_mm_2D", ";x_{Track} [mm]; TDC Phase" , (xMax-xMin)/xBin, xMin , xMax , 20 , -0.5 , 19.5);
  //new TH2D("candidateStubs_impact_mm_CW", ";x_{Track} [mm]; TDC Phase" , (xMax-xMin)/xBin, xMin , xMax , 20 , -0.5 , 19.5);
  
  new TH1D("nCandidateClusters", ";TDC Phase; Number of Events", 17 , -0.5 , 16.5);
  new TH1D("nClusters_ANA", ";TDC Phase; Number of Clusters", 17 , -0.5 , 16.5);
  new TH1D("nCandidateStubs", ";TDC Phase; Number of Events", 17 , -0.5 , 16.5);
  new TH1D("nStubs_CHIP", ";TDC Phase; Number of Events", 17 , -0.5 , 16.5);
  new TH1D("nStubs_RECO", ";TDC Phase; Number of Events", 17 , -0.5 , 16.5);
  new TH1D("nStubs_ANA", ";TDC Phase; Number of Events", 17 , -0.5 , 16.5);
  new TH1D("stubMatching", "; ; Number of Events" , 5 , -0.5 , 5.0 -0.5);
  new TH2D("stubMatching_2D", "; ; Number of Events" , 5 , -0.5 , 5.0 -0.5 ,(xMax-xMin)/xBin, xMin , xMax );
  new TH2D("stubMatching_CW", "; ; Number of Events" , 5 , -0.5 , 5.0 -0.5 , 20 , -0.5 , 19.5);
  new TH2D("stubMatching_TDC", "; ; Number of Events" , 5 , -0.5 , 5.0 -0.5 , 20 , -0.5 , 19.5);
}

TH1* Histogrammer::GetHistoByName(const char* dir, const char* hname){
  fout_->cd(dir);
  return Utility::getHist1D(hname);
}

TH1* Histogrammer::GetHistoByName(const std::string& dir, const std::string& hname) {
  fout_->cd(dir.c_str());
  return Utility::getHist1D(hname);

}

TH2* Histogrammer::Get2DHistoByName(const char* dir, const char* hname){
  fout_->cd(dir);
  return Utility::getHist2D(hname);
}

TH2* Histogrammer::Get2DHistoByName(const std::string& dir, const std::string& hname) {
  fout_->cd(dir.c_str());
  return Utility::getHist2D(hname);
}

TProfile* Histogrammer::GetProfileByName(const std::string& dir, const std::string& hname) {
  fout_->cd(dir.c_str());
  return Utility::getHistProfile(hname);
}
TProfile* Histogrammer::GetProfileByName(const char* dir, const char* hname){
  fout_->cd(dir);
  return Utility::getHistProfile(hname);
}

void Histogrammer::FillAlignmentOffsetVsZ(const char* det, const char* histo, int iz, float z, float x, float x_err){

  fout_->cd("TrackFit");
  char histname[50];
  strcpy( histname, det );
  strcat( histname, histo );
  TH1* h = (TH1*) gDirectory->Get(histname);
  h->Fill(z, x);
  h->SetBinError(iz+1, x_err);
}

void Histogrammer::closeFile() { 
  fout_->cd();
  fout_->Write();
  fout_->Close();
  isFileopen_=false;  
}
void Histogrammer::fillClusterHistograms( const char* det, std::vector<tbeam::cluster>& cvec, 
                                          const char* col) {
  fout_->cd(det);
  TString c(col);
  Utility::fillHist1D( "ncluster" + c, cvec.size() );
  for( unsigned int i =0; i<cvec.size(); i++ ) {
    Utility::fillHist1D("clusterWidth" + c,cvec[i].size);
    Utility::fillHist1D("clusterPos" + c,cvec[i].x);
    Utility::fillHistProfile("clusterWidthVsPosProf" + c,cvec[i].x,cvec[i].size);
    Utility::fillHist2D("clusterWidthVsPos2D" + c,cvec[i].x,cvec[i].size);
  } 
}

Histogrammer::~Histogrammer() {
  if(isFileopen_)  {
    std::cout << "You forgot to close the output file!!!Closing it now" << std::endl;  
    closeFile();
   }
  delete fout_; 
}
