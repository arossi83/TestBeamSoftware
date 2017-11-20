#include "Utility.h"
#include <climits>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

using std::string;

namespace Utility {
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
  
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  
    while (std::string::npos != pos || string::npos != lastPos)  {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
    
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
    
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }

  void getChannelMaskedHits( std::vector<int>& vec, const std::vector<int>& mch ) {
    std::vector<int> vtemp = vec;
    vec.clear();
    for(auto& ch : vtemp) {
      //std::cout << ch;
      if( std::find(mch.begin(), mch.end(), ch) != mch.end() )   continue;
      vec.push_back(ch);

    }
  }

  void getChannelMaskedClusters( std::vector<tbeam::cluster*>& vec, const std::vector<int>& mch ) {
    std::vector<tbeam::cluster*> vtemp = vec;
    vec.clear();
    for(auto& clus : vtemp) {
      //std::cout << ch;
      if( std::find(mch.begin(), mch.end(), clus->x) != mch.end() )   continue;
      vec.push_back(clus);
    }
  }

  void getChannelMaskedStubs( std::vector<tbeam::stub*>& vec, const std::vector<int>& mch ) {
    std::vector<tbeam::stub*> vtemp = vec;
    vec.clear();
    for(auto& stub : vtemp) {
      //std::cout << ch;
      if( std::find(mch.begin(), mch.end(), stub->x) != mch.end() )   continue;
      vec.push_back(stub);
    }
  }

  void fill2DHistofromVec( const std::vector<int>& vecC0, const std::vector<int>& vecC1,const char* h) {
    for( unsigned int i = 0; i<vecC0.size(); i++ ) {
      Utility::fillHist2D(h,vecC0.at(i),0);
    }
    for( unsigned int i = 0; i<vecC1.size(); i++ ) {
      Utility::fillHist2D(h,1015-vecC1.at(i),1);
    }
  }

  int readStubWord( std::map<std::string,std::vector<unsigned int> >& stubids, const uint32_t sWord ) {
    int ncbcSw = 0;
    if (sWord > 0) {
	for (unsigned int i = 0; i < 16; i++) {
#ifdef NOV_15
          if (i == 11 || i == 13) continue;//only for nov15
#endif
          if ((sWord >> i) & 0x1) {
            ncbcSw++;
	    if (i <= 7) stubids.at("C0").push_back(i);
            else stubids.at("C1").push_back(i-8);
	  }
	}
      }
    return ncbcSw;  
  }

  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  
  TH1* getHist1D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      std::cerr << "**** getHist1D: Histogram for <" << hname << "> not found!" << std::endl;
      return 0;
    }
    TH1 *h = 0;
    if (obj->InheritsFrom("TH1D"))
      h = dynamic_cast<TH1D*>(obj);
    else if (obj->InheritsFrom("TH1C"))
      h = dynamic_cast<TH1C*>(obj);
    //else if (obj->InheritsFrom("TH1K"))
    //  h = dynamic_cast<TH1K*>(obj);
    else if (obj->InheritsFrom("TH1S"))
      h = dynamic_cast<TH1S*>(obj);
    else if (obj->InheritsFrom("TH1I"))
      h = dynamic_cast<TH1I*>(obj);
    else
      h = dynamic_cast<TH1F*>(obj);
    
    if (!h) {
      std::cerr << "**** getHist1D: <" << hname << "> may not be a 1D Histogram" << std::endl;
      return 0;
    }
    return h;
  }
  
  TH1* getHist1D(const string& hname) {
    return getHist1D(hname.c_str());
  }
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  TH2* getHist2D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      cerr << "**** getHist2D: Histogram for <" << hname << "> not found!" << endl;
      return nullptr;
    }
    
    TH2 *h = nullptr;
    if (obj->InheritsFrom("TH2D"))
      h = dynamic_cast<TH2D*>(obj);
    else if (obj->InheritsFrom("TH2C"))
      h = dynamic_cast<TH2C*>(obj);
    else if (obj->InheritsFrom("TH2S"))
      h = dynamic_cast<TH2S*>(obj);
    else if (obj->InheritsFrom("TH2I"))
      h = dynamic_cast<TH2I*>(obj);
    else
      h = dynamic_cast<TH2F*>(obj);
    
    if (!h) {
      cerr << "**** getHist2D: <<" << hname << ">> may not be a 2D Histogram" << endl;
      return nullptr;
    }
    return h;
  }
  TH2* getHist2D(const string& hname) {
    return getHist2D(hname.c_str());
  }

  TProfile* getHistProfile(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname);
    if (!obj) {
      cerr << "**** getHistProfile: Histogram for <" << hname << "> not found!" << endl;
      return nullptr;
    }
    TProfile *h = nullptr;
    if (obj->InheritsFrom("TProfile"))
      h = dynamic_cast<TProfile*>(obj);

    if (!h) {
      cerr << "**** getHistProfile: <<" << hname << ">> may not be a 2Profile Histogram" << endl;
      return nullptr;
    }
    return h;
  }
  TProfile* getHistProfile(const string& hname) {
    return getHistProfile(hname.c_str());
  }

  // ---------------------------------------------
  // Convenience routine for track cleaning
  // ---------------------------------------------

  void removeTrackDuplicates(const tbeam::TelescopeEvent *telEv, std::vector<tbeam::Track>& tkNoOverlap) {
    //Loop over original track collection and compare all track pairs to check for duplicates
    for(unsigned int i = 0; i<telEv->xPos->size(); i++) {
      double tkX = telEv->xPos->at(i);
      double tkY = telEv->yPos->at(i);
      bool isduplicate = false;
      for (unsigned int j = i+1; j<telEv->xPos->size(); j++) {
        double tkX_j = telEv->xPos->at(j);
        double tkY_j = telEv->yPos->at(j);
        //if (fabs(tkY-tkY_j)<0.015*4 && fabs(tkX-tkX_j)<0.072*4) isduplicate = true;
        //if (fabs(tkY-tkY_j)<0.015 && fabs(tkX-tkX_j)<0.072) isduplicate = true;
        if (fabs(tkY-tkY_j)<0.001 && fabs(tkX-tkX_j)<0.001) isduplicate = true;
      }
      if (!isduplicate) {
        tbeam::Track t(i,tkX,tkY,telEv->dxdz->at(i),telEv->dydz->at(i),telEv->chi2->at(i),telEv->ndof->at(i));  
        tkNoOverlap.push_back(t);
      }      
    }
  }

  void cutTrackFei4Residuals(const tbeam::FeIFourEvent* fei4ev ,const std::vector<tbeam::Track>& tkNoOverlap, std::vector<tbeam::Track>& selectedTk, 
                             const double xResMean, const double yResMean, const double xResPitch, const double yResPitch, bool doClosestTrack) {
   
    double mindelta = 999.;
    double minresx = 999.;
    double minresy = 999.;
    int itkClosest = -1;

    for(unsigned int itk = 0; itk < tkNoOverlap.size(); itk++) {
      const tbeam::Track tTemp(tkNoOverlap.at(itk));
      double tkX = tTemp.xPos;
      double tkY = tTemp.yPos;
#ifdef OCT_16
      tkX = -1.*tTemp.yPos;
      tkY = tTemp.xPos;
#endif
      if (!doClosestTrack){
        minresx = 999.;
        minresy = 999.;
	mindelta = 999.;
      }
      for (unsigned int i = 0; i < fei4ev->col->size(); i++) {
        double xval = 8.375 - (fei4ev->row->at(i)-1)*0.05;
        double yval = 9.875 - (fei4ev->col->at(i)-1)*0.250;
#ifdef OCT_16
        yval = 9.875 - (fei4ev->col->at(i)-1)*0.250;
        xval = 8.375 - (fei4ev->row->at(i)-1)*0.05;
#endif
        double xres = xval - tkX - xResMean;
        double yres = yval - tkY - yResMean;
        if (sqrt(xres*xres+yres*yres)<mindelta){
	  mindelta = sqrt(xres*xres+yres*yres);
	  minresx = xres;
	  minresy = yres;
	  itkClosest = itk;
	}
      }
      if (!doClosestTrack && (std::fabs(minresx) < xResPitch) && (std::fabs(minresy) < yResPitch)) selectedTk.push_back(tTemp);
    }
    if (doClosestTrack && (std::fabs(minresx) < xResPitch) && (std::fabs(minresy) < yResPitch) && itkClosest!=-1){
      const tbeam::Track tClosest(tkNoOverlap.at(itkClosest));
      selectedTk.push_back(tClosest);
    }
  }

  double extrapolateTrackAtDUTwithAngles(const tbeam::Track& track, double FEI4_z, double offset, double zDUT, double theta){

    //Compute distance between DUT center and track impact at DUT along X 
    double xTkAtDUT = track.xPos + (zDUT - FEI4_z) * track.dxdz;
#ifdef OCT_16
    xTkAtDUT = track.yPos + (zDUT - FEI4_z) * track.dydz;
    xTkAtDUT *= -1.;
#endif
    xTkAtDUT = (xTkAtDUT + offset)/ (cos(theta)*(1.-track.dxdz*tan(theta)));
    return xTkAtDUT; 
  }

  std::pair<double, double> extrapolateTrackAtDUTwithAngles(const tbeam::Track& track, double FEI4_z, double offset_d0, double zDUT_d0, double deltaZ, double theta){

    //Compute distance between DUT center and track impact at DUT along X 
    double xTkAtDUT_d0 = track.xPos + (zDUT_d0 - FEI4_z) * track.dxdz;
#ifdef OCT_16
    xTkAtDUT_d0 = track.yPos + (zDUT_d0 - FEI4_z) * track.dydz;
    xTkAtDUT_d0 *= -1.;
#endif
    xTkAtDUT_d0 = (xTkAtDUT_d0 + offset_d0)/ (cos(theta)*(1.-track.dxdz*tan(theta)));

    double zDUT_d1 = zDUT_d0 + deltaZ*cos(theta);
    double xTkAtDUT_d1 = track.xPos + (zDUT_d1 - FEI4_z) * track.dxdz;
#ifdef OCT_16
    xTkAtDUT_d1 = track.yPos + (zDUT_d1 - FEI4_z) * track.dydz;
    xTkAtDUT_d1 *= -1;
#endif
    double offset_d1 = offset_d0 + sin(theta)*deltaZ;
    xTkAtDUT_d1 = (xTkAtDUT_d1 + offset_d1)/ (cos(theta)*(1.-track.dxdz*tan(theta)));

    std::pair<double, double> xTkAtDUT;
    xTkAtDUT.first = xTkAtDUT_d0;
    xTkAtDUT.second = xTkAtDUT_d1;
    return xTkAtDUT;
  }
}
