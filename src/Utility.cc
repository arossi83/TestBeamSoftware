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
          if (i == 11 || i == 13) continue;//only for nov15
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


  void removeTrackDuplicates(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *xTkNoOverlap, std::vector<double> *yTkNoOverlap){

    for(unsigned int i = 0; i<xTk->size(); i++) {
      double tkX = xTk->at(i);
      double tkY = yTk->at(i);
      bool isduplicate = false;
      for (unsigned int j = i+1; j<xTk->size(); j++) {
        double tkX_j = xTk->at(j);
        double tkY_j = yTk->at(j);
        if (fabs(tkY-tkY_j)<0.015*4 && fabs(tkX-tkX_j)<0.072*4) isduplicate = true;
      }
      if (!isduplicate){
        (*xTkNoOverlap).push_back(tkX);
        (*yTkNoOverlap).push_back(tkY);
      }
    }
  }

  void removeTrackDuplicates(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *slopeTk, std::vector<double> *xTkNoOverlap, std::vector<double> *yTkNoOverlap, std::vector<double> *slopeTkNoOverlap){

    for(unsigned int i = 0; i<xTk->size(); i++) {
      double tkX = xTk->at(i);
      double tkY = yTk->at(i);
      double tkSlope = slopeTk->at(i);
      bool isduplicate = false;
      for (unsigned int j = i+1; j<xTk->size(); j++) {
        double tkX_j = xTk->at(j);
        double tkY_j = yTk->at(j);
        if (fabs(tkY-tkY_j)<0.015*4 && fabs(tkX-tkX_j)<0.072*4) isduplicate = true;
      }
      if (!isduplicate){
        (*xTkNoOverlap).push_back(tkX);
        (*yTkNoOverlap).push_back(tkY);
	(*slopeTkNoOverlap).push_back(tkSlope);
      }
    }
  }

  void removeTrackDuplicates(const tbeam::TelescopeEvent *telEv, std::vector<tbeam::Track>& tkNoOverlap) {
    //Loop over original track collection and compare all track pairs to check for duplicates
    for(unsigned int i = 0; i<telEv->xPos->size(); i++) {
      double tkX = telEv->xPos->at(i);
      double tkY = telEv->yPos->at(i);
      bool isduplicate = false;
      for (unsigned int j = i+1; j<telEv->xPos->size(); j++) {
        double tkX_j = telEv->xPos->at(j);
        double tkY_j = telEv->yPos->at(j);
        if (fabs(tkY-tkY_j)<0.015*4 && fabs(tkX-tkX_j)<0.072*4) isduplicate = true;
      }
      if (!isduplicate) {
        tbeam::Track t(i,tkX,tkY,telEv->dxdz->at(i),telEv->dydz->at(i),telEv->chi2->at(i),telEv->ndof->at(i));  
        tkNoOverlap.push_back(t);
      }      
    }
  }

  void cutTrackFei4Residuals(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<int> *colFei4, std::vector<int> *rowFei4, std::vector<double> *xSelectedTk, std::vector<double> *ySelectedTk, double xResMean, double yResMean, double xResPitch, double yResPitch){

    for(unsigned int itk = 0; itk < xTk->size(); itk++) {
      double tkX = -1.*xTk->at(itk);//-1.*telEv()->xPos->at(itk);
      double tkY = yTk->at(itk);//telEv()->yPos->at(itk);
      double minresx = 999.;
      double minresy = 999.;
      for (unsigned int i = 0; i < colFei4->size(); i++) {
        double xval = -9.875 + (colFei4->at(i)-1)*0.250;
        double yval = -8.375 + (rowFei4->at(i)-1)*0.05;
        double xres = xval - tkX - xResMean;//fStepGaus_x->GetParameter(4);//fGausResiduals_x->GetParameter("Mean");
        double yres = yval - tkY - yResMean;//fStepGaus_y->GetParameter(4);//fGausResiduals_y->GetParameter("Mean");
        if(std::fabs(xres) < std::fabs(minresx))   minresx = xres;
        if(std::fabs(yres) < std::fabs(minresy))   minresy = yres;
      }
      //cout << "FEI4 minresx="<<minresx <<" minresy="<<minresy<< endl;
      if(std::fabs(minresx) < 1*xResPitch/2. &&
        std::fabs(minresy) < 1*yResPitch/2.) {
        //cout << "Selected track !"<<endl;
	(*xSelectedTk).push_back(xTk->at(itk));
        (*ySelectedTk).push_back(yTk->at(itk));
      }
    }

    return;
  }

  void cutTrackFei4Residuals(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *slopeTk, std::vector<int> *colFei4, std::vector<int> *rowFei4, std::vector<double> *xSelectedTk, std::vector<double> *ySelectedTk, std::vector<double> *slopeSelectedTk, double xResMean, double yResMean, double xResPitch, double yResPitch){

    for(unsigned int itk = 0; itk < xTk->size(); itk++) {
      double tkX = -1.*xTk->at(itk);//-1.*telEv()->xPos->at(itk);
      double tkY = yTk->at(itk);//telEv()->yPos->at(itk);
      double minresx = 999.;
      double minresy = 999.;
      for (unsigned int i = 0; i < colFei4->size(); i++) {
        double xval = -9.875 + (colFei4->at(i)-1)*0.250;
        double yval = -8.375 + (rowFei4->at(i)-1)*0.05;
        double xres = xval - tkX - xResMean;//fStepGaus_x->GetParameter(4);//fGausResiduals_x->GetParameter("Mean");
        double yres = yval - tkY - yResMean;//fStepGaus_y->GetParameter(4);//fGausResiduals_y->GetParameter("Mean");
        if(std::fabs(xres) < std::fabs(minresx))   minresx = xres;
        if(std::fabs(yres) < std::fabs(minresy))   minresy = yres;
      }
      //cout << "FEI4 minresx="<<minresx <<" minresy="<<minresy<< endl;
      if(std::fabs(minresx) < 1*xResPitch/2. &&
         std::fabs(minresy) < 1*yResPitch/2.) {
        //cout << "Selected track !"<<endl;
        (*xSelectedTk).push_back(xTk->at(itk));
        (*ySelectedTk).push_back(yTk->at(itk));
        (*slopeSelectedTk).push_back(slopeTk->at(itk));
      }
    }
  }

  void cutTrackFei4Residuals(const tbeam::FeIFourEvent* fei4ev ,const std::vector<tbeam::Track>& tkNoOverlap, std::vector<tbeam::Track>& selectedTk, 
                             const double xResMean, const double yResMean, const double xResPitch, const double yResPitch) {
    
    for(unsigned int itk = 0; itk < tkNoOverlap.size(); itk++) {
      const tbeam::Track tTemp(tkNoOverlap.at(itk));
      double tkX = -1.*tTemp.xPos;//-1.*telEv()->xPos->at(itk);
      double tkY = tTemp.yPos;//telEv()->yPos->at(itk);
      double minresx = 999.;
      double minresy = 999.;
      for (unsigned int i = 0; i < fei4ev->col->size(); i++) {
        double xval = -9.875 + (fei4ev->col->at(i)-1)*0.250;
        double yval = -8.375 + (fei4ev->row->at(i)-1)*0.05;
        double xres = xval - tkX - xResMean;//fStepGaus_x->GetParameter(4);//fGausResiduals_x->GetParameter("Mean");
        double yres = yval - tkY - yResMean;//fStepGaus_y->GetParameter(4);//fGausResiduals_y->GetParameter("Mean");
        if(std::fabs(xres) < std::fabs(minresx))   minresx = xres;
        if(std::fabs(yres) < std::fabs(minresy))   minresy = yres;
      }
      //cout << "FEI4 minresx="<<minresx <<" minresy="<<minresy<< endl;
      if(std::fabs(minresx) < 1*xResPitch/2. && std::fabs(minresy) < 1*yResPitch/2.) {
        selectedTk.push_back(tTemp);    
      }
    }
  }


}